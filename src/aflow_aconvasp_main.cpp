// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo
// Dane Morgan (up to 2003)
// Wahyu Setyawan (up to 2009)

#ifndef _AFLOW_PFLOW_MAIN_CPP_
#define _AFLOW_PFLOW_MAIN_CPP_

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iostream>
#include <istream>
#include <ostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <stdlib.h>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_argv.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xcombos.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xhttp.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xparser.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xscalar.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow.h"
#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "aflow_init.h"
#include "aflow_xhost.h"
#include "aflowlib/aflowlib_database.h"
#include "aflowlib/aflowlib_libraries.h"
#include "aflowlib/aflowlib_web_interface.h"
#include "flow/aflow_avasp.h"
#include "flow/aflow_bader.h"
#include "flow/aflow_ivasp.h"
#include "flow/aflow_kvasp.h"
#include "flow/aflow_pflow.h"
#include "flow/aflow_support_types.h"
#include "flow/aflow_xclasses.h"
#include "interfaces/aflow_pthreads.h"
#include "modules/APL/aflow_apl.h"  //ME20200330
#include "modules/CCE/aflow_cce.h" //RF20200203
#include "modules/COMPARE/aflow_compare_structure.h" //DX20181023
#include "modules/GFA/aflow_gfa.h" //DF20190329
#include "modules/HULL/aflow_chull.h"
#include "modules/POCC/aflow_pocc.h"  //CO20181226
#include "modules/POCC/aflow_pocc_old.h"
#include "modules/PROTOTYPES/aflow_anrl.h"
#include "modules/QCA/aflow_qca.h" //SD20220323
#include "modules/SYM/aflow_spacegroup.h"
#include "modules/SYM/aflow_symmetry.h"
#include "modules/SYM/aflow_symmetry_spacegroup.h"
#include "modules/SYM/aflow_wyckoff.h"
#include "structure/aflow_defects.h"
#include "structure/aflow_lattice.h"
#include "structure/aflow_surface.h"
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

using std::cerr;
using std::cin;
using std::cout;
using std::deque;
using std::endl;
using std::ifstream;
using std::istream;
using std::istringstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setprecision;
using std::setw;
using std::stringstream;
using std::vector;

using aurostd::xmatrix;
using aurostd::xvector;

extern double NearestNeighbor(const xstructure& a);

// ***************************************************************************

uint PflowARGs(vector<string>& argv, vector<string>& cmds, aurostd::xoption& vpflow) {  // vpflow is really XHOST.vflag_pflow
  const bool LDEBUG = (false || XHOST.DEBUG);
  // GENERAL STUFF

  vpflow.flag("PFLOW_HELP", aurostd::args2flag(argv, cmds, "--HELP|--help"));
  vpflow.flag("PROTOS", (aurostd::args2flag(argv, cmds, "--protos|--prototypes") || aurostd::args2flag(argv, cmds, "--proto")));// && (argv.size()==2));

  vpflow.args2addattachedscheme(argv, cmds, "PROTOS_ICSD", "--protos_icsd=|--prototypes_icsd=", "");

  if (!vpflow.flag("PROTOS_ICSD")) {
    vpflow.flag("PROTOS_ICSD", aurostd::args2flag(argv, cmds, "--protos_icsd|--prototypes_icsd"));
  }

  //  bool XXX=aurostd::args2flag(argv,cmds,"--xxx") && argv.at(1)=="--xxx";
  //  bool YYY=aurostd::args2flag(argv,cmds,"--yyy");

  vpflow.flag("ACE", aurostd::args2flag(argv, cmds, "--ace"));
  // DX20170818 - Added tolerance and no_scan options to Xgroups - START
  vpflow.args2addattachedscheme(argv, cmds, "AGROUP", "--sitepointgroup=|--agroup=", "");
  if (vpflow.flag("AGROUP")) {
    vpflow.flag("SYMMETRY::NO_SCAN", aurostd::args2flag(argv, cmds, "--no_scan"));
    if (aurostd::args2attachedflag(argv, "--sitepointgroup=|--agroup=")) { // DX20170803
      vpflow.args2addattachedscheme(argv, cmds, "SYMMETRY::TOLERANCE", "--sitepointgroup=|--agroup=", ""); // DX20200907 - default is system specific, leaving empty
    }
    vpflow.flag("SYMMETRY::SCREEN_ONLY", aurostd::args2flag(argv, cmds, "--screen_only")); // DX20170803
    // DX20170921 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv, cmds, "SYMMETRY::MAGNETIC", "--mag=|--magnetic=|--magmom=", ""); // DX20170803
    // DX20170921 - MAGNETIC SYMMETRY - END
    //  ME20210206 - web mode
    if (XHOST.vflag_control.flag("WWW")) {
      vpflow.flag("SYMMETRY::SCREEN_ONLY", true);
      XHOST.QUIET = true;
    }
  }
  // DX20170818 - Added tolerance and no_scan options to Xgroups - END
  vpflow.flag("AGROUP2", aurostd::args2flag(argv, cmds, "--sitepointgroup2|--agroup2"));
  vpflow.flag("AGROUP2m", aurostd::args2flag(argv, cmds, "--sitepointgroup2m|--agroup2m"));
  vpflow.flag("AFLOWIN", aurostd::args2flag(argv, cmds, "--aflowin"));
  vpflow.flag("ALPHABETIC", aurostd::args2flag(argv, cmds, "--alpha|--alphabetic"));
  vpflow.args2addattachedscheme(argv, cmds, "ALPHA_COMPOUND", "--alpha_compound=|--alpha_compounds=", "");
  vpflow.args2addattachedscheme(argv, cmds, "ALPHA_SPECIES", "--alpha_species=|--alpha_specie=", "");
  vpflow.args2addattachedscheme(argv, cmds, "ANGLES", "--angle=|--angles=", "0.0");

  vpflow.args2addattachedscheme(argv, cmds, "AFLOWLIB::ENTRY_JSON", "--aflowlib2=|--aflowlib=", "");  // SC20190812
  vpflow.args2addattachedscheme(argv, cmds, "AFLOWLIB_AUID2AURL", "--aflowlib_auid2aurl=|--auid2aurl=", "");
  vpflow.args2addattachedscheme(argv, cmds, "AFLOWLIB_AURL2AUID", "--aflowlib_aurl2auid=|--aurl2auid=", "");
  vpflow.args2addattachedscheme(argv, cmds, "AFLOWLIB_AUID2LOOP", "--aflowlib_auid2loop=|--auid2loop=", "");
  vpflow.args2addattachedscheme(argv, cmds, "AFLOWLIB_AURL2LOOP", "--aflowlib_aurl2loop=|--aurl2loop=", "");

  vpflow.flag("AFLOWSYM_PYTHON", aurostd::args2attachedflag(argv, cmds, "--aflow_sym_python|--aflowsym_python")); // DX20210202

  // DX20190206 - add AFLUX functionality to command line - START
  vpflow.args2addattachedscheme(argv, cmds, "AFLUX::SUMMONS", "--aflux=", "");  // CO20200520 - AFLUX::SUMMONS
  if (vpflow.flag("AFLUX::SUMMONS")) {  // CO20200520 - AFLUX::SUMMONS
    vpflow.flag("AFLUX::USAGE", aurostd::args2flag(argv, cmds, "--usage"));
  }
  // DX20190206 - add AFLUX functionality to command line - END
  vpflow.flag("ANALYZEDB", aurostd::args2flag(argv, cmds, "--analyze_database"));  // ME20191001

  // Commands for serializing bands and DOS data to JSON
  vpflow.args2addattachedscheme(argv, cmds, "DOSDATA2JSON", "--dosdata2json=", "./"); // EG
  if (vpflow.flag("DOSDATA2JSON")) {
    vpflow.args2addattachedscheme(argv, cmds, "DOSDATA2JSON::PARAMS", "--dos_parameters=", "");
  }  // CO20180214 removed params confusion
  vpflow.args2addattachedscheme(argv, cmds, "BANDSDATA2JSON", "--bandsdata2json=", "./"); // EG
  // End commands

  vpflow.flag("STRUCTURE2JSON", aurostd::args2flag(argv, cmds, "--structure2json|--struct2json")); // DX20190508

  // CO
  vpflow.flag("BADER", aurostd::args2flag(argv, cmds, "--bader"));
  if (vpflow.flag("BADER")) {    // CO
    // manage directory, ASSUME LOCAL IF NOT SPECIFIED
    if (XHOST.vflag_control.flag("DIRECTORY")) {
      vpflow.push_attached("BADER::DIRECTORY", XHOST.vflag_control.getattachedscheme("DIRECTORY"));
    } else {
      vpflow.push_attached("BADER::DIRECTORY", ".");
    }
    // XHOST.vflag_control.getattachedscheme("DIRECTORY"); //this is not doing anything //CO20190520
    vpflow.flag("BADER::USAGE", aurostd::args2flag(argv, cmds, "--usage"));  // usage
    vpflow.flag("BADER::CRITICAL_POINTS", aurostd::args2flag(argv, cmds, "--critical_points|--cp"));  // critical points
    vpflow.args2addattachedscheme(argv, cmds, "BADER::CALCULATE", "--calculate=|--calc=", "");    // -c
    vpflow.args2addattachedscheme(argv, cmds, "BADER::NOCALCULATE", "--nocalculate=|--nocalc=", "");  // -n
    vpflow.args2addattachedscheme(argv, cmds, "BADER::PARTITION", "--partition=|--part=", "");  // -b
    vpflow.args2addattachedscheme(argv, cmds, "BADER::REFINE_EDGE_METHOD", "--refine_edge_method=|--rem=|--r=", "");  // -r
    vpflow.args2addattachedscheme(argv, cmds, "BADER::REFERENCE", "--reference=|--ref=", "");     // -ref
    vpflow.args2addattachedscheme(argv, cmds, "BADER::VACUUM", "--vacuum=|--vac=", "");   // -vac
    vpflow.args2addattachedscheme(argv, cmds, "BADER::TERMINATE", "--terminate=|--term=", "");    // -m
    vpflow.args2addattachedscheme(argv, cmds, "BADER::PRINT_ALL", "--print_all=", "");    // -p all_atom | -p all_bader
    vpflow.args2addattachedscheme(argv, cmds, "BADER::PRINT_INDEX", "--print_index=|--print_idx=", "");     // -p sel_atom | -p sel_bader
    vpflow.args2addattachedscheme(argv, cmds, "BADER::PRINT_SELECT_ATOM", "--print_select_atom=|--print_sel_atom=", "");    // -p sel_atom
    vpflow.args2addattachedscheme(argv, cmds, "BADER::PRINT_SELECT_BADER", "--print_select_bader=|--print_sel_bader=", ""); // -p sel_bader
    vpflow.args2addattachedscheme(argv, cmds, "BADER::PRINT_SUM_ATOM", "--print_sum_atom=", "");  // -p sum_atom
    vpflow.args2addattachedscheme(argv, cmds, "BADER::PRINT_SUM_BADER", "--print_sum_bader=", "");    // -p sum_bader
    vpflow.flag("BADER::QUIET", aurostd::args2flag(argv, cmds, "--quiet|--q"));  // quiets output
    vpflow.flag("BADER::CONSOLIDATE_ATOMS2SPECIES", aurostd::args2flag(argv, cmds, "--consolidate_atoms2species|--a2s"));  // quiets output
    if (vpflow.flag("BADER::CONSOLIDATE_ATOMS2SPECIES")) {
      // we can consolidate and delete summing chgcars
      vpflow.flag("BADER::REMOVE_BADER_ATOMS", aurostd::args2flag(argv, cmds, "--remove_bader_atoms|--rba"));  // keep only consolidated files
    }
    vpflow.args2addattachedscheme(argv, cmds, "BADER::JVXL_ALL_SPECIES", "--jvxl_all_species=|--jvxl=", ""); // allows jvxl output for atoms2species
    if (vpflow.flag("BADER::JVXL_ALL_SPECIES")) {
      // when jvxl-ing, we can just remove summing chgcars and keep full consolidated chgcar, or delete all and just keep jvxl
      vpflow.flag("BADER::REMOVE_BADER_ATOMS", aurostd::args2flag(argv, cmds, "--remove_bader_atoms|--rba"));  // keep only consolidated files
      vpflow.flag("BADER::KEEP::JVXL_ONLY", aurostd::args2flag(argv, cmds, "--keep=jvxl_only|--keep_jvxl_only|--jvxl_only"));  // removes consolidated chgcar_files, keeps only jvxl
    }
  }

  vpflow.args2addattachedscheme(argv, cmds, "BANDS", "--bands=", "");

  vpflow.flag("BANDSTRUCTURE", aurostd::args2flag(argv, cmds, "--bandstructure|--bandsstructures|--bands_structures|--band_structures|--bands_structure|--band_structure|--bs"));
  vpflow.flag("BZPLOT", aurostd::args2flag(argv, cmds, "--plotbz|--bzplot"));
  vpflow.args2addattachedscheme(argv, cmds, "BZPLOTUSEKPOINTS", "--bzplotuseKPOINTS=|--bzplotdusekpoints=", "");
  vpflow.flag("BZPLOTDATA", aurostd::args2flag(argv, cmds, "--bzplotdata"));
  vpflow.args2addattachedscheme(argv, cmds, "BZPLOTDATAUSEKPOINTS", "--bzplotdatauseKPOINTS=|--bzplotdatausekpoints=", "");

  vpflow.args2addattachedscheme(argv, cmds, "FIX_BANDS", "--fix_bands=", "");

  // DX+CO START

  vpflow.args2addattachedscheme(argv, cmds, "FULLSYMMETRY", "--aflow-sym=|--AFLOW-SYM=|--AFLOWSYM=|--aflowSYM=|--aflowsym=|--full_symmetry=|--full_sym=|--fullsym=", ""); // DX20170803 Added other aliases
  if (vpflow.flag("FULLSYMMETRY")) {
    vpflow.flag("FULLSYMMETRY::NO_SCAN", aurostd::args2flag(argv, cmds, "--no_scan"));
    if (aurostd::args2attachedflag(argv, "--aflow-sym=|--AFLOW-SYM=|--AFLOWSYM=|--aflowSYM=|--aflowsym=|--full_symmetry=|--full_sym=|--fullsym=")) { // DX20170803
      vpflow.args2addattachedscheme(argv, cmds, "FULLSYMMETRY::TOLERANCE",
                                    "--aflow-sym=|--AFLOW-SYM=|--AFLOWSYM=|--aflowSYM=|--aflowsym=|--full_symmetry=|--full_sym=|--fullsym=", ""); // DX20170803 //DX20200907 - default is system specific, leaving empty
    }
    vpflow.flag("FULLSYMMETRY::SCREEN_ONLY", aurostd::args2flag(argv, cmds, "--screen_only")); // DX20170803
    // DX20170921 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv, cmds, "FULLSYMMETRY::MAGNETIC", "--mag=|--magnetic=|--magmom=", ""); // DX20170803
    // DX20170921 - MAGNETIC SYMMETRY - END
  }
  // DX+CO END

  vpflow.flag("BANDGAP_WAHYU", aurostd::args2flag(argv, cmds, "--bandgap_wahyu|--gap_wahyu"));
  vpflow.args2addattachedscheme(argv, cmds, "BANDGAP", "--bandgap=", "./");

  vpflow.flag("BANDGAPS", aurostd::args2flag(argv, cmds, "--bandgaps|--gaps"));
  //[CO20191004]vpflow.flag("BANDGAPDOS",aurostd::args2flag(argv,cmds,"--bandgap_from_dos|--bandgap_from_DOS|--bandgap_dos|--bandgapdos")); //CO20191004
  vpflow.args2addattachedscheme(argv, cmds, "BANDGAPDOS", "--bandgap_from_dos=|--bandgap_from_DOS=|--bandgap_dos=|--bandgapdos=", "./"); // CO20191004
  vpflow.flag("BANDGAPLISTDOS", aurostd::args2flag(argv, cmds, "--bandgaplist_from_dos|--bandgaplist_from_DOS"));
  vpflow.args2addattachedscheme(argv, cmds, "BZDIRECTION", "--bzdirection=|--bzdirections=|--bzd=", "");
  // DX20181102 - add transform2original option - START
  if (vpflow.flag("BZDIRECTION")) {
    if (aurostd::args2attachedflag(argv, "--bzdirection=|--bzdirections=|--bzd=")) {
      vpflow.args2addattachedscheme(argv, cmds, "BZDIRECTION::LATTICE", "--bzdirection=|--bzdirections=|--bzd=", "1");
    }
    vpflow.flag("BZDIRECTION::TRANSFORM2ORIGINAL", aurostd::args2flag(argv, cmds, "--transform2original"));
    vpflow.flag("BZDIRECTION::PRINT_TRANSFORMATION_MATRIX", aurostd::args2flag(argv, cmds, "--print_transformation_matrix"));
  }
  // DX20181102 - add transform2original option - END
  vpflow.flag("BZMAX", aurostd::args2flag(argv, cmds, "--BZmax"));

  vpflow.args2addattachedscheme(argv, cmds, "CAGES", "--cages=", "-1.0");
  // cerr << "vpflow.flag(\"CAGES\")=" << vpflow.flag("CAGES") << endl;
  // cerr << "vpflow.getattachedscheme(\"CAGES\")=" << vpflow.getattachedscheme("CAGES") << endl;

  vpflow.flag("CALCULATED_ICSD_RANDOM", (aurostd::args2flag(argv, cmds, "--calculated=icsd")) && aurostd::args2flag(argv, cmds, "--random|--rnd"));
  vpflow.args2addattachedscheme(argv, cmds, "CALCULATED", "--calculated=", "all");
  // cerr << "vpflow.flag(\"CALCULATED\")=" << vpflow.flag("CALCULATED") << endl;
  // cerr << vpflow.getattachedscheme("CALCULATED") << endl;

  vpflow.flag("CART", aurostd::args2flag(argv, cmds, "--cart|-cart|-c|--cartesian"));

  vpflow.flag("CCE_CORRECTION::USAGE", aurostd::args2flag(argv, cmds, "--cce_correction|--cce"));
  vpflow.args2addattachedscheme(argv, cmds, "CCE_CORRECTION::POSCAR_PATH", "--cce_correction=|--cce=", "");
  vpflow.flag("CCE_CORRECTION::POSCAR2CCE", aurostd::args2flag(argv, cmds, "--poscar2cce")); // ME
  vpflow.flag("CCE_CORRECTION::GET_CCE_CORRECTION", aurostd::args2flag(argv, cmds, "--get_cce_correction|--get_cce_cor|--get_cce_corrections|--get_cce_cors"));
  vpflow.flag("CCE_CORRECTION", vpflow.flag("CCE_CORRECTION::USAGE") || !vpflow.getattachedscheme("CCE_CORRECTION::POSCAR_PATH").empty());
  if (vpflow.flag("CCE_CORRECTION") && aurostd::args2flag(argv, cmds, "--usage")) {
    vpflow.flag("CCE_CORRECTION::USAGE", true);
  }
  vpflow.args2addattachedscheme(
      argv, cmds, "CCE_CORRECTION::ENTHALPIES_FORMATION_DFT",
      "--enthalpies_formation_dft=|--enthalpy_formation_dft=|--dft_formation_enthalpies=|--dft_formation_energies=|--dft_formation_enthalpy=|--dft_formation_energy=|--dftes=|--dfte=", "");
  vpflow.args2addattachedscheme(argv, cmds, "CCE_CORRECTION::FUNCTIONALS", "--functional=|--func=|--functionals=|--funcs=", "");
  vpflow.args2addattachedscheme(argv, cmds, "CCE_CORRECTION::OXIDATION_NUMBERS", "--oxidation_numbers=|--ox_nums=|--oxidation_number=|--ox_num=", "");
  vpflow.args2addattachedscheme(argv, cmds, "CCE_CORRECTION::PRINT", "--print=", "OUT"); // ME
  vpflow.flag("CCE_CORRECTION::UNIT_TEST", aurostd::args2flag(argv, cmds, "--cce_test")); // RF20200409
  vpflow.flag("CCE_CORRECTION::GET_OXIDATION_NUMBERS", aurostd::args2flag(argv, cmds, "--get_oxidation_numbers|--get_ox_nums|--get_oxidation_number|--get_ox_num|--poscar2ox_nums|--poscar2ox_num")); // RF20200725
  vpflow.flag("CCE_CORRECTION::GET_CATION_COORDINATION_NUMBERS", aurostd::args2flag(argv, cmds,
                                                                                    "--get_cation_coordination_numbers|--get_cation_coord_nums|--get_cation_coordination_number|--get_cation_coord_num|--get_"
                                                                                    "coordination_numbers_cation|--get_coordination_number_cation|--get_coordination_numbers_cations|--get_coordination_number_"
                                                                                    "cations|--get_coord_num_cation|--get_coord_nums_cation|--get_coord_nums_cations|--get_coord_num_cations|--poscar2cation_"
                                                                                    "coord_nums|--poscar2cation_coord_num")); // RF20200814
  vpflow.args2addattachedscheme(argv, cmds, "CCE_CORRECTION::DIST_TOL", "--tolerance=|dist_tol=|distance_tolerance=|dist_tolerance=|distance_tol=", ""); // RF20200819

  vpflow.flag("CHECKINTEGRITIY", aurostd::args2flag(argv, cmds, "--check_integrity|--checki"));
  vpflow.args2addattachedscheme(argv, cmds, "CHANGESUFFIX", "--suffix=", "./");

  // CO
  vpflow.args2addattachedscheme(argv, cmds, "CHGCAR2JVXL", "--chgcar2jvxl=|--c2j=", ""); // create JVXL from CHGCAR, CO
  if (vpflow.flag("CHGCAR2JVXL")) {
    vpflow.flag("CHGCAR2JVXL::USAGE", aurostd::args2flag(argv, cmds, "--usage")); // usage
    vpflow.args2addattachedscheme(argv, cmds, "CHGCAR2JVXL::OUTPUT", "--output=|--o=", "");
  }

  // CO
  vpflow.args2addattachedscheme(argv, cmds, "CHGDIFF", "--chgdiff=", "");
  if (vpflow.flag("CHGDIFF")) { // CO
    vpflow.flag("CHGDIFF::USAGE", aurostd::args2flag(argv, cmds, "--usage")); // usage
    vpflow.args2addattachedscheme(argv, cmds, "CHGDIFF:OUTPUT", "--output=|--o=", "");
  }

  vpflow.flag("CHGINT", aurostd::args2flag(argv, cmds, "--chgint") && argv.at(1) == "--chgint");
  // CO
  vpflow.args2addattachedscheme(argv, cmds, "CHGSUM", "--chgsum=", "");
  if (vpflow.flag("CHGSUM")) { // CO
    vpflow.flag("CHGSUM::USAGE", aurostd::args2flag(argv, cmds, "--usage")); // usage
    vpflow.args2addattachedscheme(argv, cmds, "CHGSUM::OUTPUT", "--output=|--o=", "");
  }

  vpflow.flag("CLEANALL", aurostd::args2flag(argv, cmds, "--cleanall|--clean_all"));

  vpflow.args2addattachedscheme(argv, cmds, "COMPARE", "--compare=", "");
  vpflow.flag("CMPSTR", aurostd::args2flag(argv, cmds, "--cmp_str") && argv.at(1) == "--cmp_str");

  vpflow.flag("CHULL::INIT", aurostd::args2flag(argv, cmds, "--convex_hull|--chull")); // initiate chull calculation
  vpflow.args2addattachedscheme(argv, cmds, "PFLOW::ALLOY", "--alloy=", ""); // define alloy
  vpflow.flag("PFLOW::LOAD_API", aurostd::args2flag(argv, cmds, "--load_API|--load_api|--loadapi|--lapi|--api")); // force load api
  vpflow.flag("PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT", aurostd::args2flag(argv, cmds, "--load_entries_entry_output|--loadentriesentryoutput|--leo")); // verbose loading entries output
  if (vpflow.flag("CHULL::INIT")) {
    // do NOT rearrange flags, order is important!
    vpflow.args2addattachedscheme(argv, cmds, "CHULL::PATH", "--destination=|--path=", ""); // determine how to get output
    vpflow.flag("CHULL::USAGE", aurostd::args2flag(argv, cmds, "--usage")); // usage
    vpflow.flag("CHULL::GNUPLOT_DOC", false); // depreciated, NO gnuplot options, all LATEX
    vpflow.flag("CHULL::SCREEN_ONLY", aurostd::args2flag(argv, cmds, "--screen_only")); // print to screen
    if (vpflow.flag("CHULL::SCREEN_ONLY")) {
      XHOST.QUIET = true;
    }
    if (XHOST.QUIET && vpflow.flag("PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT")) {
      vpflow.flag("PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT", false);
    } // keep quiet
    vpflow.args2addattachedscheme(argv, cmds, "CHULL::KEEP", "--keep=", "");
    vector<string> keep_parts;
    if (vpflow.flag("CHULL::KEEP")) {
      aurostd::string2tokens(vpflow.getattachedscheme("CHULL::KEEP"), keep_parts, ",");
    }
    vpflow.flag("CHULL::LOG", aurostd::args2flag(argv, cmds, "--keep=log|--keep_log|--keeplog|--log"));
    if (!vpflow.flag("CHULL::LOG")) { // only turn on if not already on
      for (size_t i = 0; i < keep_parts.size(); i++) {
        if (!keep_parts[i].empty()) {
          if (keep_parts[i][0] == 'L' || keep_parts[i][0] == 'l') {
            vpflow.flag("CHULL::LOG", true);
          }
        }
      }
    }
    vpflow.args2addattachedscheme(argv, cmds, "PFLOW::LOAD_LIBRARY", "--load_library=|--loadlibrary=|--ll=", ""); // get libraries to load
    vpflow.args2addattachedscheme(argv, cmds, "CHULL::NEGLECT", "--neglect=|--ban=", ""); // remove by auid from chull calculation
    vpflow.flag("CHULL::SEE_NEGLECT", aurostd::args2flag(argv, cmds, "--see_neglect|--seeneglect|--sn")); // see why compounds get neglected
    vpflow.args2addattachedscheme(argv, cmds, "CHULL::REMOVE_EXTREMA", "--remove_extreme_points=|--removeextremepoints=|--remove_extrema=|--removeextrema=|--rep=",
                                  ""); // set threshold, from bottom for enthalpy of formation, from top for entropic temperature, units of meV OR K
    vpflow.flag("CHULL::ENTROPIC_TEMPERATURE", aurostd::args2flag(argv, cmds, "--entropic_temperature|--entropictemperature|--entroptemp")); // use entropic temperature instead of enthalpy of formation
    vpflow.flag("CHULL::INCLUDE_PAW_GGA", aurostd::args2flag(argv, cmds, "--include_paw_gga|--paw_gga")); // include entries calculated with PAW_GGA
    vpflow.flag("CHULL::SKIP_STRUCTURE_COMPARISON", aurostd::args2flag(argv, cmds, "--skip_structure_comparison|--skipstructruecomparison|--skipstructcomp|--ssc")); // use entropic temperature instead of enthalpy of formation
    vpflow.flag("CHULL::SKIP_STABILITY_CRITERION_ANALYSIS", aurostd::args2flag(argv, cmds,
                                                                               "--skip_stability_criterion_analysis|--skip_stability_criterion|--skipstabilitycriterionanalysis|--skipstabilitycriterion|--skip_"
                                                                               "scriterion|--skipscriterion|--sscriterion")); // use entropic temperature instead of enthalpy of formation
    vpflow.flag("CHULL::SKIP_N+1_ENTHALPY_GAIN_ANALYSIS", aurostd::args2flag(argv, cmds,
                                                                             "--skip_n_plus_1_enthalpy_gain_analysis|--skip_n_plus_1_energy_gain_analysis|--skipnplus1enthalpygainanalysis|--"
                                                                             "skipnplus1energygainanalysis|--skip_nplus1|--skipnplus1|--snp1|--snpo")); // use entropic temperature instead of enthalpy of formation
    vpflow.flag("CHULL::INCLUDE_SKEWED_HULLS", aurostd::args2flag(argv, cmds, "--include_skewed_hulls|--include_skewed|--ish")); // use entropic temperature instead of enthalpy of formation
    vpflow.flag("CHULL::INCLUDE_UNRELIABLE_HULLS", aurostd::args2flag(argv, cmds, "--include_unreliable_hulls|--include_unreliable|--iuh")); // use entropic temperature instead of enthalpy of formation
    vpflow.flag("CHULL::INCLUDE_OUTLIERS", aurostd::args2flag(argv, cmds, "--include_outliers|--io")); // use entropic temperature instead of enthalpy of formation
    vpflow.flag("CHULL::STRICT_OUTLIER_ANALYSIS", aurostd::args2flag(argv, cmds, "--strict_outlier_analysis|--soa"));
    vpflow.flag("CHULL::INCLUDE_ILL_CONVERGED", aurostd::args2flag(argv, cmds, "--include_ill_converged|--iic")); // use entropic temperature instead of enthalpy of formation
    vpflow.flag("CHULL::LATEX_OUTPUT", aurostd::args2flag(argv, cmds, "--latex_output|--latexoutput")); // verbose latex output (cout and FileMESSAGE)
    vpflow.flag("CHULL::LATEX_INTERACTIVE", aurostd::args2flag(argv, cmds, "--latex_interactive|--latexinteractive")); // interact with latex (execute vs. execute2string)
    if (vpflow.flag("CHULL::LATEX_INTERACTIVE") && !vpflow.flag("CHULL::LATEX_OUTPUT")) {
      vpflow.flag("CHULL::LATEX_OUTPUT", true);
    } // keep verbose latex output
    vpflow.flag("CHULL::PLOT_ISO_MAX_LATENT_HEAT", aurostd::args2flag(argv, cmds, "--plot_iso_max_latent_heat|--plot_isomax|--iso_max|--isomax")); // plot iso-max latent heat lines for hull points
    // need to determine output type immediately, so we can accept/ignore output specific flags (SAFE)
    vpflow.args2addattachedscheme(argv, cmds, "CHULL::OUTPUT", "--print=|--p=|--output=|--o=", ""); // determine how to get output, could pull from XHOST, but --output was used previously (backwards compatibility)
    if (vpflow.flag("CHULL::OUTPUT")) {
      vector<string> out_forms;
      string out_form;
      aurostd::string2tokens(vpflow.getattachedscheme("CHULL::OUTPUT"), out_forms, ",");
      if (out_forms.size() > 1) {
        vpflow.flag("CHULL::MULTI_OUTPUT", true);
      } // CO20180409
      for (size_t i = 0; i < out_forms.size(); i++) {
        out_form = aurostd::toupper(out_forms[i]);
        if (out_form[0] == 'A') {
          vpflow.flag("CHULL::APOOL_OUT", true); // turn on
        } else if (out_form[0] == 'F') {
          vpflow.flag("CHULL::TEXT_DOC", true); // turn on
          vpflow.flag("CHULL::JSON_DOC", true); // turn on
          vpflow.flag("CHULL::LATEX_DOC", true); // turn on default
        } else if (out_form[0] == 'T') {
          vpflow.flag("CHULL::TEXT_DOC", true); // turn on
        } else if (out_form[0] == 'J') { // MB20190305
          if (out_form == "JUPYTER2") {
            vpflow.flag("CHULL::WRITE_JUPYTER2", true);
          } else if (out_form == "JUPYTER3") {
            vpflow.flag("CHULL::WRITE_JUPYTER3", true);
          } else if (out_form == "JUPYTER") {
            vpflow.flag("CHULL::WRITE_JUPYTER3", true);
          } else {
            vpflow.flag("CHULL::JSON_DOC", true);
          } // turn on
        } else if (out_form[0] == 'W') {
          vpflow.flag("CHULL::WEB_DOC", true); // turn on
          //[WSCHMITT20190620 - stability criterion error otherwise for Pd1]vpflow.flag("CHULL::SKIP_STRUCTURE_COMPARISON",true); //the app cannot handle more than one g-state in the visualization
          vpflow.flag("FORCE", true); // just include everything!
          XHOST.vflag_control.flag("WWW", true); // CO20201215
          // vpflow.flag("CHULL::INCLUDE_UNRELIABLE",true); //we include colors on the website
        } else if (out_form[0] == 'L' || (out_form[0] == 'P' && out_form != "PNG") || out_form == "PDF") { // Latex or Pdf
          vpflow.flag("CHULL::LATEX_DOC", true); // turn on default
        } else if (out_form == "PNG") { // PNG
          vpflow.flag("CHULL::PNG_IMAGE", true);
        } // else { //deal with later //[CO20200106 - close bracket for indenting]}
      }
      // fix CHULL::MULTI_OUTPUT, PDF + PNG not really considered "different" outputs
      if (out_forms.size() == 2 && vpflow.flag("CHULL::LATEX_DOC") && vpflow.flag("CHULL::PNG_IMAGE")) {
        vpflow.flag("CHULL::MULTI_OUTPUT", false);
      }
    } else {
      vpflow.flag("CHULL::LATEX_DOC", true);
    } // default
    vpflow.args2addattachedscheme(argv, cmds, "CHULL::DIST2HULL", "--distance_to_hull=|--distancetohull=|--distance2hull=|--dist2hull=|--d2h=", ""); // calculate distance to hull for point
    vpflow.args2addattachedscheme(argv, cmds, "CHULL::STABILITY_CRITERION", "--stability_criterion=|--stabilitycriterion=|--stable_criterion=|--scriterion=|--sc=", ""); // calculate stable criterion for point
    vpflow.args2addattachedscheme(argv, cmds, "CHULL::N+1_ENTHALPY_GAIN",
                                  "--n+1_enthalpy_gain=|--n+1_energy_gain=|--n+1enthalpygain=|--n+1energygain=|--n+1egain=|--n1egain=|--n+1_enthalpygain=|--n+1+energygain=|--n+1_egain=|--nplus1=", ""); // calculate stable criterion for point
    vpflow.args2addattachedscheme(argv, cmds, "CHULL::CALCULATE_FAKE_HULL_STABILITY_CRITERION", "--fake_hull_sc=", ""); // CO20210315
    vpflow.flag("CHULL::CALCULATE_FAKE_HULL_N+1_ENTHALPY_GAIN", aurostd::args2flag(argv, cmds, "--fake_hull_np1eg")); // SK20200325 - skip all n-dimensional points and calculate new hull
    vpflow.flag("CHULL::CALCULATE_HIGHEST_DIMENSION_ONLY", aurostd::args2flag(argv, cmds, "--calculate_highest_dimension_only|--calc_ND_only")); // CO20210407
    vpflow.args2addattachedscheme(argv, cmds, "CHULL::HULL_FORMATION_ENTHALPY", "--hull_formation_enthalpy=|--hull_enthalpy=|--hull_energy=", ""); // calculate stable criterion for point
    if (vpflow.flag("CHULL::STABILITY_CRITERION") || vpflow.flag("CHULL::N+1_ENTHALPY_GAIN") || vpflow.flag("CHULL::HULL_FORMATION_ENTHALPY")) {
      // vpflow.flag("CHULL::TEXT_DOC",false);   //turn off  //leave on, as user might request json/text format output
      // vpflow.flag("CHULL::JSON_DOC",false);   //turn off  //leave on, as user might request json/text format output
      //[CO20210201 - chull-web SS plotter]vpflow.flag("CHULL::WEB_DOC",false);    //turn off
      vpflow.flag("CHULL::LATEX_DOC", false); // turn off
    }
    if (vpflow.flag("CHULL::LATEX_DOC") || vpflow.flag("CHULL::PNG_IMAGE")) { // latex specific options
      vpflow.flag("CHULL::DOC_ONLY", aurostd::args2flag(argv, cmds, "--document_only|--documentonly|--doc_only|--doconly|--doc")); // no convex hull picture
      vpflow.flag("CHULL::NO_DOC", aurostd::args2flag(argv, cmds, "--no_document|--nodocument|--no_doc|--nodoc|--full_page_image|--fullpageimage")); // no convex hull picture
      if (vpflow.flag("CHULL::NO_DOC") && vpflow.flag("CHULL::DOC_ONLY")) {
        vpflow.flag("CHULL::DOC_ONLY", false);
      }
      vpflow.flag("CHULL::KEEP_TEX", aurostd::args2flag(argv, cmds, "--keep=tex|--keep_tex|--keeptex|--tex")); // keep .tex file to edit
      if (!vpflow.flag("CHULL::KEEP_TEX")) { // only turn on if not already on
        for (size_t i = 0; i < keep_parts.size(); i++) {
          if (!keep_parts[i].empty()) {
            if (keep_parts[i][0] == 'T' || keep_parts[i][0] == 't') {
              vpflow.flag("CHULL::KEEP_TEX", true);
            }
          }
        }
      }
      vpflow.flag("CHULL::IMAGE_ONLY", aurostd::args2flag(argv, cmds, "--image_only|--imageonly|--image|--picture_only|--pictureonly|--picture|--pic")); // image only, no doc or links, special compilation
      if (vpflow.flag("CHULL::IMAGE_ONLY") && vpflow.flag("CHULL::DOC_ONLY")) {
        vpflow.flag("CHULL::IMAGE_ONLY", false);
      } // doc takes precedence
      vpflow.flag("CHULL::LIGHT_CONTRAST", aurostd::args2flag(argv, cmds, "--light_contrast|--lightcontrast|--lc")); // lighter blue
      vpflow.flag("CHULL::LARGE_FONT", aurostd::args2flag(argv, cmds, "--large_font|--largefont|--large|--lf")); // keeps hyperlinks
      if (vpflow.flag("CHULL::PNG_IMAGE") &&
          !(vpflow.flag("CHULL::NO_DOC") || vpflow.flag("CHULL::IMAGE_ONLY"))) { // that means we get full report //!vpflow.flag("CHULL::MULTI_OUTPUT")&& removing this because png means IMAGE_ONLY exclusively, report cannot be created with png
        vpflow.flag("CHULL::IMAGE_ONLY", true);
      }
      if (vpflow.flag("CHULL::PNG_IMAGE")) {
        vpflow.args2addattachedscheme(argv, cmds, "CHULL::PNG_RESOLUTION", "--png_resolution=|--pngresolution=|--pngr=", ""); // calculate distance to hull for point
      }
    }
    vpflow.args2addattachedscheme(argv, cmds, "CHULL::DIST_FOR_EQUIVALENCE_ANALYSIS", "--distance_for_equivalence_analysis=|--dist4equi=|--d4e=", "0.0"); // calculate distance to hull for point
    if (vpflow.flag("CHULL::DIST_FOR_EQUIVALENCE_ANALYSIS")) {
      vpflow.flag("CHULL::DIST_FOR_EQUIVALENCE_ANALYSIS_ND_ONLY", aurostd::args2attachedflag(argv, cmds, "--d4e_ND|--d4e_nd"));
    }
  }

  vpflow.args2addattachedscheme(argv, cmds, "CLAT", "--clat=", "");
  vpflow.args2addattachedscheme(argv, cmds, "COMPARE", "--compare=", "");

  // DX20201220 - put all comparison functions prior to general options - START
  vpflow.flag("COMPARE_DATABASE_ENTRIES", aurostd::args2attachedflag(argv, cmds, "--compare_database_entries"));
  vpflow.flag("COMPARE_MATERIAL", aurostd::args2attachedflag(argv, cmds, "--compare_material|--compare_materials")); // DX20190424 - added plural variant
  vpflow.flag("COMPARE_STRUCTURE", aurostd::args2attachedflag(argv, cmds, "--compare_structure|--compare_structures")); // DX20190424 - added plural variant
  vpflow.flag("COMPARE_PERMUTATION", aurostd::args2attachedflag(argv, cmds,
                                                                "--compare_atom_decoration|--compare_atom_decorations|--unique_atom_decoration|--unique_atom_decorations|--compare_permutation|--compare_"
                                                                "permutations|--unique_permutation|--unique_permutations"));
  vpflow.flag("COMPARE2DATABASE", aurostd::args2attachedflag(argv, cmds, "--compare2database"));
  vpflow.flag("COMPARE2PROTOTYPES", aurostd::args2attachedflag(argv, cmds, "--compare2protos|--compare2prototypes"));

  // COMPARISON: GENERAL OPTIONS
  if (vpflow.flag("COMPARE_DATABASE_ENTRIES") || vpflow.flag("COMPARE_MATERIAL") || vpflow.flag("COMPARE_STRUCTURE") || vpflow.flag("COMPARE_PERMUTATION") || vpflow.flag("COMPARE2DATABASE") ||
      vpflow.flag("COMPARE2PROTOTYPES")) {
    // input values
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE::MISFIT_MATCH", "--misfit_match=", ""); // DX20201118
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE::MISFIT_FAMILY", "--misfit_family=", ""); // DX20201118
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE::NP", "--np=|--num_proc=", "");
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE::PAGE_SIZE", "--page_size=", ""); // ME20220426 - page size for AFLUX
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE::MAGNETIC", "--mag=|--magnetic=|--magmom=|--magmoms=|--magnetic_moment=|--magnetic_moments=", ""); // DX20170803
    // booleans
    vpflow.flag("COMPARE::USAGE", aurostd::args2flag(argv, cmds, "--usage"));
    vpflow.flag("COMPARE::OPTIMIZE_MATCH", aurostd::args2flag(argv, cmds, "--optimize|--optimize_match")); // DX20190201 - changed from fast_match to optimize_match
    vpflow.flag("COMPARE::NO_SCALE_VOLUME", aurostd::args2flag(argv, cmds, "--no_scale_volume"));
    vpflow.flag("COMPARE::KEEP_UNMATCHED", aurostd::args2flag(argv, cmds, "--keep_unmatched")); // DX20190424
    vpflow.flag("COMPARE::IGNORE_SYMMETRY", aurostd::args2flag(argv, cmds, "--ignore_symmetry")); // DX20190424
    vpflow.flag("COMPARE::IGNORE_WYCKOFF", aurostd::args2flag(argv, cmds, "--ignore_Wyckoff")); // DX20190424
    vpflow.flag("COMPARE::IGNORE_ENVIRONMENT_ANALYSIS", aurostd::args2flag(argv, cmds, "--ignore_environment|--ignore_environment_analysis|--ignore_env")); // DX20190807
    vpflow.flag("COMPARE::REMOVE_DUPLICATE_COMPOUNDS", aurostd::args2flag(argv, cmds, "--remove_duplicates|--remove_duplicate_compounds")); // DX20190201
    vpflow.flag("COMPARE::MATCH_TO_AFLOW_PROTOS", aurostd::args2flag(
                                                      argv, cmds, "--match_to_aflow_prototypes|--add_matching_aflow_prototypes|--add_matching_aflow_protos|--add_matching_prototypes|--add_matching_protos")); // DX20190724 //DX20220411 - added --match_to_aflow_prototypes
    vpflow.flag("COMPARE::ADD_AFLOW_PROTOTYPE_DESIGNATION",
                aurostd::args2flag(argv, cmds, "--add_prototype_designation|--add_aflow_prototype_designation|--add_anrl_designation|--add_prototype|--add_prototype_information")); // DX20190724
    vpflow.flag("COMPARE::UNDECORATED_COMPARISON", aurostd::args2flag(argv, cmds, "--undecorated_comparison|--undecorated|--no_atom_decoration")); // DX20191212
    vpflow.flag("COMPARE::PRIMITIVIZE", aurostd::args2flag(argv, cmds, "--primitive|--primitivize|--prim")); // DX20201006
    vpflow.flag("COMPARE::MINKOWSKI", aurostd::args2flag(argv, cmds, "--minkowski|--minkowski_reduction|--minkowski_lattice_reduction")); // DX20201006
    vpflow.flag("COMPARE::NIGGLI", aurostd::args2flag(argv, cmds, "--niggli|--niggli_reduction|--niggli_lattice_reduction")); // DX20201006
    vpflow.flag("COMPARE::SCREEN_ONLY", aurostd::args2flag(argv, cmds, "--screen_only")); // DX20170803
    vpflow.flag("COMPARE::ICSD_COMPARISON", aurostd::args2flag(argv, cmds, "--ICSD")); // DX20190201
    vpflow.flag("COMPARE::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS", aurostd::args2flag(argv, cmds, "--ignore_atom_decoration_comparison|--ignore_decoration_comparison|--ignore_decorations")); // DX20190424
    // ME20210206 - web mode
    if (XHOST.vflag_control.flag("WWW")) {
      vpflow.flag("COMPARE::SCREEN_ONLY", true);
      XHOST.QUIET = true;
    }
  }

  if (vpflow.flag("COMPARE_DATABASE_ENTRIES")) {
    // FUNCTION_SPECIFIC
    // booleans
    vpflow.flag("COMPARE_DATABASE_ENTRIES::STRUCTURE", aurostd::args2flag(argv, cmds, "--structure|--structure_comparison|--ignore_species"));
    // input values
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE_DATABASE_ENTRIES::ALLOY", "--compare_alloy=|--alloy=|--species=", "");
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE_DATABASE_ENTRIES::ARITY", "--arity=|--nspecies=", "");
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE_DATABASE_ENTRIES::CATALOG", "--catalog=|--library=", "");
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE_DATABASE_ENTRIES::PROPERTY_LIST", "--properties=|--property_list=|--property=", "");
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE_DATABASE_ENTRIES::RELAXATION_STEP", "--relaxation_step=|--relaxation=", "");
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE_DATABASE_ENTRIES::SPACE_GROUP", "--space_group=|--spacegroup=|--sg=", "");
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE_DATABASE_ENTRIES::STOICHIOMETRY", "--stoichiometry=|--stoich=", "");
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE_DATABASE_ENTRIES::WYCKOFF", "--Wyckoff=|--wyckoff=", "");
  }
  if (vpflow.flag("COMPARE_MATERIAL") || vpflow.flag("COMPARE_STRUCTURE")) {
    vpflow.flag("COMPARE_STRUCTURE::PRINT", aurostd::args2flag(argv, cmds, "--print|--print_mapping|--print_mappings"));
    // structures appended to command
    if (aurostd::args2attachedflag(argv, "--compare_material=|--compare_materials=|--compare_structure=|--compare_structures=")) { // DX20170803
      vpflow.args2addattachedscheme(argv, cmds, "COMPARE_STRUCTURE::STRUCTURE_LIST", "--compare_material=|--compare_materials=|--compare_structure=|--compare_structures=", "");
    }
    // structures from directory
    if (XHOST.vflag_control.flag("DIRECTORY")) {
      vpflow.push_attached("COMPARE_STRUCTURE::DIRECTORY", XHOST.vflag_control.getattachedscheme("DIRECTORY"));
    } else if (aurostd::args2attachedflag(argv, "--compare_material_directory|--compare_material_dir|--compare_structure_directory|--compare_structure_dir")) { // legacy
      vpflow.push_attached("COMPARE_STRUCTURE::DIRECTORY", ".");
    }
    // structures from file
    if (XHOST.vflag_control.flag("FILE")) {
      vpflow.push_attached("COMPARE_STRUCTURE::FILE", XHOST.vflag_control.getattachedscheme("FILE"));
    }
  }
  if (vpflow.flag("COMPARE_PERMUTATION")) {
    vpflow.flag("COMPARE_PERMUTATION::PRINT", aurostd::args2flag(argv, cmds, "--print_misfit|--print|--misfit"));
  }
  if (vpflow.flag("COMPARE2DATABASE")) {
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE2DATABASE::PROPERTY_LIST", "--properties=|--property_list=|--property=", "");
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE2DATABASE::RELAXATION_STEP", "--relaxation_step=|--relaxation=", "");
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE2DATABASE::CATALOG", "--catalog=|--library=", "");
  }
  if (vpflow.flag("COMPARE2PROTOTYPES")) {
    vpflow.args2addattachedscheme(argv, cmds, "COMPARE2PROTOTYPES::CATALOG", "--catalog=|--library=", "");
  }
  // DX20201220 - put all comparison functions prior to general options - END

  vpflow.flag("AFLOWMACHL::CoordCE_CSV", aurostd::args2flag(argv, cmds, "--coordce_csv")); // CO20200930
  vpflow.flag("CORNERS", aurostd::args2flag(argv, cmds, "--corner|--corners"));

  // DX20210301 - consolidated all DATA functions
  vpflow.flag("DATA", aurostd::args2attachedflag(argv, cmds, "--data"));
  vpflow.flag("DATA_CRYSTAL_POINT_GROUP", aurostd::args2attachedflag(argv, cmds, "--point_group_crystal_data|--pointgroup_crystal_data|--pgroupxtal_data|--pgroup_xtal_data")); // DX20210209
  vpflow.flag("DATA_EDATA", aurostd::args2attachedflag(argv, cmds, "--edata"));
  vpflow.flag("DATA_REAL_LATTICE", aurostd::args2attachedflag(argv, cmds, "--lattice_data|--data_lattice|--real_lattice_data|--data_real_lattice")); // DX20210209
  vpflow.flag("DATA_RECIPROCAL_LATTICE", aurostd::args2attachedflag(argv, cmds, "--reciprocal_lattice_data|--reciprocallattice_data|--klattice_data|--data_reciprocal_lattice")); // DX20210209
  vpflow.flag("DATA_SGDATA", aurostd::args2attachedflag(argv, cmds, "--sgdata|--space_group_data"));
  vpflow.flag("DATA_SUPERLATTICE", aurostd::args2attachedflag(argv, cmds, "--superlattice_data|--data_superlattice")); // DX20210209

  // EDATA/SGDATA/LATTICE DATA: GENERAL OPTIONS
  if (vpflow.flag("DATA") || vpflow.flag("DATA_CRYSTAL_POINT_GROUP") || vpflow.flag("DATA_EDATA") || vpflow.flag("DATA_REAL_LATTICE") || vpflow.flag("DATA_RECIPROCAL_LATTICE") || vpflow.flag("DATA_SGDATA") ||
      vpflow.flag("DATA_SUPERLATTICE")) {
    vpflow.flag("DATA::NO_SCAN", aurostd::args2flag(argv, cmds, "--no_scan"));
    vpflow.flag("DATA::SUPPRESS_WYCKOFF_PRINTING", aurostd::args2flag(argv, cmds, "--suppress_Wyckoff|--suppress_Wyckoff_printing|--suppress_wyckoff|--suppress_wyckoff_printing")); // DX20210211
    vpflow.args2addattachedscheme(argv, cmds, "DATA::SETTING", "--setting=", "1");
    vpflow.args2addattachedscheme(argv, cmds, "DATA::MAGNETIC", "--mag=|--magnetic=|--magmom=", "");
    if (aurostd::args2attachedflag(argv,
                                   "--data=|--lattice_data=|--data_lattice=|--real_lattice_data=|--data_real_lattice=|--point_group_crystal_data=|--pointgroup_crystal_data=|--pgroupxtal_data=|--pgroup_xtal_"
                                   "data=|--reciprocal_lattice_data=|--reciprocallattice_data=|--klattice_data=|--data_reciprocal_lattice=|--superlattice_data=|--data_superlattice=|--edata=|--sgdata=|--space_"
                                   "group_data=")) {
      vpflow.args2addattachedscheme(argv, cmds, "DATA::TOLERANCE",
                                    "--data=|--lattice_data=|--data_lattice=|--real_lattice_data=|--data_real_lattice=|--point_group_crystal_data=|--pointgroup_crystal_data=|--pgroupxtal_data=|--pgroup_xtal_"
                                    "data=|--reciprocal_lattice_data=|--reciprocallattice_data=|--klattice_data=|--data_reciprocal_lattice=|--superlattice_data=|--data_superlattice=|--edata=|--sgdata=|--space_"
                                    "group_data=",
                                    "");
    }
    vpflow.flag("DATA::USAGE", aurostd::args2flag(argv, cmds, "--usage"));
  }

  vpflow.args2addattachedscheme(argv, cmds, "DATA1", "--data1=", "");
  vpflow.flag("DATA2", aurostd::args2flag(argv, cmds, "--data2"));
  vpflow.args2addattachedscheme(argv, cmds, "DEBYE", "--debye=", "");
  vpflow.args2addattachedscheme(argv, cmds, "DIFF", "--diff=", "./");

  vpflow.args2addattachedscheme(argv, cmds, "DISP", "--disp=", "");
  vpflow.args2addattachedscheme(argv, cmds, "DIST", "--dist=", "");

  vpflow.flag("EDOS", aurostd::args2flag(argv, cmds, "--edos"));
  vpflow.flag("EFFMASS", aurostd::args2flag(argv, cmds, "--em")); // CAMILO
  vpflow.args2addattachedscheme(argv, cmds, "EWALD", "--ewald=", "");

  // DX20170818 - Added tolerance and no_scan options to Xgroups - START
  vpflow.args2addattachedscheme(argv, cmds, "EQUIVALENT", "--equivalent=|--equiv=|--inequivalent=|--inequiv=|--iatoms=|--eatoms=", "");
  if (vpflow.flag("EQUIVALENT")) {
    vpflow.flag("SYMMETRY::NO_SCAN", aurostd::args2flag(argv, cmds, "--no_scan"));
    if (aurostd::args2attachedflag(argv, "--equivalent=|--equiv=|--inequivalent=|--inequiv=|--iatoms=|--eatoms=")) { // DX20170803
      vpflow.args2addattachedscheme(argv, cmds, "SYMMETRY::TOLERANCE", "--equivalent=|--equiv=|--inequivalent=|--inequiv=|--iatoms=|--eatoms=", ""); // DX20200907 - default is system specific, leaving empty
    }
    // DX20170921 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv, cmds, "SYMMETRY::MAGNETIC", "--mag=|--magnetic=|--magmom=", ""); // DX20170803
    // DX20170921 - MAGNETIC SYMMETRY - END
  }
  // DX20170818 - Added tolerance and no_scan options to Xgroups - END
  vpflow.flag("EXTRACT_SYMMETRY", aurostd::args2flag(argv, cmds, "--extract_symmetry|--xsymmetry"));

  vpflow.args2addattachedscheme(argv, cmds, "EIGCURV", "--eigcurv=", "./"); // CAMILO

  // DX20170818 - Added tolerance and no_scan options to Xgroups - START
  vpflow.args2addattachedscheme(argv, cmds, "FGROUP", "--factorgroup=|--fgroup=", "");
  if (vpflow.flag("FGROUP")) {
    vpflow.flag("SYMMETRY::NO_SCAN", aurostd::args2flag(argv, cmds, "--no_scan"));
    if (aurostd::args2attachedflag(argv, "--factorgroup=|--fgroup=")) { // DX20170803
      vpflow.args2addattachedscheme(argv, cmds, "SYMMETRY::TOLERANCE", "--factorgroup=|--fgroup=", ""); // DX20200907 - default is system specific, leaving empty
    }
    vpflow.flag("SYMMETRY::SCREEN_ONLY", aurostd::args2flag(argv, cmds, "--screen_only")); // DX20170803
    // DX20170921 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv, cmds, "SYMMETRY::MAGNETIC", "--mag=|--magnetic=|--magmom=", ""); // DX20170803
    // DX20170921 - MAGNETIC SYMMETRY - END
    //  ME20210206 - web mode
    if (XHOST.vflag_control.flag("WWW")) {
      vpflow.flag("SYMMETRY::SCREEN_ONLY", true);
      XHOST.QUIET = true;
    }
  }
  // DX20170818 - Added tolerance and no_scan options to Xgroups - END

  vpflow.flag("FIND_CLOSED_PACKING_PLANE", aurostd::args2flag(argv, cmds, "--find_closed_packing_plane")); // CO20191110

  vpflow.flag("FRAC", aurostd::args2flag(argv, cmds, "--frac|-frac|--fractional|-fract|--fract|--direct|-direct|-f|-d"));
  vpflow.flag("FROZSL_VASPSETUP_AFLOW", aurostd::args2flag(argv, cmds, "--frozsl_vaspsetup_aflow|--frozsl_vaspsetup|--frozsl_vasp|--frozsl_setup|--phvaspsetup"));
  vpflow.flag("FROZSL_VASPSETUP_POSCAR", aurostd::args2flag(argv, cmds, "--frozsl_vaspsetup_poscar"));
  vpflow.flag("FROZSL_INPUT", aurostd::args2flag(argv, cmds, "--frozsl_input"));
  vpflow.flag("FROZSL_OUTPUT", aurostd::args2flag(argv, cmds, "--frozsl_output"));
  vpflow.flag("FROZSL_ANALYZE", aurostd::args2flag(argv, cmds, "--frozsl_analyze"));
  vpflow.flag("FROZSL_README", aurostd::args2flag(argv, cmds, "--readme=frozsl"));

  vpflow.flag("GFA::INIT", aurostd::args2flag(argv, cmds, "--gfa|--glass_forming_ability")); // DF20190329 - GFA
  if (vpflow.flag("GFA::INIT")) { // DF20190329
    vpflow.args2addattachedscheme(argv, cmds, "GFA::AE_FILE", "--atomic_environments_file=|--ae_file=|aef=", "none"); // DF20190329
    vpflow.args2addattachedscheme(argv, cmds, "GFA::FORMATION_ENTHALPY_CUTOFF", "--cutoff_formation_enthalpy=|--cutoff_enthalpy=|--cutoff_energy=|--cut=", "0.05"); // DF20190619
  } // DF20190329

  vpflow.flag("ATOMIC_ENVIRONMENT::INIT", aurostd::args2flag(argv, cmds, "--ae|--atomic_environment")); // HE20210331 - Testing
  if (vpflow.flag("ATOMIC_ENVIRONMENT::INIT")) { // HE20210331
    vpflow.args2addattachedscheme(argv, cmds, "ATOMIC_ENVIRONMENT::AUID", "--auid=", "none");
    vpflow.args2addattachedscheme(argv, cmds, "ATOMIC_ENVIRONMENT::MODE", "--mode=", "1");
    vpflow.args2addattachedscheme(argv, cmds, "ATOMIC_ENVIRONMENT::RADIUS", "--radius=", "4");

  } // HE20210331

  vpflow.flag("GENERATE_CERAMICS", aurostd::args2flag(argv, cmds, "--generate_ceramics|--gen_ceram")); // CO20200731
  if (vpflow.flag("GENERATE_CERAMICS")) { // CO20200731
    vpflow.args2addattachedscheme(argv, cmds, "GENERATE_CERAMICS::NON_METALS", "--non_metals=|--nm=", "C");
    vpflow.args2addattachedscheme(argv, cmds, "GENERATE_CERAMICS::METALS", "--metals=|--m=", "Mo,Nb,Ta,V,W");
    vpflow.args2addattachedscheme(argv, cmds, "GENERATE_CERAMICS::METAL_ARITY", "--N=", "5");
  }

  vpflow.flag("GEOMETRY", aurostd::args2flag(argv, cmds, "--geometry|--abc_angles")); // CO20190808

  vpflow.flag("GULP", aurostd::args2flag(argv, cmds, "--gulp"));
  vpflow.flag("GETTEMP", aurostd::args2flag(argv, cmds, "--getTEMP|--getTEMPS|--Init::GetTEMPs|--gettemp|--gettemps|--getemp|--getemps"));

  vpflow.args2addattachedscheme(argv, cmds, "KPOINTS", "--kpoints=|--kppra=|-kpoints=|-kppra=|-k=", "1");
  vpflow.args2addattachedscheme(argv, cmds, "FLAG::XVASP_KPOINTS_DELTA", "--delta_kpoints=|--dkpoints=|-dkpoints=|-dk=", "0.01"); // CO20171025

  vpflow.flag("KPATH", aurostd::args2flag(argv, cmds, "--kpath"));

  vpflow.flag("JOINSTRLIST", aurostd::args2flag(argv, cmds, "--join_strlist") && argv.at(1) == "--join_strlist");

  vpflow.args2addattachedscheme(argv, cmds, "HKL", "--hkl=", "");
  vpflow.args2addattachedscheme(argv, cmds, "HKL_SEARCH_TRIVIAL", "--hkl_search=", "");
  vpflow.args2addattachedscheme(argv, cmds, "HKL_SEARCH_SIMPLE", "--hkl_search_simple=", "");
  vpflow.args2addattachedscheme(argv, cmds, "HKL_SEARCH_COMPLETE", "--hkl_search_complete=", "");

  vpflow.flag("IAP::INIT", aurostd::args2flag(argv, cmds, "--iap|--IAP")); // initiate IAP writer
  if (vpflow.flag("IAP::INIT")) {
    vpflow.args2addattachedscheme(argv, cmds, "IAP::FORMAT", "--format=|--f=", "");
    vpflow.args2addattachedscheme(argv, cmds, "IAP::AFLOWLIBPATH", "--aflowlibpath=|--path=", "");
  }

  vpflow.flag("ICSD", aurostd::args2flag(argv, cmds, "--icsd"));
  vpflow.flag("ICSD_ALLLESSTHAN", aurostd::args2flag(argv, cmds, "--icsd_alllessthan"));
  vpflow.flag("ICSD_ALLMORETHAN", aurostd::args2flag(argv, cmds, "--icsd_allmorethan"));
  vpflow.flag("ICSD_BASISLT", aurostd::args2flag(argv, cmds, "--icsd_basislessthan"));
  vpflow.flag("ICSD_BASISGT", aurostd::args2flag(argv, cmds, "--icsd_basismorethan"));
  vpflow.flag("ICSD_CHEM", aurostd::args2flag(argv, cmds, "--icsd_chem"));
  vpflow.flag("ICSD_CUBIC", aurostd::args2flag(argv, cmds, "--icsd_cubic"));
  vpflow.flag("ICSD_DENSLESSTHAN", aurostd::args2flag(argv, cmds, "--icsd_denslessthan"));
  vpflow.flag("ICSD_DENSMORETHAN", aurostd::args2flag(argv, cmds, "--icsd_densmorethan"));
  vpflow.flag("ICSD_HEXAGONAL", aurostd::args2flag(argv, cmds, "--icsd_hexagonal"));
  vpflow.flag("ICSD_ID", aurostd::args2flag(argv, cmds, "--icsd_id|--icsd_ID"));
  vpflow.flag("ICSD_MAKELABEL", aurostd::args2flag(argv, cmds, "--icsd_makelabel"));
  vpflow.flag("ICSD_LESSTHAN", aurostd::args2flag(argv, cmds, "--icsd_lessthan"));
  vpflow.flag("ICSD_LISTMETALS", aurostd::args2flag(argv, cmds, "--icsd_listmetals"));
  vpflow.flag("ICSD_MONOCLINIC", aurostd::args2flag(argv, cmds, "--icsd_monoclinic"));
  vpflow.flag("ICSD_MORETHAN", aurostd::args2flag(argv, cmds, "--icsd_morethan"));
  vpflow.flag("ICSD_N_ARY", aurostd::args2flag(argv, cmds, "--icsd_n_ary"));
  vpflow.flag("ICSD_NOBROKENBASIS", aurostd::args2flag(argv, cmds, "--icsd_nobrokenbasis|--icsd_completebasis"));
  vpflow.flag("ICSD_NOPARTIALOCC", aurostd::args2flag(argv, cmds, "--icsd_nopartialocc"));
  vpflow.flag("ICSD_ORTHORHOMBIC", aurostd::args2flag(argv, cmds, "--icsd_orthorhombic|--icsd_orthorombic|--icsd_ortorhombic|--icsd_ortorombic"));
  vpflow.flag("ICSD_PROTO", aurostd::args2flag(argv, cmds, "--icsd_proto"));
  vpflow.flag("ICSD_RHOMBOHEDRAL", aurostd::args2flag(argv, cmds, "--icsd_rhombohedral|--icsd_rombohedral|--icsd_trigonal"));
  vpflow.flag("ICSD_REMOVE_AND", aurostd::args2flag(argv, cmds, "--icsd_remove_and"));
  vpflow.flag("ICSD_REMOVE_OR", aurostd::args2flag(argv, cmds, "--icsd_remove_or"));
  vpflow.flag("ICSD_REMOVEMETALS", aurostd::args2flag(argv, cmds, "--icsd_removemetals"));
  vpflow.flag("ICSD_SG", aurostd::args2flag(argv, cmds, "--icsd_sg"));
  vpflow.flag("ICSD_SGLESSTHAN", aurostd::args2flag(argv, cmds, "--icsd_sglessthan"));
  vpflow.flag("ICSD_SGMORETHAN", aurostd::args2flag(argv, cmds, "--icsd_sgmorethan"));
  vpflow.flag("ICSD_TETRAGONAL", aurostd::args2flag(argv, cmds, "--icsd_tetragonal"));
  vpflow.flag("ICSD_TRICLINIC", aurostd::args2flag(argv, cmds, "--icsd_triclinic"));
  vpflow.flag("ICSD_TRIGONAL", aurostd::args2flag(argv, cmds, "--icsd_trigonal"));
  vpflow.flag("ICSD_TRI", aurostd::args2flag(argv, cmds, "--icsd_tri|--icsd_TRI"));
  vpflow.flag("ICSD_MCL", aurostd::args2flag(argv, cmds, "--icsd_mcl|--icsd_MCL"));
  vpflow.flag("ICSD_MCLC", aurostd::args2flag(argv, cmds, "--icsd_mclc|--icsd_MCLC"));
  vpflow.flag("ICSD_ORC", aurostd::args2flag(argv, cmds, "--icsd_orc|--icsd_ORC"));
  vpflow.flag("ICSD_ORCC", aurostd::args2flag(argv, cmds, "--icsd_orcc|--icsd_ORCC"));
  vpflow.flag("ICSD_ORCF", aurostd::args2flag(argv, cmds, "--icsd_orcf|--icsd_ORCF"));
  vpflow.flag("ICSD_ORCI", aurostd::args2flag(argv, cmds, "--icsd_orci|--icsd_ORCI"));
  vpflow.flag("ICSD_TET", aurostd::args2flag(argv, cmds, "--icsd_tet|--icsd_TET"));
  vpflow.flag("ICSD_BCT", aurostd::args2flag(argv, cmds, "--icsd_bct|--icsd_BCT"));
  vpflow.flag("ICSD_RHL", aurostd::args2flag(argv, cmds, "--icsd_rhl|--icsd_RHL"));
  vpflow.flag("ICSD_HEX", aurostd::args2flag(argv, cmds, "--icsd_hex|--icsd_HEX"));
  vpflow.flag("ICSD_CUB", aurostd::args2flag(argv, cmds, "--icsd_cub|--icsd_CUB"));
  vpflow.flag("ICSD_FCC", aurostd::args2flag(argv, cmds, "--icsd_fcc|--icsd_FCC"));
  vpflow.flag("ICSD_BCC", aurostd::args2flag(argv, cmds, "--icsd_bcc|--icsd_BCC"));
  vpflow.flag("ICSD_UNIQUE", aurostd::args2flag(argv, cmds, "--icsd_unique"));
  vpflow.flag("ICSD2POSCAR", aurostd::args2flag(argv, cmds, "--icsd2poscar|--icsd2POSCAR"));
  vpflow.flag("ICSD2PROTO", aurostd::args2flag(argv, cmds, "--icsd2proto"));
  vpflow.flag("ICSD2WYCK", aurostd::args2flag(argv, cmds, "--icsd2wyck|--icsdwyck"));

  vpflow.flag("IDENTICAL", aurostd::args2flag(argv, cmds, "--identical"));
  vpflow.flag("INCELL", aurostd::args2flag(argv, cmds, "--incell"));
  vpflow.flag("INCOMPACT", aurostd::args2flag(argv, cmds, "--incompact"));
  vpflow.flag("INSPHERE", aurostd::args2flag(argv, cmds, "--insphere"));
  vpflow.args2addattachedscheme(argv, cmds, "INTPOL", "--intpol=", "");

  vpflow.flag("INWS", aurostd::args2flag(argv, cmds, "--inwignerseitz|--inws"));

  // DX20200131 - add isopointal prototype function - START
  vpflow.flag("ISOPOINTAL_PROTOTYPES", aurostd::args2attachedflag(argv, cmds, "--get_isopointal_prototypes|--isopointal_prototypes|--get_same_symmetry_prototypes"));
  if (vpflow.flag("ISOPOINTAL_PROTOTYPES")) {
    vpflow.flag("ISOPOINTAL_PROTOTYPES::USAGE", aurostd::args2flag(argv, cmds, "--usage"));
    vpflow.args2addattachedscheme(argv, cmds, "ISOPOINTAL_PROTOTYPES::CATALOG", "--catalog=|--library=", "");
  }
  // DX20200131 - add isopointal prototype function - END

  vpflow.args2addattachedscheme(argv, cmds, "INFLATE_LATTICE", "--inflate_lattice=|--ilattice=", "");
  vpflow.args2addattachedscheme(argv, cmds, "INFLATE_VOLUME", "--inflate_volume=|--ivolume=", "");

  vpflow.flag("KBAND", aurostd::args2flag(argv, cmds, "--kband"));
  vpflow.args2addattachedscheme(argv, cmds, "KILL", "--kill=", "");

  vpflow.flag("HNF", (aurostd::args2flag(argv, cmds, "--hnf|--HNF|--hfn|--HFN|--hnftol|--HNFTOL|--hfntol|--HFNTOL|--pocc_hnf"))); // CO20181226
  vpflow.flag("HNFCELL", aurostd::args2flag(argv, cmds, "--hnfcell|--hfncell|--pocc_show_unique_structures|--pocc_show_unique|--pocc_print_unique")); // CO20181226
  vpflow.flag("POCC_COUNT_TOTAL", aurostd::args2flag(argv, cmds, "--pocc_count_total_structures|--pocc_count_total")); // CO20181226
  vpflow.flag("POCC_COUNT_UNIQUE", aurostd::args2flag(argv, cmds, "--pocc_count_unique_structures|--pocc_count_unique")); // CO20181226
  vpflow.flag("POCC_COUNT_UNIQUE_FAST", aurostd::args2flag(argv, cmds, "--pocc_count_unique_fast")); // SD20230609
  vpflow.flag("POCC_SKIP_WRITING_FILES", aurostd::args2flag(argv, cmds, "--pocc_skip_writing_files")); // CO20181226
  if (vpflow.flag("HNF") || vpflow.flag("POCC_COUNT_TOTAL") || vpflow.flag("POCC_COUNT_UNIQUE") || vpflow.flag("POCC_COUNT_UNIQUE_FAST")) {
    vpflow.flag("HNFCELL", true);
  } // funnel all of these through command_line function
  if (vpflow.flag("HNFCELL")) {
    vpflow.flag("POCC_SKIP_WRITING_FILES", true);
  } // CO20190401  //no point writing files if through command_line
  if (XHOST.vflag_control.flag("WWW")) {
    vpflow.flag("POCC_SKIP_WRITING_FILES", true);
  } // CO20190401 //CO20200404 - new web flag
  vpflow.flag("MULTIENUMALL", aurostd::args2flag(argv, cmds, "--multienum|--enum"));
  vpflow.flag("MULTIENUMSORT", aurostd::args2flag(argv, cmds, "--multienumsort|--enumsort"));
  vpflow.flag("POSCAR2ENUM", (aurostd::args2flag(argv, cmds, "--poscar2multienum|--poscar2enum")));
  vpflow.flag("POSCAR2GULP", (aurostd::args2flag(argv, cmds, "--poscar2gulp")));
  vpflow.flag("POCC_INPUT", aurostd::args2flag(argv, cmds, "--pocc_input|--enum_input"));
  vpflow.args2addattachedscheme(argv, cmds, "POCC::CONVOLUTION", "--pocc_convolution=|--pocc_conv=|--poccconv=", "");

  vpflow.args2addattachedscheme(argv, cmds, "JMOL", "--jmol=", "");

  vpflow.flag("JMOLGIF", aurostd::args2flag(argv, cmds, "--jgif"));
  vpflow.args2addattachedscheme(argv, cmds, "JUSTAFTER", "--justafter=", "");
  vpflow.args2addattachedscheme(argv, cmds, "JUSTBEFORE", "--justbefore=", "");
  vpflow.args2addattachedscheme(argv, cmds, "JUSTBETWEEN", "--justbetween=", "");

  vpflow.flag("LATTICEREDUCTION", aurostd::args2flag(argv, cmds, "--latticereduction|--latreduction"));

  vpflow.args2addattachedscheme(argv, cmds, "LTCELL", "--ltcell=", "");
  vpflow.args2addattachedscheme(argv, cmds, "LTCELLFV", "--ltcellfv=", "");

  vpflow.args2addattachedscheme(argv, cmds, "LATTICE_TYPE", "--lattice_type=|--lattice=|--lattice_crystal=", ""); // DX20200820 - allow tolerance to be added
  if (aurostd::args2attachedflag(argv, "--lattice_type=|--lattice=|--lattice_crystal=")) { // DX20200820 - add tolerance
    vpflow.args2addattachedscheme(argv, cmds, "LATTICE::TOLERANCE", "--lattice_type=|--lattice=|--lattice_crystal=", ""); // DX20200907 - default is system specific, leaving empty
  }
  vpflow.args2addattachedscheme(argv, cmds, "LATTICE_LATTICE_TYPE", "--lattice_lattice_type=|--lattice_lattice=", ""); // DX20200820 - allow tolerance to be added
  if (aurostd::args2attachedflag(argv, "--lattice_lattice_type=|--lattice_lattice=")) { // DX20200820 - add tolerance
    vpflow.args2addattachedscheme(argv, cmds, "LATTICE_LATTICE::TOLERANCE", "--lattice_lattice_type=|--lattice_lattice=", ""); // DX20200907 - default is system specific, leaving empty
  }
  vpflow.flag("LATTICE_HISTOGRAM", aurostd::args2flag(argv, cmds, "--latticehistogram"));

  // DX20181023 - list prototype labels - START
  vpflow.flag("LIST_PROTOTYPE_LABELS", (aurostd::args2flag(argv, cmds, "--prototype_labels|--proto_labels"))); // && (argv.size()==2));
  if (vpflow.flag("LIST_PROTOTYPE_LABELS")) {
    vpflow.args2addattachedscheme(argv, cmds, "LIST_PROTOTYPE_LABELS::LIBRARY", "--library=|--catalog=", ""); // DX20190509 - added catalog variant
    vpflow.args2addattachedscheme(argv, cmds, "LIST_PROTOTYPE_LABELS::ARITY", "--arity=|--number_of_species=|--species_count=", "");
    vpflow.args2addattachedscheme(argv, cmds, "LIST_PROTOTYPE_LABELS::STOICHIOMETRY", "--stoichiometry=|--stoich=", "");
    vpflow.args2addattachedscheme(argv, cmds, "LIST_PROTOTYPE_LABELS::SPACE_GROUP", "--space_group_number=|--space_group=|--sg=", "");
    vpflow.args2addattachedscheme(argv, cmds, "LIST_PROTOTYPE_LABELS::WYCKOFF_STRING", "--Wyckoff_letters=|--wyckoff_letters=|--Wyckoff=|--wyckoff=", "");
    vpflow.flag("LIST_PROTOTYPE_LABELS::JSON", aurostd::args2flag(argv, cmds, "--json"));
  }
  // DX20181023 - list prototype labels - END

  vpflow.flag("MAGNETICPARAMETERS", aurostd::args2flag(argv, cmds, "--magpara") || aurostd::args2attachedflag(argv, cmds, "--magpara="));

  vpflow.args2addattachedscheme(argv, cmds, "MAXATOMS", "--maxatoms=|--max_atoms=|--atomsmax=|--atoms_max=", "");

  vpflow.flag("MAKESTRLIST", aurostd::args2flag(argv, cmds, "--make_strlist") && argv.at(1) == "--make_strlist");

  vpflow.flag("MINKOWSKI_BASIS_REDUCTION", aurostd::args2flag(argv, cmds, "--minkowski_basis_reduction|--minkowski|--mink"));
  vpflow.flag("MISCIBILITY", aurostd::args2flag(argv, cmds, "--MIX|--mix|--MISCIBILITY|--miscibility|--MISCIBILE|--miscibile"));
  vpflow.flag("MOM", aurostd::args2flag(argv, cmds, "--mom"));
  vpflow.flag("MSI", aurostd::args2flag(argv, cmds, "--msi"));

  vpflow.flag("MULTI=SH", aurostd::args2flag(argv, cmds, "--multi=sh|--multish"));

  vpflow.flag("NATOMS", aurostd::args2flag(argv, cmds, "--natoms|--numatoms"));
  vpflow.flag("NBONDXX", aurostd::args2flag(argv, cmds, "--nbondxx")); // CO20171025
  vpflow.flag("NAMES", (aurostd::args2flag(argv, cmds, "--names|--species") && (argv.at(1) == "--names" || argv.at(1) == "--species")));
  vpflow.flag("NANOPARTICLE", (aurostd::args2flag(argv, cmds, "--nanoparticle") && argv.at(1) == "--nanoparticle"));
  vpflow.flag("NDATA", aurostd::args2flag(argv, cmds, "--ndata"));
  vpflow.flag("NIGGLI", aurostd::args2flag(argv, cmds, "--niggli"));
  // vpflow.flag("NOSD",aurostd::args2flag(argv,cmds,"--nosd"));
  vpflow.flag("NNDIST", aurostd::args2flag(argv, cmds, "--nn|--nearestneighbour|--nearestneighbor"));
  vpflow.flag("NOORDERPARAMETER", aurostd::args2flag(argv, cmds, "--noorderparameter|--norderparameter|--noorder|--norder"));

  vpflow.flag("NUMNAMES", aurostd::args2flag(argv, cmds, "--numnames") && argv.at(1) == "--numnames");
  vpflow.flag("NSPECIES", aurostd::args2flag(argv, cmds, "--nspecies|--numspecies"));

  vpflow.flag("PEARSON_SYMBOL", aurostd::args2attachedflag(argv, cmds, "--pearson_symbol|--pearson|--Pearson_symbol|--Pearson")); // DX20210611 - added capitalized variants and changed to args2attachedflag
  // DX20210611 - START
  if (vpflow.flag("PEARSON_SYMBOL")) {
    vpflow.flag("PEARSON::NO_SCAN", aurostd::args2flag(argv, cmds, "--no_scan"));
    if (aurostd::args2attachedflag(argv, "--pearson_symbol=|--pearson=|--Pearson_symbol=|--Pearson=")) {
      vpflow.args2addattachedscheme(argv, cmds, "PEARSON_SYMBOL::TOLERANCE", "--pearson_symbol=|--pearson=|--Pearson_symbol=|--Pearson=", "");
    }
    if (aurostd::args2attachedflag(argv, "--tolerance_spectrum=|--tol_spectrum=")) {
      vpflow.args2addattachedscheme(argv, cmds, "PEARSON_SYMBOL::TOLERANCE_SPECTRUM", "--tolerance_spectrum=|--tol_spectrum=", "0.01,1.0,100");
    }
    vpflow.args2addattachedscheme(argv, cmds, "PEARSON_SYMBOL::MAGNETIC", "--mag=|--magnetic=|--magmom=", "");
    vpflow.flag("PEARSON_SYMBOL::USAGE", aurostd::args2flag(argv, cmds, "--usage"));
  }
  // DX20210611 - STOP
  //  vpflow.flag("POCCUPATION",aurostd::args2flag(argv,cmds,"--poccupation|--partial_occupation|--partialoccupation|--pocc"));
  //  vpflow.flag("OPARAMETER",aurostd::args2flag(argv,cmds,"--oparameter|--order_parameter|--orderparameter|--opar"));

  vpflow.flag("PDB", aurostd::args2flag(argv, cmds, "--pdb"));
  vpflow.flag("PDOS", aurostd::args2flag(argv, cmds, "--pdos") && argv.at(1) == "--pdos");
  // DX20170818 - Added tolerance and no_scan options to Xgroups - START
  vpflow.args2addattachedscheme(argv, cmds, "PGROUP", "--pointgroup=|--pgroup=", "");
  if (vpflow.flag("PGROUP")) {
    vpflow.flag("SYMMETRY::NO_SCAN", aurostd::args2flag(argv, cmds, "--no_scan"));
    if (aurostd::args2attachedflag(argv, "--pointgroup=|--pgroup=")) { // DX20170803
      vpflow.args2addattachedscheme(argv, cmds, "SYMMETRY::TOLERANCE", "--pointgroup=|--pgroup=", ""); // DX20200907 - default is system specific, leaving empty
    }
    vpflow.flag("SYMMETRY::SCREEN_ONLY", aurostd::args2flag(argv, cmds, "--screen_only")); // DX20170803
    // ME20210206 - web mode
    if (XHOST.vflag_control.flag("WWW")) {
      vpflow.flag("SYMMETRY::SCREEN_ONLY", true);
      XHOST.QUIET = true;
    }
  }
  // DX20170818 - Added tolerance and no_scan options to Xgroups - END
  // DX20170818 - Added tolerance and no_scan options to Xgroups - START
  vpflow.args2addattachedscheme(argv, cmds, "PGROUPX", "--pointgroup_crystal=|--pgroup_crystal=|--pgroup_xtal=|--pgroupx=|--pgroupX=", "");
  if (vpflow.flag("PGROUPX")) {
    vpflow.flag("SYMMETRY::NO_SCAN", aurostd::args2flag(argv, cmds, "--no_scan"));
    if (aurostd::args2attachedflag(argv, "--pointgroup_crystal=|--pgroup_crystal=|--pgroup_xtal=|--pgroupx=|--pgroupX=")) { // DX20170803
      vpflow.args2addattachedscheme(argv, cmds, "SYMMETRY::TOLERANCE", "--pointgroup_crystal=|--pgroup_crystal=|--pgroup_xtal=|--pgroupx=|--pgroupX=", ""); // DX20200907 - default is system specific, leaving empty
    }
    vpflow.flag("SYMMETRY::SCREEN_ONLY", aurostd::args2flag(argv, cmds, "--screen_only")); // DX20170803
    // DX20170921 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv, cmds, "SYMMETRY::MAGNETIC", "--mag=|--magnetic=|--magmom=", ""); // DX20170803
    // DX20170921 - MAGNETIC SYMMETRY - END
    //  ME20210206 - web mode
    if (XHOST.vflag_control.flag("WWW")) {
      vpflow.flag("SYMMETRY::SCREEN_ONLY", true);
      XHOST.QUIET = true;
    }
  }
  // DX20170818 - Added tolerance and no_scan options to Xgroups - END
  // DX20170818 - Added tolerance and no_scan options to Xgroups - START
  vpflow.args2addattachedscheme(argv, cmds, "PGROUPK", "--pointgroupklattice=|--pgroupk=", "");
  if (vpflow.flag("PGROUPK")) {
    vpflow.flag("SYMMETRY::NO_SCAN", aurostd::args2flag(argv, cmds, "--no_scan"));
    if (aurostd::args2attachedflag(argv, "--pointgroupklattice=|--pgroupk=")) { // DX20170803
      vpflow.args2addattachedscheme(argv, cmds, "SYMMETRY::TOLERANCE", "--pointgroupklattice=|--pgroupk=", ""); // DX20200907 - default is system specific, leaving empty
    }
    vpflow.flag("SYMMETRY::SCREEN_ONLY", aurostd::args2flag(argv, cmds, "--screen_only")); // DX20170803
    // ME20210206 - web mode
    if (XHOST.vflag_control.flag("WWW")) {
      vpflow.flag("SYMMETRY::SCREEN_ONLY", true);
      XHOST.QUIET = true;
    }
  }
  // DX20170818 - Added tolerance and no_scan options to Xgroups - END
  // DX20200206 - add Patterson symmetry - START
  vpflow.args2addattachedscheme(argv, cmds, "PGROUPK_PATTERSON", "--pointgroupk_Patterson=|--pgroupk_Patterson=", "");
  if (vpflow.flag("PGROUPK_PATTERSON")) {
    vpflow.flag("SYMMETRY::NO_SCAN", aurostd::args2flag(argv, cmds, "--no_scan"));
    if (aurostd::args2attachedflag(argv, "--pointgroupk_Patterson=|--pgroupk_Patterson=")) { // DX20170803
      vpflow.args2addattachedscheme(argv, cmds, "SYMMETRY::TOLERANCE", "--pointgroupk_Patterson=|--pgroupk_Patterson=", ""); // DX20200907 - default is system specific, leaving empty
    }
    vpflow.flag("SYMMETRY::SCREEN_ONLY", aurostd::args2flag(argv, cmds, "--screen_only")); // DX20170803
    // ME20210206 - web mode
    if (XHOST.vflag_control.flag("WWW")) {
      vpflow.flag("SYMMETRY::SCREEN_ONLY", true);
      XHOST.QUIET = true;
    }
  }
  // DX20200206 - add Patterson symmetry - END
  // DX20171205 - Added pgroupk_xtal - START
  vpflow.args2addattachedscheme(argv, cmds, "PGROUPK_XTAL", "--pointgroupkcrystal=|--pgroupk_xtal=", "");
  if (vpflow.flag("PGROUPK_XTAL")) {
    vpflow.flag("SYMMETRY::NO_SCAN", aurostd::args2flag(argv, cmds, "--no_scan"));
    if (aurostd::args2attachedflag(argv, "--pointgroupkcrystal=|--pgroupk_xtal=")) { // DX20170803
      vpflow.args2addattachedscheme(argv, cmds, "SYMMETRY::TOLERANCE", "--pointgroupkcrystal=|--pgroupk_xtal=", ""); // DX20200907 - default is system specific, leaving empty
    }
    vpflow.flag("SYMMETRY::SCREEN_ONLY", aurostd::args2flag(argv, cmds, "--screen_only")); // DX20170803
    // ME20210206 - web mode
    if (XHOST.vflag_control.flag("WWW")) {
      vpflow.flag("SYMMETRY::SCREEN_ONLY", true);
      XHOST.QUIET = true;
    }
  }
  // DX20171205 - Added pgroupk_xtal - END
  vpflow.flag("PLANEDENS", aurostd::args2flag(argv, cmds, "--planedens") && argv.at(1) == "--planedens");
  vpflow.args2addattachedscheme(argv, cmds, "PLATON", "--platon=", "");

  // ME20200313 - Added guards for when plots commands are used without = sign
  if (aurostd::args2flag(argv, cmds, "--plotband|--plotbands|--plot_band|--plot_bands")) {
    vpflow.flag("PLOT_BAND", true);
    vpflow.addattachedscheme("PLOT_BAND", "./", true);
  } else {
    vpflow.args2addattachedscheme(argv, cmds, "PLOT_BAND", "--plotband=|--plotbands=|--plot_band=|--plot_bands=", "./"); // ME20190614
  }
  vpflow.args2addattachedscheme(argv, cmds, "PLOT_BANDSPINSPLIT", "--plotband_spinsplit=|--plot_band_spin_split=|--plotbands_spinsplit=|--plot_bands_spin_split=", "./");
  vpflow.args2addattachedscheme(argv, cmds, "PLOT_BAND2", "--plotband2=|--plot_band2=|--plotbands2=|--plot_bands2=", "./");
  if (aurostd::args2flag(argv, cmds, "--plotbanddos|--plotbandsdos|--plot_band_dos|--plot_bands_dos|--plotdosband|--plotdosbands|--plot_dos_band|--plot_dos_bands")) { // CO20210712
    vpflow.flag("PLOT_BANDDOS", true);
    vpflow.addattachedscheme("PLOT_BANDDOS", "./", true);
  } else {
    vpflow.args2addattachedscheme(argv, cmds, "PLOT_BANDDOS", "--plotbanddos=|--plotbandsdos=|--plot_band_dos=|--plot_bands_dos=|--plotdosband=|--plotdosbands=|--plot_dos_band=|--plot_dos_bands=", "./"); // ME20190614 //CO20210712
  }
  if (aurostd::args2flag(argv, cmds, "--plotdos|--plot_dos")) {
    vpflow.flag("PLOT_DOS", true);
    vpflow.addattachedscheme("PLOT_DOS", "./", true);
  } else {
    vpflow.args2addattachedscheme(argv, cmds, "PLOT_DOS", "--plotdos=|--plot_dos=", "./");
  }
  vpflow.flag("PLOT_ALL_ATOMS", aurostd::args2flag(argv, cmds, "--plot_all_atoms|--plotallatoms")); // CO20191010
  vpflow.args2addattachedscheme(argv, cmds, "PLOT_DOSWEB", "--plotdosweb=|--plot_dos_web=", "./");
  if (aurostd::args2flag(argv, cmds, "--plotpedos|--plotpdos|--plot_pedos|--plot_pdos")) {
    vpflow.flag("PLOT_PDOS", true);
    vpflow.addattachedscheme("PLOT_PDOS", "./", true);
  } else {
    vpflow.args2addattachedscheme(argv, cmds, "PLOT_PDOS", "--plotpedos=|--plotpdos=|--plot_pedos=|--plot_pdos=", "./"); // ME20190614
  }
  if (aurostd::args2flag(argv, cmds, "--plotpedosall|--plotpdosall|--plot_pedos_all|--plot_pdos_all")) {
    vpflow.flag("PLOT_PDOSALL", true);
    vpflow.addattachedscheme("PLOT_PDOSALL", "./", true);
  } else {
    vpflow.args2addattachedscheme(argv, cmds, "PLOT_PDOSALL", "--plotpedosall=|--plotpdosall=|--plot_pedos_all=|--plot_pdos_all=", "./"); // ME20190614
  }
  vpflow.args2addattachedscheme(argv, cmds, "PLOT_PEDOSALL_AFLOWLIB", "--plotpedos_nonequivalent=", "./");
  // ME20190614 BEGIN
  if (aurostd::args2flag(argv, cmds, "--plotthermo|--plot_thermo")) {
    vpflow.flag("PLOT_THERMO", true);
    vpflow.addattachedscheme("PLOT_THERMO", "./", true);
  } else {
    vpflow.args2addattachedscheme(argv, cmds, "PLOT_THERMO", "--plotthermo=|--plot_thermo=", "./");
  }
  if (aurostd::args2flag(argv, cmds, "--plottcond|--plotthermalconductivity|--plot_tcond|--plot_thermal_conductivity")) {
    vpflow.flag("PLOT_TCOND", true);
    vpflow.addattachedscheme("PLOT_TCOND", "./", true);
  } else {
    vpflow.args2addattachedscheme(argv, cmds, "PLOT_TCOND", "--plottcond=|--plotthermalconductivity=|--plot_tcond=|--plot_thermal_conductivity=", "./");
  }
  if (aurostd::args2flag(argv, cmds, "--plotphdos|--plot_phdos")) {
    vpflow.flag("PLOT_PHDOS", true);
    vpflow.addattachedscheme("PLOT_PHDOS", "./", true);
  } else {
    vpflow.args2addattachedscheme(argv, cmds, "PLOT_PHDOS", "--plotphdos=|--plot_phdos=", "./");
  }
  if (aurostd::args2flag(argv, cmds, "--plotphdisp|--plotphonondispersion|--plotphononsdispersion|--pphdis|--plot_phdisp|--plot_phonon_dispersion|--plot_phonons_dispersion")) {
    vpflow.flag("PLOT_PHDISP", true);
    vpflow.addattachedscheme("PLOT_PHDISP", "./", true);
  } else {
    vpflow.args2addattachedscheme(argv, cmds, "PLOT_PHDISP", "--plotphdisp=|--plotphonondispersion=|--plotphononsdispersion=|--pphdis=|--plot_phdisp=|--plot_phonon_dispersion=|--plot_phonons_dispersion=", "./");
  }
  if (aurostd::args2flag(argv, cmds, "--plotphdispdos|--plot_phdispdos")) {
    vpflow.flag("PLOT_PHDISPDOS", true);
    vpflow.addattachedscheme("PLOT_PHDISPDOS", "./", true);
  } else {
    vpflow.args2addattachedscheme(argv, cmds, "PLOT_PHDISPDOS", "--plotphdispdos=|--plot_phdispdos=", "./");
  }
  // Additional DOS/band structure options
  vpflow.flag("PLOTTER::NOSHIFT", aurostd::args2flag(argv, cmds, "--noshift"));
  vpflow.flag("PLOTTER::NOWATERMARK", aurostd::args2flag(argv, cmds, "--nowatermark"));
  vpflow.args2addattachedscheme(argv, cmds, "PLOTTER::PROJECTION", "--projection=", "ORBITALS");
  vpflow.args2addattachedscheme(argv, cmds, "PLOTTER::UNIT", "--unit=", "EV");
  vpflow.args2addattachedscheme(argv, cmds, "PLOTTER::PRINT", "--print=", "pdf");
  vpflow.args2addattachedscheme(argv, cmds, "PLOTTER::TITLE", "--title=", "");
  vpflow.args2addattachedscheme(argv, cmds, "PLOTTER::OUTFILE", "--outfile=", ""); // ME20200313 - user defined output file name
  // ME20190614 END
  // AS20200909 BEGIN
  if (aurostd::args2flag(argv, cmds, "--plotthermoqha")) {
    vpflow.flag("PLOT_THERMO_QHA", true);
    vpflow.addattachedscheme("PLOT_THERMO_QHA", "./", true);
  } else {
    vpflow.args2addattachedscheme(argv, cmds, "PLOT_THERMO_QHA", "--plotthermoqha=", "./");
  }
  vpflow.args2addattachedscheme(argv, cmds, "PLOTTER::EOSMODEL", "--eosmodel=", "SJ"); // AS20210705
  // AS20200909 END
  // AS20210701 BEGIN
  if (aurostd::args2flag(argv, cmds, "--plotgrdisp|--plotgruneisendispersion|--plotgrueneisendispersion")) {
    vpflow.flag("PLOT_GRUENEISEN_DISPERSION", true);
    vpflow.addattachedscheme("PLOT_GRUENEISEN_DISPERSION", "./", true);
  } else {
    vpflow.args2addattachedscheme(argv, cmds, "PLOT_GRUENEISEN_DISPERSION", "--plotgrdisp=", "./");
  }
  // AS20210701 END

  vpflow.flag("POCC", aurostd::args2flag(argv, cmds, "--pocc") && argv.at(1) == "--pocc");
  vpflow.args2addattachedscheme(argv, cmds, "POCC_DOS", "--pocc_dos=", "./");
  vpflow.args2addattachedscheme(argv, cmds, "POCC_MAG", "--pocc_mag=", "./");
  vpflow.args2addattachedscheme(argv, cmds, "POCC_BANDGAP", "--pocc_bandgap=", "./");
  vpflow.args2addattachedscheme(argv, cmds, "POCC_MINIMUM_CONFIGURATION", "--pocc_minimum_configuration|--pocc_min", "./"); // CO20190803

  // ME20181113
  //   vpflow.args2addattachedscheme(argv,cmds,"MODULE","--module=",""); //ME20190112

  // AFLOW modules
  vpflow.args2addattachedscheme(argv, cmds, "AFLOW::MODULE", "--module=", ""); // CO20180214

  vpflow.flag("RENDER", aurostd::args2flag(argv, cmds, "--render"));
  if (vpflow.flag("RENDER")) {
    vpflow.args2addattachedscheme(argv, cmds, "RENDER_OUT", "--directory=|-d=", "");
  }
  vpflow.flag("POSCAR", aurostd::args2flag(argv, cmds, "--poscar"));
  vpflow.flag("POSCAR2AFLOWIN", aurostd::args2flag(argv, cmds, "--poscar2aflowin|--poscaraflowin|--poscar2aflow|--poscaraflow"));
  vpflow.flag("POSCAR2WYCKOFF", aurostd::args2flag(argv, cmds, "--poscar2wyckoff"));
  // CO
  vpflow.args2addattachedscheme(argv, cmds, "PREPARE_CHGCAR_4_JMOL", "--prepare_chgcar_4_jmol=|--prep4jmol=", ""); // can be multiple
  if (vpflow.flag("PREPARE_CHGCAR_4_JMOL")) {
    vpflow.flag("PREPARE_CHGCAR_4_JMOL::USAGE", aurostd::args2flag(argv, cmds, "--usage")); // usage
    vpflow.args2addattachedscheme(argv, cmds, "PREPARE_CHGCAR_4_JMOL::OUTCAR", "--outcar=", ""); // singular
    vpflow.flag("PREPARE_CHGCAR_4_JMOL::ZIP", aurostd::args2flag(argv, cmds, "--zip"));
  }
  vpflow.flag("PRIM", aurostd::args2flag(argv, cmds, "--prim|--prim0"));
  vpflow.flag("PRIM1", aurostd::args2flag(argv, cmds, "--prim1"));
  vpflow.flag("PRIM2", aurostd::args2flag(argv, cmds, "--prim2"));
  vpflow.flag("PRIM3", aurostd::args2flag(argv, cmds, "--prim3"));

  vpflow.args2addattachedscheme(argv, cmds, "PROTO", "--proto=|--proto_icsd=", ""); // --proto=123:A:B:C --proto_icsd=Gd1Mn2Si2_ICSD_54947
  vpflow.args2addattachedscheme(argv, cmds, "PARAMS", "--params=|--parameters=", ""); // --proto=123:A:B:C --proto_icsd=Gd1Mn2Si2_ICSD_54947
  vpflow.args2addattachedscheme(argv, cmds, "POCC_PARAMS", "--pocc_params=", ""); // --pocc_params=S0-1xC_S1-0.5xE-0.5xF_S2-0.3333xA-0.3333xB-0.3333xD
  vpflow.args2addattachedscheme(argv, cmds, "POCC_TOL", "--pocc_tol=", ""); // --pocc_params=S0-1xC_S1-0.5xE-0.5xF_S2-0.3333xA-0.3333xB-0.3333xD
  if (vpflow.flag("POCC_PARAMS")) {
    if (LDEBUG) {
      cerr << "PflowARGs(): BEFORE pocc_params=" << vpflow.getattachedscheme("POCC_PARAMS") << endl;
    }
    vpflow.push_attached("POCC_PARAMS", vpflow.getattachedscheme("POCC_PARAMS")); // CO20190126 - do later //pflow::FIX_POCC_PARAMS()
    if (LDEBUG) {
      cerr << "PflowARGs(): AFTER  pocc_params=" << vpflow.getattachedscheme("POCC_PARAMS") << endl;
    }
  }
  vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW", "--aflow_proto=|--aflow_proto_icsd=", ""); // --aflow_proto=123:A:B:C  --aflow_proto_icsd=Gd1Mn2Si2_ICSD_54947
  if (vpflow.flag("PROTO") || vpflow.flag("PROTO_AFLOW")) { // CO20180622 - set these flags if structure is downloaded from web
    vpflow.flag("PROTO::VASP", aurostd::args2flag(argv, cmds, "--vasp") || (!aurostd::args2flag(argv, cmds, "--qe") && !aurostd::args2flag(argv, cmds, "--abinit") && !aurostd::args2flag(argv, cmds, "--aims") &&
                                                                            !aurostd::args2flag(argv, cmds, "--cif") && !aurostd::args2flag(argv, cmds, "--elk") && aurostd::args2flag(argv, cmds, "--lmp"))); // DX20190131 //DX20200313 - added elk /SD20240111 - added lmp
    vpflow.flag("PROTO::ITC", aurostd::args2flag(argv, cmds, "--itc")); // CO20220613
    vpflow.flag("PROTO::QE", aurostd::args2flag(argv, cmds, "--qe"));
    vpflow.flag("PROTO::ABCCAR", aurostd::args2flag(argv, cmds, "--abccar")); // DX20190123 - add abccar output
    vpflow.flag("PROTO::ABINIT", aurostd::args2flag(argv, cmds, "--abinit"));
    vpflow.flag("PROTO::AIMS", aurostd::args2flag(argv, cmds, "--aims"));
    vpflow.flag("PROTO::CIF", aurostd::args2flag(argv, cmds, "--cif")); // DX20190123 - add cif output
    vpflow.flag("PROTO::ELK", aurostd::args2flag(argv, cmds, "--elk")); // DX20200313 - add elk output
    vpflow.flag("PROTO::LMP", aurostd::args2flag(argv, cmds, "--lmp")); // SD20240111 - add lmp output
    vpflow.flag("PROTO::HEX", aurostd::args2flag(argv, cmds, "--hex"));
    vpflow.flag("PROTO::RHL", aurostd::args2flag(argv, cmds, "--rhl"));
    vpflow.flag("PROTO::ADD_EQUATIONS", aurostd::args2flag(argv, cmds, "--add_equations")); // DX20180615 - add equation info
    vpflow.flag("PROTO::EQUATIONS_ONLY", aurostd::args2flag(argv, cmds, "--equations_only")); // DX20180615 - add equation info
    vpflow.flag("PROTO::SYMBOLS_ONLY", aurostd::args2flag(argv, cmds, "--parameter_symbols_only")); // HE20250524 - add parameter info
    vpflow.flag("PROTO::WEBPAGE", aurostd::args2flag(argv, cmds, "--webpage")); // HE20250524 - allow different rules for the webpage generation
    vpflow.flag("PROTO::STRICT", aurostd::args2flag(argv, cmds, "--strict")); // HE20250527 - enforce atom count check
    vpflow.flag("PROTO::USE_ANRL_LATTICE_PARAM", aurostd::args2flag(argv, cmds, "--use_anrl_lattice_param|--use_anrl_lattice_parameter")); // DX20190227 - add original lattice parameter option
    vpflow.flag("BADER", false); // CO20181226 don't run bader
    vpflow.flag("ITC", false); // PRIORITIES  //CO20220613
    vpflow.flag("QE", false); // PRIORITIES
    vpflow.flag("ABCCAR", false); // PRIORITIES //DX20190123 - add ABCCAR
    vpflow.flag("ABINIT", false); // PRIORITIES
    vpflow.flag("AIMS", false); // PRIORITIES
    vpflow.flag("CIF", false); // PRIORITIES //DX20190123 - add CIF
    vpflow.flag("ELK", false); // PRIORITIES //DX20200313 - add ELK
    vpflow.flag("LMP", false); // PRIORITIES //SD20240111 - add LMP
  }
  //[CO20180627 - MOVED UP]vpflow.args2addattachedscheme(argv,cmds,"PROTO_AFLOW","--aflow_proto=|--aflow_proto_icsd=","");      // --aflow_proto=123:A:B:C  --aflow_proto_icsd=Gd1Mn2Si2_ICSD_54947
  if (vpflow.flag("PROTO_AFLOW")) {
    string vlist;
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::PRESSURE", "--pressure=|--PRESSURE=|--pstress=|--PSTRESS=", "");
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::POTENTIAL", "--potential=|--pp=|--potentials=", "");
    if (vpflow.flag("PROTO_AFLOW::POTENTIAL")) {
      vlist += "--potential=" + vpflow.getattachedscheme("PROTO_AFLOW::POTENTIAL") + " ";
    }

    // AFLOW modules
    if (vpflow.flag("AFLOW::MODULE")) {
      vlist += "--module=" + vpflow.getattachedscheme("AFLOW::MODULE") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::APL_SUPERCELL", "--apl_supercell=", ""); // CO20180214
    if (vpflow.flag("PROTO_AFLOW::APL_SUPERCELL")) {
      vlist += "--apl_supercell=" + vpflow.getattachedscheme("PROTO_AFLOW::APL_SUPERCELL") + " ";
    }

    vpflow.flag("PROTO_AFLOW::BADER", aurostd::args2flag(argv, cmds, "--bader|--run_bader|--BADER|--RUN_BADER")); // CO20180214
    if (vpflow.flag("PROTO_AFLOW::BADER")) {
      vlist += "--bader "; // CO20180214
    }
    vpflow.flag("PROTO_AFLOW::SPIN_REMOVE_RELAX_1", aurostd::args2flag(argv, cmds, "--spin_remove_relax_1|--SPIN_REMOVE_RELAX_1")); // CO20180214
    if (vpflow.flag("PROTO_AFLOW::SPIN_REMOVE_RELAX_1")) {
      vlist += "--spin_remove_relax_1 "; // CO20180214
    }
    vpflow.flag("PROTO_AFLOW::SPIN_REMOVE_RELAX_2", aurostd::args2flag(argv, cmds, "--spin_remove_relax_2|--SPIN_REMOVE_RELAX_2")); // CO20180214
    if (vpflow.flag("PROTO_AFLOW::SPIN_REMOVE_RELAX_2")) {
      vlist += "--spin_remove_relax_2 "; // CO20180214
    }

    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::POTIM", "--potim=|--POTIM=", "");
    if (vpflow.flag("PROTO_AFLOW::POTIM")) {
      vlist += "--potim=" + vpflow.getattachedscheme("PROTO_AFLOW::POTIM") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::RELAX_TYPE", "--relax_type=|--RELAX_TYPE=", ""); // CO20180214
    if (vpflow.flag("PROTO_AFLOW::RELAX_TYPE")) {
      vlist += "--relax_type=" + vpflow.getattachedscheme("PROTO_AFLOW::RELAX_TYPE") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::RELAX_MODE", "--relax_mode=|--RELAX_MODE=", "");
    if (vpflow.flag("PROTO_AFLOW::RELAX_MODE")) {
      vlist += "--relax_mode=" + vpflow.getattachedscheme("PROTO_AFLOW::RELAX_MODE") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::RELAX_COUNT", "--relax_count=|--RELAX_COUNT=", ""); // CO20181226
    if (vpflow.flag("PROTO_AFLOW::RELAX_COUNT")) {
      vlist += "--relax_count=" + vpflow.getattachedscheme("PROTO_AFLOW::RELAX_COUNT") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::RUN_RELAX_STATIC", "--run_relax_static|--RUN_RELAX_STATIC|--relax_static", ""); // CO20181226
    if (vpflow.flag("PROTO_AFLOW::RUN_RELAX_STATIC")) {
      vlist += "--run_relax_static ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::RUN_RELAX_STATIC_BANDS", "--run_relax_static_bands|--RUN_RELAX_STATIC_BANDS|--relax_static_bands|--bands|--band", ""); // CO20181226
    if (vpflow.flag("PROTO_AFLOW::RUN_RELAX_STATIC_BANDS")) {
      vlist += "--run_relax_static_bands ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::RUN_RELAX_STATIC_DIELECTRIC", "--run_relax_static_dielectric|--RUN_RELAX_STATIC_DIELECTRIC|--relax_static_dielectric|--dielectric", ""); // CO20181226
    if (vpflow.flag("PROTO_AFLOW::RUN_RELAX_STATIC_DIELECTRIC")) {
      vlist += "--run_relax_static_dielectric ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::RUN_STATIC", "--run_static|--RUN_STATIC|--static", ""); // CO20181226
    if (vpflow.flag("PROTO_AFLOW::RUN_STATIC")) {
      vlist += "--run_static ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::RUN_STATIC_BANDS", "--run_static_bands|--RUN_STATIC_BANDS|--static_bands", ""); // CO20181226
    if (vpflow.flag("PROTO_AFLOW::RUN_STATIC_BANDS")) {
      vlist += "--run_static_bands ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::RUN_STATIC_DIELECTRIC", "--run_static_dielectric|--RUN_STATIC_DIELECTRIC|--static_dielectric", ""); // CO20181226
    if (vpflow.flag("PROTO_AFLOW::RUN_STATIC_DIELECTRIC")) {
      vlist += "--run_static_dielectric ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::PRECISION", "--prec=|--PREC=|--precision=|--PRECISION=", "");
    if (vpflow.flag("PROTO_AFLOW::PRECISION")) {
      vlist += "--precision=" + vpflow.getattachedscheme("PROTO_AFLOW::PRECISION") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::ALGORITHM", "--algo=|--ALGO=|--algorithm=|--ALGORITHM=", "");
    if (vpflow.flag("PROTO_AFLOW::ALGORITHM")) {
      vlist += "--algorithm=" + vpflow.getattachedscheme("PROTO_AFLOW::ALGORITHM") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::METAGGA", "--metagga=|--METAGGA=", "");
    if (vpflow.flag("PROTO_AFLOW::METAGGA")) {
      vlist += "--metagga=" + vpflow.getattachedscheme("PROTO_AFLOW::METAGGA") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::IVDW", "--ivdw=|--IVDW=", "");
    if (vpflow.flag("PROTO_AFLOW::IVDW")) {
      vlist += "--ivdw=" + vpflow.getattachedscheme("PROTO_AFLOW::IVDW") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::TYPE", "--type=|--TYPE=", "");
    if (vpflow.flag("PROTO_AFLOW::TYPE")) {
      vlist += "--type=" + vpflow.getattachedscheme("PROTO_AFLOW::TYPE") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::CONVERT_UNIT_CELL", "--convert_unit_cell=|--CONVERT_UNIT_CELL=", "");
    if (vpflow.flag("PROTO_AFLOW::CONVERT_UNIT_CELL")) {
      vlist += "--convert_unit_cell=" + vpflow.getattachedscheme("PROTO_AFLOW::CONVERT_UNIT_CELL") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::VOLUME_PLUS_EQUAL", "--volume_plus_equal=|--VOLUME_PLUS_EQUAL=", "");
    if (vpflow.flag("PROTO_AFLOW::VOLUME_PLUS_EQUAL")) {
      vlist += "--volume_plus_equal=" + vpflow.getattachedscheme("PROTO_AFLOW::VOLUME_PLUS_EQUAL") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::VOLUME_MULTIPLY_EQUAL", "--volume_multiply_equal=|--VOLUME_MULTIPLY_EQUAL=", "");
    if (vpflow.flag("PROTO_AFLOW::VOLUME_MULTIPLY_EQUAL")) {
      vlist += "--volume_multiply_equal=" + vpflow.getattachedscheme("PROTO_AFLOW::VOLUME_MULTIPLY_EQUAL") + " ";
    }
    vpflow.flag("PROTO_AFLOW::VOLUME_PRESERVED", aurostd::args2flag(argv, cmds, "--volume_preserved|--no_volume_adjustment")); // CO20180214
    if (vpflow.flag("PROTO_AFLOW::VOLUME_PRESERVED")) {
      vlist += "--volume_preserved "; // CO20180214
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::EDIFFG", "--ediffg=|--EDIFFG=", "");
    if (vpflow.flag("PROTO_AFLOW::EDIFFG")) {
      vlist += "--ediffg=" + vpflow.getattachedscheme("PROTO_AFLOW::EDIFFG") + " ";
    }

    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::KSCHEME", "--kscheme=|--KSCHEME=", ""); // CO20181226
    if (vpflow.flag("PROTO_AFLOW::KSCHEME")) {
      vlist += "--kscheme=" + vpflow.getattachedscheme("PROTO_AFLOW::KSCHEME") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::KSCHEME_STATIC", "--kscheme_static=|--KSCHEME_STATIC=", ""); // CO20181226
    if (vpflow.flag("PROTO_AFLOW::KSCHEME_STATIC")) {
      vlist += "--kscheme_static=" + vpflow.getattachedscheme("PROTO_AFLOW::KSCHEME_STATIC") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::KSCHEME_DIELECTRIC", "--kscheme_dielectric=|--KSCHEME_DIELECTRIC=", "");
    if (vpflow.flag("PROTO_AFLOW::KSCHEME_DIELECTRIC")) {
      vlist += "--kscheme_dielectric=" + vpflow.getattachedscheme("PROTO_AFLOW::KSCHEME_DIELECTRIC") + " ";
    }

    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::KPPRA", "--kppra=|--KPPRA=", "");
    if (vpflow.flag("PROTO_AFLOW::KPPRA") && vpflow.flag("KPOINTS")) {
      vpflow.flag("KPOINTS", false);
    } // CO20200223 - --kppra collision
    if (vpflow.flag("PROTO_AFLOW::KPPRA")) {
      vlist += "--kppra=" + vpflow.getattachedscheme("PROTO_AFLOW::KPPRA") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::KPPRA_STATIC", "--kppra_static=|--KPPRA_STATIC=", "");
    if (vpflow.flag("PROTO_AFLOW::KPPRA_STATIC")) {
      vlist += "--kppra_static=" + vpflow.getattachedscheme("PROTO_AFLOW::KPPRA_STATIC") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::KPPRA_DIELECTRIC", "--kppra_dielectric=|--KPPRA_DIELECTRIC=", "");
    if (vpflow.flag("PROTO_AFLOW::KPPRA_DIELECTRIC")) {
      vlist += "--kppra_dielectric=" + vpflow.getattachedscheme("PROTO_AFLOW::KPPRA_DIELECTRIC") + " ";
    }
    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::BANDS_GRID", "--bands_grid=|--BANDS_GRID=", "");
    if (vpflow.flag("PROTO_AFLOW::BANDS_GRID")) {
      vlist += "--bands_grid=" + vpflow.getattachedscheme("PROTO_AFLOW::BANDS_GRID") + " ";
    }

    vpflow.args2addattachedscheme(argv, cmds, "PROTO_AFLOW::ENMAX_MULTIPLY", "--enmax_multiply=|--ENMAX_MULTIPLY=", "");
    if (vpflow.flag("PROTO_AFLOW::ENMAX_MULTIPLY")) {
      vlist += "--enmax_multiply=" + vpflow.getattachedscheme("PROTO_AFLOW::ENMAX_MULTIPLY") + " ";
    }

    vpflow.flag("PROTO_AFLOW::USAGE", aurostd::args2flag(argv, cmds, "--usage"));
    vpflow.flag("PROTO_AFLOW::POTENTIAL_TYPE", aurostd::args2flag(argv, cmds, "--potential_type|--pottype|--potentials_type|--potstype|--pott")); // CO20191110
    if (vpflow.flag("PROTO_AFLOW::POTENTIAL_TYPE")) {
      vlist += "--potential_type ";
    }
    vpflow.flag("PROTO_AFLOW::POTENTIAL_COMPLETE", aurostd::args2flag(argv, cmds, "--potential_complete|--potcomplete|--potentials_complete|--potscomplete|--potc|--potential_type_date|--potentialtypedate")); // CO20191110
    if (vpflow.flag("PROTO_AFLOW::POTENTIAL_COMPLETE")) {
      vlist += "--potential_complete ";
    }
    vpflow.flag("PROTO_AFLOW::MISSING", aurostd::args2flag(argv, cmds, "--missing"));
    if (vpflow.flag("PROTO_AFLOW::MISSING")) {
      vlist += "--missing ";
    }
    vpflow.flag("PROTO_AFLOW::NOAUTOPP", aurostd::args2flag(argv, cmds, "--noautopp"));
    if (vpflow.flag("PROTO_AFLOW::NOAUTOPP")) {
      vlist += "--noautopp ";
    }
    vpflow.flag("PROTO_AFLOW::LDAU", aurostd::args2flag(argv, cmds, "--ldau|--ldau2"));
    if (vpflow.flag("PROTO_AFLOW::LDAU")) {
      vlist += "--ldau ";
    }
    vpflow.flag("PROTO_AFLOW::NOLDAU", aurostd::args2flag(argv, cmds, "--noldau|--noldau2"));
    if (vpflow.flag("PROTO_AFLOW::NOLDAU")) {
      vlist += "--noldau ";
    }
    vpflow.flag("PROTO_AFLOW::NEGLECT_NOMIX", aurostd::args2flag(argv, cmds, "--neglect_nomix|--neglectnomix"));
    if (vpflow.flag("PROTO_AFLOW::NEGLECT_NOMIX")) {
      vlist += "--neglect_nomix ";
    }
    vpflow.flag("PROTO_AFLOW::STDOUT", aurostd::args2flag(argv, cmds, "--stdout"));
    if (vpflow.flag("PROTO_AFLOW::STDOUT")) {
      vlist += "--stdout ";
    }
    vpflow.flag("PROTO_AFLOW::QE", aurostd::args2flag(argv, cmds, "--qe"));
    if (vpflow.flag("PROTO_AFLOW::QE")) {
      vlist += "--qe ";
    }
    vpflow.flag("PROTO_AFLOW::ABCCAR", aurostd::args2flag(argv, cmds, "--abccar")); // DX20190131                                                                     //DX20190123 - added CIF
    if (vpflow.flag("PROTO_AFLOW::ABCCAR")) {
      vlist += "--abccar "; // DX20190131
    }
    vpflow.flag("PROTO_AFLOW::ABINIT", aurostd::args2flag(argv, cmds, "--abinit"));
    if (vpflow.flag("PROTO_AFLOW::ABINIT")) {
      vlist += "--abinit ";
    }
    vpflow.flag("PROTO_AFLOW::AIMS", aurostd::args2flag(argv, cmds, "--aims"));
    if (vpflow.flag("PROTO_AFLOW::AIMS")) {
      vlist += "--aims ";
    }
    vpflow.flag("PROTO_AFLOW::CIF", aurostd::args2flag(argv, cmds, "--cif")); // DX20190131                                                                           //DX20190123 - added CIF
    if (vpflow.flag("PROTO_AFLOW::CIF")) {
      vlist += "--cif "; // DX20190131
    }
    vpflow.flag("PROTO_AFLOW::ELK", aurostd::args2flag(argv, cmds, "--elk")); // DX20200313 - added ELK
    if (vpflow.flag("PROTO_AFLOW::ELK")) {
      vlist += "--elk ";
    }
    vpflow.flag("PROTO_AFLOW::LMP", aurostd::args2flag(argv, cmds, "--lmp")); // SD20240111 - added LMP
    if (vpflow.flag("PROTO_AFLOW::LMP")) {
      vlist += "--lmp ";
    }
    vpflow.flag("PROTO_AFLOW::HEX", aurostd::args2flag(argv, cmds, "--hex"));
    if (vpflow.flag("PROTO_AFLOW::HEX")) {
      vlist += "--hex ";
    }
    vpflow.flag("PROTO_AFLOW::RHL", aurostd::args2flag(argv, cmds, "--rhl"));
    if (vpflow.flag("PROTO_AFLOW::RHL")) {
      vlist += "--rhl ";
    }
    vpflow.flag("PROTO_AFLOW::VASP", aurostd::args2flag(argv, cmds, "--vasp"));
    if (vpflow.flag("PROTO_AFLOW::VASP")) {
      vlist += "--vasp ";
    }
    vpflow.flag("PROTO_AFLOW::ITC", aurostd::args2flag(argv, cmds, "--itc")); // CO20220613
    if (vpflow.flag("PROTO_AFLOW::ITC")) {
      vlist += "--itc ";
    }
    vpflow.flag("PROTO_AFLOW::HTQC", aurostd::args2flag(argv, "--htqc"));
    if (vpflow.flag("PROTO_AFLOW::HTQC")) {
      vlist += "--htqc ";
    }
    vpflow.flag("PROTO_AFLOW::LIST", aurostd::args2flag(argv, cmds, "--list"));
    //     if(vpflow.flag("PROTO_AFLOW::LIST")) vlist+="--listDEBUG ";
    if (vpflow.flag("PROTO_AFLOW::LIST")) {
      vpflow.push_attached("PROTO_AFLOW::LIST_VCMD", vlist);
    }
    vpflow.flag("KPOINTS", false); // PRIORITIES
    vpflow.flag("ITC", false); // PRIORITIES  //CO20220613
    vpflow.flag("QE", false); // PRIORITIES
    vpflow.flag("ABCCAR", false); // PRIORITIES //DX20190123 - add ABCCAR
    vpflow.flag("ABINIT", false); // PRIORITIES
    vpflow.flag("AIMS", false); // PRIORITIES
    vpflow.flag("CIF", false); // PRIORITIES //DX20190123 - add CIF
    vpflow.flag("ELK", false); // PRIORITIES //DX20200313 - add ELK
    vpflow.flag("LMP", false); // PRIORITIES //SD20240111 - add LMP
  }

  vpflow.flag("PROTOINITIAL", aurostd::args2flag(argv, cmds, "--initial"));

  // PSEUDOPOTENTIAL CHECK
  vpflow.args2addattachedscheme(argv, cmds, "PSEUDOPOTENTIALS_CHECK", "--pseudopotentials_check=|--pp_check=|--ppk=", "");
  if (vpflow.flag("PSEUDOPOTENTIALS_CHECK")) {
    vpflow.flag("PSEUDOPOTENTIALS_CHECK::USAGE", aurostd::args2flag(argv, cmds, "--usage"));
  }

  // ME20211103
  if (aurostd::args2flag(argv, cmds, "--python_modules|--create_python_modules")) {
    vpflow.flag("PYTHON_MODULES", true);
    vpflow.addattachedscheme("PYTHON_MODULES", "", true);
  } else {
    vpflow.args2addattachedscheme(argv, cmds, "PYTHON_MODULES", "--python_modules=|--create_python_modules=", "");
  }

  // MOVE ON
  vpflow.flag("ITC", aurostd::args2flag(argv, cmds, "--itc") && !vpflow.flag("PROTO_AFLOW") && !vpflow.flag("PROTO")); // CO20220613
  vpflow.flag("QE", aurostd::args2flag(argv, cmds, "--qe") && !vpflow.flag("PROTO_AFLOW") && !vpflow.flag("PROTO"));
  vpflow.flag("ABCCAR", aurostd::args2flag(argv, cmds, "--abccar") && !vpflow.flag("PROTO_AFLOW") && !vpflow.flag("PROTO")); // DX20190123 - moved ABCCAR to here
  vpflow.flag("ABINIT", aurostd::args2flag(argv, cmds, "--abinit") && !vpflow.flag("PROTO_AFLOW") && !vpflow.flag("PROTO"));
  vpflow.flag("AIMS", aurostd::args2flag(argv, cmds, "--aims") && !vpflow.flag("PROTO_AFLOW") && !vpflow.flag("PROTO"));
  vpflow.flag("ELK", aurostd::args2flag(argv, cmds, "--elk") && !vpflow.flag("PROTO_AFLOW") && !vpflow.flag("PROTO")); // DX20200313
  vpflow.flag("ATAT", aurostd::args2flag(argv, cmds, "--atat") && !vpflow.flag("PROTO_AFLOW") && !vpflow.flag("PROTO")); // SD20220123
  vpflow.flag("LMP", aurostd::args2flag(argv, cmds, "--lmp") && !vpflow.flag("PROTO_AFLOW") && !vpflow.flag("PROTO")); // SD20240111

  // DX20190123 - CIF to work with PROTO - START
  // CO20220715 - rewriting a bit for clarity
  if (!(vpflow.flag("PROTO_AFLOW") || vpflow.flag("PROTO"))) { // DX20190123 - to differentiate between proto and istream  //DX20190123 - ensure not a prototype run
    vpflow.flag("CIF", aurostd::args2flag(argv, cmds, "--cif"));
    vpflow.args2addattachedscheme(argv, cmds, "CIF", "--cif=|--CIF=", "");
    if (vpflow.flag("CIF")) {
      vpflow.flag("CIF::USAGE", aurostd::args2flag(argv, cmds, "--usage|--USAGE"));
      vpflow.flag("CIF::NO_SYMMETRY", aurostd::args2flag(argv, cmds, "--no_symmetry|--no_sym"));
      vpflow.flag("CIF::NO_SCAN", aurostd::args2flag(argv, cmds, "--no_scan")); // DX20210610
      if (aurostd::args2attachedflag(argv, "--cif=|--CIF=")) { // DX20170803
        vpflow.args2addattachedscheme(argv, cmds, "CIF::TOLERANCE", "--cif=|--CIF=", ""); // DX20200907 - default is system specific, leaving empty
      }
      if (aurostd::args2attachedflag(argv, "--setting=")) {
        vpflow.args2addattachedscheme(argv, cmds, "CIF::SETTING", "--setting=", "1");
      }
    }
  }
  // DX20190123 - CIF to work with PROTO - END
  vpflow.args2addattachedscheme(argv, cmds, "QSUB", "--qsub=|--qstart=|--bsub=|--sbatch=", "");
  vpflow.args2addattachedscheme(argv, cmds, "QDEL", "--qdel=|--scancel=|--bkill=", "");

  vpflow.flag("QMVASP", aurostd::args2flag(argv, cmds, "--qmvasp"));

  vpflow.flag("PFLOW::QUEUE_STATUS", aurostd::args2flag(argv, cmds, "--queue_status|--queue|--q")); // CO20200526

  vpflow.flag("QCA::INIT", aurostd::args2flag(argv, cmds, "--quasi_chem_approx|--qca")); // SD20220323 - initiate quasi-chemical approx calculation
  if (vpflow.flag("QCA::INIT")) {
    vpflow.flag("QCA::USAGE", aurostd::args2flag(argv, cmds, "--usage"));
    vpflow.flag("QCA::BINODAL", aurostd::args2flag(argv, cmds, "--binodal"));
    vpflow.flag("QCA::USE_SG", aurostd::args2flag(argv, cmds, "--use_sg"));
    vpflow.flag("QCA::SCREEN_ONLY", aurostd::args2flag(argv, cmds, "--screen_only"));
    vpflow.flag("QCA::IMAGE_ONLY", aurostd::args2flag(argv, cmds, "--image_only|--image"));
    vpflow.args2addattachedscheme(argv, cmds, "QCA::PRINT", "--print=|--p=|--output=|--o=", "");
    vpflow.args2addattachedscheme(argv, cmds, "QCA::DIRECTORY", "--directory=", "./");
    vpflow.args2addattachedscheme(argv, cmds, "QCA::PLATTICE", "--plattice=|--plat=", "");
    vpflow.args2addattachedscheme(argv, cmds, "QCA::ELEMENTS", "--elements=|--elem=", "");
    vpflow.args2addattachedscheme(argv, cmds, "QCA::MAX_NUM_ATOMS", "--max_num_atoms=|--mna=", "");
    vpflow.args2addattachedscheme(argv, cmds, "QCA::CV_CUTOFF", "--cv_cutoff=|--cv_cut=", "");
    vpflow.args2addattachedscheme(argv, cmds, "QCA::CONC_CURVE_RANGE", "--conc_curve_range=|--conc_curve=", "");
    vpflow.args2addattachedscheme(argv, cmds, "QCA::CONC_NPTS", "--conc_npts=", "");
    vpflow.args2addattachedscheme(argv, cmds, "QCA::TEMP_RANGE", "--temp_range=|--temp=", "");
    vpflow.args2addattachedscheme(argv, cmds, "QCA::TEMP_NPTS", "--temp_npts=", "");
    vpflow.args2addattachedscheme(argv, cmds, "QCA::AFLOWLIB_DIRECTORY", "--aflowlib_directory=|--aflowlib_dir=", "");
    vpflow.args2addattachedscheme(argv, cmds, "QCA::AFLOW_MAX_NUM_ATOMS", "--aflow_max_num_atoms=", "");
  }

  vpflow.flag("SOLIQUIDY::INIT", aurostd::args2flag(argv, cmds, "--soliquidy")); ////HE20240528
  if (vpflow.flag("SOLIQUIDY::INIT")) {
    vpflow.args2addattachedscheme(argv, cmds, "SOLIQUIDY::OUTPUT", "--out=|--outfolder=|--output=|-o=", "");
    vpflow.args2addattachedscheme(argv, cmds, "SOLIQUIDY::WORKLIST", "--worklist=|--work=", "");
    vpflow.args2addattachedscheme(argv, cmds, "SOLIQUIDY::AUID", "--auid=", "");
    vpflow.flag("SOLIQUIDY::X3D", aurostd::args2flag(argv, cmds, "--x3d|--3d|--render"));
  }

  vpflow.args2addattachedscheme(argv, cmds, "RASMOL", "--rasmol=", "");

  vpflow.flag("RAYTRACE", (aurostd::args2flag(argv, cmds, "--raytrace") && argv.at(1) == "--raytrace"));
  vpflow.flag("RBANAL", aurostd::args2flag(argv, cmds, "--rbanal") && argv.at(1) == "--rbanal");
  vpflow.flag("RBDIST", aurostd::args2flag(argv, cmds, "--rbdist") && argv.at(1) == "--rbdist");
  vpflow.flag("RDF", aurostd::args2flag(argv, cmds, "--rdf")); // CO20220711
  if (!vpflow.flag("RDF")) {
    vpflow.args2addattachedscheme(argv, cmds, "RDF", "--rdf=", "");
  } // CO20220711 - also look for '=' variant, must be XOR with flag (without '=')
  vpflow.flag("RDF::RAW_COUNTS", aurostd::args2flag(argv, cmds, "--raw_counts")); // CO20220627
  vpflow.args2addattachedscheme(argv, cmds, "RDFCMP", "--rdfcmp=", "");

  vpflow.flag("RMATOM", aurostd::args2flag(argv, cmds, "--rm_atom") && argv.at(1) == "--rm_atom");
  vpflow.flag("RMCOPIES", aurostd::args2flag(argv, cmds, "--rm_copies") && argv.at(1) == "--rm_copies");
  vpflow.flag("RSM", aurostd::args2flag(argv, cmds, "--rsm"));

  vpflow.args2addattachedscheme(argv, cmds, "SCALE", "--scale=", "0.0");

  vpflow.flag("SD", (aurostd::args2flag(argv, cmds, "--sd") && argv.at(1) == "--sd"));
  vpflow.flag("SETCM", (aurostd::args2flag(argv, cmds, "--setcm") && argv.at(1) == "--setcm"));
  vpflow.flag("SETORIGIN", (aurostd::args2flag(argv, cmds, "--setorigin") && argv.at(1) == "--setorigin"));
  vpflow.flag("SEWALD", aurostd::args2flag(argv, cmds, "--sewald") && argv.at(1) == "--sewald");
  vpflow.args2addattachedscheme(argv, cmds, "SHELL", "--shell=", "");
  vpflow.args2addattachedscheme(argv, cmds, "SHIFT", "--shift=", "");
  // DX20170818 - Added tolerance and no_scan options to Xgroups - START
  vpflow.args2addattachedscheme(argv, cmds, "SGROUP", "--spacegroup=|--sgroup=", "");
  if (vpflow.flag("SGROUP")) {
    vpflow.flag("SYMMETRY::NO_SCAN", aurostd::args2flag(argv, cmds, "--no_scan"));
    if (aurostd::args2attachedflag(argv, "--spacegroup=|--sgroup=")) { // DX20170803
      vpflow.args2addattachedscheme(argv, cmds, "SYMMETRY::TOLERANCE", "--spacegroup=|--sgroup=", ""); // DX20200907 - default is system specific, leaving empty
    }
    vpflow.args2addattachedscheme(argv, cmds, "SYMMETRY::SGROUP_RADIUS", "--radius=", ""); // DX20170803
    vpflow.flag("SYMMETRY::SCREEN_ONLY", aurostd::args2flag(argv, cmds, "--screen_only")); // DX20170803
    // DX20170921 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv, cmds, "SYMMETRY::MAGNETIC", "--mag=|--magnetic=|--magmom=", ""); // DX20170803
    // DX20170921 - MAGNETIC SYMMETRY - END
    //  ME20210206 - web mode
    if (XHOST.vflag_control.flag("WWW")) {
      vpflow.flag("SYMMETRY::SCREEN_ONLY", true);
      XHOST.QUIET = true;
    }
  }
  // DX20170818 - Added tolerance and no_scan options to Xgroups - END
  //  vpflow.flag("SPLINE",aurostd::args2flag(argv,cmds,"--spline") && argv.at(1)=="--spline");

  if (!vpflow.flag("LIST_PROTOTYPE_LABELS")) { // DX20210708 - protect against other commands
    vpflow.args2addattachedscheme(argv, cmds, "SG::AFLOW", "--aflowSG=|--space_group=|--sg=", ""); // DX20210611 - added aliases
    vpflow.args2addattachedscheme(argv, cmds, "SG::AFLOW_LABEL", "--aflowSG_label=|--space_group_label=|--sg_label=", ""); // DX20210611 - added aliases
    vpflow.args2addattachedscheme(argv, cmds, "SG::AFLOW_NUMBER", "--aflowSG_number=|--space_group_number=|--sg_number=", ""); // DX20210611 - added aliases
  }
  // DX20170926 - Create flags for SG functions - START
  if (vpflow.flag("SG::AFLOW") || vpflow.flag("SG::AFLOW_LABEL") || vpflow.flag("SG::AFLOW_NUMBER")) {
    vpflow.flag("SG::NO_SCAN", aurostd::args2flag(argv, cmds, "--no_scan"));
    if (aurostd::args2attachedflag(argv, "--aflowSG=|--aflowSG_label=|--aflowSG_number=|--space_group=|--space_group_label=|--space_group_number=|--sg=|--sg_label=|--sg_number=")) { // DX20210611 - added aliases
      vpflow.args2addattachedscheme(argv, cmds, "SG::TOLERANCE", "--aflowSG=|--aflowSG_label=|--aflowSG_number=|--space_group=|--space_group_label=|--space_group_number=|--sg=|--sg_label=|--sg_number=", ""); // DX20200907 - default is system specific, leaving empty //DX20210611 - added aliases
    }
    if (aurostd::args2attachedflag(argv, "--tolerance_spectrum=|--tol_spectrum=")) {
      vpflow.args2addattachedscheme(argv, cmds, "SG::TOLERANCE_SPECTRUM", "--tolerance_spectrum=|--tol_spectrum=", "0.01,1.0,100"); // DX20200907 - added default
    }
    // DX20170921 - MAGNETIC SYMMETRY - START
    vpflow.args2addattachedscheme(argv, cmds, "SG::MAGNETIC", "--mag=|--magnetic=|--magmom=", ""); // DX20170803
    // DX20170921 - MAGNETIC SYMMETRY - END
  }
  // DX20170926 - Create flags for SG functions - END

  vpflow.args2addattachedscheme(argv, cmds, "SG::PLATON", "--platonSG=", "");
  vpflow.args2addattachedscheme(argv, cmds, "SG::PLATON_LABEL", "--platonSG_label=", "");
  vpflow.args2addattachedscheme(argv, cmds, "SG::PLATON_NUMBER", "--platonSG_number=", "");
  // DX20170926 - Create flags for SG functions - START
  if (vpflow.flag("SG::PLATON") || vpflow.flag("SG::PLATON_LABEL") || vpflow.flag("SG::PLATON_NUMBER")) {
    if (aurostd::args2attachedflag(argv, "--platonSG=|--platonSG_label=|--platonSG_number=")) {
      vpflow.args2addattachedscheme(argv, cmds, "SG::TOLERANCE", "--platonSG=|--platonSG_label=|--platonSG_number=", "1");
    }
  }
  // DX20170926 - Create flags for SG functions - END

  vpflow.args2addattachedscheme(argv, cmds, "SG::FINDSYM", "--findsymSG=", "");
  vpflow.args2addattachedscheme(argv, cmds, "SG::FINDSYM_LABEL", "--findsymSG_label=", "");
  vpflow.args2addattachedscheme(argv, cmds, "SG::FINDSYM_NUMBER", "--findsymSG_number=", "");
  vpflow.args2addattachedscheme(argv, cmds, "SG::FINDSYM_PRINT", "--findsym_print=", "");
  vpflow.args2addattachedscheme(argv, cmds, "SG::FINDSYM_EXEC", "--findsym=", "");
  // DX20170926 - Create flags for SG functions - START
  if (vpflow.flag("SG::FINDSYM") || vpflow.flag("SG::FINDSYM_LABEL") || vpflow.flag("SG::FINDSYM_NUMBER") || vpflow.flag("SG::FINDSYM_PRINT") || vpflow.flag("SG::FINDSYM_EXEC")) {
    if (aurostd::args2attachedflag(argv, "--findsymSG=|--findsymSG_label=|--findsymSG_number=|--findsym_print=|--findsym=")) {
      vpflow.args2addattachedscheme(argv, cmds, "SG::TOLERANCE", "--findsymSG=|--findsymSG_label=|--findsymSG_number=|--findsym_print=|--findsym=", "1");
    }
  }
  // DX20170926 - Create flags for SG functions - END

  vpflow.args2addattachedscheme(argv, cmds, "SLAB", "--slab=", "");

  vpflow.flag("SOF", aurostd::args2flag(argv, cmds, "--sof"));
  vpflow.flag("SPECIES", aurostd::args2flag(argv, cmds, "--species"));
  vpflow.flag("STATDIEL", aurostd::args2flag(argv, cmds, "--statdiel")); // CAMILO

  vpflow.flag("GENERALIZED_STACKING_FAULT_ENERGY", aurostd::args2flag(argv, cmds, "--stacking_fault")); // CO20190321
  if (vpflow.flag("GENERALIZED_STACKING_FAULT_ENERGY")) {
    vpflow.args2addattachedscheme(argv, cmds, "GENERALIZED_STACKING_FAULT_ENERGY::SHEAR_DIRECTION", "--shear_direction=|--shear=", ""); // CO20190321
    vpflow.args2addattachedscheme(argv, cmds, "GENERALIZED_STACKING_FAULT_ENERGY::STEP_SIZE", "--step_size=", ""); // CO20190321
    vpflow.args2addattachedscheme(argv, cmds, "GENERALIZED_STACKING_FAULT_ENERGY::STEPS", "--steps=", ""); // CO20190321
    vpflow.args2addattachedscheme(argv, cmds, "GENERALIZED_STACKING_FAULT_ENERGY::FIXED_LAYERS", "--fixed_layers=", ""); // CO20190321
    vpflow.flag("GENERALIZED_STACKING_FAULT_ENERGY::SPIN_OFF", aurostd::args2flag(argv, cmds, "--spin_off")); // CO20190321
    vpflow.flag("GENERALIZED_STACKING_FAULT_ENERGY::PARTIAL_DISSOCIATION", aurostd::args2flag(argv, cmds, "--partial_dissociation")); // CO20190321
  }

  vpflow.flag("CLEAVAGE_ENERGY", aurostd::args2flag(argv, cmds, "--cleavage_energy")); // CO20190321
  if (vpflow.flag("CLEAVAGE_ENERGY")) {
    vpflow.args2addattachedscheme(argv, cmds, "CLEAVAGE_ENERGY::RELAXATION_LAYERS", "--relaxation_layers=", ""); // CO20190321
    vpflow.flag("CLEAVAGE_ENERGY::SPIN_OFF", aurostd::args2flag(argv, cmds, "--spin_off")); // CO20190321
  }
  vpflow.args2addattachedscheme(argv, cmds, "CREATE_SLAB::PLANE_INTEREST", "--plane_interest=|--plane=", ""); // CO20190321
  vpflow.args2addattachedscheme(argv, cmds, "CREATE_SLAB::TOTAL_LAYERS", "--total_layers=", ""); // CO20190321
  vpflow.args2addattachedscheme(argv, cmds, "CREATE_SLAB::VACUUM", "--vacuum=", ""); // CO20190321

  vpflow.flag("STDCONVCELL", aurostd::args2flag(argv, cmds, "--sc|--standard_conventional|--std_conv|--sconv"));
  vpflow.flag("STDPRIMCELL", aurostd::args2flag(argv, cmds, "--sp|--standard_primitive|--std_prim|--sprim"));
  // DX20190128 - add structure2ANRL - START
  vpflow.flag("STRUCTURE2ANRL", aurostd::args2attachedflag(argv, cmds, "--prototype") && !vpflow.flag("LIST_PROTOTYPE_LABELS")); // DX20190416 - changed from --anrl to --prototype //DX20190509 - added check for prototype labels
  if (vpflow.flag("STRUCTURE2ANRL")) {
    vpflow.args2addattachedscheme(argv, cmds, "STRUCTURE2ANRL::SETTING", "--setting=", "");
    vpflow.args2addattachedscheme(argv, cmds, "STRUCTURE2ANRL::TOLERANCE", "--tolerance=", ""); // DX20191028
    vpflow.args2addattachedscheme(argv, cmds, "STRUCTURE2ANRL::TOLERANCE_SPECTRUM", "--tolerance_spectrum=", ""); // DX20200820
    vpflow.flag("STRUCTURE2ANRL::FORCE_WYCKOFF", aurostd::args2flag(argv, cmds, "--force_Wyckoff|--force_wyckoff|--force_Wyckoff_order|--force_wyckoff_order")); // DX20191028
    vpflow.flag("STRUCTURE2ANRL::PRINT_ELEMENT_NAMES", aurostd::args2flag(argv, cmds, "--print_element_names|--print_element_name|--print_element")); // DX20210622
    vpflow.flag("STRUCTURE2ANRL::PRINT_ATOMIC_NUMBERS", aurostd::args2flag(argv, cmds, "--print_atomic_numbers|--print_atomic_number")); // DX20210622
  }
  // DX20190128 - add structure2ANRL - END

  vpflow.flag("SUMPDOS", (aurostd::args2flag(argv, cmds, "--sumpdos") && argv.at(1) == "--sumpdos"));

  vpflow.args2addattachedscheme(argv, cmds, "SUPERCELL", "--supercell=", "");
  vpflow.args2addattachedscheme(argv, cmds, "SUPERCELLSTRLIST", "--supercell_strlist=", "");

  vpflow.flag("SWAP", aurostd::args2flag(argv, cmds, "--swap"));

  vpflow.flag("TERDATA", aurostd::args2flag(argv, cmds, "--terdata") || aurostd::args2attachedflag(argv, cmds, "--terdata="));
  vpflow.flag("TERDATA_EXIST", aurostd::args2flag(argv, cmds, "--terdata_exist"));

  vpflow.flag("UFFENERGY", aurostd::args2flag(argv, cmds, "--uffenergy|--ue"));

  // DATABASE
  vpflow.flag("PATCHDB", aurostd::args2flag(argv, cmds, "--patch_database")); // ME20200829
  vpflow.flag("UPDATEDB", aurostd::args2flag(argv, cmds, "--update_database")); // ME20191001
  vpflow.flag("REBUILDDB", aurostd::args2flag(argv, cmds, "--rebuild_database")); // ME20191001
  vpflow.args2addattachedscheme(argv, cmds, "DBPATCHFILES", "--patchfiles=", "");

  // DX20180710 - we do not want to run if the flag was used in proto - vpflow.flag("VASP",aurostd::args2flag(argv,cmds,"--vasp"));
  vpflow.flag("VASP", aurostd::args2flag(argv, cmds, "--vasp|--vasp4") && !vpflow.flag("PROTO_AFLOW") && !vpflow.flag("PROTO")); // DX20180710 - check if used in proto //CO20210119 - added vasp4
  vpflow.flag("VASP5", aurostd::args2flag(argv, cmds, "--vasp5") && !vpflow.flag("PROTO_AFLOW") && !vpflow.flag("PROTO")); // DX20180710 - check if used in proto //CO20210119 - added vasp5

  // ME20200330
  vpflow.flag("VISUALIZE_PHONONS", aurostd::args2flag(argv, cmds, "--visualize_phonons"));
  if (vpflow.flag("VISUALIZE_PHONONS")) {
    vpflow.args2addattachedscheme(argv, cmds, "ADISP::AMPLITUDE", "--amplitude=", "");
    vpflow.args2addattachedscheme(argv, cmds, "ADISP::BRANCHES", "--branches=", "");
    vpflow.args2addattachedscheme(argv, cmds, "ADISP::FORMAT", "--format=", "");
    vpflow.args2addattachedscheme(argv, cmds, "ADISP::PERIODS", "--periods=", "");
    vpflow.args2addattachedscheme(argv, cmds, "ADISP::QPOINTS", "--qpoints=|--qpoint=|--q=", "");
    vpflow.args2addattachedscheme(argv, cmds, "ADISP::STEPS", "--steps=", "");
    vpflow.args2addattachedscheme(argv, cmds, "ADISP::SUPERCELL", "--scell=", "");
    vpflow.push_attached("ADISP::DIRECTORY", XHOST.vflag_control.getattachedscheme("DIRECTORY"));
  }

  vpflow.args2addattachedscheme(argv, cmds, "VOLUME::EQUAL", "--volume=", "");
  vpflow.args2addattachedscheme(argv, cmds, "VOLUME::MULTIPLY_EQUAL", "--volume*=", "");
  vpflow.args2addattachedscheme(argv, cmds, "VOLUME::PLUS_EQUAL", "--volume+=", "");

  //[CO20200404 - refer to XHOST]vpflow.flag("WWW",aurostd::args2flag(argv,cmds,"--web|--www|--http"));
  vpflow.flag("WYCKOFF", aurostd::args2flag(argv, cmds, "--wyckoff|--wy"));

  vpflow.args2addattachedscheme(argv, cmds, "XRAY", "--xray=", "");
  vpflow.args2addattachedscheme(argv, cmds, "XRAY_PEAKS", "--xray_peaks=", ""); // CO20190520
  vpflow.args2addattachedscheme(argv, cmds, "PLOT_XRAY", "--plot_xray=", ""); // CO20190520
  vpflow.flag("PLOT_XRAY::FORCE_GENERIC_TITLE", aurostd::args2flag(argv, cmds, "--force_generic_title|--force_title|--title")); // CO20190629
  vpflow.args2addattachedscheme(argv, cmds, "PLOT_XRAY_FILE", "--plot_xray_file=", ""); // CO20190520
  vpflow.args2addattachedscheme(argv, cmds, "XYZ", "--xyz=", "");
  vpflow.flag("XYZWS", aurostd::args2flag(argv, cmds, "--xyzwignerseitz|--xyzws"));

  //  vpflow.flag("XXX",aurostd::args2flag(argv,cmds,"--xxx"));
  vpflow.flag("XFIXX", (aurostd::args2flag(argv, cmds, "--xfixX|--xfixx") && argv.size() == 4));

  vpflow.flag("XTALFINDER_PYTHON", aurostd::args2attachedflag(argv, cmds, "--aflow_xtalfinder_python|--xtalfinder_python")); // DX20201228

  vpflow.flag("XXX", (aurostd::args2flag(argv, cmds, "--xxx")));
  if (vpflow.flag("XXX")) {
    cout << "XXX" << endl;
  }

  vpflow.args2addattachedscheme(argv, cmds, "LIB2RAW", "--lib2raw=", "");
  if (vpflow.flag("LIB2RAW")) {
    vpflow.flag("LIB2RAW_LOCAL", aurostd::args2flag(argv, cmds, "--local"));
  }
  vpflow.args2addattachedscheme(argv, cmds, "LIB2LIB", "--lib2lib=", ""); // CT20181212
  if (vpflow.flag("LIB2LIB")) {
    vpflow.flag("LIB2LIB_LOCAL", aurostd::args2flag(argv, cmds, "--local"));
  } // CT20181212
  vpflow.flag("FORCE", aurostd::args2flag(argv, cmds, "--force"));

  vpflow.flag("XPLUG::INIT", aurostd::args2flag(argv, cmds, "--xplug"));
  if (vpflow.flag("XPLUG::INIT")) {
    vpflow.flag("XPLUG::CLEAN", aurostd::args2flag(argv, cmds, "--clean|--C"));
    vpflow.args2addattachedscheme(argv, cmds, "XPLUG::DIRECTORY", "--directory=|--D=", "./");
    vpflow.args2addattachedscheme(argv, cmds, "XPLUG::LAST_MOD_DAYS", "--last_mod=|--lm=|--mod_last=|--ml=", "3.0");
    vpflow.args2addattachedscheme(argv, cmds, "XPLUG::ZIP_SIZE_GB", "--zip_size=|--zs=|--size_zip=|--sz=", "9");
    vpflow.args2addattachedscheme(argv, cmds, "XPLUG::PREFIX", "--prefix=", "");
    vpflow.args2addattachedscheme(argv, cmds, "XPLUG::RELATIVE", "--relative=", "");
  }

  //  cerr << "vpflow.flag(\"LIB2RAW\")=" << vpflow.flag("LIB2RAW") << endl;
  //  cerr << "vpflow.getattachedscheme(\"LIB2RAW\")=" << vpflow.getattachedscheme("LIB2RAW") << endl;

  vpflow.args2addattachedscheme(argv, cmds, "XRD_DIST", "--xrd_dist=|--XRD_DIST=", "");

  // XELEMENTS STUFF
  vpflow.args2addattachedscheme(argv, cmds, "XELEMENT", "--xelement=|--XELEMENT=|--element=|--ELEMENT=", "");

  vpflow.args2addattachedscheme(argv, cmds, "ZVAL", "--zval=|--ZVAL=", "");
  vpflow.args2addattachedscheme(argv, cmds, "ZVAL::CELL", "--zval_cell=|--ZVAL_CELL=|--zvalcell=|--ZVALCELL=", "");
  vpflow.args2addattachedscheme(argv, cmds, "ZVAL::ATOM", "--zval_atom=|--ZVAL_ATOM=|--zvalatom=|--ZVALATOM=", "");
  vpflow.args2addattachedscheme(argv, cmds, "POMASS", "--pomass=|--POMASS=", "");
  vpflow.args2addattachedscheme(argv, cmds, "POMASS::CELL", "--pomass_cell=|--POMASS_CELL=", "");
  vpflow.args2addattachedscheme(argv, cmds, "POMASS::ATOM", "--pomass_atom=|--POMASS_ATOM=|--pomassatom=|--POMASSATOM=", "");

  // Richard's symmetry functions (RHT)
  // DX START
  vpflow.flag("ORTHODEFECT_RHT", aurostd::args2flag(argv, cmds, "--OrthoDefect")); // RHT
  vpflow.flag("REVERSE_SPACEGROUP_RHT", aurostd::args2flag(argv, cmds, "--revsg")); // RHT
  vpflow.flag("PRIMITIVE_LATTICE_RHT", aurostd::args2flag(argv, cmds, "--primr | --fastprimitivecell | --fprim")); // RHT
  // DX20210610 - START
  vpflow.args2addattachedscheme(argv, cmds, "WYCKOFF_POSITIONS", "--Wyckoff=|--Wyckoff_positions=|--wyckoff=|--wyckoff_positions=|--wyccar=", "");
  if (vpflow.flag("WYCKOFF_POSITIONS")) {
    // DX20180807 - added more wyccar flags (--usage, --no_scan, setting, --magmom) - START
    vpflow.flag("WYCKOFF_POSITIONS::USAGE", aurostd::args2flag(argv, cmds, "--usage|--USAGE"));
    vpflow.flag("WYCKOFF_POSITIONS::NO_SCAN", aurostd::args2flag(argv, cmds, "--no_scan"));
    vpflow.args2addattachedscheme(argv, cmds, "WYCKOFF_POSITIONS::PRINT_WYCCAR",
                                  "--wyccar=", ""); // DX20210525 - treat wyccar as a special way of printing //DX20210708 - needs to be args2addattachedscheme to account for possible tolerance input
    vpflow.flag("WYCKOFF_POSITIONS::PRINT_LETTERS_ONLY", aurostd::args2flag(argv, cmds, "--letters_only"));
    vpflow.flag("WYCKOFF_POSITIONS::PRINT_SITE_SYMMETRIES_ONLY", aurostd::args2flag(argv, cmds, "--site_symmetries_only"));
    vpflow.flag("WYCKOFF_POSITIONS::PRINT_MULTIPLICITIES_ONLY", aurostd::args2flag(argv, cmds, "--multiplicities_only"));
    if (aurostd::args2attachedflag(argv, "--Wyckoff=|--Wyckoff_positions=|--wyckoff=|--wyckoff_positions=|--wyccar=")) {
      vpflow.args2addattachedscheme(argv, cmds, "WYCKOFF_POSITIONS::TOLERANCE", "--Wyckoff=|--Wyckoff_positions=|--wyckoff=|--wyckoff_positions=|--wyccar=", ""); // DX20200907 - default is system specific, leaving empty
    }
    if (aurostd::args2attachedflag(argv, "--setting=")) {
      vpflow.args2addattachedscheme(argv, cmds, "WYCKOFF_POSITIONS::SETTING", "--setting=", "1");
    }
    vpflow.args2addattachedscheme(argv, cmds, "WYCKOFF_POSITIONS::MAGNETIC", "--mag=|--magnetic=|--magmom=", "");
    // DX20180807 - added more wyccar flags (--usage, --no_scan, setting, --magmom) - END
  }
  // DX20210610 - END
  //  end Richard's symmetry (RHT)
  //  WE MIGHT NEED TO PUT THEM AROUND IN ALPHABETIC ORDER, keep the //RHT
  // DX END

  // *************************************
  // cluster expansion method
  vpflow.args2addattachedscheme(argv, cmds, "CE::CLUSTEREXPANSION", "--cluster-expansion=|--ce=", "");

  // special Quasirandom Structure (SQS)
  vpflow.args2addattachedscheme(argv, cmds, "CE::SQS", "--special-quasirandom-structure=|--sqs=", "");
  // get all clusters
  vpflow.args2addattachedscheme(argv, cmds, "CE::CLUSTERS", "--cluster=|--clusters=", "");
  // get all superlattices
  vpflow.args2addattachedscheme(argv, cmds, "CE::SUPERLATTICE", "--superlattice=", "");

  // *************************************
  // effective mass
  vpflow.flag("EFFECTIVEMASS", aurostd::args2flag(argv, cmds, "--effective-mass|--em")); // && (argv.size() == 3));

  // QHA

  if (LDEBUG) {
    cout << "PflowARGs: xscheme=" << vpflow.xscheme << endl;
  }
  if (LDEBUG) {
    cout << "PflowARGs: vxscheme.size()=" << vpflow.vxscheme.size() << endl;
  }
  if (LDEBUG) {
    cout << "PflowARGs: argv.size()=" << argv.size() << endl;
  }

  return vpflow.vxscheme.size();
}

namespace pflow {
  int main(vector<string>& argv, vector<string>& cmds) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    // cerr << "C username=" << XHOST.user << endl;
    // cerr << "C groupname=" << XHOST.group << endl;
    // cerr << "C XHOST.ostrPID=" << XHOST.ostrPID.str() << endl;
    // cerr << "C XHOST.ostrTID=" << XHOST.ostrTID.str() << endl; //CO20200502 - threadID
    // cerr.flush();
    // cerr << argv.size() << endl;
    aurostd::xoption vpflow;
    const ostringstream aus;
    // GENERAL **************************************************
    const string progname = argv.at(0);
    //  std::vector<string> cmds;
    // INFORMATION **************************************************

    // cerr << "************************************************************" << endl;
    // cerr << "* AFLOW IN mode PFLOW                                   *" << endl;
    // cerr << "************************************************************" << endl;

    _aflags aflags;
    aflags.Directory = "./";
    aurostd::args2flag(argv, cmds, "--np="); // put them in cmds
    aurostd::args2flag(argv, cmds, "--npmax"); // put them in cmds

    vpflow = XHOST.vflag_pflow;

    if (vpflow.flag("PFLOW_HELP") && argv.size() == 2) {
      cout << "**************************************************************************************************" << endl;
      cout << "  aflow --readme=pflow | --readme=processor | --readme=aconvasp | --readme_aconvasp" << endl;
      cout << "     Returns the HELP information for the \"processing machinery\"" << endl;
      // cout << AFLOW_AConvaspHelp();
      cout << aflow::Banner("BANNER_BIG");
      return 1;
    }
    if (vpflow.flag("PFLOW_HELP") && argv.size() == 3) {
      helpIndividualOption(argv);
    }
    if (vpflow.flag("PROTOS")) {
      cout << aflowlib::PrototypesHelp();
      cout << aflow::Banner("BANNER_BIG");
      return 1;
    }

    if (vpflow.flag("FIX_BANDS")) {
      pflow::FIXBANDS(aflags, vpflow.getattachedscheme("FIX_BANDS"));
      return 0;
    }

    const string EXTRACT_KPOINTS = aurostd::args2string(argv, cmds, "--extract_kpoints|--xkpoints", "nan");
    const string EXTRACT_INCAR = aurostd::args2string(argv, cmds, "--extract_incar|--xincar", "nan");
    const string EXTRACT_POSCAR = aurostd::args2string(argv, cmds, "--extract_poscar|--xposcar", "nan");
    const string EXTRACT_POTCAR = aurostd::args2string(argv, cmds, "--extract_potcar|--xpotcar", "nan");
    const string EXTRACT_PARTCAR = aurostd::args2string(argv, cmds, "--extract_partcar|--xpartcar", "nan"); // CO20181226

    // aflow style operations
    aflags.AFLOW_PERFORM_CLEAN = XHOST.vflag_aflow.flag("CLEAN");
    vpflow.flag("CLEAN", aflags.AFLOW_PERFORM_CLEAN);
    aflags.AFLOW_PERFORM_DIRECTORY = XHOST.vflag_control.flag("DIRECTORY");
    aflags.AFLOW_PERFORM_FILE = XHOST.vflag_control.flag("FILE");
    if (vpflow.flag("CLEAN")) {
      if (!aflags.AFLOW_PERFORM_DIRECTORY) {
        cerr << "AFLOW: to use --clean, you must specify the directory" << endl;
        return 0;
      } else {
        aurostd::xoption opts_clean; // CO20210716
        opts_clean.flag("SAVE_CONTCAR", aurostd::args2flag(argv, cmds, "--contcar_save|--save_contcar")); // CO20210716 - saves contcar no matter what
        opts_clean.flag("SAVE_CONTCAR_OUTCAR_COMPLETE", aurostd::args2flag(argv, cmds, "--contcar_save_outcar_complete|--save_contcar_outcar_complete")); // CO20210716 - saves contcar only if outcar is complete
        KBIN::Clean(aflags, opts_clean);
      }
    }

    // if(XXX) pflow::XXX(argv,cin);
    //  if(pflow::CheckCommands(argv,cmds)==false) return 0;

    xstructure a;
    //  a.iomode=IOAFLOW_AUTO;
    a.iomode = IOVASP_AUTO;
    const xvector<double> _emptyv(1);
    const xmatrix<double> _emptym(1, 1);

    // check whether a function is run or not
    // if not, give an error message
    bool _PROGRAMRUN = false;
    // *********************************************************************
    if (argv.size() == 1 && !_PROGRAMRUN) {
      cout << aflow::Banner("BANNER_TINY") << endl;
      //    cout << pflow::Intro_pflow("aconvasp");
      cout << pflow::Intro_pflow("aflow");
      cout << aflow::Banner("BANNER_TINY") << endl;
      _PROGRAMRUN = true;
    }
    // CO
    if (!argv.empty() && !_PROGRAMRUN) {
      if (vpflow.flag("BADER")) {
        cout << bader_functions::BaderCalc(vpflow);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("AFLOWMACHL::CoordCE_CSV")) {
        aflowMachL::writeCoordCECSV();
        _PROGRAMRUN = true;
      } // CO20200930
      if (vpflow.flag("CHGCAR2JVXL")) {
        cout << pflow::CHGCAR2JVXL(vpflow);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("CHGDIFF")) {
        cout << pflow::CHGDIFF(vpflow);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("CHGSUM")) {
        cout << pflow::CHGSUM(vpflow);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("CLEAVAGE_ENERGY")) {
        pflow::CleavageEnergyCalculation(vpflow, cin);
        _PROGRAMRUN = true;
      } // CO20190520
      // DX20190425 START
      if (vpflow.flag("COMPARE_DATABASE_ENTRIES")) {
        cout << compare::compareDatabaseEntries(vpflow);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("COMPARE_MATERIAL") || vpflow.flag("COMPARE_STRUCTURE")) {
        if (XHOST.vflag_control.flag("DIRECTORY") || XHOST.vflag_control.flag("FILE") || vpflow.flag("COMPARE_STRUCTURE::STRUCTURE_LIST")) {
          cout << compare::compareInputStructures(vpflow);
          _PROGRAMRUN = true;
        } else { // ME20210429 - Use std::cin otherwise
          cout << compare::compareInputStructures(vpflow, cin);
          _PROGRAMRUN = true;
        }
      }
      // DX20190425 END
      if (vpflow.flag("COMPARE_PERMUTATION")) {
        cout << compare::compareAtomDecorations(cin, vpflow);
        _PROGRAMRUN = true;
      } // DX20190201
      if (vpflow.flag("GFA::INIT")) {
        pflow::GLASS_FORMING_ABILITY(vpflow);
        _PROGRAMRUN = true;
      } // DF20190329 - GFA
      if (vpflow.flag("ATOMIC_ENVIRONMENT::INIT")) {
        pflow::ATOMIC_ENVIRONMENT(vpflow);
        _PROGRAMRUN = true;
      } // HE20210331 - Testing
      if (vpflow.flag("FIND_CLOSED_PACKING_PLANE")) {
        pflow::findClosedPackingPlane(cin);
        _PROGRAMRUN = true;
      } // CO20190808
      if (vpflow.flag("GENERATE_CERAMICS")) {
        cout << pflow::GENERATE_CERAMICS_PRINT(vpflow) << endl;
        _PROGRAMRUN = true;
      } // CO20200731
      // DX+CO START
      if (vpflow.flag("FULLSYMMETRY")) {
        pflow::CalculateFullSymmetry(cin, vpflow, cout);
        _PROGRAMRUN = true;
      }
      // DX+CO END
      if (vpflow.flag("STRUCTURE2ANRL")) {
        cout << anrl::structure2anrl(cin, vpflow);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("COMPARE2PROTOTYPES")) {
        cout << compare::printMatchingPrototypes(cin, vpflow);
        _PROGRAMRUN = true;
      } // DX20190314
      if (vpflow.flag("COMPARE2DATABASE")) {
        cout << compare::printCompare2Database(cin, vpflow);
        _PROGRAMRUN = true;
      } // DX20190201
      if (vpflow.flag("GENERALIZED_STACKING_FAULT_ENERGY")) {
        pflow::GeneralizedStackingFaultEnergyCalculation(vpflow, cin);
        _PROGRAMRUN = true;
      } // CO20190520
      if (vpflow.flag("PREPARE_CHGCAR_4_JMOL")) {
        cout << bader_functions::prepare_CHGCAR_4_Jmol(vpflow);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PFLOW::QUEUE_STATUS")) {
        cout << pflow::getQueueStatus(vpflow);
        _PROGRAMRUN = true;
      } // CO20200526
      // DX
      if (vpflow.flag("SOLIQUIDY::INIT")) {
        soliquidy::run_cmd(vpflow, cin);
        _PROGRAMRUN = true;
      } // HE20240528
    }
    // *********************************************************************
    if (argv.size() == 2 && !_PROGRAMRUN) {
      // put cartesian or fractional
      if (vpflow.flag("BZPLOT")) {
        LATTICE::BZPLOTDATA("", cin, 1);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("BZPLOTDATA")) {
        LATTICE::BZPLOTDATA("", cin, 0);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("ICSD_MAKELABEL")) {
        pflow::ICSD(argv, cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("KPATH")) {
        pflow::KPATH(cin, aurostd::args2attachedutype<double>(argv, "--grid=", -1), XHOST.vflag_control.flag("WWW"));
        _PROGRAMRUN = true;
      } // CO20200329 - default value -1 so we can decide grid automatically  //CO20200404 - new web flag
      if (vpflow.flag("NANOPARTICLE")) {
        cout << pflow::NANOPARTICLE(cin, xvector<double>(0));
        _PROGRAMRUN = true;
      }
      // if(vpflow.flag("PROTO_GUS_CPP")) {pflow::PROTO_GUS_CPP(argv); _PROGRAMRUN=true;}
      // QHA
    }
    // *********************************************************************
    if (argv.size() >= 2 && !_PROGRAMRUN) {
      // CO20220613 - must go first, as many manipulations are now possible at once
      // do not micromanage the user, if they want to run nonsensical options simultaneously, let them
      if (vpflow.flag("ABCCAR") || vpflow.flag("ABINIT") || vpflow.flag("AIMS") || vpflow.flag("ATAT") || vpflow.flag("CIF") || vpflow.flag("ELK") || vpflow.flag("LMP") || vpflow.flag("QE") ||
          vpflow.flag("VASP") || vpflow.flag("VASP5") || vpflow.flag("ITC") || vpflow.flag("INCELL") || vpflow.flag("INCOMPACT") || vpflow.flag("INWS") || vpflow.flag("MINKOWSKI_BASIS_REDUCTION") ||
          vpflow.flag("NIGGLI") || vpflow.flag("STDCONVCELL") || vpflow.flag("STDPRIMCELL") || vpflow.flag("CART") || vpflow.flag("FRAC")) {
        // use functions as much as possible to avoid repetition of the workflow
        xstructure xstr(cin);
        // change iomode - will also ReScale(1.0)
        if (vpflow.flag("ABCCAR")) {
          xstr.xstructure2abccar();
        }
        if (vpflow.flag("ABINIT")) {
          xstr.xstructure2abinit();
        }
        if (vpflow.flag("AIMS")) {
          xstr.xstructure2aims();
        }
        if (vpflow.flag("ATAT")) {
          xstr.xstructure2atat();
        }
        if (vpflow.flag("CIF")) {
          pflow::getCIFOptions(xstr, vpflow);
          xstr.xstructure2cif();
        } // DX20250319 - CIF options are required here (since the code was moved)
        if (vpflow.flag("ELK")) {
          xstr.xstructure2elk();
        }
        if (vpflow.flag("LMP")) {
          xstr.xstructure2lmp();
        }
        if (vpflow.flag("QE")) {
          xstr.xstructure2qe();
        }
        if (vpflow.flag("VASP") || vpflow.flag("VASP5")) {
          xstr.xstructure2vasp();
          if (vpflow.flag("VASP5")) {
            a.is_vasp4_poscar_format = false;
            a.is_vasp5_poscar_format = true;
          }
        }
        if (vpflow.flag("ITC")) {
          xstr.xstructure2itc();
        }
        // perform structure/lattice normalization
        if (vpflow.flag("INCELL")) {
          xstr.BringInCell();
        }
        if (vpflow.flag("INCOMPACT")) {
          xstr.BringInCompact();
        }
        if (vpflow.flag("INWS")) {
          xstr.BringInWignerSeitz();
        }
        if (vpflow.flag("MINKOWSKI_BASIS_REDUCTION")) {
          xstr.MinkowskiBasisReduction();
        }
        if (vpflow.flag("NIGGLI")) {
          xstr.NiggliUnitCellForm();
        }
        if (vpflow.flag("STDCONVCELL")) {
          xstr.Standard_Conventional_UnitCellForm();
        }
        if (vpflow.flag("STDPRIMCELL")) {
          xstr.Standard_Primitive_UnitCellForm();
        }
        // change coord_flag
        if (vpflow.flag("CART")) {
          xstr.SetCoordinates(_COORDS_CARTESIAN_);
        }
        if (vpflow.flag("FRAC")) {
          xstr.SetCoordinates(_COORDS_FRACTIONAL_);
        }
        cout << xstr;
        _PROGRAMRUN = true;
      }

      // A
      if (vpflow.flag("ACE")) {
        pflow::ACE(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("AFLOWIN")) {
        cout << pflow::AFLOWIN(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("AGROUP")) {
        pflow::SYMMETRY_GROUPS(aflags, cin, vpflow, cout);
        _PROGRAMRUN = true;
      } // DX20170818
      if (vpflow.flag("AGROUP2")) {
        pflow::AGROUP2(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("AGROUP2m")) {
        pflow::AGROUP2m(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("ALPHABETIC")) {
        cout << pflow::ALPHABETIC(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("ANGLES")) {
        pflow::ANGLES(vpflow.getattachedscheme("ANGLES"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("ALPHA_COMPOUND")) {
        cout << pflow::ALPHACompound(vpflow.getattachedscheme("ALPHA_COMPOUND"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("ALPHA_SPECIES")) {
        cout << pflow::ALPHASpecies(vpflow.getattachedscheme("ALPHA_SPECIES"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("AFLOWLIB::ENTRY_JSON")) {
        aflowlib::WEB_Aflowlib_Entry(vpflow.getattachedscheme("AFLOWLIB::ENTRY_JSON"), cout);
        _PROGRAMRUN = true;
      } // SC20190812
      if (vpflow.flag("AFLOWLIB_AUID2AURL")) {
        cout << aflowlib::AflowlibLocator(vpflow.getattachedscheme("AFLOWLIB_AUID2AURL"), "AFLOWLIB_AUID2AURL");
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("AFLOWLIB_AURL2AUID")) {
        cout << aflowlib::AflowlibLocator(vpflow.getattachedscheme("AFLOWLIB_AURL2AUID"), "AFLOWLIB_AURL2AUID");
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("AFLOWLIB_AUID2LOOP")) {
        cout << aflowlib::AflowlibLocator(vpflow.getattachedscheme("AFLOWLIB_AUID2LOOP"), "AFLOWLIB_AUID2LOOP");
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("AFLOWLIB_AURL2LOOP")) {
        cout << aflowlib::AflowlibLocator(vpflow.getattachedscheme("AFLOWLIB_AURL2LOOP"), "AFLOWLIB_AURL2LOOP");
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("AFLUX::SUMMONS")) {
        cout << aflowlib::AFLUXCall(vpflow) << endl;
        _PROGRAMRUN = true;
      } // DX20190206 - add AFLUX command line functionality //CO20200520 - AFLUX::SUMMONS
      if (vpflow.flag("AFLOWSYM_PYTHON")) {
        SYM::writePythonScript(cout);
        _PROGRAMRUN = true;
      } // DX20210202
      if (vpflow.flag("AFLUX")) {
        cout << aflowlib::AFLUXCall(vpflow) << endl;
        _PROGRAMRUN = true;
      } // DX20190206 - add AFLUX command line functionality
      //[CO20220614 - moved up]if(vpflow.flag("ATAT")) {cout << input2ATATxstr(cin); _PROGRAMRUN=true;} //SD20220123
      // B
      if (vpflow.flag("BANDGAP_WAHYU")) {
        AConvaspBandgap(argv);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("BANDGAP")) {
        pflow::BANDGAP(vpflow, cout);
        _PROGRAMRUN = true;
      } // CAMILO  //CO20171006
      if (vpflow.flag("BZPLOTDATAUSEKPOINTS")) {
        LATTICE::BZPLOTDATA(vpflow.getattachedscheme("BZPLOTDATAUSEKPOINTS"), cin, 10);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("BZPLOTUSEKPOINTS")) {
        LATTICE::BZPLOTDATA(vpflow.getattachedscheme("BZPLOTUSEKPOINTS"), cin, 11);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("BANDSTRUCTURE")) {
        pflow::BANDSTRUCTURE(aflags);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("BANDGAPS")) {
        AConvaspBandgaps(cin, cout);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("BANDGAPDOS")) {
        pflow::BANDGAP_DOS(vpflow, cout);
        _PROGRAMRUN = true;
      } // CO20191004
      if (vpflow.flag("BANDGAPDOS_WAHYU")) {
        AConvaspBandgapFromDOS(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("BANDGAPLISTDOS")) {
        AConvaspBandgapListFromDOS(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("BZDIRECTION")) {
        if (vpflow.flag("BZDIRECTION::LATTICE")) {
          cerr << "[" << vpflow.getattachedscheme("BZDIRECTION") << "]" << endl;
          cout << pflow::BZDirectionsLATTICE(vpflow.getattachedscheme("BZDIRECTION::LATTICE"));
        } else {
          cout << pflow::BZDirectionsSTRUCTURE(cin, vpflow); // DX20181101
        }
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("BZMAX")) {
        pflow::BZMAX(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("BANDS")) {
        pflow::BANDS(vpflow.getattachedscheme("BANDS"), cin);
        _PROGRAMRUN = true;
      }
      // C
      if (vpflow.flag("CAGES") && !AFLOW_PTHREADS::FLAG) {
        pflow::CAGES(aflags, vpflow.getattachedscheme("CAGES"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("CAGES") && AFLOW_PTHREADS::FLAG) {
        pflow::CAGES(aflags, vpflow.getattachedscheme("CAGES"), cin);
        _PROGRAMRUN = true;
      }
      //[CO20220614 - moved up]if(vpflow.flag("CART")) {cout << pflow::CART(cin); _PROGRAMRUN=true;}
      if (vpflow.flag("CCE_CORRECTION")) {
        cce::run(vpflow);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("CCE_CORRECTION::POSCAR2CCE")) {
        cce::run(vpflow, std::cin);
        _PROGRAMRUN = true;
      } // ME20200508  //CO20201105
      if (vpflow.flag("CCE_CORRECTION::GET_CCE_CORRECTION")) {
        cce::run(vpflow, std::cin);
        _PROGRAMRUN = true;
      } // RF20200916  //CO20201105
      if (vpflow.flag("CCE_CORRECTION::GET_OXIDATION_NUMBERS")) {
        cce::run(vpflow, std::cin);
        _PROGRAMRUN = true;
      } // RF20200725 //CO20201105
      if (vpflow.flag("CCE_CORRECTION::GET_CATION_COORDINATION_NUMBERS")) {
        cce::run(vpflow, std::cin);
        _PROGRAMRUN = true;
      } // RF20200814 //CO20201105
      if (vpflow.flag("CHECKINTEGRITIY")) {
        pflow::CheckIntegritiy();
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("CHANGESUFFIX")) {
        pflow::ChangeSuffix(vpflow.getattachedscheme("CHANGESUFFIX"));
        _PROGRAMRUN = true;
      } // KY20131222
      //[CO20220614 - moved up]if(vpflow.flag("CIF") && !vpflow.flag("PROTO_AFLOW") && !vpflow.flag("PROTO")) {pflow::CIF(cin,vpflow); _PROGRAMRUN=true;} //DX20180806 - added vpflow
      if (vpflow.flag("CLEANALL")) {
        pflow::CLEANALL(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("CORNERS")) {
        cout << pflow::CORNERS(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("CALCULATED_ICSD_RANDOM")) {
        cout << aflowlib::CALCULATED_ICSD_RANDOM();
        _PROGRAMRUN = true;
        return 0;
      }
      if (vpflow.flag("CALCULATED")) {
        cout << aflowlib::CALCULATED();
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("CLAT")) {
        pflow::CLAT(vpflow.getattachedscheme("CLAT"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("CHULL::INIT")) {
        chull::convexHull(vpflow);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("COMPARE")) {
        pflow::COMPARE(vpflow.getattachedscheme("COMPARE"));
        _PROGRAMRUN = true;
      }
      // D
      if (vpflow.flag("DATA")) {
        pflow::DATA(cin, vpflow, "DATA", cout);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("DATA_CRYSTAL_POINT_GROUP")) {
        pflow::DATA(cin, vpflow, "CRYSTAL_POINT_GROUP", cout);
        _PROGRAMRUN = true;
      } // DX20210209
      if (vpflow.flag("DATA_EDATA")) {
        pflow::DATA(cin, vpflow, "EDATA", cout);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("DATA_REAL_LATTICE")) {
        pflow::DATA(cin, vpflow, "REAL_LATTICE", cout);
        _PROGRAMRUN = true;
      } // DX20210209
      if (vpflow.flag("DATA_RECIPROCAL_LATTICE")) {
        pflow::DATA(cin, vpflow, "RECIPROCAL_LATTICE", cout);
        _PROGRAMRUN = true;
      } // DX20210209
      if (vpflow.flag("DATA_SGDATA")) {
        pflow::DATA(cin, vpflow, "SGDATA", cout);
        _PROGRAMRUN = true;
      } // DX20170818 //DX20210301 - consolidated into pflow::DATA()
      if (vpflow.flag("DATA_SUPERLATTICE")) {
        pflow::DATA(cin, vpflow, "SUPERLATTICE", cout);
        _PROGRAMRUN = true;
      } // DX20210209
      if (vpflow.flag("DEBYE")) {
        pflow::DEBYE(vpflow.getattachedscheme("DEBYE"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("DIFF")) {
        pocc::DIFF(vpflow.getattachedscheme("DIFF"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("DISP")) {
        pflow::DISP(vpflow.getattachedscheme("DISP"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("DIST")) {
        pflow::DIST(vpflow.getattachedscheme("DIST"), cin);
        _PROGRAMRUN = true;
      }
      // if(DYNADIEL) {pflow::DYNADIEL(argv) ; _PROGRAMRUN=true ;} // CAMILO
      //  E
      if (vpflow.flag("EDOS")) {
        pflow::EDOS(argv);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("EFFMASS")) {
        pflow::EFFMASS(argv, cout);
        _PROGRAMRUN = true;
      } // CAMILO
      //  if(vpflow.flag("EFFECTIVEMASS")) {pflow::EffectiveMass(argv,aurostd::args2string(argv,"--em","./"),cout); _PROGRAMRUN=true;}
      if (vpflow.flag("EIGCURV")) {
        pflow::EIGCURV(vpflow.getattachedscheme("EIGCURV"), cout);
        _PROGRAMRUN = true;
      } // CAMILO
      //[CO20220614 - moved up]if(vpflow.flag("ELK")) {cout << input2ELKxstr(cin); _PROGRAMRUN=true;} //DX20200313
      if (vpflow.flag("EQUIVALENT")) {
        cout << pflow::EQUIVALENT(aflags, cin, vpflow);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("EWALD")) {
        pflow::EWALD(vpflow.getattachedscheme("EWALD"), cin);
        _PROGRAMRUN = true;
      }
      // F
      if (vpflow.flag("FROZSL_VASPSETUP_AFLOW")) {
        cout << pflow::FROZSL_VASPSETUP(argv, 0);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("FROZSL_VASPSETUP_POSCAR")) {
        cout << pflow::FROZSL_VASPSETUP(argv, 1);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("FROZSL_ANALYZE")) {
        cout << pflow::FROZSL_ANALYZE(cin);
        _PROGRAMRUN = true;
      }

      if (vpflow.flag("FROZSL_README")) {
        cout << aurostd::EmbData::get_content("README_AFLOW_FROZSL.TXT", "README") << endl;
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("FROZSL_INPUT")) {
        cout << pflow::FROZSL_INPUT();
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("FROZSL_OUTPUT")) {
        cout << pflow::FROZSL_OUTPUT();
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("FGROUP")) {
        pflow::SYMMETRY_GROUPS(aflags, cin, vpflow, cout);
        _PROGRAMRUN = true;
      } // DX20170818
      //[CO20220614 - moved up]if(vpflow.flag("FRAC")) {cout << pflow::FRAC(cin); _PROGRAMRUN=true;}
      // G
      if (vpflow.flag("GETTEMP")) {
        AFLOW_getTEMP(argv);
        _PROGRAMRUN = true;
        // cout << "TEMPERATURE " << Message(__AFLOW_FILE__) << endl;return 0; _PROGRAMRUN=true;
      } // CO20200106 - patching for auto-indenting
      if (vpflow.flag("GEOMETRY")) {
        cout << pflow::GEOMETRY(cin) << endl;
        _PROGRAMRUN = true;
      } // CO20191110
      if (vpflow.flag("GULP")) {
        pflow::GULP(cin);
        _PROGRAMRUN = true;
      }
      // H
      if (vpflow.flag("HKL")) {
        pflow::HKL(vpflow.getattachedscheme("HKL"), aflags, cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("HKL_SEARCH_TRIVIAL")) {
        pflow::HKLSearch(vpflow.getattachedscheme("HKL_SEARCH_TRIVIAL"), aflags, cin, "HKL_SEARCH_TRIVIAL");
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("HKL_SEARCH_SIMPLE")) {
        pflow::HKLSearch(vpflow.getattachedscheme("HKL_SEARCH_SIMPLE"), aflags, cin, "HKL_SEARCH_SIMPLE");
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("HKL_SEARCH_COMPLETE")) {
        pflow::HKLSearch(vpflow.getattachedscheme("HKL_SEARCH_COMPLETE"), aflags, cin, "HKL_SEARCH_COMPLETE");
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("HNFCELL")) {
        pflow::POCC_COMMAND_LINE(vpflow, cin, cout);
        _PROGRAMRUN = true;
      } // CO20181226
      // I
      if (vpflow.flag("IAP::INIT")) {
        aflowMachL::WriteFileIAPCFG(vpflow);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("ICSD") || vpflow.flag("ICSD_CHEM") || vpflow.flag("ICSD_PROTO") || vpflow.flag("ICSD_ID") || vpflow.flag("ICSD_LESSTHAN") || vpflow.flag("ICSD_MORETHAN") ||
          vpflow.flag("ICSD_DENSLESSTHAN") || vpflow.flag("ICSD_DENSMORETHAN") || vpflow.flag("ICSD_SG") || vpflow.flag("ICSD_SGLESSTHAN") || vpflow.flag("ICSD_SGMORETHAN") || vpflow.flag("ICSD_TRICLINIC") ||
          vpflow.flag("ICSD_MONOCLINIC") || vpflow.flag("ICSD_ORTHORHOMBIC") || vpflow.flag("ICSD_TETRAGONAL") || vpflow.flag("ICSD_RHOMBOHEDRAL") || vpflow.flag("ICSD_TRIGONAL") || vpflow.flag("ICSD_HEXAGONAL") ||
          vpflow.flag("ICSD_CUBIC") || vpflow.flag("ICSD_UNIQUE") || vpflow.flag("ICSD_BASISLT") || vpflow.flag("ICSD_BASISGT") || vpflow.flag("ICSD_NOBROKENBASIS") || vpflow.flag("ICSD_NOPARTIALOCC") ||
          vpflow.flag("ICSD_N_ARY") || vpflow.flag("ICSD_REMOVE_AND") || vpflow.flag("ICSD_REMOVE_OR") || vpflow.flag("ICSD_REMOVEMETALS") || vpflow.flag("ICSD_ALLLESSTHAN") || vpflow.flag("ICSD_ALLMORETHAN") ||
          vpflow.flag("ICSD_TRI") || vpflow.flag("ICSD_MCL") || vpflow.flag("ICSD_MCLC") || vpflow.flag("ICSD_ORC") || vpflow.flag("ICSD_ORCC") || vpflow.flag("ICSD_ORCF") || vpflow.flag("ICSD_ORCI") ||
          vpflow.flag("ICSD_TET") || vpflow.flag("ICSD_BCT") || vpflow.flag("ICSD_RHL") || vpflow.flag("ICSD_HEX") || vpflow.flag("ICSD_CUB") || vpflow.flag("ICSD_FCC") || vpflow.flag("ICSD_BCC")) {
        pflow::ICSD(argv, cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("ICSD_LISTMETALS")) {
        pflow::ICSD_ListMetals();
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("ICSD2POSCAR")) {
        pflow::ICSD_2POSCAR(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("ICSD2PROTO")) {
        pflow::ICSD_2PROTO(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("ICSD2WYCK")) {
        pflow::ICSD_2WYCK(cin, vpflow.flag("SOF"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("IDENTICAL")) {
        cout << pflow::IDENTICAL(cin);
        _PROGRAMRUN = true;
      }
      //[CO20220614 - moved up]if(vpflow.flag("INCELL")) {cout << pflow::INCELL(cin); _PROGRAMRUN=true;}
      //[CO20220614 - moved up]if(vpflow.flag("INCOMPACT")) {cout << pflow::INCOMPACT(cin); _PROGRAMRUN=true;}
      if (vpflow.flag("INTPOL")) {
        pflow::INTPOL(vpflow.getattachedscheme("INTPOL"));
        _PROGRAMRUN = true;
      }
      //[CO20220614 - moved up]if(vpflow.flag("INWS")) {cout << pflow::INWS(cin); _PROGRAMRUN=true;}
      if (vpflow.flag("INFLATE_LATTICE")) {
        cout << pflow::INFLATE_LATTICE(vpflow.getattachedscheme("INFLATE_LATTICE"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("INFLATE_VOLUME")) {
        cout << pflow::INFLATE_VOLUME(vpflow.getattachedscheme("INFLATE_VOLUME"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("ISOPOINTAL_PROTOTYPES")) {
        cout << compare::isopointalPrototypes(cin, vpflow) << endl;
        _PROGRAMRUN = true;
      } // DX20200131
      // J
      if (vpflow.flag("JUSTAFTER")) {
        sflow::JUST(vpflow.getattachedscheme("JUSTAFTER"), cin, "JUSTAFTER");
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("JUSTBEFORE")) {
        sflow::JUST(vpflow.getattachedscheme("JUSTBEFORE"), cin, "JUSTBEFORE");
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("JUSTBETWEEN")) {
        sflow::JUST(vpflow.getattachedscheme("JUSTBETWEEN"), cin, "JUSTBETWEEN");
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("JMOL")) {
        pflow::JMOL(vpflow.getattachedscheme("JMOL"), cin);
        _PROGRAMRUN = true;
      }
      // K
      if (vpflow.flag("KBAND")) {
        pflow::KBAND(argv);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("KPOINTS")) {
        pflow::KPOINTS(vpflow.getattachedscheme("KPOINTS"), cin, cout);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("FLAG::XVASP_KPOINTS_DELTA")) {
        pflow::KPOINTS_DELTA(vpflow, cin, cout);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("KILL")) {
        sflow::KILL(vpflow.getattachedscheme("KILL"));
        _PROGRAMRUN = true;
      }
      // L
      if (vpflow.flag("LATTICEREDUCTION")) {
        cout << pflow::LATTICEREDUCTION(cin);
        _PROGRAMRUN = true;
      }
      // DX20200820 [OBSOELTE] if(vpflow.flag("LATTICE_TYPE")) {cout << pflow::LATTICE_TYPE(cin); _PROGRAMRUN=true;}
      // DX20200820 [OBSOELTE] if(vpflow.flag("LATTICE_LATTICE_TYPE")) {cout << pflow::LATTICE_LATTICE_TYPE(cin); _PROGRAMRUN=true;}
      if (vpflow.flag("LATTICE_TYPE")) {
        cout << pflow::LATTICE_TYPE(cin, vpflow);
        _PROGRAMRUN = true;
      } // DX20200820 - added vpflow
      if (vpflow.flag("LATTICE_LATTICE_TYPE")) {
        cout << pflow::LATTICE_LATTICE_TYPE(cin, vpflow);
        _PROGRAMRUN = true;
      } // DX20200820 - added vpflow
      if (vpflow.flag("LATTICE_HISTOGRAM")) {
        CheckLatticeHistogram();
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("LIST_PROTOTYPE_LABELS")) {
        cout << pflow::listPrototypeLabels(vpflow) << endl;
        _PROGRAMRUN = true;
      } // DX20190201
      if (vpflow.flag("LIB2RAW")) {
        XHOST.sensors_allowed = false;
        aflowlib::LIB2RAW(vpflow.getattachedscheme("LIB2RAW"), vpflow.flag("FORCE"), vpflow.flag("LIB2RAW_LOCAL"));
        XHOST.sensors_allowed = true;
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("LIB2LIB")) {
        XHOST.sensors_allowed = false;
        aflowlib::LIB2LIB(vpflow.getattachedscheme("LIB2LIB"), vpflow.flag("FORCE"), vpflow.flag("LIB2LIB_LOCAL"));
        XHOST.sensors_allowed = true;
        _PROGRAMRUN = true;
      } // CT20181212
      if (vpflow.flag("LTCELL")) {
        cout << pflow::LTCELL(vpflow.getattachedscheme("LTCELL"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("LTCELLFV")) {
        cout << pflow::LTCELL(vpflow.getattachedscheme("LTCELLFV"), cin);
        _PROGRAMRUN = true;
      }
      // M
      if (vpflow.flag("MULTI=SH")) {
        AFLOW_PTHREADS::MULTI_sh(argv);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("MAGNETICPARAMETERS")) {
        pflow::MagneticParameters(aurostd::args2attachedstring(argv, "--magpara=", "./"), cout);
        _PROGRAMRUN = true;
      }
      //[CO20220614 - moved up]if(vpflow.flag("MINKOWSKI_BASIS_REDUCTION")) {cout << pflow::MINKOWSKIBASISREDUCTION(cin); _PROGRAMRUN=true;}
      if (vpflow.flag("MSI")) {
        pflow::MSI(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("MOM")) {
        pflow::MOM(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("MULTIENUMALL")) {
        pocc::MultienumPrintAllXstr(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("MULTIENUMSORT")) {
        pocc::MultienumPrintSortedXstr(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("MAXATOMS")) {
        cout << pflow::ATOMSMAX(vpflow.getattachedscheme("MAXATOMS"), cin);
        _PROGRAMRUN = true;
      }
      // N
      if (vpflow.flag("NAMES") && argv.size() >= 2) {
        cout << pflow::NAMES(argv, cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("NUMNAMES") && argv.size() >= 2) {
        cout << pflow::NUMNAMES(argv, cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("NATOMS")) {
        cout << pflow::NATOMS(cin) << endl;
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("NBONDXX")) {
        cout << pflow::NBONDXX(cin);
        _PROGRAMRUN = true;
      } // CO20171025
      if (vpflow.flag("NDATA")) {
        pflow::NDATA(cin);
        _PROGRAMRUN = true;
      }
      //[CO20220614 - moved up]if(vpflow.flag("NIGGLI")) {cout << pflow::NIGGLI(cin); _PROGRAMRUN=true;}
      if (vpflow.flag("NNDIST")) {
        cout << pflow::NNDIST(cin) << endl;
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("NOORDERPARAMETER")) {
        cout << pflow::NOORDERPARAMETER(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("NSPECIES")) {
        cout << pflow::NSPECIES(cin) << endl;
        _PROGRAMRUN = true;
      }
      // O
      // if(vpflow.flag("OPARAMETER")) {pflow::OPARAMETER(argv,cin); _PROGRAMRUN=true;}
      // P

      // Serializers for DOS and bands
      if (vpflow.flag("DOSDATA2JSON")) {
        estructure::DOSDATA_JSON(vpflow, cout);
        _PROGRAMRUN = true;
      } // EG
      if (vpflow.flag("BANDSDATA2JSON")) {
        estructure::BANDSDATA_JSON(vpflow, cout);
        _PROGRAMRUN = true;
      } // EG
      // End serializers

      if (vpflow.flag("PLOT_BANDSPINSPLIT")) {
        estructure::PLOT_BAND_SPINSPLIT(vpflow.getattachedscheme("PLOT_BANDSPINSPLIT"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PLOT_BAND2")) {
        estructure::PLOT_BAND2(vpflow.getattachedscheme("PLOT_BAND2"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PLOT_DOSWEB")) {
        estructure::PLOT_DOSWEB(vpflow.getattachedscheme("PLOT_DOSWEB"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PLOT_PEDOS")) {
        estructure::PLOT_PEDOS(vpflow.getattachedscheme("PLOT_PEDOS"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PLOT_PEDOSALL")) {
        estructure::PLOT_PEDOSALL(vpflow.getattachedscheme("PLOT_PEDOSALL"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PLOT_PEDOSALL_AFLOWLIB")) {
        estructure::PLOT_PEDOSALL_AFLOWLIB(vpflow.getattachedscheme("PLOT_PEDOSALL_AFLOWLIB"), aflags);
        _PROGRAMRUN = true;
      }
      // ME20190614 BEGIN
      if (vpflow.flag("PLOT_BAND")) {
        aurostd::xoption pltopts = plotter::getPlotOptionsEStructure(vpflow, "PLOT_BAND");
        plotter::PLOT_BAND(pltopts);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PLOT_DOS")) {
        aurostd::xoption pltopts = plotter::getPlotOptionsEStructure(vpflow, "PLOT_DOS");
        plotter::PLOT_DOS(pltopts);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PLOT_BANDDOS")) {
        aurostd::xoption pltopts = plotter::getPlotOptionsEStructure(vpflow, "PLOT_BANDDOS");
        plotter::PLOT_BANDDOS(pltopts);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PLOT_PDOS")) {
        aurostd::xoption pltopts = plotter::getPlotOptionsEStructure(vpflow, "PLOT_PDOS", true);
        plotter::PLOT_PDOS(pltopts);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PLOT_PDOSALL")) {
        aurostd::xoption pltopts = plotter::getPlotOptionsEStructure(vpflow, "PLOT_PDOSALL", false);
        pltopts.push_attached("DATASET", "-1");
        plotter::PLOT_PDOS(pltopts);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PLOT_PHDOS")) {
        aurostd::xoption pltopts = plotter::getPlotOptionsPhonons(vpflow, "PLOT_PHDOS");
        plotter::PLOT_PHDOS(pltopts);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PLOT_PHDISP")) {
        aurostd::xoption pltopts = plotter::getPlotOptionsPhonons(vpflow, "PLOT_PHDISP");
        plotter::PLOT_PHDISP(pltopts);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PLOT_PHDISPDOS")) {
        aurostd::xoption pltopts = plotter::getPlotOptionsPhonons(vpflow, "PLOT_PHDISPDOS");
        plotter::PLOT_PHDISPDOS(pltopts);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PLOT_THERMO")) {
        aurostd::xoption plotopts = plotter::getPlotOptions(vpflow, "PLOT_THERMO");
        plotter::PLOT_THERMO(plotopts);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PLOT_TCOND")) {
        aurostd::xoption plotopts = plotter::getPlotOptions(vpflow, "PLOT_TCOND");
        plotter::PLOT_TCOND(plotopts);
        _PROGRAMRUN = true;
      }
      // ME20190614 END
      // AS20200909 BEGIN
      if (vpflow.flag("PLOT_THERMO_QHA")) {
        aurostd::xoption plotopts = plotter::getPlotOptionsQHAthermo(vpflow, "PLOT_THERMO_QHA");
        plotter::PLOT_THERMO_QHA(plotopts);
        _PROGRAMRUN = true;
      }
      // AS20200909 END
      // AS20210701 BEGIN
      if (vpflow.flag("PLOT_GRUENEISEN_DISPERSION")) {
        aurostd::xoption plotopts = plotter::getPlotOptionsPhonons(vpflow, "PLOT_GRUENEISEN_DISPERSION");
        plotter::PLOT_GRUENEISEN_DISPERSION(plotopts);
        _PROGRAMRUN = true;
      }
      // AS20210701 END
      if (vpflow.flag("PROTOS_ICSD")) {
        cout << aflowlib::PrototypesIcsdHelp(vpflow.getattachedscheme("PROTOS_ICSD"));
        cout << aflow::Banner("BANNER_BIG");
        return 1;
      }
      // if(POCCUPATION) {pflow::POCCUPATION(argv,cin); _PROGRAMRUN=true;}
      if (vpflow.flag("POCC_DOS")) {
        pocc::POCC_DOS(cout, vpflow.getattachedscheme("POCC_DOS"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("POCC_MAG")) {
        pocc::POCC_MAG(vpflow.getattachedscheme("POCC_MAG"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("POCC_BANDGAP")) {
        pocc::POCC_BANDGAP(vpflow.getattachedscheme("POCC_BANDGAP"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("POCC_MINIMUM_CONFIGURATION")) {
        cout << pocc::POCC_MINIMUM_CONFIGURATION(vpflow) << endl;
        _PROGRAMRUN = true;
      } // CO20191110
      //
      if (vpflow.flag("PROTO")) {
        if ((vpflow.flag("PROTO::SYMBOLS_ONLY")) || (vpflow.flag("PROTO::EQUATIONS_ONLY"))) {
          pflow::PROTO_LIBRARIES(vpflow);
          _PROGRAMRUN = true;
        } else {
          cout << pflow::PROTO_LIBRARIES(vpflow);
          _PROGRAMRUN = true;
        }
      }
      if (vpflow.flag("PROTO_AFLOW")) {
        pflow::PROTO_AFLOW(vpflow, false);
        _PROGRAMRUN = true;
      } // non reversed
      if (vpflow.flag("PLATON") && argv.size() >= 2) {
        cout << pflow::PLATON(vpflow.getattachedscheme("PLATON"), cin);
        _PROGRAMRUN = true;
      } // << endl;
      if (vpflow.flag("POMASS")) {
        pflow::ZVAL("POMASS," + vpflow.getattachedscheme("POMASS"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("POMASS::CELL")) {
        pflow::ZVAL("POMASS::CELL," + vpflow.getattachedscheme("POMASS::CELL"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("POMASS::ATOM")) {
        pflow::ZVAL("POMASS::ATOM," + vpflow.getattachedscheme("POMASS::ATOM"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PEARSON_SYMBOL")) {
        cout << pflow::PEARSON_SYMBOL(cin, vpflow);
        _PROGRAMRUN = true;
      } // DX20210611 - addded vpflow
      if (vpflow.flag("PDB")) {
        pflow::PDB(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PGROUP")) {
        pflow::SYMMETRY_GROUPS(aflags, cin, vpflow, cout);
        _PROGRAMRUN = true;
      } // DX20170818
      if (vpflow.flag("PGROUPX")) {
        pflow::SYMMETRY_GROUPS(aflags, cin, vpflow, cout);
        _PROGRAMRUN = true;
      } // DX20170818
      if (vpflow.flag("PGROUPK")) {
        pflow::SYMMETRY_GROUPS(aflags, cin, vpflow, cout);
        _PROGRAMRUN = true;
      } // DX20170818
      if (vpflow.flag("PGROUPK_PATTERSON")) {
        pflow::SYMMETRY_GROUPS(aflags, cin, vpflow, cout);
        _PROGRAMRUN = true;
      } // DX20200206
      if (vpflow.flag("PGROUPK_XTAL")) {
        pflow::SYMMETRY_GROUPS(aflags, cin, vpflow, cout);
        _PROGRAMRUN = true;
      } // DX20171205
      if (vpflow.flag("POSCAR")) {
        cout << pflow::POSCAR(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("RENDER")) {
        pflow::RENDER(cin, vpflow.getattachedscheme("RENDER_OUT"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("POSCAR2AFLOWIN")) {
        cout << pflow::POSCAR2AFLOWIN(cin, vpflow.getattachedscheme("AFLOW_PROTO::MODULE"));
        _PROGRAMRUN = true;
      } // Modified - ME20181113//ME20190112
      if (vpflow.flag("POSCAR2ENUM")) {
        pocc::POSCAR2ENUM(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("POSCAR2GULP")) {
        pocc::POSCAR2GULP(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("POCC_INPUT")) {
        pflow::POCC_INPUT();
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("POCC::CONVOLUTION")) {
        pocc::POCC_Convolution(vpflow);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("POSCAR2WYCKOFF")) {
        pflow::POSCAR2WYCKOFF(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PRIM")) {
        cout << pflow::PRIM(cin, 0);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PRIM1")) {
        cout << pflow::PRIM(cin, 1);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PRIM2")) {
        cout << pflow::PRIM(cin, 2);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PRIM3")) {
        cout << pflow::PRIM(cin, 3);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PSEUDOPOTENTIALS_CHECK")) {
        pflow::PSEUDOPOTENTIALS_CHECK(vpflow, vpflow.getattachedscheme("PSEUDOPOTENTIALS_CHECK"), cout);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PYTHON_MODULES")) {
        pflow::PYTHON_MODULES(vpflow.getattachedscheme("PYTHON_MODULES"));
        _PROGRAMRUN = true;
      } // ME20211103

      // Q
      //[CO20220614 - moved up]if(vpflow.flag("QE")) {cout << input2QExstr(cin); _PROGRAMRUN=true;}
      if (vpflow.flag("QCA::INIT")) {
        qca::quasiChemicalApprox(vpflow);
        _PROGRAMRUN = true;
      } // SD20220323
      if (vpflow.flag("QDEL")) {
        sflow::QDEL(vpflow.getattachedscheme("QDEL"));
        _PROGRAMRUN = true;
      } // NEW
      if (vpflow.flag("QMVASP")) {
        pflow::QMVASP(vpflow);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("QSUB")) {
        sflow::QSUB(vpflow.getattachedscheme("QSUB"));
        _PROGRAMRUN = true;
      } // NEW
      // R
      if (vpflow.flag("RDF")) {
        pflow::RDF(vpflow.getattachedscheme("RDF"), cin, vpflow.flag("RDF::RAW_COUNTS"));
        _PROGRAMRUN = true;
      } // CO20220627
      if (vpflow.flag("RDFCMP")) {
        pflow::RDFCMP(vpflow.getattachedscheme("RDFCMP"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("RMCOPIES")) {
        cout << pflow::RMCOPIES(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("RSM")) {
        pflow::RSM(argv, cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("RASMOL")) {
        pflow::RASMOL(vpflow.getattachedscheme("RASMOL"), cin);
        _PROGRAMRUN = true;
      }
      // ME20191001 START
      // ME20200829 - Added patch functionality
      if (vpflow.flag("REBUILDDB") || vpflow.flag("UPDATEDB") || vpflow.flag("PATCHDB")) {
        int return_code = 199; // Placeholder: return code in the 100s do not run the analysis
        aflowlib::AflowDB db(DEFAULT_AFLOW_DB_FILE, DEFAULT_AFLOW_DB_DATA_PATH, DEFAULT_AFLOW_DB_LOCK_FILE, XHOST.vschema, XHOST.vschema_internal);
        const string patchfiles = vpflow.getattachedscheme("DBPATCHFILES");
        // Hierarchy: rebuild > update > patch
        if (vpflow.flag("REBUILDDB") || vpflow.flag("UPDATEDB")) {
          return_code = db.rebuildDatabase(patchfiles, vpflow.flag("REBUILDDB"));
        } else if (vpflow.flag("PATCHDB")) {
          // false: do not check timestamps (always patch)
          return_code = db.patchDatabase(patchfiles, false);
        }
        if ((return_code < 100) || (return_code >= 200)) {
          db.analyzeDatabase(DEFAULT_AFLOW_DB_STATS_FILE);
        }
        _PROGRAMRUN = true;
        return return_code;
      }
      if (vpflow.flag("ANALYZEDB")) {
        aflowlib::AflowDB db(DEFAULT_AFLOW_DB_FILE, XHOST.vschema, XHOST.vschema_internal);
        db.analyzeDatabase(DEFAULT_AFLOW_DB_STATS_FILE);
        _PROGRAMRUN = true;
      }
      // ME20191001 END

      if (vpflow.flag("RMATOM")) {
        cout << pflow::RMATOM(cin, aurostd::args2utype(argv, "--rm_atom", 0));
        _PROGRAMRUN = true;
      }
      // S
      if (vpflow.flag("SHELL")) {
        pflow::SHELL(vpflow.getattachedscheme("SHELL"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("SHIFT")) {
        cout << pflow::SHIFT(vpflow.getattachedscheme("SHIFT"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("SLAB")) {
        cout << slab::MAKE_SLAB(vpflow.getattachedscheme("SLAB"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("STATDIEL")) {
        pflow::STATDIEL(argv);
        _PROGRAMRUN = true;
      } // CAMILO

      if (vpflow.flag("SG::AFLOW") && argv.size() >= 2) {
        cout << pflow::SG(vpflow, cin, "AFLOW", "ALL") << endl;
        _PROGRAMRUN = true;
      } // DX20170926
      if (vpflow.flag("SG::AFLOW_LABEL") && argv.size() >= 2) {
        cout << pflow::SG(vpflow, cin, "AFLOW", "LABEL") << endl;
        _PROGRAMRUN = true;
      } // DX20170926
      if (vpflow.flag("SG::AFLOW_NUMBER") && argv.size() >= 2) {
        cout << pflow::SG(vpflow, cin, "AFLOW", "NUMBER") << endl;
        _PROGRAMRUN = true;
      } // DX20170926
      if (vpflow.flag("SG::PLATON") && argv.size() >= 2) {
        cout << pflow::SG(vpflow, cin, "PLATON", "ALL") << endl;
        _PROGRAMRUN = true;
      } // DX20170926
      if (vpflow.flag("SG::PLATON_LABEL") && argv.size() >= 2) {
        cout << pflow::SG(vpflow, cin, "PLATON", "LABEL") << endl;
        _PROGRAMRUN = true;
      } // DX20170926
      if (vpflow.flag("SG::PLATON_NUMBER") && argv.size() >= 2) {
        cout << pflow::SG(vpflow, cin, "PLATON", "NUMBER") << endl;
        _PROGRAMRUN = true;
      } // DX20170926
      if (vpflow.flag("SG::FINDSYM") && argv.size() >= 2) {
        cout << pflow::SG(vpflow, cin, "FINDSYM", "ALL") << endl;
        _PROGRAMRUN = true;
      } // DX20170926
      if (vpflow.flag("SG::FINDSYM_LABEL") && argv.size() >= 2) {
        cout << pflow::SG(vpflow, cin, "FINDSYM", "LABEL") << endl;
        _PROGRAMRUN = true;
      } // DX20170926
      if (vpflow.flag("SG::FINDSYM_NUMBER") && argv.size() >= 2) {
        cout << pflow::SG(vpflow, cin, "FINDSYM", "NUMBER") << endl;
        _PROGRAMRUN = true;
      } // DX20170926

      if (vpflow.flag("SG::FINDSYM_PRINT")) {
        pflow::FINDSYM(vpflow, 0, cin);
        _PROGRAMRUN = true;
      } // DX20170926
      if (vpflow.flag("SG::FINDSYM_EXEC")) {
        pflow::FINDSYM(vpflow, 1, cin);
        _PROGRAMRUN = true;
      } // DX20170926

      if (vpflow.flag("SUPERCELL")) {
        cout << pflow::SUPERCELL(vpflow.getattachedscheme("SUPERCELL"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("SUPERCELLSTRLIST")) {
        pflow::SUPERCELLSTRLIST(vpflow.getattachedscheme("SUPERCELLSTRLIST"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("SD") && argv.size() >= 2) {
        cout << pflow::SD(argv, cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("SG")) {
        pflow::SG(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("SGROUP")) {
        pflow::SYMMETRY_GROUPS(aflags, cin, vpflow, cout);
        _PROGRAMRUN = true;
      } // DX20170818
      if (vpflow.flag("SPECIES")) {
        cout << pflow::SPECIES(cin);
        _PROGRAMRUN = true;
      }
      //[CO20220614 - moved up]if(vpflow.flag("STDCONVCELL")) {cout << GetStandardConventional(xstructure(cin,IOAFLOW_AUTO)); _PROGRAMRUN=true;}
      //[CO20220614 - moved up]if(vpflow.flag("STDPRIMCELL")) {cout << GetStandardPrimitive(xstructure(cin,IOAFLOW_AUTO)); _PROGRAMRUN=true;}
      if (vpflow.flag("SCALE")) {
        cout << pflow::SCALE(vpflow.getattachedscheme("SCALE"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("STRUCTURE2JSON")) {
        const xstructure xstr(cin, IOAFLOW_AUTO);
        cout << xstructure2json(xstr).toString() << endl;
        _PROGRAMRUN = true;
      } // DX20190508

      // T
      // U
      if (vpflow.flag("UFFENERGY")) {
        pocc::UFFENERGY(cin);
        _PROGRAMRUN = true;
      }
      // V
      //[CO20220614 - moved up]if(vpflow.flag("VASP")||vpflow.flag("VASP5")) {cout << input2VASPxstr(cin,vpflow.flag("VASP5")); _PROGRAMRUN=true;} //added bool for vasp5
      if (vpflow.flag("VISUALIZE_PHONONS")) {
        apl::createAtomicDisplacementSceneFile(vpflow);
        _PROGRAMRUN = true;
      } // ME20200330
      if (vpflow.flag("VOLUME::EQUAL")) {
        cout << pflow::VOLUME("VOLUME::EQUAL," + vpflow.getattachedscheme("VOLUME::EQUAL"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("VOLUME::MULTIPLY_EQUAL")) {
        cout << pflow::VOLUME("VOLUME::MULTIPLY_EQUAL," + vpflow.getattachedscheme("VOLUME::MULTIPLY_EQUAL"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("VOLUME::PLUS_EQUAL")) {
        cout << pflow::VOLUME("VOLUME::PLUS_EQUAL," + vpflow.getattachedscheme("VOLUME::PLUS_EQUAL"), cin);
        _PROGRAMRUN = true;
      }
      // X
      if (vpflow.flag("XPLUG::INIT")) {
        aflowlib::XPLUG(vpflow);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("XYZ")) {
        pflow::XYZ(vpflow.getattachedscheme("XYZ"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("XYZWS")) {
        pflow::XYZWS(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("XRAY")) {
        pflow::XRAY(vpflow.getattachedscheme("XRAY"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("XRAY_PEAKS")) {
        pflow::XRAY_PEAKS(vpflow, cin);
        _PROGRAMRUN = true;
      } // CO20190520
      if (vpflow.flag("PLOT_XRAY")) {
        pflow::PLOT_XRAY(vpflow, cin);
        _PROGRAMRUN = true;
      } // CO20190520
      if (vpflow.flag("PLOT_XRAY_FILE")) {
        pflow::PLOT_XRAY(vpflow);
        _PROGRAMRUN = true;
      } // CO20190520
      if (vpflow.flag("XRD_DIST")) {
        pflow::GetAtomicPlaneDist(vpflow.getattachedscheme("XRD_DIST"), cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("XELEMENT")) {
        pflow::XelementPrint(vpflow.getattachedscheme("XELEMENT"), cout);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("XTALFINDER_PYTHON")) {
        compare::writePythonScript(cout);
        _PROGRAMRUN = true;
      }
      // Y
      // Z
      if (vpflow.flag("ZVAL")) {
        pflow::ZVAL("ZVAL," + vpflow.getattachedscheme("ZVAL"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("ZVAL::CELL")) {
        pflow::ZVAL("ZVAL::CELL," + vpflow.getattachedscheme("ZVAL::CELL"));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("ZVAL::ATOM")) {
        pflow::ZVAL("ZVAL::ATOM," + vpflow.getattachedscheme("ZVAL::ATOM"));
        _PROGRAMRUN = true;
      }

      // Richard's Functions:
      if (vpflow.flag("REVERSE_SPACEGROUP_RHT")) {
        cout << SYM::ReverseSpaceGroup(argv) << endl;
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("ORTHODEFECT_RHT")) {
        SYM::OrthoDefect(cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PRIMITIVE_LATTICE_RHT")) {
        xstructure str(cin);
        str.GetPrimitiveCell();
        cout << str << endl;
        _PROGRAMRUN = true;
      }
      // DX START
      // WYCCAR FUNCTIONS (AUTO and MANUAL)

      if (vpflow.flag("WYCKOFF_POSITIONS")) {
        cout << pflow::WyckoffPositions(vpflow, cin);
        _PROGRAMRUN = true;
      } // DX20180807 - put into pflow function
      // End Richard's Functions
      // DX END
    }
    // *********************************************************************
    if (argv.size() == 3 && !_PROGRAMRUN) {
      if (vpflow.flag("CHGINT")) {
        pflow::CHGINT(argv);
        _PROGRAMRUN = true;
      }
      if (EXTRACT_KPOINTS != "nan") {
        cout << pflow::EXTRACT_xcar(aflags, argv, "KPOINTS", EXTRACT_KPOINTS);
        _PROGRAMRUN = true;
      }
      if (EXTRACT_INCAR != "nan") {
        cout << pflow::EXTRACT_xcar(aflags, argv, "INCAR", EXTRACT_INCAR);
        _PROGRAMRUN = true;
      }
      if (EXTRACT_POSCAR != "nan") {
        cout << pflow::EXTRACT_xcar(aflags, argv, "POSCAR", EXTRACT_POSCAR);
        _PROGRAMRUN = true;
      }
      if (EXTRACT_POTCAR != "nan") {
        cout << pflow::EXTRACT_xcar(aflags, argv, "POTCAR", EXTRACT_POTCAR);
        _PROGRAMRUN = true;
      }
      if (EXTRACT_PARTCAR != "nan") {
        cout << pflow::EXTRACT_xcar(aflags, argv, "PARTCAR", EXTRACT_PARTCAR);
        _PROGRAMRUN = true;
      } // CO20181226
      if (vpflow.flag("EXTRACT_SYMMETRY")) {
        cout << pflow::EXTRACT_Symmetry(aflags, argv);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("KPATH")) {
        pflow::KPATH(cin, aurostd::args2attachedutype<double>(argv, "--grid=", -1), XHOST.vflag_control.flag("WWW"));
        _PROGRAMRUN = true;
      } // CO20200329 - default value -1 so we can decide grid automatically //CO20200404 - new web flag
      if (vpflow.flag("INSPHERE")) {
        pflow::XYZINSPHERE(cin, aurostd::args2utype(argv, "--insphere", 0.0));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("NANOPARTICLE")) {
        cout << pflow::NANOPARTICLE(cin, aurostd::args2xvectorutype<double>(argv, "--nanoparticle", argv.size() - 2));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("POCC")) {
        pflow::POCC(argv);
        _PROGRAMRUN = true;
      }

      if (vpflow.flag("RAYTRACE")) {
        pflow::RAYTRACE(argv);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("SETORIGIN")) {
        cout << pflow::SETORIGIN(cin, (uint) aurostd::string2utype<int>(argv.at(2)));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("SEWALD")) {
        pflow::SEWALD(argv, cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("SGROUP")) {
        pflow::SGROUP(aflags, cin, aurostd::args2utype(argv, "--sgroup|--spacegroup", KBIN_SYMMETRY_SGROUP_RADIUS_DEFAULT));
        _PROGRAMRUN = true;
      }

      if (vpflow.flag("WYCKOFF")) {
        cout << pflow::WYCKOFF(argv, cin);
        _PROGRAMRUN = true;
      }
    }
    // *********************************************************************
    if (argv.size() >= 3 && !_PROGRAMRUN) {
      if (vpflow.flag("MISCIBILITY")) {
        cout << pflow::MISCIBILITY(argv);
        _PROGRAMRUN = true;
      }
    }
    // *********************************************************************
    if (argv.size() == 4 && !_PROGRAMRUN) {
      if (vpflow.flag("JOINSTRLIST")) {
        pflow::JOINSTRLIST(argv);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("MAKESTRLIST")) {
        pflow::MAKESTRLIST(argv);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("NANOPARTICLE")) {
        cout << pflow::NANOPARTICLE(cin, aurostd::args2xvectorutype<double>(argv, "--nanoparticle", argv.size() - 2));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PDOS")) {
        pflow::PDOS(argv);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PLANEDENS")) {
        pflow::PLANEDENS(argv);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("PRIM") && AFLOW_PTHREADS::FLAG) {
        cout << pflow::PRIM(cin, 0);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("RBANAL")) {
        pflow::RBANAL(argv);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("SUMPDOS")) {
        pflow::SUMPDOS(argv);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("SWAP")) {
        cout << pflow::xstrSWAP(argv, cin);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("XFIXX")) {
        aflowlib::XFIX_LIBRARY_ALL("LIB2", argv);
        _PROGRAMRUN = true;
      }
    }
    // *********************************************************************
    if (argv.size() == 5 && !_PROGRAMRUN) {
      if (vpflow.flag("RBDIST")) {
        pflow::RBDIST(argv);
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("SETCM")) {
        cout << pflow::SETCM(cin, aurostd::args2xvectorutype<double>(argv, "--setcm", 3));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("SETORIGIN")) {
        cout << pflow::SETORIGIN(cin, aurostd::args2xvectorutype<double>(argv, "--setorigin", 3));
        _PROGRAMRUN = true;
      }
      if (vpflow.flag("CMPSTR")) {
        pflow::CMPSTR(argv);
        _PROGRAMRUN = true;
      }
      // superlattice
    }
    // *********************************************************************
    // ----------------------------------------------- ERRORS

    if (!_PROGRAMRUN) {
      cerr << "aflow: the number of arguments is not correct" << endl;
      cerr << "Try \'aflow --help\' for more information" << endl;
      for (size_t i = 0; i < argv.size(); i++) {
        // cerr << "argv.at(" << i << ")=" << argv.at(i) << endl; // return 0;
        cerr << argv[i] << " ";
      }
      cerr << endl;
    }
    // *********************************************************************
    // return 1;
    return 0; // CO
  }
} // namespace pflow

// ***************************************************************************
// pflow::CheckCommands
// ***************************************************************************
namespace pflow {
  bool CheckCommands(vector<string> argv, const vector<string>& cmds) {
    string _cmd;
    vector<string> tokens;
    bool found = false;
    // check identities
    for (int i = argv.size() - 1; i >= 1 && !found; i--) {
      _cmd = aurostd::RemoveWhiteSpaces(string(argv.at(i)));
      for (size_t j = 0; j < cmds.size() && !found; j++) { // cerr << _cmd << " " << cmds[j] << endl;
        if (_cmd == cmds[j]) {
          found = true;
        }
      }
    }
    // check with =
    for (int i = argv.size() - 1; i >= 1 && !found; i--) {
      aurostd::RemoveWhiteSpaces(string(argv.at(i)));
      aurostd::string2tokens(aurostd::RemoveWhiteSpaces(string(argv.at(i))), tokens, "=");
      _cmd = tokens.at(0) + "=";
      for (size_t j = 0; j < cmds.size() && !found; j++) { // cerr << _cmd << " " << cmds[j] << endl;
        if (_cmd == cmds[j]) {
          found = true;
        }
      }
    }
    // not found
    if (!found) {
      cerr << aflow::Banner("BANNER_TINY") << endl;
      cerr << "ERROR - pflow::CheckCommands: command not found: " << _cmd << endl;
      return false;
    }
    return true;
  }
} // namespace pflow

// ***************************************************************************
// pflow::Intro_pflow
// ***************************************************************************
// patched by CO20200106 to avoid long string construction (patching for indents)
namespace pflow {
  string Intro_pflow(string x) {
    stringstream strstream;
    const string tab = "  ";
    string xspaces;
    for (size_t i = 0; i < x.size(); i++) {
      xspaces += " ";
    } // spaces size of x
    // intro(strstream);
    strstream << "******* BEGIN POSTPROCESSING MODE ******************************************************************" << endl;
    strstream << tab << x << " --help [-h] option_name" << endl;
    strstream << endl;
    strstream << tab << x << " --abccar < POSCAR | WYCCAR" << endl;
    strstream << tab << x << " --abinit < POSCAR" << endl;
    strstream << tab << x << " --aims < POSCAR" << endl;
    strstream << tab << x << " --ace < POSCAR" << endl;
    strstream << tab << x << " --use_aflow.in=XXX" << endl;
    strstream << tab << x << " --aflowin < POSCAR" << endl;
    strstream << tab << x << " --aflowSG[_label,_number][=<tolerance_value>|=tight|=loose] < POSCAR" << endl;
    strstream << tab << x
              << " --aflow-sym|--AFLOW-SYM|--AFLOWSYM|--aflowSYM|--aflowsym|--full_symmetry|--full_sym|--fullsym[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=txt|json] [--screen_only] "
                 "[--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR"
              << endl;
    strstream << tab << x << " --poscar2aflowin < POSCAR" << endl;
    strstream << tab << x << " --angle=cutoff < POSCAR" << endl;
    strstream << tab << x << " --agroup|--sitepointgroup[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=txt|json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR" << endl;
    strstream << tab << x << " --agroup2 < POSCAR" << endl;
    strstream << tab << x << " --agroup2m < POSCAR" << endl;
    strstream << tab << x << " --alphabetic < POSCAR" << endl;
    strstream << tab << x << " --alpha_compound=string1,string2...." << endl;
    strstream << tab << x << " --alpha_species=string1,string2..." << endl;
    strstream << tab << x << " --aflowlib=entry [--print=txt|json]" << endl;
    strstream << tab << x << " --aflowlib_auid2aurl=auid1,auid2....|--auid2aurl=..." << endl;
    strstream << tab << x << " --aflowlib_aurl2auid=aurl1,aurl2.... [ --aurl2auid=..." << endl;
    strstream << tab << x << " --aflowlib_auid2loop=auid1,auid2....|--auid2loop=..." << endl;
    strstream << tab << x << " --aflowlib_aurl2loop=aurl1,aurl2.... [ --aurl2loop=..." << endl;
    strstream << tab << x << " --atat < POSCAR" << endl;
    strstream << tab << x << " [options] --bader -D DIRECTORY" << endl;
    strstream << tab << xspaces << " " << "options are:  --usage" << endl;
    strstream << tab << xspaces << " " << "              --critical_points|--cp" << endl;
    strstream << tab << xspaces << " " << "              --calculate=|--calc=bader|voronoi" << endl;
    strstream << tab << xspaces << " " << "              --nocalculate=|--nocalc=bader|voronoi" << endl;
    strstream << tab << xspaces << " " << "              --partition=|--part=neargrid|ongrid" << endl;
    strstream << tab << xspaces << " " << "              --refine_edge_method=|--rem=-1|-2|-3" << endl;
    strstream << tab << xspaces << " " << "              --reference=|--ref=REF_FILE" << endl;
    strstream << tab << xspaces << " " << "              --vacuum=|--vac=off|auto|DENSITY_THRESHOLD" << endl;
    strstream << tab << xspaces << " " << "              --terminate=|--term=known|max" << endl;
    strstream << tab << xspaces << " " << "              --print_all=atom|bader|both" << endl;
    strstream << tab << xspaces << " " << "              --print_index=|--print_idx=atom|bader|both" << endl;
    strstream << tab << xspaces << " " << "              --print_select_atom=|--print_sel_atom=[LIST OR RANGE]" << endl;
    strstream << tab << xspaces << " " << "              --print_select_bader=|--print_sel_bader=[LIST OR RANGE]" << endl;
    strstream << tab << xspaces << " " << "              --print_sum_atom=[LIST OR RANGE]" << endl;
    strstream << tab << xspaces << " " << "              --print_sum_bader=[LIST OR RANGE]" << endl;
    strstream << tab << xspaces << " " << "              --quiet|--q" << endl;
    strstream << tab << xspaces << " " << "              --consolidate_atoms2species|--a2s" << endl;
    strstream << tab << xspaces << " " << "              --remove_bader_atoms|--rba" << endl;
    strstream << tab << xspaces << " " << "              --jvxl_all_species=|--jvxl=CUTOFF1,CUTOFF2...[::DOWNSAMPLE1,DOWNSAMPLE2,..]|CUTOFF1[,DOWNSAMPLE1:CUTOFF2,DOWNSAMPLE2:...]|CUTOFF[,DOWNSAMPLE]" << endl;
    strstream << tab << xspaces << " " << "              --keep=jvxl_only|--jvxl_only" << endl;
    strstream << tab << x << " --bands=PROOUT < POSCAR [AFLOW: NEED VERIFICATION]" << endl;
    strstream << tab << x << " --bandgap[=bands_directory1[,bands_directory2,...]]" << endl;
    strstream << tab << x << " --bandgaps < file_containing_bands_directories" << endl;
    strstream << tab << x << " --bandgapdos[=bands_directory1[,bands_directory2,...]]" << endl;
    strstream << tab << x << " --bandgaplistdos < DOSCAR_list" << endl;
    strstream << tab << x << " --bandstructure|--bs" << endl;
    strstream << tab << x << " --bzdirections|--bzd < POSCAR" << endl;
    strstream << tab << x << " --bzdirections=|--bzd=LATTICE" << endl;
    strstream << tab << x << " --BZmax < POSCAR" << endl;
    strstream << tab << x << " --bzplot|--plotbz < POSCAR" << endl;
    strstream << tab << x << " --bzplotuseKPOINTS=KPOINTS < POSCAR" << endl;
    strstream << tab << x << " --bzplotdata < POSCAR" << endl;
    strstream << tab << x << " --bzplotdatauseKPOINTS=KPOINTS < POSCAR" << endl;
    strstream << tab << x << " --cart [-c] < POSCAR" << endl;
    strstream << tab << x
              << " [--cce (prints user instructions and exits)] --cce=POSCAR_FILE_PATH [--oxidation_numbers=ox_num_1,ox_num_2,...] [--enthalpies_formation_dft=form_enthalpy_1,form_enthalpy_2,...] "
                 "[--functionals=functional_1,functional_2,...]"
              << endl;
    strstream << tab << x << " [options] --chgcar2jvxl=|--c2j=CHGCAR11[,CHGCAR2,...]::CUTOFF1,CUTOFF2...[::DOWNSAMPLE1,DOWNSAMPLE2,...]|CHGCAR1,CUTOFF1[,DOWNSAMPLE1:CHGCAR2,CUTOFF2[,DOWNSAMPLE2:...]]|CHGCAR,CUTOFF[,DOWNSAMPLE]"
              << endl;
    strstream << tab << xspaces << " " << "options are:  --usage" << endl;
    strstream << tab << xspaces << " " << "              --output=|--o=OUTPUT_FILE" << endl;
    strstream << tab << x << " [options]--chgdiff=CHGCAR1,CHGCAR2" << endl;
    strstream << tab << xspaces << " " << "options are:  --usage" << endl;
    strstream << tab << xspaces << " " << "              --output=|--o=CHGCAR_OUT" << endl;
    strstream << tab << x << " --chgsum=CHGCAR1,CHGCAR2,..." << endl;
    strstream << tab << xspaces << " " << "options are:  --usage" << endl;
    strstream << tab << xspaces << " " << "              --output=|--o=CHGCAR_OUT" << endl;
    strstream << tab << x << " --check_integrity|--checki" << endl;
    strstream << tab << x << " --cif < POSCAR" << endl;
    strstream << tab << x << " --clean -D DIRECTORY" << endl;
    strstream << tab << x << " --clean_all < LIST_DIRECTORIES" << endl;
    strstream << tab << x << " --compare2database < POSCAR" << endl; // DX20210611
    strstream << tab << x << " --compare2prototypes < POSCAR" << endl; // DX20210611
    strstream << tab << x << " --compare_database_entries [--alloy=AlgAlMn...] [--nspecies=<number>] [--stoichiometry=1,2,3,...] [--space_group=225,186,227,...]" << endl; // DX20210611
    strstream << tab << x << " --compare_materials=POSCAR1,POSCAR2,...| -D <path> | -F=<filename> [--np=xx (default 1)]" << endl; // DX20210611
    strstream << tab << x << " --compare_structures=POSCAR1,POSCAR2,...| -D <path> | -F=<filename> [--np=xx (default 1)]" << endl; // DX20210611
    strstream << tab << x << " --convex_hull=|--chull --alloy=MnPdPt[,AlCuZn,...] [--np=1] [chull_options] [--destination=[DIRECTORY]]" << endl;
    strstream << tab << xspaces << " " << "chull_options are:" << endl;
    strstream << endl;
    strstream << tab << xspaces << " " << "GENERAL OPTIONS:" << endl;
    strstream << tab << xspaces << " " << "              --usage" << endl;
    strstream << tab << xspaces << " " << "              --output=|--o=|--print=|--p=latex|pdf|png|json|text|jupyter|jupyter2|jupyter3" << endl;
    strstream << tab << xspaces << " " << "              --screen_only" << endl;
    strstream << tab << xspaces << " " << "              --keep=tex|--keep_tex|--keeptex|--tex" << endl;
    strstream << tab << xspaces << " " << "              --keep=log|--keep_log|--keeplog|--log" << endl;
    strstream << tab << xspaces << " " << "              --keep=tex,log" << endl;
    strstream << endl;
    strstream << tab << xspaces << " " << "LOADING OPTIONS:" << endl;
    strstream << tab << xspaces << " " << "              --load_library=|--loadlibrary=|--ll=icsd|lib1|lib2|lib3" << endl;
    strstream << tab << xspaces << " " << "              --load_API|--load_api|--loadapi|--lapi|--api" << endl;
    strstream << tab << xspaces << " " << "              --load_entries_entry_output|--loadentriesentryoutput|--leo" << endl;
    strstream << tab << xspaces << " " << "              --neglect=|--ban=aflow:bb0d45ab555bc208);aflow:fb9eaa58604ce774" << endl;
    strstream << tab << xspaces << " " << "              --see_neglect|--seeneglect|--sn" << endl;
    strstream << tab << xspaces << " " << "              --remove_extreme_points=|--removeextremepoints=|--remove_extrema=|--removeextrema=|--rep=-1000" << endl;
    strstream << tab << xspaces << " " << "              --entropic_temperature|--entropictemperature|--entroptemp" << endl;
    strstream << tab << xspaces << " " << "              --include_paw_gga|--paw_gga" << endl;
    strstream << endl;
    strstream << tab << xspaces << " " << "ANALYSIS OPTIONS:" << endl;
    strstream << tab << xspaces << " " << "              --distance_to_hull=|--distancetohull=|--distance2hull=|--dist2hull=|--d2h=aflow:bb0d45ab555bc208,aflow:fb9eaa58604ce774" << endl;
    strstream << tab << xspaces << " " << "              --stability_criterion=|--stabilitycriterion=|--stable_criterion=|--scriterion=|--sc=aflow:bb0d45ab555bc208,aflow:fb9eaa58604ce774" << endl;
    strstream << tab << xspaces << " "
              << "              "
                 "--n+1_enthalpy_gain=|--=|--n+1enthalpygain=|--n+1energygain=|--n+1egain=|--n1egain=|--n+1_enthalpygain=|--n+1+energygain=|--n+1_egain=|--nplus1=aflow:bb0d45ab555bc208,aflow:fb9eaa58604ce774"
              << endl;
    strstream << tab << xspaces << " " << "              --hull_formation_enthalpy=|--hull_enthalpy=|--hull_energy=0.25,0.25" << endl;
    strstream << tab << xspaces << " " << "              --skip_structure_comparison|--skipstructruecomparison|--skipstructcomp|--ssc" << endl;
    strstream << tab << xspaces << " " << "              --skip_stability_criterion_analysis|--skipstabilitycriterionanalysis|--skipscriterion|--sscriterion" << endl;
    strstream << tab << xspaces << " "
              << "              --skip_n_plus_1_enthalpy_gain_analysis|--skip_n_plus_1_energy_gain_analysis|--skipnplus1enthalpygainanalysis|--skipnplus1energygainanalysis|--skipnplus1|--snp1|--snpo" << endl;
    strstream << tab << xspaces << " " << "              --include_skewed_hulls|--include_skewed|--ish" << endl;
    strstream << tab << xspaces << " " << "              --include_unreliable_hulls|--include_unreliable|--iuh" << endl;
    strstream << tab << xspaces << " " << "              --include_outliers|--io" << endl;
    strstream << tab << xspaces << " " << "              --strict_outlier_analysis|--soa" << endl;
    strstream << tab << xspaces << " " << "              --include_ill_converged|--iic" << endl;
    strstream << tab << xspaces << " " << "              --force" << endl;
    strstream << endl;
    strstream << tab << xspaces << " " << "LATEX/PDF/PNG OPTIONS:" << endl;
    strstream << tab << xspaces << " " << "              --image_only|--imageonly|--image" << endl;
    strstream << tab << xspaces << " " << "              --no_document|--nodocument|--no_doc|--nodoc|--full_page_image|--fullpageimage" << endl;
    strstream << tab << xspaces << " " << "              --document_only|--documentonly|--doc_only|--doconly|--doc" << endl;
    strstream << tab << xspaces << " " << "              --latex_output|--latexoutput" << endl;
    strstream << tab << xspaces << " " << "              --latex_interactive|--latexinteractive" << endl;
    strstream << tab << xspaces << " " << "              --light_contrast|--lightcontrast|--lc" << endl;
    strstream << tab << xspaces << " " << "              --large_font|--largefont|--large|--lf" << endl;
    strstream << tab << xspaces << " " << "              --png_resolution=|--pngresolution=|--pngr=300" << endl;
    strstream << tab << xspaces << " " << "              --plot_iso_max_latent_heat|--plot_isomax|--iso_max|--isomax" << endl;
    strstream << tab << x << " --corners|--corner < POSCAR" << endl;
    strstream << tab << x << " --data[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=txt|json] < POSCAR" << endl;
    strstream << tab << x << " --data1=rcut < POSCAR" << endl;
    strstream << tab << x << " --data2 < POSCAR" << endl;
    strstream << tab << x << " --debye=THERMO[.bz2|.gz|.xz]" << endl;
    strstream << tab << x << " --diff=POSCAR1,POSCAR2" << endl;
    strstream << tab << x << " --disp=cutoff < POSCAR" << endl;
    strstream << tab << x << " --dist=cutoff < POSCAR" << endl;
    strstream << tab << x << " --delta_kpoints=number [or --dkpoints=number|-dkpoints=number|-dk=number] < POSCAR" << endl;
    strstream << tab << x << " --edata[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=txt|json] < POSCAR" << endl;
    strstream << tab << x << " --effective-mass|--em DIRECTORY" << endl;
    strstream << tab << x << " --eigcurv=DIRECTORY(with bands)" << endl;
    strstream << tab << x << " --element=Z|name|symbol[,property[,property]....]" << endl;
    strstream << tab << x << " --enum|--multienum < POSCAR" << endl;
    strstream << tab << x << " --enumsort|--multienumsort < POSCAR" << endl;
    strstream << tab << x << " --equivalent|--equiv|--inequivalent|--inequiv|--iatoms|--eatoms[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=txt|json] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR"
              << endl;
    strstream << tab << x << " --ewald[=eta] < POSCAR" << endl;
    strstream << tab << x << " --findsym[=tolerance_relative: default " << DEFAULT_FINDSYM_TOL << "] < POSCAR" << endl;
    strstream << tab << x << " --findsym_print[=tolerance_relative: default " << DEFAULT_FINDSYM_TOL << "] < POSCAR" << endl;
    strstream << tab << x << " --findsymSG[_label,_number][=tolerance_relative: default " << DEFAULT_FINDSYM_TOL << "] < POSCAR" << endl;
    strstream << tab << x << " --frac [-f,-d,--fract,--fractional,--direct] < POSCAR" << endl;
    strstream << tab << x << " --getTEMP [--runstat|--runbar|--refresh=X|--warning_beep=T|--warning_halt=T]" << endl;
    strstream << tab << x << " --gfa --alloy=CaCu [--ae_file=All_atomic_environments_read.dat] [--cutoff_energy=0.05]" << endl;
    strstream << tab << x << " --generate_ceramics --nm=C --m=Mo,Nb,Ta,V,W -N=5" << endl; // CO20200731
    strstream << tab << x << " --get_isopointal_prototype < POSCAR" << endl; // DX20210611
    strstream << tab << x << " --gulp < POSCAR" << endl;
    strstream << tab << x << " --identical < POSCAR" << endl;
    strstream << tab << x << " --incell < POSCAR" << endl;
    strstream << tab << x << " --incompact < POSCAR" << endl;
    strstream << tab << x << " --insphere radius < POSCAR" << endl;
    strstream << tab << x << " --inwignerseitz [--inws] < POSCAR" << endl;
    strstream << tab << x << " --inflate_lattice=coefficient|--ilattice coefficient=coefficient < POSCAR" << endl;
    strstream << tab << x << " --inflate_volume=coefficient|--ivolume=coefficient < POSCAR" << endl;
    strstream << tab << x << " --jgif < POSCAR" << endl;
    strstream << tab << x << " --jmol[=n1[,n2[,n3[,color[,true/false]]]]] < POSCAR" << endl;
    strstream << tab << x << " --pocc_hnf < POSCAR" << endl;
    strstream << tab << x << " --pocc_count_total < PARTCAR" << endl;
    strstream << tab << x << " --pocc_count_unique < PARTCAR" << endl;
    strstream << tab << x << " --pocc_show_unique < PARTCAR" << endl;
    strstream << tab << x << " --kpath [--grid=XX] < POSCAR" << endl;
    strstream << tab << x << " --kpoints=KDENS [or --kppra=KDENS|-k=KDENS] < POSCAR" << endl;
    strstream << tab << x << " --maxatoms=N|--max_atoms=N|--atoms_max=N|--atomsmax=N < POSCAR" << endl;
    strstream << tab << x << " --msi < POSCAR" << endl;
    strstream << tab << x << " --latticehistogram < POSCAR" << endl;
    strstream << tab << x << " --latticereduction < POSCAR" << endl;
    strstream << tab << x << " --lattice_type|--lattice|--lattice_crystal < POSCAR" << endl;
    strstream << tab << x << " --lattice_lattice_type|--lattice_lattice < POSCAR" << endl;
    strstream << tab << x << " --latticehistogram < POSCAR" << endl;
    strstream << tab << x << " --use_LOCK=XXX" << endl;
    strstream << tab << x << " --ltcell=a11,a12,a13,a21,a22,a23,a31,a32,a33 < POSCAR" << endl;
    strstream << tab << x << " --ltcell=a11,a22,a33 < POSCAR" << endl;
    strstream << tab << x << " --ltcell=file < POSCAR" << endl;
    strstream << tab << x << " --ltcellfv=v1,v2,v3,phi < POSCAR" << endl;
    strstream << tab << x << " --magpara=directory|--magpara" << endl;
    strstream << tab << x << " --minkowski_basis_reduction|--minkowski|--mink < POSCAR" << endl;
    strstream << tab << x << " --miscibility|--mix string" << endl;
    strstream << tab << x << " --mom < POSCAR" << endl;
    strstream << tab << x << " --natoms|--numatoms < POSCAR" << endl;
    strstream << tab << x << " --nbondxx < POSCAR" << endl;
    strstream << tab << x << " --names|--species A1 A2 ... < POSCAR" << endl;
    strstream << tab << x << " --nanoparticle radius distance < POSCAR" << endl;
    strstream << tab << x << " --ndata < POSCAR" << endl;
    strstream << tab << x << " --niggli < POSCAR" << endl;
    strstream << tab << x << " --nn < POSCAR" << endl;
    strstream << tab << x << " --noorderparameter < POSCAR" << endl;
    strstream << tab << x << " --nosd < POSCAR" << endl;
    strstream << tab << x << " --numnames A1 A2 ... < POSCAR" << endl;
    strstream << tab << x << " --nspecies|--numspecies< POSCAR" << endl;
    strstream << tab << x << " --OrthoDefect < POSCAR" << endl;
    strstream << tab << x << " --pdb < POSCAR" << endl;
    strstream << tab << x << " --rsm < POSCAR" << endl;
    strstream << tab << x << " --rsm --z atom_num1 atom_num2 atom_num3 < POSCAR" << endl;
    strstream << tab << x << " --rsm --z atom_symbol1 atom_symbol2 atom_symbol3 < POSCAR" << endl;
    strstream << tab << x << " --pearson_symbol|--pearson < POSCAR" << endl;
    strstream << tab << x << " --platon[=EQUAL|EXACT][,ang,d1,d2,d3] < POSCAR|platonSG" << endl;
    strstream << tab << x << " --platonSG[_label,_number][=EQUAL|EXACT][,ang,d1,d2,d3] < POSCAR" << endl;
    strstream << tab << x << " --plotband|--plotbands[=directory[,Emin[,Emax]]]] [--keep=gpl] [--print=pdf|gif|eps|jpg|png] [--title=] [--outfile=]" << endl;
    strstream << tab << x << " --plotband_spinsplit[=directory[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]]]" << endl;
    strstream << tab << x << " --plotbanddos|--plotbandsdos[=directory[,Emin[,Emax[,DOSSCALE]]]]] [--keep=gpl] [--noshift] [--print=pdf|eps|gif|jpg|png] [--projection=atoms|lm|none|orbitals] [--title=] [--outfile=]" << endl;
    strstream << tab << x << " --plotdos[=directory[,Emin[,Emax[,DOSSCALE]]]]] [--print=pdf|eps|gif|jpg|png] [--keep=gpl] [--noshift] [--projection=atoms|lm|none|orbitals] [--title=] [--outfile=]" << endl;
    strstream << tab << x << " --plotdosweb[=directory[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]]" << endl;
    strstream << tab << x << " --plotpdos|--plotpedos[=directory[,atom[,Emin[,Emax[,DOSSCALE]]]]] [--keep=gpl] [--noshift] [--print=pdf|eps|gif|jpg|png] [--projection=lm|none|orbitals] [--title=] [--outfile=]" << endl;
    strstream << tab << x << " --plotpdosall|--plotpedosall[=directory[,Emin[,Emax[,DOSSCALE]]]] [--keep=gpl] [--noshift] [--print=pdf|eps|gif|jpg|png] [--projection=lm|none|orbitals] [--title=] [--outfile=]" << endl;
    strstream << tab << x << " --plotpedosall_nonquivalent[=directory[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]]" << endl;
    strstream << tab << x << " --plotphdisp|--plotphonondispersion|--pphdis[=directory,[Emin,[Emax]]] [--keep=gpl] [--print=pdf|eps|gif|jpg|png] [--title=] [--unit=THz|Hz|eV|meV|rcm|cm-1] [--outfile=]" << endl;
    strstream << tab << x << " --plotphdos[=directory,[Emin,[Emax[,DOSSCALE]]]] [--keep=gpl] [--print=pdf|eps|gif|jpg|png] [--title=] [--unit=THz|Hz|eV|meV|rcm|cm-1] [--outfile=]" << endl;
    strstream << tab << x << " --plotphdispdos[=directory,[Emin,[Emax[,DOSSCALE]]]] [--keep=gpl] [--print=pdf|eps|gif|jpg|png] [--title=] [--unit=THz|Hz|eV|meV|rcm|cm-1] [--outfile=]" << endl;
    strstream << tab << x << " --plotthermo[=directory[,Tmin[,Tmax]]] [--keep=gpl] [--print=pdf|eps|gif|jpg|png] [--title=] [--outfile=]" << endl;
    strstream << tab << x << " --plotcond|--plothermalconductivity[=directory[,Tmin[,Tmax]]] [--keep=gpl] [--print=pdf|eps|gif|jpg|png] [--title=] [--outfile=]" << endl;
    strstream << tab << x << " --plotthermoqha[=directory[,Tmin[,Tmax]]] [--keep=gpl] [--print=pdf|eps|gif|jpg|png] [--title=] [--eosmodel=] [--outfile=]" << endl; // AS202210705
    strstream << tab << x << " --plotgrdisp|--plotgrueneisendispersion|--plotgruneisendispersion[=directory,[Emin,[Emax]]] [--keep=gpl] [--print=pdf|eps|gif|jpg|png] [--title=]  [--outfile=]" << endl; // AS20210705
    strstream << tab << x << " --pomass[=directory]" << endl;
    strstream << tab << x << " --pomass_atom[=directory]" << endl;
    strstream << tab << x << " --pomass_cell[=directory]" << endl;
    strstream << tab << x << " --pocc_input" << endl;
    strstream << tab << x << " --pocc_dos[=directory[,T[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]]]" << endl;
    strstream << tab << x << " --pocc_mag[=directory[,T[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]]]" << endl;
    strstream << tab << x << " --pocc_bandgap[=directory[,T[,DOS_Emin[,DOS_Emax[,DOSSCALE]]]]]" << endl;
    strstream << tab << x << " --pocc_minimum_configuration[=directory]" << endl;
    strstream << tab << x << " --poscar < ABCCAR | WYCCAR" << endl;
    strstream << tab << x << " --poscar2aflowin < POSCAR" << endl;
    strstream << tab << x << " --poscar2wyckoff < POSCAR" << endl;
    strstream << tab << x << " --poscar2gulp < POSCAR" << endl;
    strstream << tab << x << " [options] --prepare_chgcar_4_jmol=|--prep4jmol=CHGCAR1[,CHGACAR2,...]" << endl;
    strstream << tab << xspaces << " " << "options are:  --usage" << endl;
    strstream << tab << xspaces << " " << "              --outcar=OUTCAR" << endl;
    strstream << tab << xspaces << " " << "              --zip" << endl;
    strstream << tab << x << " --prim < POSCAR" << endl;
    strstream << tab << x << " --prim2 < POSCAR" << endl;
    strstream << tab << x << " --primr|--fastprimitivecell|--fprim < POSCAR" << endl;
    strstream << tab << x << " --prototype < POSCAR" << endl;
    strstream << tab << x
              << " --pseudopotentials_check=[POTCAR|OUTCAR]["
                 "|.bz2|.gz|.xz] | --pp_check= | --ppk="
              << endl;
    strstream << tab << x << " --python_modules[=module1,module2] | --create_python_modules=[module1,module2] [-D directory]" << endl;
    strstream << tab << x << " --qe < POSCAR" << endl;
    strstream << tab << x << " --qmvasp [--static] [-D directory]" << endl;
    strstream << tab << x << " --quasi_chem_approx|--qca --plattice=|--plat=fcc --elements=|--elem=Au,Pt[,Zn] [qca_options] [--directory=[DIRECTORY]]" << endl;
    strstream << tab << xspaces << " " << "options are:" << endl;
    strstream << endl;
    strstream << tab << xspaces << " " << "GENERAL OPTIONS:" << endl;
    strstream << tab << xspaces << " " << "              --usage" << endl;
    strstream << tab << xspaces << " " << "              --screen_only" << endl;
    strstream << tab << xspaces << " " << "              --image_only|--image" << endl;
    strstream << tab << xspaces << " " << "              --aflowlib_directory=|--aflowlib_dir=..." << endl;
    strstream << tab << xspaces << " " << "              --print=|--p=|--output=|--o=txt" << endl;
    strstream << endl;
    strstream << tab << xspaces << " " << "BINODAL OPTIONS:" << endl;
    strstream << tab << xspaces << " " << "              --binodal" << endl;
    strstream << tab << xspaces << " " << "              --use_sg" << endl;
    strstream << tab << xspaces << " " << "              --aflow_max_num_atoms=4" << endl;
    strstream << tab << xspaces << " " << "              --max_num_atoms=|--mna=8" << endl;
    strstream << tab << xspaces << " " << "              --cv_cutoff=|--cv_cut=0.05" << endl;
    strstream << tab << xspaces << " " << "              --conc_curve_range=|--conc_curve=0,1,1,0" << endl;
    strstream << tab << xspaces << " " << "              --conc_npts=20" << endl;
    strstream << tab << xspaces << " " << "              --temp_range=|--temp=300,5000" << endl;
    strstream << tab << xspaces << " " << "              --temp_npts=150" << endl;
    strstream << tab << x << " --rasmol[=n1[,n2[,n3]]] < POSCAR" << endl;
    strstream << tab << x << " --revsg [#] [n] [l] [m]" << endl;
    strstream << tab << x << " --rm_atom iatom < POSCAR" << endl;
    strstream << tab << x << " --rm_copies < POSCAR" << endl;
    strstream << tab << x << " --rdf[=rmax[,nbins[,sigma[,window_gaussian]]]] [--raw_counts] < POSCAR" << endl; // CO20220627
    strstream << tab << x << " --scale=s < POSCAR" << endl;
    strstream << tab << x << " --sd A1 A2 ... < POSCAR" << endl;
    strstream << tab << x << " --setcm cm1 cm2 cm3 < POSCAR" << endl;
    strstream << tab << x << " --setorigin r1 r2 r3 | atom# < POSCAR" << endl;
    strstream << tab << x << " --sewald eta < POSCAR" << endl;
    strstream << tab << x << " --sgdata|--space_group_data[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=txt|json] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR" << endl;
    strstream << tab << x << " --shell=ns,r1,r2,name,dens < POSCAR" << endl;
    strstream << tab << x << " --shift=Sx,Sy,Sz[,cCdD] < POSCAR" << endl;
    strstream << tab << x << " --sitepointgroup|--agroup < POSCAR" << endl;
    strstream << tab << x << " --slab=h,k,l[,#filled_layers[,#vacuum layers]] < POSCAR" << endl;
    strstream << tab << x << " --spacegroup radius < POSCAR" << endl;
    strstream << tab << x << " --species < POSCAR" << endl;
    strstream << tab << x << " --statdiel OUTCAR*" << endl;
    strstream << tab << x << " --std_conv|--standard_conventional|--sconv|--sc < POSCAR" << endl;
    strstream << tab << x << " --std_prim|--standard_primitive|--sprim|--sp < POSCAR" << endl;
    strstream << tab << x << " --suffix=[directory,]\"from2to\" where from/to=[n=none;r1=relax1;r2=relax2;r3=relax3;s=static;b=bands] or without abbreviations" << endl;
    strstream << tab << x << " --supercell=a11,a12,a13,a21,a22,a23,a31,a32,a33 < POSCAR" << endl;
    strstream << tab << x << " --supercell=a11,a22,a33 < POSCAR" << endl;
    strstream << tab << x << " --supercell=file < POSCAR" << endl;
    strstream << tab << x << " --supercell_strlist=a11,a12,a13,a21,a22,a23,a31,a32,a33,strlist" << endl;
    strstream << tab << x << " --supercell_strlist=a11,a22,a33,strlist" << endl;
    strstream << tab << x << " --supercell_strlist=file,strlist" << endl;
    strstream << tab << x << " --swap specie0 specie1 < POSCAR" << endl;
    strstream << tab << x << " --uffenergy|--ue  < POSCAR" << endl;
    strstream << tab << x << " --unique_atom_decorations < POSCAR" << endl; // DX20210611
    strstream << tab << x << " --vasp < GEOM" << endl;
    strstream << tab << x << " --itc < GEOM" << endl; // CO20220613
    strstream << tab << x << " --volume=x|--volume*=x|--volume+=x < POSCAR" << endl;
    strstream << tab << x << " --wyccar [TOL] < POSCAR" << endl;
    strstream << tab << x << " --wyccman|--WyccarManual|--wm [TOL] [ITERATION] < POSCAR" << endl;
    strstream << tab << x << " --xray=lambda < POSCAR" << endl;
    strstream << tab << x << " --xray_peaks=lambda < POSCAR" << endl;
    strstream << tab << x << " --xrd_dist=h,k,l < POSCAR" << endl;
    strstream << tab << x << " --xyz[=n1[,n2[,n3]]] < POSCAR" << endl;
    strstream << tab << x << " --xyzwignerseitz [--xyzws] < POSCAR" << endl;
    strstream << tab << x << " --zval[=directory]" << endl;
    strstream << tab << x << " --zval_atom[=directory]" << endl;
    strstream << tab << x << " --zval_cell[=directory]" << endl;
    strstream << endl;
    strstream << tab << x << " --pointgroup|--pgroup[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=txt|json] [--screen_only] < POSCAR" << endl;
    strstream << tab << x
              << " --pointgroup_crystal|--pgroup_crystal|--pgroup_xtal|--pgroupx|--pgroupX[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=txt|json] [--screen_only] "
                 "[--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR"
              << endl;
    strstream << tab << x << " --pgl < POSCAR" << endl;
    strstream << tab << x << " --pgx < POSCAR" << endl;
    strstream << tab << x << " --factorgroup|--fgroup[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=txt|=json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR" << endl;
    strstream << tab << x << " --sitepointgroup|--agroup < POSCAR" << endl;
    strstream << tab << x << " --spacegroup|--sgroup[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=txt|=json] [--screen_only] [--radius] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR" << endl;
    strstream << tab << x << " --pointgroupklattice|--pgroupk[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=txt|json] [--screen_only] < POSCAR" << endl;
    strstream << tab << x << " --pointgroupkcrystal|--pgroupk_xtal[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=txt|json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR" << endl;
    strstream << tab << x << " --pointgroupk_Patterson|--pgroupk_Patterson[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=txt|json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR"
              << endl;
    strstream << endl;
    strstream << " Prototypes HTQC" << endl;
    strstream << tab << x << " --prototypes|--protos" << endl;
    strstream << tab << x << " [options] --proto=label*[:speciesA*[:speciesB*]..[:volumeA*[:volumeB*]..|:volume]] [--params=... [--hex]]" << endl;
    strstream << endl;
    strstream << " Prototypes ICSD (from the local machine or aflowlib servers)" << endl;
    strstream << tab << x << " [--server=nnnnnnnnnn] --prototypes_icsd[=N]|--protos_icsd[=N]" << endl;
    strstream << tab << x << " [--server=nnnnnnnnnn] [--vasp|--itc|--qe|--abinit|--aims] --proto_icsd=label" << endl; // CO20220613
    strstream << endl;
    strstream << " Order Parameter [EXPERIMENTA]" << endl;
    strstream << tab << x << " --order_parameter [--order_parameter_sum XXX] [--Li8Mn16O32] < POSCAR [EXPERIMENTA]" << endl;
    strstream << endl;
    strstream << " Partial Occupation [EXPERIMENTA]" << endl;
    strstream << tab << x << " --partial_occupation|--partialoccupation|--pocc < POSCAR [EXPERIMENTA]" << endl;
    strstream << endl;
    strstream << " AFLOW Operations" << endl;
    strstream << tab << x << " [options] --aflow_proto=label*:speciesA*[:speciesB*][:volumeA*[:volumeB*]|:volume]" << endl;
    strstream << tab << x << " [options] --aflow_proto_icsd=label* potential_type (ICSD capable)" << endl;
    strstream << tab << xspaces << " " << "options are:  --potential=pot_LDA|pot_GGA|potpaw_LDA|potpaw_GGA|potpaw_PBE|potpaw_LDA_KIN|potpaw_PBE_KIN" << endl;
    strstream << tab << xspaces << " " << "              --potential_complete" << endl;
    strstream << tab << xspaces << " " << "              --module=[APL|QHA|AAPL]" << endl;
    strstream << tab << xspaces << " " << "              --apl_supercell=NxNxN" << endl;
    strstream << tab << xspaces << " " << "              --usage" << endl;
    strstream << tab << xspaces << " " << "              --missing" << endl;
    strstream << tab << xspaces << " " << "              --noautopp" << endl;
    strstream << tab << xspaces << " " << "              --bader (default: DEFAULT_VASP_FORCE_OPTION_BADER)" << endl;
    strstream << tab << xspaces << " " << "              --spin_remove_relax_1 (default: DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1)" << endl;
    strstream << tab << xspaces << " " << "              --spin_remove_relax_2 (default: DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2)" << endl;
    strstream << tab << xspaces << " " << "              --kscheme=[M|G] (default: DEFAULT_KSCHEME in .aflow.rc) --kscheme_static=[M|G] (default: DEFAULT_KSCHEME_STATIC in .aflow.rc)" << endl;
    strstream << tab << xspaces << " " << "              --kppra=NNNN (default: DEFAULT_KPPRA in .aflow.rc) --kppra_static=NNNN (default: DEFAULT_KPPRA_STATIC in .aflow.rc) --bands_grid=NNNN  (default: DEFAULT_BANDS_GRID in .aflow.rc )"
              << endl;
    strstream << tab << xspaces << " " << "              --enmax_multiply=NNNN (default: VASP_PREC_ENMAX_XXXX in .aflow.rc)" << endl;
    strstream << tab << xspaces << " " << "              --pressure=0,1,2 (kB) (default:0.0)" << endl;
    strstream << tab << xspaces << " " << "              --potim=XXX (default 0.05) (VASP)" << endl;
    strstream << tab << xspaces << " " << "              --relax_type=[ALL|IONS|CELL_SHAPE|CELL_VOLUME|IONS_CELL_VOLUME|IONS_CELL_SHAPE]" << endl;
    strstream << tab << xspaces << " " << "              --relax_mode=[ENERGY|FORCES|ENERGY_FORCES|FORCES_ENERGY] (default: DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME in .aflow.rc) (VASP)" << endl;
    strstream << tab << xspaces << " " << "              --relax_count=XX (default: DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT in .aflow.rc) (VASP)" << endl;
    strstream << tab << xspaces << " " << "              --run_relax_static" << endl;
    strstream << tab << xspaces << " " << "              --run_relax_static_bands" << endl;
    strstream << tab << xspaces << " " << "              --run_static" << endl;
    strstream << tab << xspaces << " " << "              --run_static_bands" << endl;
    strstream << tab << xspaces << " " << "              --precision=[(LOW|MEDIUM|NORMAL|HIGH|ACCURATE), PRESERVED] (default: DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME in .aflow.rc) (VASP)" << endl;
    strstream << tab << xspaces << " " << "              --algorithm=[(NORMAL|VERYFAST|FAST|ALL|DAMPED), PRESERVED] (default: DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME in .aflow.rc) (VASP)" << endl;
    strstream << tab << xspaces << " " << "              --metagga=[TPSS|RTPSS|M06L|MBJL|SCAN|MS0|MS1|MS2|NONE] (default: DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME in .aflow.rc) (VASP)" << endl;
    strstream << tab << xspaces << " " << "              --ivdw=[number_for_VASP_see_manual_for_IVDW|0] (default: DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME in .aflow.rc) (VASP)" << endl;
    strstream << tab << xspaces << " " << "              --type=[METAL|INSULATOR|SEMICONDUCTOR|DEFAULT] (default: DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME in .aflow.rc (VASP)" << endl;
    strstream << tab << xspaces << " " << "              --convert_unit_cell=[SPRIM, SCONV, NIGGLI, MINK, INCELL, COMPACT, WS, CART, FRAC, PRES]" << endl;
    strstream << tab << xspaces << " " << "              --volume_plus_equal=XXX" << endl;
    strstream << tab << xspaces << " " << "              --volume_multiply_equal=XXX" << endl;
    strstream << tab << xspaces << " " << "              --volume_preserved" << endl;
    strstream << tab << xspaces << " " << "              --ediffg=XXX  (default: DEFAULT_VASP_PREC_EDIFFG in .aflow.rc)" << endl;
    strstream << tab << xspaces << " " << "              --ldau2" << endl;
    strstream << tab << xspaces << " " << "              --noldau2" << endl;
    strstream << tab << xspaces << " " << "              --neglect_nomix" << endl;
    strstream << tab << xspaces << " " << "              --stdout" << endl;
    strstream << tab << xspaces << " " << "              --qe" << endl;
    strstream << tab << xspaces << " " << "              --abinit" << endl;
    strstream << tab << xspaces << " " << "              --aims" << endl;
    strstream << tab << xspaces << " " << "              --params=... { check aflow --readme=anrl }" << endl;
    strstream << tab << xspaces << " " << "              --hex        { check aflow --readme=anrl }" << endl;
    strstream << tab << xspaces << " " << "              --list" << endl;
    strstream << tab << x << " --extract_kpoints|--xkpoints " << _AFLOWIN_ << "" << endl;
    strstream << tab << x << " --extract_incar|--xincar " << _AFLOWIN_ << "" << endl;
    strstream << tab << x << " --extract_poscar|--xposcar " << _AFLOWIN_ << "" << endl;
    strstream << tab << x << " --extract_potcar|--xpotcar " << _AFLOWIN_ << "" << endl;
    strstream << endl;
    strstream << " CAGES SEARCH" << endl;
    strstream << tab << x << " --cages[=roughness] < POSCAR" << endl;
    strstream << endl;
    strstream << " HKL and HKL searches for surfaces" << endl;
    strstream << tab << x << " --miller=h,k,l[,nlayer[,elayer]] < POSCAR" << endl;
    strstream << tab << x << " --hkl=h,k,l[,bond] < POSCAR" << endl;
    strstream << tab << x << " --hkl_search[=khlmax[,bond[,step]]] < POSCAR" << endl;
    strstream << tab << x << " --hkl_search_simple[=cutoff[,bond[,khlmax[,step]]]] < POSCAR	" << endl;
    strstream << tab << x << " --hkl_search_complete[=cutoff[,bond[,khlmax[,step]]]] < POSCAR" << endl;
    strstream << endl;
    strstream << tab << x << " --chgint CHGCAR" << endl;
    strstream << tab << x << " --clat=a,b,c,alpha,beta,gamma" << endl;
    strstream << tab << x << " --cmp_str POSCAR1 POSCAR2 rcut" << endl;
    strstream << tab << x << " --compare=a,b,c,d,e,f,g,h,k,j,i,l" << endl;
    strstream << tab << x << " --intpol=file1,file2,nimages,nearest_image_flag" << endl;
    strstream << tab << x << " --join_strlist strlist1 strlist2" << endl;
    strstream << tab << x << " --make_strlist OUTCAR XDATCAR" << endl;
    strstream << tab << x << " --pdos pdos.in PROOUT [AFLOW: NEED VERIFICATION]" << endl;
    strstream << tab << x << " --planedens dens2d.in CHGCAR" << endl;
    strstream << tab << x << " --pocc PROOUT  [AFLOW: NEED VERIFICATION]" << endl;
    strstream << tab << x << " --raytrace rtfile" << endl;
    strstream << tab << x << " --rbanal nim nearest_image_flag" << endl;
    strstream << tab << x << " --rbdist POSCAR1 POSCAR2 n|N|e|E" << endl;
    strstream << tab << x << " --rdfcmp=rmax,nbins,sigma,nshmax,POSCAR1,POSCAR2" << endl;
    strstream << tab << x << " --spline npt < file" << endl;
    strstream << tab << x << " --sumpdos pdos.in PROOUT [AFLOW: NEED VERIFICATION]" << endl;
    strstream << endl;
    strstream << " ICSD MODE" << endl;
    strstream << tab << "ICSD input can be ternary.icsd, binary.icsd or other collections" << endl;
    strstream << tab << "of databases like sub-databases generated with and/or as follow." << endl;
    strstream << tab << x << " --icsd symbol/Z symbol/Z symbol/Z < input.icsd" << endl;
    strstream << tab << x << " --icsd_alllessthan symbol/Z" << endl;
    strstream << tab << x << " --icsd_allmorethan symbol/Z" << endl;
    strstream << tab << x << " --icsd_basislessthan #basis < input.icsd" << endl;
    strstream << tab << x << " --icsd_basismorethan #basis < input.icsd" << endl;
    strstream << tab << x << " --icsd_chem MgB4 < input.icsd" << endl;
    strstream << tab << x << " --icsd_cubic < input.icsd" << endl;
    strstream << tab << x << " --icsd_triclinic < input.icsd" << endl;
    strstream << tab << x << " --icsd_monoclinic < input.icsd" << endl;
    strstream << tab << x << " --icsd_orthorhombic < input.icsd" << endl;
    strstream << tab << x << " --icsd_tetragonal < input.icsd" << endl;
    strstream << tab << x << " --icsd_rhombohedral < input.icsd" << endl;
    strstream << tab << x << " --icsd_trigonal < input.icsd" << endl;
    strstream << tab << x << " --icsd_hexagonal < input.icsd" << endl;
    strstream << tab << x << " --icsd_cubic --icsd_orthorhombic < input.icsd" << endl;
    strstream << tab << x << " --icsd_tri < input.icsd" << endl;
    strstream << tab << x << " --icsd_mcl < input.icsd" << endl;
    strstream << tab << x << " --icsd_mclc < input.icsd" << endl;
    strstream << tab << x << " --icsd_orc < input.icsd" << endl;
    strstream << tab << x << " --icsd_orcc < input.icsd" << endl;
    strstream << tab << x << " --icsd_orcf < input.icsd" << endl;
    strstream << tab << x << " --icsd_orci < input.icsd" << endl;
    strstream << tab << x << " --icsd_tet < input.icsd" << endl;
    strstream << tab << x << " --icsd_bct < input.icsd" << endl;
    strstream << tab << x << " --icsd_rhl < input.icsd" << endl;
    strstream << tab << x << " --icsd_hex < input.icsd" << endl;
    strstream << tab << x << " --icsd_cub < input.icsd" << endl;
    strstream << tab << x << " --icsd_fcc < input.icsd" << endl;
    strstream << tab << x << " --icsd_bcc < input.icsd" << endl;
    strstream << tab << x << " --icsd_denslessthan X.X < input.icsd" << endl;
    strstream << tab << x << " --icsd_densmorethan Y.Y < input.icsd" << endl;
    strstream << tab << x << " --icsd_denslessthan X.X --icsd_densmorethan Y.Y < input.icsd" << endl;
    strstream << tab << x << " --icsd_id #ICSD_ID < input.icsd" << endl;
    strstream << tab << x << " --icsd_makelabel < input" << endl;
    strstream << tab << x << " --icsd_lessthan symbol/Z < input.icsd" << endl;
    strstream << tab << x << " --icsd_morethan symbol/Z < input.icsd" << endl;
    strstream << tab << x << " --icsd_lessthan symbol/Z --icsd_morethan symbol/Z < input.icsd" << endl;
    strstream << tab << x << " --icsd_listmetals" << endl;
    strstream << tab << x << " --icsd_nobrokenbasis < input.icsd" << endl;
    strstream << tab << x << " --icsd_nopartialocc < input.icsd" << endl;
    strstream << tab << x << " --icsd_n_ary #species < input.icsd" << endl;
    strstream << tab << x << " --icsd_proto #Nspecies1 #Nspecies2 #Nspecies3 ... < input.icsd" << endl;
    strstream << tab << x << " --icsd_remove_all symbol/Z symbol/Z" << endl;
    strstream << tab << x << " --icsd_remove_or symbol/Z symbol/Z" << endl;
    strstream << tab << x << " --icsd_removemetals" << endl;
    strstream << tab << x << " --icsd_sg #SG < input.icsd" << endl;
    strstream << tab << x << " --icsd_sglessthan #SG < input.icsd" << endl;
    strstream << tab << x << " --icsd_sgmorethan #SG < input.icsd" << endl;
    strstream << tab << x << " --icsd_sgmorethan #SG --icsd_sglessthan #SG < input.icsd" << endl;
    strstream << tab << x << " --icsd_unique < input.icsd" << endl;
    strstream << tab << x << " --icsd2aflowin < input.icsd" << endl;
    strstream << tab << x << " --icsd2poscar < input.icsd" << endl;
    strstream << tab << x << " --icsd2proto < input.icsd" << endl;
    strstream << tab << x << " --icsd2wyck < input.icsd" << endl;
    strstream << tab << x << " --icsd2wyck --sof < input.icsd" << endl;
    strstream << tab << x << " --icsdproto2aflowin < input.proto" << endl;
    strstream << endl;
    strstream << " FROZSL" << endl;
    strstream << tab << "scripting for using the FROZSL phonon framework" << endl;
    strstream << tab << x << " --frozsl_vaspsetup_aflow|--frozsl_vaspsetup < FROZSL.output" << endl;
    strstream << tab << x << " --frozsl_vaspsetup_aflow --file" << endl;
    strstream << tab << x << " --frozsl_vaspsetup_poscar < FROZSL.output" << endl;
    strstream << tab << x << " --frozsl_vaspsetup_poscar --file" << endl;
    strstream << tab << x << " --frozsl_analyze < aflow.frozsl.out" << endl;
    strstream << tab << x << " --frozsl_input" << endl;
    strstream << tab << x << " --frozsl_output" << endl;
    strstream << tab << x << " --readme=frozsl" << endl;
    strstream << endl;
    strstream << " ICSD/LIBN MODE (only for duke.edu computers)" << endl;
    strstream << tab << "scripting for the ICSD/LIBN libraries" << endl;
    strstream << tab << x << " [--force] --lib2raw=FCC/La1Se1_ICSD_27104" << endl;
    strstream << tab << x << " [--force] --lib2raw=AgCdZr/T0001.A2BC" << endl;
    strstream << tab << x << " [--force] --lib2raw=all[,directory]" << endl;
    strstream << tab << x << " [--force] --lib2scrub=[all,lib0,lib1,...,lib9,icsd]" << endl;
    strstream << endl;
    strstream << tab << x << " --xfixX system structure" << endl;
    strstream << endl;
    strstream << " CALCULATED DATABASE" << endl;
    strstream << tab << x << " --calculated" << endl;
    strstream << tab << x << " --calculated=icsd --random" << endl;
    strstream << tab << x << " --calculated[[=]all|icsd|lib1|lib2|lib3]" << endl;
    strstream << endl;
    strstream << " CLUSTER EXPANSION (very primitive, only fcc and bcc)" << endl;
    strstream << tab << x << " --cluster-expansion=...|--ce=structure_type,A,B,EA,EB" << endl;
    strstream << tab << x << " --superlattice=structure_type,n_min,n_max < POSCAR" << endl;
    strstream << tab << x << " --superlattice=VASP,structure_type,A,B < superlattice_name (CHECK) ???" << endl;
    strstream << tab << x << " --cluster=structure_type,n_min,n_max,m_min,m_max" << endl;
    strstream << tab << x << " --special-quasirandom-structure=..|--sqs=structure_type,atom_num,neighbor_num,sl_num_min,sl_num_max,A,B" << endl;
    strstream << tab << x << " --special-quasirandom-structure=..|--sqs=structure_type n1 n2 < POSCAR (CHECK)" << endl;
    strstream << endl;
    strstream << " TERNARY CONVEX HULL (only for duke.edu computers)" << endl;
    strstream << tab << x << " --terdata=A:B:C  [--fonts=XX|--keep=eps|--print=jpg|--print=gif|--print=png]" << endl;
    strstream << tab << x << " --terdata_exist list" << endl;
    strstream << endl;
    strstream << "******* END POSTPROCESSING MODE ********************************************************************" << endl;
    strstream << endl;
    // strstream << "Note: all the routines, except **, are tested to conform to convasp output" << endl;
    return strstream.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::ABCCAR
// ***************************************************************************
namespace pflow {
  xstructure ABCCAR(istream& input) {
    xstructure str(input, IOAFLOW_AUTO);
    str.iomode = IOVASP_ABCCAR;
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::ACE
// ***************************************************************************
namespace pflow {
  void ACE(istream& input) {
    //  xstructure a(input,IOAFLOW_AUTO);
    // pflow::PrintACE(a,cout);
    pflow::PrintACE(xstructure(input, IOAFLOW_AUTO), cout);
  }
} // namespace pflow

// DX20170927 - add spin info to xstructure - START
//  ***************************************************************************
//  pflow::AddSpinToXstructure (collinear version)
//  ***************************************************************************
namespace pflow {
  bool AddSpinToXstructure(xstructure& a, vector<double>& vmag) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (vmag.size() != a.atoms.size()) {
      cerr << XPID << "pflow::AddSpinToXstructure (collinear): ERROR: Number of magnetic moments (" << vmag.size() << ") does not match the number of atoms (" << a.atoms.size() << ")." << endl;
      return false;
    }
    // Only collinear for now 20170927
    for (size_t i = 0; i < a.atoms.size(); i++) {
      a.atoms[i].spin = vmag[i];
      a.atoms[i].spin_is_given = true;
      if (LDEBUG) {
        cerr << XPID << "pflow::AddSpinToXstructure (collinear): atom " << i << " magnetic moment: " << a.atoms[i].spin << endl;
      }
    }
    return true;
  }
} // namespace pflow
// DX20170927 - add spin info to xstructure - START

// DX20171205 - add spin info to xstructure - START
//  ***************************************************************************
//  pflow::AddSpinToXstructure (non-collinear version)
//  ***************************************************************************
namespace pflow {
  bool AddSpinToXstructure(xstructure& a, vector<xvector<double>>& vmag_noncoll) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (vmag_noncoll.size() != a.atoms.size()) {
      cerr << XPID << "pflow::AddSpinToXstructure (non-collinear): ERROR: Number of magnetic moments (" << vmag_noncoll.size() << ") does not match the number of atoms (" << a.atoms.size() << ")." << endl;
      return false;
    }
    // Only collinear for now 20170927
    for (size_t i = 0; i < a.atoms.size(); i++) {
      a.atoms[i].noncoll_spin = vmag_noncoll[i];
      a.atoms[i].noncoll_spin_is_given = true;
      if (LDEBUG) {
        cerr << XPID << "pflow::AddSpinToXstructure (non-collinear): atom " << i << " magnetic moment: " << a.atoms[i].noncoll_spin << endl;
      }
    }
    return true;
  }
} // namespace pflow
// DX20171205 - add spin info to xstructure - START

// DX20190801 - consolidated to a function - START
//  ***************************************************************************
//  pflow::ProcessAndAddSpinToXstructure
//  ***************************************************************************
namespace pflow {
  void ProcessAndAddSpinToXstructure(xstructure& a, const string& magmom_info) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    stringstream message;

    // ---------------------------------------------------------------------------
    // check if spin information is empty
    if (magmom_info.empty()) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " Spin information is empty; no spin information added to structure." << endl;
      }
      return;
    }

    const uint num_atoms = a.atoms.size(); // DX20191107 - int to uint
    bool is_noncoll = false;
    vector<xvector<double>> vmag_noncoll;
    bool is_coll = false;
    vector<double> vmag;

    // ---------------------------------------------------------------------------
    // check for non-collinear spin
    if (GetNonCollinearMagneticInfo(num_atoms, magmom_info, vmag_noncoll)) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " Non-collinear spin system detected." << endl;
      }
      is_noncoll = true;
      if (!AddSpinToXstructure(a, vmag_noncoll)) {
        message << "(non-collinear): Number of magnetic moments (" << vmag_noncoll.size() << ") does not match the number of atoms (" << a.atoms.size() << ")." << endl;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // DX20171205 - added non-collinear
      }
    }
    // ---------------------------------------------------------------------------
    // check for non-collinear spin
    if (!is_noncoll) {
      if (GetCollinearMagneticInfo(num_atoms, magmom_info, vmag)) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " Collinear spin system detected." << endl;
        }
        is_coll = true;
        if (!AddSpinToXstructure(a, vmag)) {
          message << "(collinear): Number of magnetic moments (" << vmag.size() << ") does not match the number of atoms (" << a.atoms.size() << ")." << endl;
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // DX20171205 - added non-collinear
        }
      }
    }
    // ---------------------------------------------------------------------------
    // could not detect; input error
    if (!is_noncoll && !is_coll) {
      message << "Could not detect collinear or non-collinear spin(s). Check spin input.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // DX20171205 - added non-collinear
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::AFLOWIN
// ***************************************************************************
namespace pflow {
  string AFLOWIN(istream& input) {
    stringstream oss;
    const xstructure str(input, IOAFLOW_AUTO);
    oss << AFLOWIN_SEPARATION_LINE << endl;
    oss << _VASP_POSCAR_MODE_EXPLICIT_START_ << endl;
    oss << str << "";
    oss << _VASP_POSCAR_MODE_EXPLICIT_STOP_ << endl;
    oss << AFLOWIN_SEPARATION_LINE << endl;
    return oss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::POSCAR2AFLOWIN
// ***************************************************************************
namespace pflow {
  string POSCAR2AFLOWIN(istream& input, const string& module) {
    const stringstream oss;
    xstructure str(input, IOAFLOW_AUTO);
    _xvasp xvasp;
    AVASP_DefaultValuesBinary_AFLOWIN(xvasp);
    xvasp.AVASP_prototype_mode = LIBRARY_MODE_PROTOTYPE;
    xvasp.AVASP_flag_PRECISION_scheme = "H";
    str.iomode = IOVASP_POSCAR; // ME20220113 - set to POSCAR or else the input format will be used
    xvasp.str = str;
    if (xvasp.str.atoms.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, "pflow::POSCAR2AFLOWIN():", "POSCAR has no atoms", _INPUT_ILLEGAL_);
    } // CO20200102
    if (xvasp.str.atoms[0].name_is_given == false) {
      throw aurostd::xerror(__AFLOW_FILE__, "pflow::POSCAR2AFLOWIN():", "POSCAR is missing species information", _INPUT_ILLEGAL_);
    } // CO20200102
    if (!module.empty()) {
      xvasp.aopts.push_attached("AFLOWIN_FLAG::MODULE", module); // ME20181113
    }
    KBIN::setModules(xvasp); // ME20181110
    stringstream aflowin;
    AVASP_MakeSingleAFLOWIN(xvasp, aflowin, false, -1);
    return oss.str();
  }
} // namespace pflow

namespace pflow {
  bool SYMMETRY_GROUPS(_aflags& aflags, istream& input, aurostd::xoption& vpflow, ostream& oss) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    _kflags kflags; // DX20170815 - Add in consistency checks
    const xstructure _a(input, IOAFLOW_AUTO);
    bool osswrite = true;

    char mode = '\0';
    string aliases;
    string sym_specific_options;
    string options;

    // DX20201228 - print python script
    if (XHOST.vflag_control.flag("PRINT_MODE::PYTHON")) {
      SYM::writePythonScript(oss);
      return true;
    }

    // AGROUP
    if (vpflow.flag("AGROUP")) {
      mode = _AGROUP_;
      aliases = "--sitepointgroup|--agroup";
      sym_specific_options = "[--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]";
      options = vpflow.getattachedscheme("AGROUP");
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_FGROUP = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL = false; // DX20170815 - Add in consistency checks //DX20171205 - Added pgroupk_xtal
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON = false; // DX20200129
      kflags.KBIN_SYMMETRY_CALCULATE_SGROUP = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_IATOMS = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_AGROUP = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_PGROUP_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE = false;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE = false;
      kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE = false; // DX20200129
      kflags.KBIN_SYMMETRY_SGROUP_WRITE = false;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE = true;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE = true;
    }
    // FGROUP
    else if (vpflow.flag("FGROUP")) {
      mode = _FGROUP_;
      aliases = "--factorgroup|--fgroup";
      sym_specific_options = "[--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]";
      options = vpflow.getattachedscheme("FGROUP");
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_FGROUP = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL = false; // DX20170815 - Add in consistency checks //DX20171205 - Added pgroupk_xtal
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON = false; // DX20200129
      kflags.KBIN_SYMMETRY_CALCULATE_SGROUP = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_IATOMS = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_AGROUP = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_PGROUP_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE = false;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE = false;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE = false;
      kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE = false; // DX20200129
      kflags.KBIN_SYMMETRY_SGROUP_WRITE = false;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE = false;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE = false;
    }
    // PGROUP
    else if (vpflow.flag("PGROUP")) {
      mode = _PGROUP_;
      aliases = "--pointgroup|--pgroup";
      options = vpflow.getattachedscheme("PGROUP");
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_FGROUP = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL = false; // DX20170815 - Add in consistency checks //DX20171205 - Added pgroupk_xtal
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON = false; // DX20200129
      kflags.KBIN_SYMMETRY_CALCULATE_SGROUP = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_IATOMS = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_AGROUP = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_PGROUP_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE = false;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE = false;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE = false;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE = false;
      kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE = false; // DX20200129
      kflags.KBIN_SYMMETRY_SGROUP_WRITE = false;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE = false;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE = false;
    }
    // PGROUPK
    else if (vpflow.flag("PGROUPK")) {
      mode = _PGROUPK_;
      aliases = "--pointgroupklattice|--pgroupk";
      options = vpflow.getattachedscheme("PGROUPK");
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_FGROUP = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL = false; // DX20170815 - Add in consistency checks //DX20171205 - Added pgroupk_xtal
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON = false; // DX20200129
      kflags.KBIN_SYMMETRY_CALCULATE_SGROUP = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_IATOMS = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_AGROUP = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_PGROUP_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE = true;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE = false;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE = false;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE = false;
      kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE = false; // DX20200129
      kflags.KBIN_SYMMETRY_SGROUP_WRITE = false;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE = false;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE = false;
    }
    // PGROUPX
    else if (vpflow.flag("PGROUPX")) {
      mode = _PGROUP_XTAL_;
      aliases = "--pointgroup_crystal|--pgroup_crystal|--pgroup_xtal|--pgroupx|--pgroupX";
      sym_specific_options = "[--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]";
      options = vpflow.getattachedscheme("PGROUPX");
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_FGROUP = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL = false; // DX20170815 - Add in consistency checks //DX20171205 - Added pgroupk_xtal
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON = false; // DX20200129
      kflags.KBIN_SYMMETRY_CALCULATE_SGROUP = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_IATOMS = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_AGROUP = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_PGROUP_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE = false;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE = false;
      kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE = false; // DX20200129
      kflags.KBIN_SYMMETRY_SGROUP_WRITE = false;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE = false;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE = false;
    }
    // PGROUPK_XTAL
    else if (vpflow.flag("PGROUPK_XTAL")) {
      mode = _PGROUPK_XTAL_;
      aliases = "--pointgroupkcrystal|--pgroupk_xtal";
      sym_specific_options = "[--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]";
      options = vpflow.getattachedscheme("PGROUPK_XTAL");
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_FGROUP = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL = true; // DX20170815 - Add in consistency checks //DX20171205 - Added pgroupk_xtal
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON = false; // DX20200129
      kflags.KBIN_SYMMETRY_CALCULATE_SGROUP = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_IATOMS = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_AGROUP = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_PGROUP_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE = false;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE = false; // DX20200129
      kflags.KBIN_SYMMETRY_SGROUP_WRITE = false;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE = false;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE = false;
    }
    // DX20200206 - add Patterson symmetry - START
    //  PGROUPK_PATTERSON
    else if (vpflow.flag("PGROUPK_PATTERSON")) {
      mode = _PGROUPK_PATTERSON_;
      aliases = "--pointgroupk_Patterson|--pgroupk_Patterson";
      sym_specific_options = "[--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]";
      options = vpflow.getattachedscheme("PGROUPK_PATTERSON");
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP = true;
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK = false;
      kflags.KBIN_SYMMETRY_CALCULATE_FGROUP = true;
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL = true;
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL = true;
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON = true;
      kflags.KBIN_SYMMETRY_CALCULATE_SGROUP = false;
      kflags.KBIN_SYMMETRY_CALCULATE_IATOMS = false;
      kflags.KBIN_SYMMETRY_CALCULATE_AGROUP = false;
      kflags.KBIN_SYMMETRY_PGROUP_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE = false;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE = true;
      kflags.KBIN_SYMMETRY_SGROUP_WRITE = false;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE = false;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE = false;
    }
    // DX20200206 - add Patterson symmetry - END
    //  SGROUP
    else if (vpflow.flag("SGROUP")) {
      mode = _SGROUP_;
      aliases = "--spacegroup|--sgroup";
      sym_specific_options = "[--radius] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]";
      options = vpflow.getattachedscheme("SGROUP");
      kflags.KBIN_SYMMETRY_SGROUP_RADIUS = KBIN_SYMMETRY_SGROUP_RADIUS_DEFAULT;
      if (vpflow.flag("SYMMETRY::SGROUP_RADIUS")) {
        kflags.KBIN_SYMMETRY_SGROUP_RADIUS = aurostd::string2utype<double>(vpflow.getattachedscheme("SYMMETRY::SGROUP_RADIUS"));
      }
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_FGROUP = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL = false; // DX20170815 - Add in consistency checks //DX20171205 - Added pgroupk_xtal
      kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON = false; // DX20200129
      kflags.KBIN_SYMMETRY_CALCULATE_SGROUP = true; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_IATOMS = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_CALCULATE_AGROUP = false; // DX20170815 - Add in consistency checks
      kflags.KBIN_SYMMETRY_PGROUP_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE = false;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE = true;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE = false;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE = false;
      kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE = false; // DX20200129
      kflags.KBIN_SYMMETRY_SGROUP_WRITE = true;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE = false;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE = false;
    }
    // EQUIVALENT / IATOMS (performed in "EQUIVALENT" function)

    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() == 1) {
      if (tokens[0] == "usage" || tokens[0] == "USAGE") {
        init::MessageOption(options, __AFLOW_FUNC__,
                            aurostd::liststring2string("aflow " + aliases + "[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=txt|json] [--screen_only] " + sym_specific_options +
                                                       " < POSCAR  default: tolerance=(minimum_interatomic_distance)/100.0, print=txt")); // DX20200724 - removed return
      }
    }
    // DX20170804 - need to rescale, so we make a fast copy and calculate
    xstructure a(_a);
    a.ReScale(1.0);
    string directory = aurostd::getPWD();
    if (XHOST.vflag_control.flag("DIRECTORY")) {
      directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");
    }
    aflags.Directory = directory;
    const string print_directory = " [dir=" + a.directory + "]";
    if (!aurostd::FileExist(directory)) {
      oss << "ERROR: Unable to locate " << directory << "." << endl;
      oss << "Exiting." << endl;
      oss << endl;
      // return oss.str();
      return false;
    }

    // DX20170921 - MAGNETIC SYMMETRY - START
    if (vpflow.flag("SYMMETRY::MAGNETIC") && (mode == _AGROUP_ || mode == _FGROUP_ || mode == _PGROUP_XTAL_ || mode == _PGROUPK_XTAL_ || mode == _SGROUP_)) {
      const string magmom_info = vpflow.getattachedscheme("SYMMETRY::MAGNETIC");
      ProcessAndAddSpinToXstructure(a, magmom_info); // DX20191108 - condensed into a single function
    }
    // DX20170921 - MAGNETIC SYMMETRY - END

    // Set tolerance
    double tolerance = pflow::getSymmetryTolerance(a, vpflow.getattachedscheme("SYMMETRY::TOLERANCE")); // DX20200820 - consolidated setting tolerance into a function
    // DX20170803 - Add format flag - START
    string format = "txt";
    if (XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      format = "txt";
    } else if (XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
      format = "json";
    } else { // default is txt
      format = "txt";
    }
    bool print = false;
    if (vpflow.flag("SYMMETRY::SCREEN_ONLY")) {
      print = true;
      kflags.KBIN_SYMMETRY_PGROUP_WRITE = false;
      kflags.KBIN_SYMMETRY_PGROUPK_WRITE = false;
      kflags.KBIN_SYMMETRY_FGROUP_WRITE = false;
      kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE = false;
      kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE = false;
      kflags.KBIN_SYMMETRY_SGROUP_WRITE = false;
      kflags.KBIN_SYMMETRY_IATOMS_WRITE = false;
      kflags.KBIN_SYMMETRY_AGROUP_WRITE = false;
      osswrite = false;
    }
    // Perform full scan
    const bool force_perform = true; // if no_scan fails, still return true at default tolerance (even though it cannot be validated)
    if (vpflow.flag("SYMMETRY::NO_SCAN")) {
      a.sym_eps_no_scan = true; // DX20210406
    }

    bool tocompress = true;
    ofstream FileMESSAGE("/dev/null");

    // while(symmetry_commensurate==false)
    if (print == false) // DX20170803 - PRINT
    { // CO20200106 - patching for auto-indenting
      for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_OUT + XHOST.vext[iext]) ||
            aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT + XHOST.vext[iext]) || // DX20200129
            aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_OUT + XHOST.vext[iext]) || // DX20180118 - added pgroupk_xtal
            aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_OUT + XHOST.vext[iext]) ||
            aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_OUT + XHOST.vext[iext])) {
          tocompress = true;
        }
      }
      for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_OUT)) {
          aurostd::RemoveFile(directory, std::regex("/" + DEFAULT_AFLOW_PGROUP_OUT + ".*"));
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_OUT)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUP_XTAL_OUT + ".*"));
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT + ".*"));
        } // DX20200129
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_OUT)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_OUT + ".*"));
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_OUT)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_XTAL_OUT + ".*"));
        } // DX20180118 - added pgroupk_xtal
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_OUT)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_FGROUP_OUT + ".*"));
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_OUT)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_AGROUP_OUT + ".*"));
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_OUT)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_SGROUP_OUT + ".*"));
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_OUT)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_IATOMS_OUT + ".*"));
        }
      }
      for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_JSON + XHOST.vext[iext]) ||
            aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON + XHOST.vext[iext]) || // DX20200129
            aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_JSON + XHOST.vext[iext]) || // DX20180118 - added pgroupk_xtal
            aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_JSON + XHOST.vext[iext]) ||
            aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_JSON + XHOST.vext[iext])) {
          tocompress = true;
        }
      }
      for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_JSON)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUP_JSON + ".*"));
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_JSON)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUP_XTAL_JSON + ".*"));
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON + ".*"));
        } // DX20200129
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_JSON)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_JSON + ".*"));
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_JSON)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_XTAL_JSON + ".*"));
        } // DX20180118 - added pgroupk_xtal
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_JSON)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_FGROUP_JSON + ".*"));
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_JSON)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_AGROUP_JSON + ".*"));
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_JSON)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_SGROUP_JSON + ".*"));
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_JSON)) {
          aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_IATOMS_JSON + ".*"));
        }
      }
    }
    if (!pflow::PerformFullSymmetry(a, tolerance, a.sym_eps_no_scan, force_perform, FileMESSAGE, aflags, kflags, osswrite, oss, format)) {
      return false;
    }
    // DX20170803 - Print to symmetry operators to screen - START
    if (print == true) {
      if (XHOST.vflag_control.flag("WWW")) {
        KBIN_SymmetryToScreenWeb(a, oss, mode);
      } else {
        KBIN_SymmetryToScreen(a, format, oss, mode);
      }
    }
    // DX20170803 - Print to symmetry operators to screen - END
    //  BZIP if necessary
    if (tocompress) {
      if (format == "txt") { // DX20170803 - FORMAT
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUP_OUT);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_OUT);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT); // DX20200129
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_OUT);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_OUT); // DX20180118 - added pgroupk_xtal
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_FGROUP_OUT);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_SGROUP_OUT);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_AGROUP_OUT);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_IATOMS_OUT);
        }
      } else if (format == "json") { // DX20170803 - FORMAT
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUP_JSON);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_JSON);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON); // DX20200129
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_JSON);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_JSON); // DX20180118 - added pgroupk_xtal
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_FGROUP_JSON);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_SGROUP_JSON);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_AGROUP_JSON);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_IATOMS_JSON);
        }
      }
    }
    // return oss.str(); //string("AFLOW "+string(AFLOW_VERSION)+" Symmetry Fixed in "+directory+"  (need aflow>=2948)\n");
    return true; // string("AFLOW "+string(AFLOW_VERSION)+" Symmetry Fixed in "+directory+"  (need aflow>=2948)\n");
  }
} // namespace pflow

// ***************************************************************************
// pflow::AGROUP
// ***************************************************************************
namespace pflow {
  void AGROUP(_aflags& aflags, istream& input) {
    //  cout << aflow::Banner("BANNER_TINY") << endl;
    aflags.QUIET = true;
    xstructure a(input, IOAFLOW_AUTO);
    const bool WRITE = true;
    ofstream File("/dev/null");
    // DX20170815 - Add in consistency checks bool verbose=true;
    // DX20170815 - Add in consistency checks SYM::CalculatePointGroup(File,a,aflags,WRITE,verbose,cout);
    // DX20170815 - Add in consistency checks SYM::CalculateFactorGroup(File,a,aflags,WRITE,verbose,cout);
    // DX20170815 - Add in consistency checks SYM::CalculateSpaceGroup(File,a,aflags,false,verbose,cout);
    // DX20170815 - Add in consistency checks SYM::CalculateInequivalentAtoms(File,a,aflags,WRITE,verbose,cout);
    // DX20170815 - Add in consistency checks SYM::CalculateSitePointGroup(File,a,aflags,WRITE,true,cout);
    _kflags kflags; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_FGROUP = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL = false; // DX20170815 - Add in consistency checks //DX20171205 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON = false; // DX20200129
    kflags.KBIN_SYMMETRY_CALCULATE_SGROUP = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_IATOMS = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_AGROUP = true; // DX20170815 - Add in consistency checks
    pflow::PerformFullSymmetry(a, File, aflags, kflags, WRITE, cout); // DX20170815 - Add in consistency checks
  }
} // namespace pflow

// ***************************************************************************
// pflow::AGROUP2
// ***************************************************************************
namespace pflow {
  void AGROUP2(istream& input) {
    xstructure a(input, IOAFLOW_AUTO);
    const int Nbasis = a.atoms.size();
    const vector<vector<uint>> Bsitesym(Nbasis);
    const bool ComMidss = false; // calculate common site symmetry of A,B,and middle_AB?
    SYM::CalculateSitePointGroup2(a, ComMidss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::AGROUP2m
// ***************************************************************************
namespace pflow {
  void AGROUP2m(istream& input) {
    xstructure a(input, IOAFLOW_AUTO);
    const int Nbasis = a.atoms.size();
    const vector<vector<uint>> Bsitesym(Nbasis);
    const bool ComMidss = true; // calculate common site symmetry of A,B,and middle_AB?
    SYM::CalculateSitePointGroup2(a, ComMidss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::ALPHABETIC
// ***************************************************************************
namespace pflow {
  xstructure ALPHABETIC(istream& input) {
    xstructure a(input, IOAFLOW_AUTO);
    //  cerr << "ALPHA=" << a.SpeciesGetAlphabetic() << endl;
    a.SpeciesPutAlphabetic();
    return a;
  }
} // namespace pflow

// ***************************************************************************
// pflow::ALPHACompound
// ***************************************************************************
namespace pflow {
  string ALPHACompound(const string& options) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.empty()) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --alpha_compound=string1,string2,....");
    }
    // move on
    stringstream output;
    for (size_t i = 0; i < tokens.size(); i++) {
      string system = string(tokens[i]);
      vector<string> vsystem;
      vector<double> vnumber;
      XATOM_AlphabetizationCompound(system, vsystem, vnumber);
      output << system << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return output.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::ALPHASpecies
// ***************************************************************************
namespace pflow {
  string ALPHASpecies(const string& options) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.empty()) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --alpha_species=string1,string2,....");
    }
    // move on
    stringstream output;
    for (size_t i = 0; i < tokens.size(); i++) {
      string system = string(tokens[i]);
      vector<string> vsystem;
      vector<double> vnumber;
      XATOM_AlphabetizationSpecies(system, vsystem, vnumber);
      output << system << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return output.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::ANGLES
// ***************************************************************************
namespace pflow {
  void ANGLES(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 1) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --angle=cutoff < POSCAR");
    }
    // move on
    double cutoff = 0.0; // some defaults
    if (!tokens.empty()) {
      cutoff = aurostd::string2utype<double>(tokens[0]);
    }

    const xstructure a(input, IOAFLOW_AUTO);
    cout << aflow::Banner("BANNER_TINY") << endl;
    pflow::PrintAngles(a, cutoff, cout);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::ATOMSMAX
// ***************************************************************************
namespace pflow {
  string ATOMSMAX(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 1) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --maxatoms=N | --max_atoms=N | --atoms_max=N | --atomsmax=N < POSCAR");
    }
    // move on
    uint N = 0; // some defaults
    if (!tokens.empty()) {
      N = aurostd::string2utype<uint>(tokens[0]);
    }

    const xstructure a(input, IOAFLOW_AUTO);
    stringstream oss;
    if (a.atoms.size() > N) {
      oss << "MAX ATOMS SIZE = " << N << endl;
      oss << "your structure has " << a.atoms.size() << " atoms" << endl;
    } else {
      oss << a;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return oss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::BANDS
// ***************************************************************************
namespace pflow {
  void BANDS(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 1) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --bands=PROOUT < POSCAR");
    }
    // move on
    string filename; // some defaults
    if (!tokens.empty()) {
      filename = (tokens[0]);
    }

    //  cout << aflow::Banner("BANNER_TINY") << endl;
    const xstructure str(input, IOAFLOW_AUTO);
    pflow::projdata prd;
    prd.PROOUTinfile = filename;
    pflow::ReadInProj(prd);
    prd.lat = pflow::GetScaledLat(str);
    pflow::PrintBands(prd);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::BANDGAP // CAMILO
// ***************************************************************************
namespace pflow {
  void BANDGAP(aurostd::xoption& vpflow, ostream& oss) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    const string input = aurostd::RemoveWhiteSpaces(vpflow.getattachedscheme("BANDGAP"));
    if (input.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "input empty", _INPUT_MISSING_);
    }
    vector<string> dirs;
    aurostd::string2tokens(input, dirs, ",");
    for (size_t i = 0; i < dirs.size(); i++) {
      if (!PrintBandGap(dirs[i], oss)) {
        oss << __AFLOW_FUNC__ << " " + dirs[i] + " failed" << endl;
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
  void BANDGAP_DOS(aurostd::xoption& vpflow, ostream& oss) { // CO20191004
    const bool LDEBUG = (false || XHOST.DEBUG);
    const string input = aurostd::RemoveWhiteSpaces(vpflow.getattachedscheme("BANDGAPDOS"));
    if (input.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "input empty", _INPUT_MISSING_);
    }
    vector<string> dirs;
    aurostd::string2tokens(input, dirs, ",");
    for (size_t i = 0; i < dirs.size(); i++) {
      if (!PrintBandGap_DOS(dirs[i], oss)) {
        oss << __AFLOW_FUNC__ << " " + dirs[i] + " failed" << endl;
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::BANDSTRUCTURE
// ***************************************************************************
namespace pflow {
  void BANDSTRUCTURE(_aflags& aflags) {
    //  cout << aflow::Banner("BANNER_TINY") << endl;
    aflags.QUIET = true;
    string directory_LIB;
    string directory_RAW;
    directory_LIB = aflags.Directory;
    directory_RAW = directory_LIB + "./BANDS";
    aurostd::RemoveDirectory(directory_RAW);
    aurostd::DirectoryMake(directory_RAW);
    vector<string> vspecies;
    vector<string> vfiles;
    aflowlib::_aflowlib_entry data;
    aflowlib::GetSpeciesDirectory(directory_LIB, vspecies);
    for (size_t i = 0; i < vspecies.size(); i++) {
      data.species += vspecies[i];
      if (i < vspecies.size() - 1) {
        data.species += ",";
      }
    }
    aflowlib::LIB2RAW_Loop_Bands(directory_LIB, directory_RAW, vfiles, data, "pflow::BANDSTRUCTURE");
  }
} // namespace pflow

// ***************************************************************************
// pflow::BZDirectionsLATTICE
// ***************************************************************************
namespace pflow {
  string BZDirectionsLATTICE(const string& options) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 1) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --bzdirections= | --bzd=LATTICE");
    } else { // move on
      if (tokens.at(0) != "TRI1a" && tokens.at(0) != "TRI1b" && tokens.at(0) != "TRI2a" && tokens.at(0) != "TRI2b" && tokens.at(0) != "MCL" && tokens.at(0) != "MCLC1" && tokens.at(0) != "MCLC2" &&
          tokens.at(0) != "MCLC3" && tokens.at(0) != "MCLC4" && tokens.at(0) != "MCLC5" && tokens.at(0) != "ORC" && tokens.at(0) != "ORCC" && tokens.at(0) != "ORCF1" && tokens.at(0) != "ORCF2" &&
          tokens.at(0) != "ORCF3" && tokens.at(0) != "ORCI" && tokens.at(0) != "TET" && tokens.at(0) != "BCT1" && tokens.at(0) != "BCT2" && tokens.at(0) != "RHL1" && tokens.at(0) != "RHL2" &&
          tokens.at(0) != "HEX" && tokens.at(0) != "CUB" && tokens.at(0) != "FCC" && tokens.at(0) != "BCC") {
        cout << "LATTICES = (the ones with \"\") " << endl;
        cout << "1. TRI order: kalpha,kbeta,kgamma  > 90 (kgamma<kalpha, kgamma<kbeta) " << endl;
        cout << "   or kalpha,kbeta,kgamma  < 90 (kgamma>kalpha, kgamma>kbeta) " << endl;
        cout << "   special case when kgamma=90 " << endl;
        cout << "   \"TRI1a\" kalpha>90 kbeta>90 kgamma>90 " << endl;
        cout << "   \"TRI1b\" kalpha<90 kbeta<90 kgamma<90 " << endl;
        cout << "   \"TRI2a\" kalpha>90 kbeta>90 kgamma=90 " << endl;
        cout << "   \"TRI2b\" kalpha<90 kbeta<90 kgamma=90 " << endl;
        cout << "2. \"MCL\" unique (order b<=c) " << endl;
        cout << "3. MCLC (order alpha<90) " << endl;
        cout << "   \"MCLC1\"  kgamma>90 " << endl;
        cout << "   \"MCLC2\"  kgamma=90 " << endl;
        cout << "   \"MCLC3\"  kgamma<90, b*cos(alpha)/c + (b*sin(alpha)/a)^2 < 1 " << endl;
        cout << "   \"MCLC4\"  kgamma<90, b*cos(alpha)/c + (b*sin(alpha)/a)^2 = 1 " << endl;
        cout << "   \"MCLC5\"  kgamma<90, b*cos(alpha)/c + (b*sin(alpha)/a)^2 > 1 " << endl;
        cout << "4. \"ORC\" unique (order a<b<c) " << endl;
        cout << "5. \"ORCC\" unique (order a<b) " << endl;
        cout << "6. ORCF (order a<b<c) " << endl;
        cout << R"(   "ORCF1" "ORCF_invb2+invc2<inva2"  for 1/a^2 > 1/b^2 + 1/c^2 )" << endl;
        cout << R"(   "ORCF2" "ORCF_inva2<invb2+invc2"  for 1/a^2 < 1/b^2 + 1/c^2 )" << endl;
        cout << "   \"ORCF3\"                           for 1/a^2 = 1/b^2 + 1/c^2 " << endl;
        cout << "7. \"ORCI\" unique (order a<b<c) " << endl;
        cout << "8. \"TET\" unique (order a a c) " << endl;
        cout << "9. BCT (order a a c) " << endl;
        cout << R"(   "BCT1" "BCT_c<a" for c<a )" << endl;
        cout << R"(   "BCT2" "BCT_c>a" for c>a )" << endl;
        cout << "10. \"RHL1\" alpha<90 " << endl;
        cout << "    \"RHL2\" alpha>90 " << endl;
        cout << "11. \"HEX\" unique (order 60 90 90) " << endl;
        cout << "12. \"CUB\" unique " << endl;
        cout << "13. \"FCC\" unique (order 60 60 60) " << endl;
        cout << "14. \"BCC\" unique " << endl;
        init::ErrorOption(options, __AFLOW_FUNC__, "aflow --bzdirections= | --bzd=LATTICE");
      }
    }

    const double grid = 20.0;
    const xmatrix<double> rlattice(3, 3);
    bool foundBZ;
    rlattice[1][1] = 1.0;
    rlattice[2][2] = 1.0;
    rlattice[3][3] = 1.0;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return LATTICE::KPOINTS_Directions(tokens.at(0), rlattice, grid, IOVASP_AUTO, foundBZ);
  }
} // namespace pflow

// ***************************************************************************
// pflow::BZDirectionsSTRUCTURE
// ***************************************************************************
namespace pflow {
  string BZDirectionsSTRUCTURE(istream& input, aurostd::xoption& vpflow) {
    xstructure a(input, IOAFLOW_AUTO);
    bool foundBZ;
    const double grid = 20.0;
    xstructure str_sp;
    xstructure str_sc; // DX20181213 - need primitive lattice information for kpoint directions
    a.GetLatticeType(str_sp, str_sc); // DX20181213 - need primitive lattice information for kpoint directions

    // if not transforming kpoints, set transformation matrix to identity
    if (!vpflow.flag("BZDIRECTION::TRANSFORM2ORIGINAL")) {
      a.transform_coordinates_original2new = aurostd::eye<double>(); // CO20190520
    } else if (a.volume_changed_original2new == true) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "The volume changed between the input and final structure; cannot transform kpoints", _RUNTIME_ERROR_); // CO20200624
    }
    cerr << "LATTICE INFORMATION" << endl;
    cerr << " Real space lattice primitive           = " << a.bravais_lattice_type << endl;
    cerr << " Real space lattice variation           = " << a.bravais_lattice_variation_type << endl; // WSETYAWAN mod
    //  cerr << " Real space conventional lattice        = " << a.bravais_conventional_lattice_type << endl;
    cerr << " Real space Pearson symbol              = " << a.pearson_symbol << endl;
    cerr << " Reciprocal lattice primitive           = " << a.reciprocal_lattice_type << endl;
    cerr << " Reciprocal lattice variation           = " << a.reciprocal_lattice_variation_type << endl; // WSETYAWAN mod
    if (vpflow.flag("BZDIRECTION::TRANSFORM2ORIGINAL") && vpflow.flag("BZDIRECTION::PRINT_TRANSFORMATION_MATRIX")) {
      cerr << "TRANSFORMATION MATRIX (P)" << endl;
      cerr << roundoff(aurostd::inverse(a.transform_coordinates_original2new), 1e-8) << endl;
      cerr << "TRANSFORMATION MATRIX (Q=P^-1)" << endl;
      cerr << roundoff(a.transform_coordinates_original2new, 1e-8) << endl;
    }
    //  cerr << " Reciprocal conventional lattice        = " << a.reciprocal_conventional_lattice_type << endl;
    cerr << "KPOINTS FILE" << endl;
    // DX20181101 return LATTICE::KPOINTS_Directions(a.bravais_lattice_variation_type,a.lattice,grid,a.iomode,foundBZ);
    return LATTICE::KPOINTS_Directions(a.bravais_lattice_variation_type, str_sp.lattice, a.transform_coordinates_original2new, grid, a.iomode, foundBZ); // DX20181213 - pass primitive lattice, not the original
  }
} // namespace pflow

// ***************************************************************************
// pflow::CAGES
// ***************************************************************************
namespace pflow {
  void CAGES(_aflags& aflags, const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (!tokens.empty() && tokens.size() != 1) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow [--np=NP] --cages[=roughness] < POSCAR");
    }
    double roughness = -1.0;
    if (!tokens.empty()) {
      roughness = aurostd::string2utype<double>(tokens.at(0));
    }
    const xstructure a(input, IOAFLOW_AUTO);
    vector<acage> cagesirreducible;
    vector<acage> cagesreducible;
    vector<acage> cages4;
    vector<acage> cages3;
    vector<acage> cages2;
    // ME20210521
    const bool ofwrite = !XHOST.vflag_control.flag("WWW");
    GetCages(a, aflags, cagesirreducible, cagesreducible, cages4, cages3, cages2, roughness, ofwrite, true, cout);
    //  cout << "REDUCIBLE_SIZE " << cagesreducible.size() << endl;
    // cout << "IRREDUCIBLE_SIZE " << cagesirreducible.size() << endl;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::CART
// ***************************************************************************
namespace pflow {
  xstructure CART(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    xstructure b(a);
    b.SetCoordinates(_COORDS_CARTESIAN_);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::ChangeSuffix(const string& options)
// ***************************************************************************
namespace pflow {
  void ChangeSuffix(const string& _operation) {
    // aflow --suffix=[directory,]"from2to"
    //  Change the suffixes of VASP files. Easy conversion between AFLOW format and VASP format. (KY20131222)
    //	Mnemonic: from2to with from/to =[n=none;r1=relax1;r2=relax2;r3=relax3;s=static;b=bands] or without abbreviations.
    vector<string> vfile{"AECCAR0", "AECCAR1", "AECCAR2", "INCAR", "POSCAR", "POTCAR", "KPOINTS", "CONTCAR", "OUTCAR", "EIGENVAL",    "DOSCAR",        "IBZKPT",
                         "OSZICAR", "PCDAT",   "XDATCAR", "CHG",   "CHGCAR", "PROCAR", "ELFCAR",  "WAVECAR", "LOCPOT", "vasprun.xml", DEFAULT_VASP_OUT};
    vector<string> tokens;
    aurostd::string2tokens(_operation, tokens, ",");
    string directory;
    string operation;
    if (tokens.empty()) {
      directory = "";
      operation = _operation;
    }
    if (tokens.size() == 1) {
      directory = "";
      operation = tokens.at(0);
    }
    if (tokens.size() == 2) {
      directory = tokens.at(0) + "/";
      operation = tokens.at(1);
    }
    if (!operation.empty()) {
      for (size_t i = 0; i < vfile.size(); i++) {
        string from = "";
        string to = "";
        string f = aurostd::CleanFileName(directory + vfile[i]);
        for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
          if (aurostd::FileExist(f + XHOST.vext[iext])) {
            aurostd::execute(XHOST.command(XHOST.vzip.at(iext)) + " -dqf " + f + XHOST.vext[iext]);
          }
          if (aurostd::FileExist(f + ".relax1" + XHOST.vext[iext])) {
            aurostd::execute(XHOST.command(XHOST.vzip.at(iext)) + " -dqf " + f + ".relax1" + XHOST.vext[iext]);
          }
          if (aurostd::FileExist(f + ".relax2" + XHOST.vext[iext])) {
            aurostd::execute(XHOST.command(XHOST.vzip.at(iext)) + " -dqf " + f + ".relax2" + XHOST.vext[iext]);
          }
          if (aurostd::FileExist(f + ".relax3" + XHOST.vext[iext])) {
            aurostd::execute(XHOST.command(XHOST.vzip.at(iext)) + " -dqf " + f + ".relax3" + XHOST.vext[iext]);
          }
          if (aurostd::FileExist(f + ".static" + XHOST.vext[iext])) {
            aurostd::execute(XHOST.command(XHOST.vzip.at(iext)) + " -dqf " + f + ".static" + XHOST.vext[iext]);
          }
          if (aurostd::FileExist(f + ".bands" + XHOST.vext[iext])) {
            aurostd::execute(XHOST.command(XHOST.vzip.at(iext)) + " -dqf " + f + ".bands" + XHOST.vext[iext]);
          }
        }
        if (aurostd::substring2bool(operation, "n2") || aurostd::substring2bool(operation, "none2") || aurostd::substring2bool(operation, "normal2")) {
          from = "";
        }
        if (aurostd::substring2bool(operation, "r12") || aurostd::substring2bool(operation, "relax12")) {
          from = ".relax1";
        }
        if (aurostd::substring2bool(operation, "r22") || aurostd::substring2bool(operation, "relax22")) {
          from = ".relax2";
        }
        if (aurostd::substring2bool(operation, "r32") || aurostd::substring2bool(operation, "relax32")) {
          from = ".relax3";
        }
        if (aurostd::substring2bool(operation, "s2") || aurostd::substring2bool(operation, "static2")) {
          from = ".static";
        }
        if (aurostd::substring2bool(operation, "b2") || aurostd::substring2bool(operation, "bands2")) {
          from = ".bands";
        }
        if (aurostd::substring2bool(operation, "2n") || aurostd::substring2bool(operation, "2none") || aurostd::substring2bool(operation, "2normal")) {
          to = "";
        }
        if (aurostd::substring2bool(operation, "2r1") || aurostd::substring2bool(operation, "2relax1")) {
          to = ".relax1";
        }
        if (aurostd::substring2bool(operation, "2r2") || aurostd::substring2bool(operation, "2relax2")) {
          to = ".relax2";
        }
        if (aurostd::substring2bool(operation, "2r3") || aurostd::substring2bool(operation, "2relax3")) {
          to = ".relax3";
        }
        if (aurostd::substring2bool(operation, "2s") || aurostd::substring2bool(operation, "2static")) {
          to = ".static";
        }
        if (aurostd::substring2bool(operation, "2b") || aurostd::substring2bool(operation, "2bands")) {
          to = ".bands";
        }

        if (!from.empty() || !to.empty()) {
          if (aurostd::FileExist(f + from)) {
            cout << XPID << "pflow::ChangeSuffix: mv " << f << from << " " << f << to << endl;
            aurostd::file2file(f + from, f + to);
          }
        }
      }
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::CheckIntegritiy
// ***************************************************************************
// Shidong Wang 2011
namespace pflow {
  void CheckIntegritiy() {
    // Check whether function `isequal', `isdifferent', `identical' etc
    // give correct answer or not
    // issue of g++ optimization (-O -O1 -O2 -O3).

    // xvector operators

    const int test_xv_size = 5;
    const int _max_int = 100;
    const int _tol_int = 0;
    const float _tol_float = 1.0e-6;
    const double _tol_double = 1.0e-6;

    xvector<int> test_int_xv_a(test_xv_size);
    xvector<int> test_int_xv_b(test_xv_size);
    xvector<float> test_float_xv_a(test_xv_size);
    xvector<float> test_float_xv_b(test_xv_size);
    xvector<double> test_double_xv_a(test_xv_size);
    xvector<double> test_double_xv_b(test_xv_size);

    // initialze random seed
    srand(time(nullptr));

    // initializae test xvectors
    for (int i = 1; i < test_xv_size + 1; i++) {
      test_int_xv_a[i] = rand() % _max_int;
      test_int_xv_b[i] = rand() % _max_int;
      test_float_xv_a[i] = PI * float((rand() % _max_int));
      test_float_xv_b[i] = PI * float((rand() % _max_int));
      test_double_xv_a[i] = EULERSNUMBER * float((rand() % _max_int));
      test_double_xv_b[i] = EULERSNUMBER * float((rand() % _max_int));
    }

    const int line_width = 60;
    for (int i = 0; i < line_width; i++) {
      cout << "*";
    }
    cout << endl;
    cout << "Shidong Wang - 2011 " << endl;
    for (int i = 0; i < line_width; i++) {
      cout << "*";
    }
    cout << endl;
    cout << "xvector\n";

    // template<class utype> bool
    // identical(const xvector<utype>&,const xvector<utype>&,const utype&) __xprototype;
    // cout << "Indentical Integer with given tolerance " << _tol_int << " : \n"

    if (identical(test_int_xv_a, test_int_xv_a, _tol_int) && (!identical(test_int_xv_a, test_int_xv_b, _tol_int))) {
      cout << ">>> Testing Integer INDENTICAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer INDENTICAL with given tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xv_a << endl
           << " and \n"
           << test_int_xv_a << endl
           << " are identical? " << std::boolalpha << identical(test_int_xv_a, test_int_xv_a, _tol_int) << endl;
      cout << "2) \n" << test_int_xv_a << endl << " and \n" << test_int_xv_b << endl << " are not identical? " << std::boolalpha << (!identical(test_int_xv_a, test_int_xv_b, _tol_int)) << endl;
    }

    if (identical(test_float_xv_a, test_float_xv_a, _tol_float) && (!identical(test_float_xv_a, test_float_xv_b, _tol_float))) {
      cout << ">>> Testing Float INDENTICAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float INDENTICAL with given tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xv_a << endl
           << " and \n"
           << test_float_xv_a << endl
           << " are identical? " << std::boolalpha << identical(test_float_xv_a, test_float_xv_a, _tol_float) << endl;
      cout << "2) \n" << test_float_xv_a << endl << " and \n" << test_float_xv_b << endl << " are not identical? " << std::boolalpha << (!identical(test_float_xv_a, test_float_xv_b, _tol_float)) << endl;
    }

    if (identical(test_double_xv_a, test_double_xv_a, _tol_double) && (!identical(test_double_xv_a, test_double_xv_b, _tol_double))) {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_double_xv_a << endl
           << " and \n"
           << test_double_xv_a << endl
           << " are identical? " << std::boolalpha << identical(test_double_xv_a, test_double_xv_a, _tol_double) << endl;
      cout << "2) \n" << test_double_xv_a << endl << " and \n" << test_double_xv_b << endl << " are not identical? " << std::boolalpha << (!identical(test_double_xv_a, test_double_xv_b, _tol_double)) << endl;
    }

    // template<class utype> bool
    // identical(const xvector<utype>&,const xvector<utype>&) __xprototype;
    if (identical(test_int_xv_a, test_int_xv_a) && (!identical(test_int_xv_a, test_int_xv_b))) {
      cout << ">>> Testing Integer INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer INDENTICAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xv_a << endl
           << " and \n"
           << test_int_xv_a << endl
           << " are identical? " << std::boolalpha << identical(test_int_xv_a, test_int_xv_a) << endl;
      cout << "2) \n" << test_int_xv_a << endl << " and \n" << test_int_xv_b << endl << " are not identical? " << std::boolalpha << (!identical(test_int_xv_a, test_int_xv_b)) << endl;
    }

    if (identical(test_float_xv_a, test_float_xv_a) && (!identical(test_float_xv_a, test_float_xv_b))) {
      cout << ">>> Testing Float INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float INDENTICAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xv_a << endl
           << " and \n"
           << test_float_xv_a << endl
           << " are identical? " << std::boolalpha << identical(test_float_xv_a, test_float_xv_a) << endl;
      cout << "2) \n" << test_float_xv_a << endl << " and \n" << test_float_xv_b << endl << " are not identical? " << std::boolalpha << (!identical(test_float_xv_a, test_float_xv_b)) << endl;
    }

    if (identical(test_double_xv_a, test_double_xv_a) && (!identical(test_double_xv_a, test_double_xv_b))) {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_double_xv_a << endl
           << " and \n"
           << test_double_xv_a << endl
           << " are identical? " << std::boolalpha << identical(test_double_xv_a, test_double_xv_a) << endl;
      cout << "2) \n" << test_double_xv_a << endl << " and \n" << test_double_xv_b << endl << " are not identical? " << std::boolalpha << (!identical(test_double_xv_a, test_double_xv_b)) << endl;
    }

    // template<class utype> bool
    // operator==(const xvector<utype>&,const xvector<utype>&) __xprototype;
    if ((test_int_xv_a == test_int_xv_a) && !(test_int_xv_a == test_int_xv_b)) {
      cout << ">>> Testing Integer == with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer == with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xv_a << endl
           << " and \n"
           << test_int_xv_a << endl
           << " are identical? " << std::boolalpha << identical(test_int_xv_a, test_int_xv_a) << endl;
      cout << "2) \n" << test_int_xv_a << endl << " and \n" << test_int_xv_b << endl << " are not identical? " << std::boolalpha << (!identical(test_int_xv_a, test_int_xv_b)) << endl;
    }

    if ((test_float_xv_a == test_float_xv_a) && (!(test_float_xv_a == test_float_xv_b))) {
      cout << ">>> Testing Float == with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float == with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xv_a << endl
           << " and \n"
           << test_float_xv_a << endl
           << " are identical? " << std::boolalpha << (test_float_xv_a == test_float_xv_a) << endl;
      cout << "2) \n" << test_float_xv_a << endl << " and \n" << test_float_xv_b << endl << " are not identical? " << std::boolalpha << (!(test_float_xv_a == test_float_xv_b)) << endl;
    }

    if (identical(test_double_xv_a, test_double_xv_a) && (!identical(test_double_xv_a, test_double_xv_b))) {
      cout << ">>> Testing Double == with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double == with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_double_xv_a << endl
           << " and \n"
           << test_double_xv_a << endl
           << " are identical? " << std::boolalpha << (test_double_xv_a == test_double_xv_a) << endl;
      cout << "2) \n" << test_double_xv_a << endl << " and \n" << test_double_xv_b << endl << " are not identical? " << std::boolalpha << (!(test_double_xv_a == test_double_xv_b)) << endl;
    }

    // template<class utype> bool
    // isdifferent(const xvector<utype>&,const xvector<utype>&,const utype&) __xprototype;
    if (!isdifferent(test_int_xv_a, test_int_xv_a, _tol_int) && isdifferent(test_int_xv_a, test_int_xv_b, _tol_int)) {
      cout << ">>> Testing Integer ISDIFFERENT with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISDIFFERENT with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xv_a << endl
           << " and \n"
           << test_float_xv_a << endl
           << " are different? " << std::boolalpha << isdifferent(test_int_xv_a, test_int_xv_a, _tol_int) << endl;
      cout << "2) \n" << test_float_xv_a << endl << " and \n" << test_float_xv_b << endl << " are not different? " << std::boolalpha << (!isdifferent(test_int_xv_a, test_int_xv_b, _tol_int)) << endl;
    }

    if (!isdifferent(test_float_xv_a, test_float_xv_a, _tol_float) && (isdifferent(test_float_xv_a, test_float_xv_b, _tol_float))) {
      cout << ">>> Testing Float ISDIFFERENT with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISDIFFERENT with given tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xv_a << endl
           << " and \n"
           << test_float_xv_a << endl
           << " are different? " << std::boolalpha << isdifferent(test_float_xv_a, test_float_xv_a, _tol_float) << endl;
      cout << "2) \n" << test_float_xv_a << endl << " and \n" << test_float_xv_b << endl << " are not different? " << std::boolalpha << (isdifferent(test_float_xv_a, test_float_xv_b, _tol_float)) << endl;
    }

    if (!isdifferent(test_double_xv_a, test_double_xv_a, _tol_double) && (isdifferent(test_double_xv_a, test_double_xv_b, _tol_double))) {
      cout << ">>> Testing Double ISDIFFERENT with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISDIFFERENT with given tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_double_xv_a << endl
           << " and \n"
           << test_double_xv_a << endl
           << " are different? " << std::boolalpha << isdifferent(test_double_xv_a, test_double_xv_a, _tol_double) << endl;
      cout << "2) \n" << test_double_xv_a << endl << " and \n" << test_double_xv_b << endl << " are not different? " << std::boolalpha << (!isdifferent(test_double_xv_a, test_double_xv_b, _tol_double)) << endl;
    }

    // template<class utype> bool
    // isdifferent(const xvector<utype>&,const xvector<utype>&) __xprototype;
    if (!isdifferent(test_int_xv_a, test_int_xv_a) && isdifferent(test_int_xv_a, test_int_xv_b)) {
      cout << ">>> Testing Integer ISDIFFERENT xvector with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISDIFFERENT xvector with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xv_a << endl
           << " and \n"
           << test_int_xv_a << endl
           << " are different? " << std::boolalpha << isdifferent(test_int_xv_a, test_int_xv_a) << endl;
      cout << "2) \n" << test_int_xv_a << endl << " and \n" << test_int_xv_b << endl << " are not different? " << std::boolalpha << (!isdifferent(test_int_xv_a, test_int_xv_b)) << endl;
    }

    if (!isdifferent(test_float_xv_a, test_float_xv_a) && (isdifferent(test_float_xv_a, test_float_xv_b))) {
      cout << ">>> Testing Float ISDIFFERENT with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISDIFFERENT with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xv_a << endl
           << " and \n"
           << test_float_xv_a << endl
           << " are different? " << std::boolalpha << isdifferent(test_float_xv_a, test_float_xv_a) << endl;
      cout << "2) \n" << test_int_xv_a << endl << " and \n" << test_int_xv_b << endl << " are not different? " << std::boolalpha << (isdifferent(test_float_xv_a, test_float_xv_b)) << endl;
    }

    if (!isdifferent(test_double_xv_a, test_double_xv_a) && (isdifferent(test_double_xv_a, test_double_xv_b))) {
      cout << ">>> Testing Double ISDIFFERENT with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISDIFFERENT with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_double_xv_a << endl
           << " and \n"
           << test_double_xv_a << endl
           << " are different? " << std::boolalpha << isdifferent(test_double_xv_a, test_double_xv_a) << endl;
      cout << "2) \n" << test_int_xv_a << endl << " and \n" << test_int_xv_b << endl << " are not different? " << std::boolalpha << (!isdifferent(test_double_xv_a, test_double_xv_b)) << endl;
    }

    // template<class utype> bool
    // aurostd::isequal(const xvector<utype>&,const xvector<utype>&,const utype&) __xprototype;
    if (aurostd::isequal(test_int_xv_a, test_int_xv_a, _tol_int) && (!aurostd::isequal(test_int_xv_a, test_int_xv_b, _tol_int))) {
      cout << ">>> Testing Integer ISEQUAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISEQUAL with given tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xv_a << endl
           << " and \n"
           << test_int_xv_a << endl
           << " are identical? " << std::boolalpha << aurostd::isequal(test_int_xv_a, test_int_xv_a, _tol_int) << endl;
      cout << "2) \n" << test_int_xv_a << endl << " and \n" << test_int_xv_b << endl << " are not identical? " << std::boolalpha << (!aurostd::isequal(test_int_xv_a, test_int_xv_b, _tol_int)) << endl;
    }

    if (aurostd::isequal(test_float_xv_a, test_float_xv_a, _tol_float) && (!aurostd::isequal(test_float_xv_a, test_float_xv_b, _tol_float))) {
      cout << ">>> Testing Float ISEQUAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISEQUAL with given tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xv_a << endl
           << " and \n"
           << test_float_xv_a << endl
           << " are identical? " << std::boolalpha << aurostd::isequal(test_float_xv_a, test_float_xv_a, _tol_float) << endl;
      cout << "2) \n" << test_float_xv_a << endl << " and \n" << test_float_xv_b << endl << " are not identical? " << std::boolalpha << (!aurostd::isequal(test_float_xv_a, test_float_xv_b, _tol_float)) << endl;
    }

    if (aurostd::isequal(test_double_xv_a, test_double_xv_a, _tol_double) && (!aurostd::isequal(test_double_xv_a, test_double_xv_b, _tol_double))) {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_double_xv_a << endl
           << " and \n"
           << test_double_xv_a << endl
           << " are identical? " << std::boolalpha << aurostd::isequal(test_double_xv_a, test_double_xv_a, _tol_double) << endl;
      cout << "2) \n"
           << test_double_xv_a << endl
           << " and \n"
           << test_double_xv_b << endl
           << " are not identical? " << std::boolalpha << (!aurostd::isequal(test_double_xv_a, test_double_xv_b, _tol_double)) << endl;
    }

    // template<class utype> bool
    // isequal(const xvector<utype>&,const xvector<utype>&) __xprototype;
    if (aurostd::isequal(test_int_xv_a, test_int_xv_a) && (!aurostd::isequal(test_int_xv_a, test_int_xv_b))) {
      cout << ">>> Testing Integer ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISEQUAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xv_a << endl
           << " and \n"
           << test_int_xv_a << endl
           << " are identical? " << std::boolalpha << aurostd::isequal(test_int_xv_a, test_int_xv_a) << endl;
      cout << "2) \n" << test_int_xv_a << endl << " and \n" << test_int_xv_b << endl << " are not identical? " << std::boolalpha << (!aurostd::isequal(test_int_xv_a, test_int_xv_b)) << endl;
    }

    if (aurostd::isequal(test_float_xv_a, test_float_xv_a) && (!aurostd::isequal(test_float_xv_a, test_float_xv_b))) {
      cout << ">>> Testing Float ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISEQUAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xv_a << endl
           << " and \n"
           << test_float_xv_a << endl
           << " are identical? " << std::boolalpha << aurostd::isequal(test_float_xv_a, test_float_xv_a) << endl;
      cout << "2) \n" << test_float_xv_a << endl << " and \n" << test_float_xv_b << endl << " are not identical? " << std::boolalpha << (!aurostd::isequal(test_float_xv_a, test_float_xv_b)) << endl;
    }

    if (aurostd::isequal(test_double_xv_a, test_double_xv_a) && (!aurostd::isequal(test_double_xv_a, test_double_xv_b))) {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_double_xv_a << endl
           << " and \n"
           << test_double_xv_a << endl
           << " are identical? " << std::boolalpha << aurostd::isequal(test_double_xv_a, test_double_xv_a) << endl;
      cout << "2) \n" << test_double_xv_a << endl << " and \n" << test_double_xv_b << endl << " are not identical? " << std::boolalpha << (!aurostd::isequal(test_double_xv_a, test_double_xv_b)) << endl;
    }

    // template<class utype> bool
    // operator!=(const xvector<utype>&,const xvector<utype>&) __xprototype;
    if (!(test_int_xv_a != test_int_xv_a) && ((test_int_xv_a != test_int_xv_b))) {
      cout << ">>> Testing Integer != with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer != with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xv_a << endl
           << " and \n"
           << test_int_xv_a << endl
           << " are different? " << std::boolalpha << (test_int_xv_a != test_int_xv_a) << endl;
      cout << "2) \n" << test_int_xv_a << endl << " and \n" << test_int_xv_b << endl << " are not different? " << std::boolalpha << ((test_int_xv_a != test_int_xv_b)) << endl;
    }

    if (!(test_float_xv_a != test_float_xv_a) && ((test_float_xv_a != test_float_xv_b))) {
      cout << ">>> Testing Float != with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float != with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xv_a << endl
           << " and \n"
           << test_int_xv_a << endl
           << " are different? " << std::boolalpha << (test_float_xv_a != test_float_xv_a) << endl;
      cout << "2) \n" << test_int_xv_a << endl << " and \n" << test_int_xv_b << endl << " are not different? " << std::boolalpha << ((test_float_xv_a != test_float_xv_b)) << endl;
    }

    if (!(test_double_xv_a != test_double_xv_a) && ((test_double_xv_a != test_double_xv_b))) {
      cout << ">>> Testing Double != with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double != with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xv_a << endl
           << " and \n"
           << test_int_xv_a << endl
           << " are different? " << std::boolalpha << (test_double_xv_a != test_double_xv_a) << endl;
      cout << "2) \n" << test_int_xv_a << endl << " and \n" << test_int_xv_b << endl << " are not different? " << std::boolalpha << ((test_double_xv_a != test_double_xv_b)) << endl;
    }

    const int test_xm_size = 3;
    xmatrix<int> test_int_xm_a(test_xm_size, test_xm_size);
    xmatrix<int> test_int_xm_b(test_xm_size, test_xm_size);
    xmatrix<float> test_float_xm_a(test_xm_size, test_xm_size);
    xmatrix<float> test_float_xm_b(test_xm_size, test_xm_size);
    xmatrix<double> test_double_xm_a(test_xm_size, test_xm_size);
    xmatrix<double> test_double_xm_b(test_xm_size, test_xm_size);

    // initialze random seed
    srand(time(nullptr));

    // initializae test xmatrices
    for (int i = 1; i < test_xm_size + 1; i++) {
      for (int j = 1; j < test_xm_size + 1; j++) {
        test_int_xm_a[i][j] = rand() % _max_int;
        test_int_xm_b[i][j] = rand() % _max_int;
        test_float_xm_a[i][j] = PI * float((rand() % _max_int));
        test_float_xm_b[i][j] = PI * float((rand() % _max_int));
        test_double_xm_a[i][j] = EULERSNUMBER * double((rand() % _max_int));
        test_double_xm_b[i][j] = EULERSNUMBER * double((rand() % _max_int));
      }
    }

    for (int i = 0; i < line_width; i++) {
      cout << "*";
    }
    cout << endl;

    cout << "xmatrix\n";

    // template<class utype> bool
    // identical(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;
    // cout << "Indentical Integer with given tolerance " << _tol_int << " : \n"

    if (identical(test_int_xm_a, test_int_xm_a, _tol_int) && (!identical(test_int_xm_a, test_int_xm_b, _tol_int))) {
      cout << ">>> Testing Integer INDENTICAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer INDENTICAL with given tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xm_a << endl
           << " and \n"
           << test_int_xm_a << endl
           << " are identical? " << std::boolalpha << identical(test_int_xm_a, test_int_xm_a, _tol_int) << endl;
      cout << "2) \n" << test_int_xm_a << endl << " and \n" << test_int_xm_b << endl << " are not identical? " << std::boolalpha << (!identical(test_int_xm_a, test_int_xm_b, _tol_int)) << endl;
    }

    if (identical(test_float_xm_a, test_float_xm_a, _tol_float) && (!identical(test_float_xm_a, test_float_xm_b, _tol_float))) {
      cout << ">>> Testing Float INDENTICAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float INDENTICAL with given tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xm_a << endl
           << " and \n"
           << test_float_xm_a << endl
           << " are identical? " << std::boolalpha << identical(test_float_xm_a, test_float_xm_a, _tol_float) << endl;
      cout << "2) \n" << test_float_xm_a << endl << " and \n" << test_float_xm_b << endl << " are not identical? " << std::boolalpha << (!identical(test_float_xm_a, test_float_xm_b, _tol_float)) << endl;
    }

    if (identical(test_double_xm_a, test_double_xm_a, _tol_double) && (!identical(test_double_xm_a, test_double_xm_b, _tol_double))) {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_double_xm_a << endl
           << " and \n"
           << test_double_xm_a << endl
           << " are identical? " << std::boolalpha << identical(test_double_xm_a, test_double_xm_a, _tol_double) << endl;
      cout << "2) \n" << test_double_xm_a << endl << " and \n" << test_double_xm_b << endl << " are not identical? " << std::boolalpha << (!identical(test_double_xm_a, test_double_xm_b, _tol_double)) << endl;
    }

    // template<class utype> bool
    // identical(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;
    if (identical(test_int_xm_a, test_int_xm_a) && (!identical(test_int_xm_a, test_int_xm_b))) {
      cout << ">>> Testing Integer INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer INDENTICAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xm_a << endl
           << " and \n"
           << test_int_xm_a << endl
           << " are identical? " << std::boolalpha << identical(test_int_xm_a, test_int_xm_a) << endl;
      cout << "2) \n" << test_int_xm_a << endl << " and \n" << test_int_xm_b << endl << " are not identical? " << std::boolalpha << (!identical(test_int_xm_a, test_int_xm_b)) << endl;
    }

    if (identical(test_float_xm_a, test_float_xm_a) && (!identical(test_float_xm_a, test_float_xm_b))) {
      cout << ">>> Testing Float INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float INDENTICAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xm_a << endl
           << " and \n"
           << test_float_xm_a << endl
           << " are identical? " << std::boolalpha << identical(test_float_xm_a, test_float_xm_a) << endl;
      cout << "2) \n" << test_float_xm_a << endl << " and \n" << test_float_xm_b << endl << " are not identical? " << std::boolalpha << (!identical(test_float_xm_a, test_float_xm_b)) << endl;
    }

    if (identical(test_double_xm_a, test_double_xm_a) && (!identical(test_double_xm_a, test_double_xm_b))) {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double INDENTICAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_double_xm_a << endl
           << " and \n"
           << test_double_xm_a << endl
           << " are identical? " << std::boolalpha << identical(test_double_xm_a, test_double_xm_a) << endl;
      cout << "2) \n" << test_double_xm_a << endl << " and \n" << test_double_xm_b << endl << " are not identical? " << std::boolalpha << (!identical(test_double_xm_a, test_double_xm_b)) << endl;
    }

    // template<class utype> bool
    // operator==(const xvector<utype>&,const xvector<utype>&) __xprototype;
    if ((test_int_xm_a == test_int_xm_a) && !(test_int_xm_a == test_int_xm_b)) {
      cout << ">>> Testing Integer == with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer == with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xm_a << endl
           << " and \n"
           << test_int_xm_a << endl
           << " are identical? " << std::boolalpha << identical(test_int_xm_a, test_int_xm_a) << endl;
      cout << "2) \n" << test_int_xm_a << endl << " and \n" << test_int_xm_b << endl << " are not identical? " << std::boolalpha << (!identical(test_int_xm_a, test_int_xm_b)) << endl;
    }

    if ((test_float_xm_a == test_float_xm_a) && (!(test_float_xm_a == test_float_xm_b))) {
      cout << ">>> Testing Float == with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float == with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xm_a << endl
           << " and \n"
           << test_float_xm_a << endl
           << " are identical? " << std::boolalpha << (test_float_xm_a == test_float_xm_a) << endl;
      cout << "2) \n" << test_float_xm_a << endl << " and \n" << test_float_xm_b << endl << " are not identical? " << std::boolalpha << (!(test_float_xm_a == test_float_xm_b)) << endl;
    }

    if (identical(test_double_xm_a, test_double_xm_a) && (!identical(test_double_xm_a, test_double_xm_b))) {
      cout << ">>> Testing Double == with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double == with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_double_xm_a << endl
           << " and \n"
           << test_double_xm_a << endl
           << " are identical? " << std::boolalpha << (test_double_xm_a == test_double_xm_a) << endl;
      cout << "2) \n" << test_double_xm_a << endl << " and \n" << test_double_xm_b << endl << " are not identical? " << std::boolalpha << (!(test_double_xm_a == test_double_xm_b)) << endl;
    }

    // template<class utype> bool
    // isdifferent(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;
    if (!isdifferent(test_int_xm_a, test_int_xm_a, _tol_int) && isdifferent(test_int_xm_a, test_int_xm_b, _tol_int)) {
      cout << ">>> Testing Integer ISDIFFERENT with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISDIFFERENT with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xm_a << endl
           << " and \n"
           << test_int_xm_a << endl
           << " are different? " << std::boolalpha << isdifferent(test_int_xm_a, test_int_xm_a, _tol_int) << endl;
      cout << "2) \n" << test_int_xm_a << endl << " and \n" << test_int_xm_b << endl << " are not different? " << std::boolalpha << (!isdifferent(test_int_xm_a, test_int_xm_b, _tol_int)) << endl;
    }

    if (!isdifferent(test_float_xm_a, test_float_xm_a, _tol_float) && (isdifferent(test_float_xm_a, test_float_xm_b, _tol_float))) {
      cout << ">>> Testing Float ISDIFFERENT with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISDIFFERENT with given tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xm_a << endl
           << " and \n"
           << test_float_xm_a << endl
           << " are different? " << std::boolalpha << isdifferent(test_float_xm_a, test_float_xm_a, _tol_float) << endl;
      cout << "2) \n" << test_float_xm_a << endl << " and \n" << test_float_xm_b << endl << " are not different? " << std::boolalpha << (isdifferent(test_float_xm_a, test_float_xm_b, _tol_float)) << endl;
    }

    if (!isdifferent(test_double_xm_a, test_double_xm_a, _tol_double) && (isdifferent(test_double_xm_a, test_double_xm_b, _tol_double))) {
      cout << ">>> Testing Double ISDIFFERENT with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISDIFFERENT with given tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_double_xm_a << endl
           << " and \n"
           << test_double_xm_a << endl
           << " are different? " << std::boolalpha << isdifferent(test_double_xm_a, test_double_xm_a, _tol_double) << endl;
      cout << "2) \n" << test_double_xm_a << endl << " and \n" << test_double_xm_b << endl << " are not different? " << std::boolalpha << (!isdifferent(test_double_xm_a, test_double_xm_b, _tol_double)) << endl;
    }

    // template<class utype> bool
    // isdifferent(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;
    if (!isdifferent(test_int_xm_a, test_int_xm_a) && isdifferent(test_int_xm_a, test_int_xm_b)) {
      cout << ">>> Testing Integer ISDIFFERENT xvector with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISDIFFERENT xvector with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xm_a << endl
           << " and \n"
           << test_int_xm_a << endl
           << " are different? " << std::boolalpha << isdifferent(test_int_xm_a, test_int_xm_a) << endl;
      cout << "2) \n" << test_int_xm_a << endl << " and \n" << test_int_xm_b << endl << " are not different? " << std::boolalpha << (!isdifferent(test_int_xm_a, test_int_xm_b)) << endl;
    }

    if (!isdifferent(test_float_xm_a, test_float_xm_a) && (isdifferent(test_float_xm_a, test_float_xm_b))) {
      cout << ">>> Testing Float ISDIFFERENT with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISDIFFERENT with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xm_a << endl
           << " and \n"
           << test_float_xm_a << endl
           << " are different? " << std::boolalpha << isdifferent(test_float_xm_a, test_float_xm_a) << endl;
      cout << "2) \n" << test_float_xm_a << endl << " and \n" << test_float_xm_b << endl << " are not different? " << std::boolalpha << (isdifferent(test_float_xm_a, test_float_xm_b)) << endl;
    }

    if (!isdifferent(test_double_xm_a, test_double_xm_a) && (isdifferent(test_double_xm_a, test_double_xm_b))) {
      cout << ">>> Testing Double ISDIFFERENT with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISDIFFERENT with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_double_xm_a << endl
           << " and \n"
           << test_double_xm_a << endl
           << " are different? " << std::boolalpha << isdifferent(test_double_xm_a, test_double_xm_a) << endl;
      cout << "2) \n" << test_double_xm_a << endl << " and \n" << test_double_xm_b << endl << " are not different? " << std::boolalpha << (!isdifferent(test_double_xm_a, test_double_xm_b)) << endl;
    }

    // template<class utype> bool
    // isequal(const xmatrix<utype>&,const xmatrix<utype>&,const utype&) __xprototype;
    if (aurostd::isequal(test_int_xm_a, test_int_xm_a, _tol_int) && (!aurostd::isequal(test_int_xm_a, test_int_xm_b, _tol_int))) {
      cout << ">>> Testing Integer ISEQUAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISEQUAL with given tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xm_a << endl
           << " and \n"
           << test_int_xm_a << endl
           << " are identical? " << std::boolalpha << aurostd::isequal(test_int_xm_a, test_int_xm_a, _tol_int) << endl;
      cout << "2) \n" << test_int_xm_a << endl << " and \n" << test_int_xm_b << endl << " are not identical? " << std::boolalpha << (!aurostd::isequal(test_int_xm_a, test_int_xm_b, _tol_int)) << endl;
    }

    if (aurostd::isequal(test_float_xm_a, test_float_xm_a, _tol_float) && (!aurostd::isequal(test_float_xm_a, test_float_xm_b, _tol_float))) {
      cout << ">>> Testing Float ISEQUAL with given tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISEQUAL with given tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xm_a << endl
           << " and \n"
           << test_float_xm_a << endl
           << " are identical? " << std::boolalpha << aurostd::isequal(test_float_xm_a, test_float_xm_a, _tol_float) << endl;
      cout << "2) \n" << test_float_xm_a << endl << " and \n" << test_float_xm_b << endl << " are not identical? " << std::boolalpha << (!aurostd::isequal(test_float_xm_a, test_float_xm_b, _tol_float)) << endl;
    }

    if (aurostd::isequal(test_double_xm_a, test_double_xm_a, _tol_double) && (!aurostd::isequal(test_double_xm_a, test_double_xm_b, _tol_double))) {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_double_xm_a << endl
           << " and \n"
           << test_double_xm_a << endl
           << " are identical? " << std::boolalpha << aurostd::isequal(test_double_xm_a, test_double_xm_a, _tol_double) << endl;
      cout << "2) \n"
           << test_double_xm_a << endl
           << " and \n"
           << test_double_xm_b << endl
           << " are not identical? " << std::boolalpha << (!aurostd::isequal(test_double_xm_a, test_double_xm_b, _tol_double)) << endl;
    }

    // template<class utype> bool
    // isequal(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;
    if (aurostd::isequal(test_int_xm_a, test_int_xm_a) && (!aurostd::isequal(test_int_xm_a, test_int_xm_b))) {
      cout << ">>> Testing Integer ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer ISEQUAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xm_a << endl
           << " and \n"
           << test_int_xm_a << endl
           << " are identical? " << std::boolalpha << aurostd::isequal(test_int_xm_a, test_int_xm_a) << endl;
      cout << "2) \n" << test_int_xm_a << endl << " and \n" << test_int_xm_b << endl << " are not identical? " << std::boolalpha << (!aurostd::isequal(test_int_xm_a, test_int_xm_b)) << endl;
    }

    if (aurostd::isequal(test_float_xm_a, test_float_xm_a) && (!aurostd::isequal(test_float_xm_a, test_float_xm_b))) {
      cout << ">>> Testing Float ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float ISEQUAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xm_a << endl
           << " and \n"
           << test_float_xm_a << endl
           << " are identical? " << std::boolalpha << aurostd::isequal(test_float_xm_a, test_float_xm_a) << endl;
      cout << "2) \n" << test_float_xm_a << endl << " and \n" << test_float_xm_b << endl << " are not identical? " << std::boolalpha << (!aurostd::isequal(test_float_xm_a, test_float_xm_b)) << endl;
    }

    if (aurostd::isequal(test_double_xm_a, test_double_xm_a) && (!aurostd::isequal(test_double_xm_a, test_double_xm_b))) {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double ISEQUAL with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_double_xm_a << endl
           << " and \n"
           << test_double_xm_a << endl
           << " are identical? " << std::boolalpha << aurostd::isequal(test_double_xm_a, test_double_xm_a) << endl;
      cout << "2) \n" << test_double_xm_a << endl << " and \n" << test_double_xm_b << endl << " are not identical? " << std::boolalpha << (!aurostd::isequal(test_double_xm_a, test_double_xm_b)) << endl;
    }

    // template<class utype> bool
    // operator!=(const xmatrix<utype>&,const xmatrix<utype>&) __xprototype;
    if (!(test_int_xm_a != test_int_xm_a) && ((test_int_xm_a != test_int_xm_b))) {
      cout << ">>> Testing Integer != with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Integer != with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_int_xm_a << endl
           << " and \n"
           << test_int_xm_a << endl
           << " are different? " << std::boolalpha << (test_int_xm_a != test_int_xm_a) << endl;
      cout << "2) \n" << test_int_xm_a << endl << " and \n" << test_int_xm_b << endl << " are not different? " << std::boolalpha << ((test_int_xm_a != test_int_xm_b)) << endl;
    }

    if (!(test_float_xm_a != test_float_xm_a) && ((test_float_xm_a != test_float_xm_b))) {
      cout << ">>> Testing Float != with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Float != with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_float_xm_a << endl
           << " and \n"
           << test_float_xm_a << endl
           << " are different? " << std::boolalpha << (test_float_xm_a != test_float_xm_a) << endl;
      cout << "2) \n" << test_float_xm_a << endl << " and \n" << test_float_xm_b << endl << " are not different? " << std::boolalpha << ((test_float_xm_a != test_float_xm_b)) << endl;
    }

    if (!(test_double_xm_a != test_double_xm_a) && ((test_double_xm_a != test_double_xm_b))) {
      cout << ">>> Testing Double != with default tolerance:  PASS\n";
    } else {
      cout << ">>> Testing Double != with default tolerance:  ***FAIL***!\n"
           << " 1) \n"
           << test_double_xm_a << endl
           << " and \n"
           << test_double_xm_a << endl
           << " are different? " << std::boolalpha << (test_double_xm_a != test_double_xm_a) << endl;
      cout << "2) \n" << test_double_xm_a << endl << " and \n" << test_double_xm_b << endl << " are not different? " << std::boolalpha << ((test_double_xm_a != test_double_xm_b)) << endl;
    }
  }
} // namespace pflow
// ***************************************************************************
// pflow::CHGDIFF
// ***************************************************************************
namespace pflow {
  string CHGDIFF(aurostd::xoption vpflow) {
    // handles flags for CHGDIFF
    string chgcar1_file;
    string chgcar2_file;
    string output_file;
    ostringstream oss;

    // DEBUG
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "BEGIN FLAGS" << endl;
    }

    const string usage_usage = "aflow --chgdiff=CHGCAR1,CHGCAR2 [chgdiff_options]";
    const string usage_options = aurostd::liststring2string("options = --usage", "          --output=|-o=CHGCAR_OUT");

    // output usage
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "CHECK USAGE" << endl;
    }
    if (vpflow.flag("CHGDIFF::USAGE")) {
      init::MessageOption(vpflow.getattachedscheme("CHGDIFF"), "pflow::CHGDIFF()", aurostd::liststring2string(usage_usage, usage_options));
      return oss.str();
    }

    // no input
    vector<string> input_tokens;
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "CHECK INPUT" << endl;
    }
    if (!vpflow.flag("CHGDIFF")) {
      oss << endl;
      oss << __AFLOW_FUNC__ << "ERROR: No input given." << endl;
      oss << __AFLOW_FUNC__ << "Exiting." << endl;
      oss << endl;
    }

    // check inputs
    string misc_option;
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "vpflow.getattachedscheme(\"CHGDIFF\")=" << vpflow.getattachedscheme("CHGDIFF") << endl;
    }
    misc_option = vpflow.getattachedscheme("CHGDIFF");
    aurostd::string2tokens(misc_option, input_tokens, ",");
    if (input_tokens.size() != 2) {
      oss << endl;
      oss << __AFLOW_FUNC__ << "ERROR: Incorrect input arguments. List two CHGCARs separated by commas." << endl;
      oss << __AFLOW_FUNC__ << "Exiting." << endl;
      oss << endl;
    }

    chgcar1_file = input_tokens.at(0);
    chgcar2_file = input_tokens.at(1);

    // get output, if specified (or standardize)
    if (!vpflow.flag("CHGDIFF::OUTPUT")) {
      output_file = "aflow_CHGDIFF.out";
    } else {
      output_file = vpflow.getattachedscheme("CHGDIFF::OUTPUT");
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "CHGCAR_OUT=" << output_file << endl;
    }

    CHGDIFF(chgcar1_file, chgcar2_file, output_file, oss);
    return oss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGDIFF
// ***************************************************************************
namespace pflow {
  bool CHGDIFF(const string& chgcar1_file, const string& chgcar2_file, const string& output_file, ostream& oss) {
    // RETURNS CHGCAR_OUT=CHGCAR_INPUT_1-CHGCAR_INPUT_2
    // Read in CHGCAR or AECCAR files
    // format_dim is as follows: numcolumns chg_tot, npts and numcolumns for
    // augmentation occupancies, numcolumns chg_diff, npts and numcolumns for
    // augmentation occupancies
    // read chgcars

    stringstream chgcar1_ss;
    stringstream chgcar2_ss;
    stringstream chgcar1_header;
    stringstream chgcar2_header;
    const double TOL = 1e-5;
    xstructure structure1;
    xstructure structure2;
    vector<int> ngrid1(3);
    vector<int> ngrid2(3);
    vector<int> format_dim1;
    vector<int> format_dim2;
    vector<double> chg_tot1;
    vector<double> chg_tot2;
    vector<double> chg_diff1;
    vector<double> chg_diff2;

    // DEBUG
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "BEGIN FUNCTION" << endl;
    }

    // get input 1
    if (!aurostd::FileExist(chgcar1_file)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << "ERROR: " << chgcar1_file << " does not exist." << endl;
      oss << __AFLOW_FUNC__ << "Exiting." << endl;
      oss << endl;
      return false;
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "CHGCAR1=" << chgcar1_file << endl;
    }
    aurostd::compressfile2stringstream(chgcar1_file, chgcar1_ss);

    // get input 2
    if (!aurostd::FileExist(chgcar2_file)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << "ERROR: " << chgcar2_file << " does not exist." << endl;
      oss << __AFLOW_FUNC__ << "Exiting." << endl;
      oss << endl;
      return false;
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "CHGCAR2=" << chgcar2_file << endl;
    }
    aurostd::compressfile2stringstream(chgcar2_file, chgcar2_ss);

    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "CHECK FORMAT OF " << chgcar1_file << endl;
    }
    if (!pflow::ReadCHGCAR(structure1, chgcar1_header, ngrid1, format_dim1, chg_tot1, chg_diff1, chgcar1_ss, oss)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << "ERROR: Input " << chgcar1_file << " format not recognized." << endl;
      oss << endl;
      return false; // CHGCAR format not recognized
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "SUCCESSFULLY READ " << chgcar1_file << endl;
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "CHECK FORMAT OF " << chgcar2_file << endl;
    }
    if (!pflow::ReadCHGCAR(structure2, chgcar2_header, ngrid2, format_dim2, chg_tot2, chg_diff2, chgcar2_ss, oss)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << "ERROR: Input " << chgcar2_file << " format not recognized." << endl;
      oss << endl;
      return false; // CHGCAR format not recognized
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "SUCCESSFULLY READ " << chgcar2_file << endl;
    }

    // check formats
    aurostd::matrix<double> lat1 = pflow::GetLat(structure1); // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> lat2 = pflow::GetLat(structure2); // CO20200404 pflow::matrix()->aurostd::matrix()
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "CHECK IF FORMATS OF CHGCARS MATCH" << endl;
    }
    if (format_dim1 != format_dim2) {
      oss << endl;
      oss << __AFLOW_FUNC__ << "WARNING: Format for input " << chgcar1_file << " does not match that of input " << chgcar2_file << "." << endl;
      oss << __AFLOW_FUNC__ << "WARNING: Using format from input " << chgcar1_file << "." << endl;
      oss << endl;
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "CHECK IF GRIDS OF CHGCARS MATCH" << endl;
    }
    if (ngrid1 != ngrid2) {
      oss << endl;
      oss << __AFLOW_FUNC__ << "ERROR: Grids for two CHGCAR's do not match. " << endl;
      oss << __AFLOW_FUNC__ << "ERROR: ngrid of " << chgcar1_file << ": " << ngrid1.at(0) << " " << ngrid1.at(1) << " " << ngrid1.at(2) << endl;
      oss << __AFLOW_FUNC__ << "ERROR: ngrid of " << chgcar2_file << ": " << ngrid2.at(0) << " " << ngrid2.at(1) << " " << ngrid2.at(2) << endl;
      oss << __AFLOW_FUNC__ << "ERROR: This will give nonsense. " << endl;
      oss << endl;
      return false;
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "CHECK IF LATTICE PARAMETERS OF CHGCARS MATCH" << endl;
    }
    if (pflow::norm(pflow::VVdiff(lat1[0], lat2[0])) > TOL || pflow::norm(pflow::VVdiff(lat1[1], lat2[1])) > TOL || pflow::norm(pflow::VVdiff(lat1[2], lat2[2])) > TOL) {
      oss << endl;
      oss << __AFLOW_FUNC__ << "WARNING: Lattice parameters for two CHGCARs do not match. " << endl;
      oss << __AFLOW_FUNC__ << "WARNING: lattice of " << chgcar1_file << ": " << endl;
      oss << __AFLOW_FUNC__ << matrix2xmatrix(lat1) << endl;
      oss << __AFLOW_FUNC__ << "WARNING: lattice of " << chgcar2_file << ": " << endl;
      oss << matrix2xmatrix(lat2) << endl;
      oss << __AFLOW_FUNC__ << "WARNING: This could give nonsense if there is much difference. " << endl;
      oss << __AFLOW_FUNC__ << "WARNING: Output will use lattice of " << chgcar1_file << "." << endl;
      oss << endl;
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "OVERALL FORMAT OF CHGCARS LOOKS OK" << endl;
    }

    // Get difference
    const vector<double> chg_tot_2m1 = pflow::VVdiff(chg_tot1, chg_tot2);
    const vector<double> chg_diff_2m1 = pflow::VVdiff(chg_diff1, chg_diff2);
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "PRINTING CHGDIFF TO " << output_file << endl;
    }
    pflow::PrintCHGCAR(structure1, chgcar1_header, ngrid1, format_dim1, chg_tot_2m1, chg_diff_2m1, output_file, oss);
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "DONE" << endl;
    }
    oss << __AFLOW_FUNC__ << output_file << " generated." << endl;
    return true;
  }
} // namespace pflow

// DX+CO START
//  ***************************************************************************
//  pflow::CHGINT
//  ***************************************************************************
namespace pflow {
  void CHGINT(vector<string> argv) {
    ifstream chgfile(argv.at(2).c_str());
    xstructure str;
    vector<int> ngrid(3);
    // Read in charge
    vector<double> chg_tot;
    vector<double> chg_diff;
    pflow::ReadChg(str, ngrid, chg_tot, chg_diff, chgfile);
    // Integrate charge
    vector<aurostd::matrix<double>> rad_chg_int; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> vor_chg_int; // CO20200404 pflow::matrix()->aurostd::matrix()
    pflow::GetChgInt(rad_chg_int, vor_chg_int, str, ngrid, chg_tot, chg_diff);
    // Print results
    pflow::PrintChgInt(rad_chg_int, vor_chg_int, cout);
  }
} // namespace pflow
// DX+CO END

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  string CHGSUM(aurostd::xoption vpflow) {
    // handles flags for CHGSUM

    vector<string> chgcar_files;
    ostringstream oss;
    string output_file;

    // DEBUG
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "BEGIN FLAGS" << endl;
    }

    const string usage_usage = "aflow --chgsum=CHGCAR1,CHGCAR2,... [chgsum_options]";
    const string usage_options = aurostd::liststring2string("options = --usage", "          --output=|-o=CHGCAR_OUT");

    // output usage
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "CHECK USAGE" << endl;
    }
    if (vpflow.flag("CHGSUM::USAGE")) {
      init::MessageOption(vpflow.getattachedscheme("CHGSUM"), "pflow::CHGSUM", aurostd::liststring2string(usage_usage, usage_options));
      return oss.str();
    }

    // no input
    vector<string> input_tokens;
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "CHECK INPUT" << endl;
    }
    if (!vpflow.flag("CHGSUM")) {
      oss << endl;
      oss << __AFLOW_FUNC__ << "ERROR: No input given." << endl;
      oss << __AFLOW_FUNC__ << "Exiting." << endl;
      oss << endl;
      oss << endl;
      return oss.str();
    }

    // check inputs
    string misc_option;
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "vpflow.getattachedscheme(\"CHGSUM\")=" << vpflow.getattachedscheme("CHGSUM") << endl;
    }
    misc_option = vpflow.getattachedscheme("CHGSUM");
    aurostd::string2tokens(misc_option, input_tokens, ",");
    if (input_tokens.size() < 2) {
      oss << endl;
      oss << __AFLOW_FUNC__ << "ERROR: Incorrect input arguments. List at least two CHGCARs separated by commas." << endl;
      oss << __AFLOW_FUNC__ << "Exiting." << endl;
      oss << endl;
      return oss.str();
    } else {
      // get inputs
      for (size_t i = 0; i < input_tokens.size(); i++) {
        chgcar_files.push_back(input_tokens[i]);
      }
    }

    // get output, if specified (or standardize)
    if (!vpflow.flag("CHGSUM::OUTPUT")) {
      output_file = "aflow_CHGSUM.out";
    } else {
      output_file = vpflow.getattachedscheme("CHGSUM::OUTPUT");
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "CHGCAR_OUT=" << output_file << endl;
    }
    CHGSUM(chgcar_files, output_file, oss);
    return oss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  bool CHGSUM(const string& chgcar_in1, const string& chgcar_in2, ostream& oss) {
    // 2 INPUTS, NO OUTPUT
    const string output_file = "aflow_chgsum.out";
    return CHGSUM(chgcar_in1, chgcar_in2, output_file, oss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  bool CHGSUM(string& species_header, const string& chgcar_in1, const string& chgcar_in2, const string& output_file, ostream& oss) {
    // 2 INPUTS WITH SPECIES_HEADER
    vector<string> chgcar_files;
    chgcar_files.push_back(chgcar_in1);
    chgcar_files.push_back(chgcar_in2);
    return CHGSUM(species_header, chgcar_files, output_file, oss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  bool CHGSUM(const string& chgcar_in1, const string& chgcar_in2, const string& output_file, ostream& oss) {
    // 2 INPUTS, NO SPECIES_HEADER
    string species_header;
    return CHGSUM(species_header, chgcar_in1, chgcar_in2, output_file, oss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  bool CHGSUM(const vector<string>& chgcar_files, ostream& oss) {
    // VECTOR INPUT, NO SPECIES_HEADER OR OUTPUT
    string species_header;
    return CHGSUM(species_header, chgcar_files, oss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  bool CHGSUM(const vector<string>& chgcar_files, const string& output_file, ostream& oss) {
    // VECTOR INPUT, NO SPECIES_HEADER
    string species_header;
    return CHGSUM(species_header, chgcar_files, output_file, oss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  bool CHGSUM(string& species_header, const vector<string>& chgcar_files, ostream& oss) {
    // VECTOR INPUT, NO OUTPUT
    const string output_file = "aflow_chgsum.out";
    return CHGSUM(species_header, chgcar_files, output_file, oss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::CHGSUM
// ***************************************************************************
namespace pflow {
  bool CHGSUM(string& species_header, const vector<string>& chgcar_files, const string& output_file, ostream& oss) {
    // RETURNS CHGCAR_OUT=\sum CHGCAR_INPUTs
    // Read in CHGCAR or AECCAR files
    // format_dim is as follows: numcolumns chg_tot, npts and numcolumns for
    // augmentation occupancies, numcolumns chg_diff, npts and numcolumns for
    // augmentation occupancies

    const double TOL = 1e-5;
    xstructure structure1;
    xstructure structure2;
    stringstream chgcar_ss;
    stringstream chgcar1_header;
    stringstream chgcar2_header;

    // DEBUG
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "BEGIN FUNCTION" << endl;
    }

    for (size_t i = 0; i < chgcar_files.size(); i++) {
      if (!aurostd::FileExist(chgcar_files[i])) {
        oss << endl;
        oss << __AFLOW_FUNC__ << "ERROR: " << chgcar_files[i] << " does not exist." << endl;
        oss << __AFLOW_FUNC__ << "Exiting." << endl;
        oss << endl;
        return false;
      }
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << "CHGCAR_IN_" << i + 1 << "=" << chgcar_files[i] << endl;
      }
    }

    // read first chgcar
    vector<int> ngrid1(3);
    vector<int> ngrid2(3);
    vector<int> format_dim1;
    vector<int> format_dim2;
    vector<double> chg_tot1;
    vector<double> chg_tot2;
    vector<double> chg_diff1;
    vector<double> chg_diff2;
    aurostd::matrix<double> lat1;
    aurostd::matrix<double> lat2; // CO20200404 pflow::matrix()->aurostd::matrix()

    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "CHECK " << chgcar_files.at(0) << endl;
    }
    aurostd::compressfile2stringstream(chgcar_files.at(0), chgcar_ss);
    if (!pflow::ReadCHGCAR(structure1, chgcar1_header, ngrid1, format_dim1, chg_tot1, chg_diff1, chgcar_ss, oss)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << "ERROR: Input " << chgcar_files.at(0) << " format not recognized." << endl;
      oss << __AFLOW_FUNC__ << "Exiting." << endl;
      oss << endl;
      return false; // CHGCAR format not recognized
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "SUCCESSFULLY READ " << chgcar_files.at(0) << endl;
    }

    // for later checks
    lat1 = pflow::GetLat(structure1);

    // scroll through other structures
    for (size_t i = 1; i < chgcar_files.size(); i++) {
      chgcar_ss.str("");
      //            structure2.~xstructure();
      chgcar2_header.str("");
      ngrid2.clear();
      format_dim2.clear();
      chg_tot2.clear();
      chg_diff2.clear();

      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << "CHECK " << chgcar_files.at(i) << endl;
      }
      aurostd::compressfile2stringstream(chgcar_files.at(i), chgcar_ss);
      if (!pflow::ReadCHGCAR(structure2, chgcar2_header, ngrid2, format_dim2, chg_tot2, chg_diff2, chgcar_ss, oss)) {
        oss << endl;
        oss << __AFLOW_FUNC__ << "ERROR: Input " << chgcar_files.at(i) << " format not recognized." << endl;
        oss << __AFLOW_FUNC__ << "Exiting." << endl;
        oss << endl;
        return false; // CHGCAR format not recognized
      }
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << "SUCCESSFULLY READ CHGCAR_INPUT_" << i << endl;
      }
      lat2 = pflow::GetLat(structure2);

      // check formats
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << "CHECK IF FORMATS OF " << chgcar_files.at(0) << " and " << chgcar_files.at(i) << " MATCH" << endl;
      }
      if (format_dim1 != format_dim2) {
        oss << endl;
        oss << __AFLOW_FUNC__ << "WARNING: Format for " << chgcar_files.at(0) << " does not match that of " << chgcar_files.at(i) << "." << endl;
        oss << __AFLOW_FUNC__ << "WARNING: Using format from " << chgcar_files.at(0) << "." << endl;
        oss << endl;
      }
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << "CHECK IF GRIDS OF " << chgcar_files.at(0) << " and " << chgcar_files.at(i) << " MATCH" << endl;
      }
      if (ngrid1 != ngrid2) {
        oss << endl;
        oss << __AFLOW_FUNC__ << "ERROR: Grid for " << chgcar_files.at(0) << " does not match that of " << chgcar_files.at(i) << "." << endl;
        oss << __AFLOW_FUNC__ << "ERROR: ngrid of " << chgcar_files.at(0) << ": " << ngrid1.at(0) << " " << ngrid1.at(1) << " " << ngrid1.at(2) << endl;
        oss << __AFLOW_FUNC__ << "ERROR: ngrid of " << chgcar_files.at(i) << ": " << ngrid2.at(0) << " " << ngrid2.at(1) << " " << ngrid2.at(2) << endl;
        oss << __AFLOW_FUNC__ << "ERROR: This will give nonsense." << endl;
        oss << __AFLOW_FUNC__ << "Exiting." << endl;
        oss << endl;
        return false;
      }
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << "CHECK IF LATTICE PARAMETERS OF " << chgcar_files.at(0) << " and " << chgcar_files.at(i) << " MATCH" << endl;
      }
      if (pflow::norm(pflow::VVdiff(lat1[0], lat2[0])) > TOL || pflow::norm(pflow::VVdiff(lat1[1], lat2[1])) > TOL || pflow::norm(pflow::VVdiff(lat1[2], lat2[2])) > TOL) {
        oss << endl;
        oss << __AFLOW_FUNC__ << "WARNING: Lattice parameters for " << chgcar_files.at(0) << " and " << chgcar_files.at(i) << " do not match. " << endl;
        oss << __AFLOW_FUNC__ << "WARNING: lattice of " << chgcar_files.at(0) << ": " << endl;
        oss << __AFLOW_FUNC__ << matrix2xmatrix(lat1) << endl;
        oss << __AFLOW_FUNC__ << "WARNING: lattice " << chgcar_files.at(i) << ": " << endl;
        oss << matrix2xmatrix(lat2) << endl;
        oss << __AFLOW_FUNC__ << "WARNING: This could give nonsense if there is much difference." << endl;
        oss << __AFLOW_FUNC__ << "WARNING: Output will use lattice of " << chgcar_files.at(0) << "." << endl;
        oss << endl;
      }
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << "OVERALL FORMAT OF " << chgcar_files.at(i) << " LOOKS OK" << endl;
      }

      // Get sum
      chg_tot1 = pflow::VVsum(chg_tot1, chg_tot2);
      chg_diff1 = pflow::VVsum(chg_diff1, chg_diff2);
    }

    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "PRINTING CHGSUM TO " << output_file << endl;
    }

    // EDIT CHGCAR_HEADER1 FORMATTING FOR BADER
    if (!species_header.empty()) {
      bader_functions::adjust_header(species_header, chgcar1_header);
    } // CO20180627

    // print chgcar
    pflow::PrintCHGCAR(structure1, chgcar1_header, ngrid1, format_dim1, chg_tot1, chg_diff1, output_file, oss);
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << "DONE" << endl;
    }
    oss << __AFLOW_FUNC__ << output_file << " generated." << endl;
    return true;
  }
} // namespace pflow

// ***************************************************************************
// pflow::CIF
// ***************************************************************************
namespace pflow {
  void CIF(istream& input, aurostd::xoption& vpflow) {
    if (vpflow.flag("CIF::USAGE")) {
      init::MessageOption("--usage", "pflow::CIF",
                          aurostd::liststring2string("aflow --cif[=<tolerance_value>|=tight|=loose] [--setting=1|2|aflow] [--no_symmetry] < POSCAR  default: tolerance=(minimum_interatomic_distance)/100.0, "
                                                     "setting=1"));
      return;
    }
    xstructure a(input, IOAFLOW_AUTO);
    a.ReScale(1.0); // DX20190123 - rescale to 1.0
    bool no_symmetry = false;
    if (vpflow.flag("CIF::NO_SYMMETRY")) {
      no_symmetry = true;
    }
    // DX20180803 - consider space group - pflow::PrintCIF(cout,a);
    if (no_symmetry) {
      pflow::PrintCIF(cout, a, 1); // DX20180803 - add space group information
    } else {
      // put tolerance flag check in loop to save time if we aren't calculating, defaultTolerance can be expensive
      double tolerance = pflow::getSymmetryTolerance(a, vpflow.getattachedscheme("CIF::TOLERANCE")); // DX20200820 - consolidated setting tolerance into a function
      const uint setting = pflow::getSpaceGroupSetting(vpflow.getattachedscheme("CIF::SETTING")); // DX20210421 - consolidated space group setting into function
      if (vpflow.flag("DATA::NO_SCAN")) {
        a.sym_eps_no_scan = true;
      } // DX20210611
      a.spacegroupnumber = a.SpaceGroup_ITC(tolerance, -1, setting, a.sym_eps_no_scan); // DX20210611 - added no_scan
      a.lattice = a.standard_lattice_ITC; // DX20180904 - need to update the lattice; may have rotated
      pflow::PrintCIF(cout, a, a.spacegroupnumber, setting); // DX20180803 - add space group information
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::getCIFSettings
// ***************************************************************************
namespace pflow {
  /// @brief Extract command-line options to generate a CIF
  ///
  /// @param a      Object containing crystal structure information
  /// @param vpflow Object containing command line options
  void getCIFOptions(xstructure& a, aurostd::xoption& vpflow) {
    if (vpflow.flag("CIF::USAGE")) {
      init::MessageOption("--usage", "pflow::CIF",
                          aurostd::liststring2string("aflow --cif[=<tolerance_value>|=tight|=loose] [--setting=1|2|aflow] [--no_symmetry] < POSCAR  default: tolerance=(minimum_interatomic_distance)/100.0, "
                                                     "setting=1"));
      return;
    }
    if (vpflow.flag("CIF::NO_SYMMETRY")) {
      a.spacegroupnumber = 1; // force to CIF to print as P1 structure
    }
    // put tolerance flag check in loop to save time if we aren't calculating, defaultTolerance can be expensive
    a.sym_eps = pflow::getSymmetryTolerance(a, vpflow.getattachedscheme("CIF::TOLERANCE"));
    a.setting_ITC = pflow::getSpaceGroupSetting(vpflow.getattachedscheme("CIF::SETTING"));
    if (vpflow.flag("DATA::NO_SCAN")) {
      a.sym_eps_no_scan = true;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::CLAT
// ***************************************************************************
namespace pflow {
  void CLAT(const string& options) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 6) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --clat=a,b,c,alpha,beta,gamma");
    }
    const xvector<double> data(6);
    if (!tokens.empty()) {
      data[1] = aurostd::string2utype<double>(tokens[0]);
    }
    if (tokens.size() >= 2) {
      data[2] = aurostd::string2utype<double>(tokens[1]);
    }
    if (tokens.size() >= 3) {
      data[3] = aurostd::string2utype<double>(tokens[2]);
    }
    if (tokens.size() >= 4) {
      data[1] = aurostd::string2utype<double>(tokens[3]);
    }
    if (tokens.size() >= 5) {
      data[2] = aurostd::string2utype<double>(tokens[4]);
    }
    if (tokens.size() >= 6) {
      data[3] = aurostd::string2utype<double>(tokens[5]);
    }
    cout << aflow::Banner("BANNER_TINY") << endl;
    pflow::PrintClat(data, cout);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::CLEANALL
// ***************************************************************************
namespace pflow {
  void CLEANALL(istream& input) {
    vector<string> vinput;
    aurostd::stream2vectorstring(input, vinput);
    for (size_t iinput = 0; iinput < vinput.size(); iinput++) {
      cout << vinput[iinput] << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::CMPSTR
// ***************************************************************************
namespace pflow {
  void CMPSTR(vector<string> argv) {
    cout << aflow::Banner("BANNER_TINY") << endl;
    // Read in input file.
    ifstream infile1(argv.at(2).c_str());
    ifstream infile2(argv.at(3).c_str());
    const double rcut = atof(argv.at(4).c_str());
    xstructure str1;
    xstructure str2;
    infile1 >> str1;
    infile2 >> str2;
    pflow::PrintCmpStr(str1, str2, rcut, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::COMPARE
// ***************************************************************************
namespace pflow {
  void COMPARE(const string& options) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 12) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --compare=a,b,c,d,e,f,g,h,k,j,i,l");
    }
    const xvector<double> aa(12);
    for (int i = 1; i <= 12; i++) {
      aa(i) = aurostd::string2utype<double>(tokens.at(i - 1));
    }
    double bb1 = 0.0;
    double bb2 = 0.0;
    double bb3 = 0.0;
    double bb4 = 0.0;
    double bb5 = 0.0;
    double bb6 = 0.0;
    cout << "PERCENTAGES" << "  ";
    cout << (bb1 = 100 * std::abs(aa(1) - aa(7)) / ((aa(1) + aa(7)) / 2.0)) << "  ";
    cout << (bb2 = 100 * std::abs(aa(2) - aa(8)) / ((aa(2) + aa(8)) / 2.0)) << "  ";
    cout << (bb3 = 100 * std::abs(aa(3) - aa(9)) / ((aa(3) + aa(9)) / 2.0)) << "  ";
    cout << (bb4 = 100 * std::abs(aa(4) - aa(10)) / ((aa(4) + aa(10)) / 2.0)) << "  ";
    cout << (bb5 = 100 * std::abs(aa(5) - aa(11)) / ((aa(5) + aa(11)) / 2.0)) << "  ";
    cout << (bb6 = 100 * std::abs(aa(6) - aa(12)) / ((aa(6) + aa(12)) / 2.0)) << "  ";
    cout << endl;
    cout << "MIN " << min(bb1, min(bb2, min(bb3, min(bb4, min(bb5, bb6))))) << endl;
    cout << "MAX " << max(bb1, max(bb2, max(bb3, max(bb4, max(bb5, bb6))))) << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::CORNERS
// ***************************************************************************
// Stefano Curtarolo (Dec-2009)
namespace pflow {
  xstructure CORNERS(istream& input) {
    xstructure a(input, IOAFLOW_AUTO);
    a.AddCorners();
    return a;
  }
} // namespace pflow

// ***************************************************************************
// pflow::CPRIM
// ***************************************************************************
// Stefano Curtarolo (Dec-2009)
namespace pflow {
  xstructure CPRIM(istream& input) {
    cerr << XPID << "pflow::CPRIM: THIS IS A DEBUG FUNCTION FOR CODING PURPOSES" << endl;
    const xstructure str_in(input, IOAFLOW_AUTO);
    xstructure str_sp;
    xstructure str_sc;
    // str_in.SetCoordinates(_COORDS_CARTESIAN_);
    const bool full_sym = false; // DX20170829 - Speed increase
    LATTICE::Standard_Lattice_Structure(str_in, str_sp, str_sc, full_sym); // DX20170829 - Speed increase
    cerr << "ORIGINAL" << endl;
    cerr << Getabc_angles(str_in.lattice, DEGREES) << endl;
    cerr << str_in;
    cerr << "STANDARD_PRIMITIVE" << endl;
    cerr << Getabc_angles(str_sp.lattice, DEGREES) << endl;
    cerr << str_sp;
    cerr << "STANDARD_CONVENTIONAL" << endl;
    cerr << Getabc_angles(str_sc.lattice, DEGREES) << endl;
    cerr << str_sc;
    return str_sp;
  }
} // namespace pflow

// ***************************************************************************
// pflow::DATA()
// ***************************************************************************
// Determines the crystallographic data (lattice data or extended)
// Stefano Curtarolo
// Updated by David Hicks (DX) //DX20210302
// Added separated real, reciprocal, and superlattice analyses
// Moved SGDATA() inside
namespace pflow {
  bool DATA(istream& input, const aurostd::xoption& vpflow, const string& smode, ostream& oss) {
    const bool LDEBUG = (false || XHOST.DEBUG);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN: smode=" << smode << endl;
    }

    // ---------------------------------------------------------------------------
    // usage
    if (vpflow.flag("DATA::USAGE")) {
      // main command
      string usage;
      if (smode == "DATA") {
        usage = "aflow --data";
      } else if (smode == "CRYSTAL_POINT_GROUP") {
        usage = "aflow --point_group_crystal_data|--pointgroup_crystal_data|--pgroupxtal_data|--pgroup_xtal_data[=<tolerance_value>|=tight|=loose]";
      } else if (smode == "EDATA") {
        usage = "aflow --edata[=<tolerance_value>|=tight|=loose]";
      } else if (smode == "REAL_LATTICE") {
        usage = "aflow --lattice_data|--data_lattice|--real_lattice_data|--data_real_lattice[=<tolerance_value>|=tight|=loose]";
      } else if (smode == "RECIPROCAL_LATTICE") {
        usage = "aflow --reciprocal_lattice_data|--reciprocallattice_data|--klattice_data|--data_reciprocal_lattice[=<tolerance_value>|=tight|=loose]";
      } else if (smode == "SGDATA") {
        usage = "aflow --sgdata|--space_group_data[=<tolerance_value>|=tight|=loose]";
      } else if (smode == "SUPERLATTICE") {
        usage = "aflow --superlattice_data|--data_superlattice[=<tolerance_value>|=tight|=loose]";
      }

      // options
      string options = "options: [--no_scan] [--print=txt|json]";
      if (smode == "EDATA" || smode == "SGDATA") {
        options += " [--setting=1|2|aflow] [--suppress_Wyckoff|--suppress_Wyckoff_printing|--suppress_wyckoff|--suppress_wyckoff_printing] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]]";
      }
      vector<string> voptions;
      aurostd::string2tokens(options, voptions, " ");
      voptions.insert(voptions.begin(), usage);

      init::MessageOption("--usage", "pflow::DATA()", voptions);
      return true; // CO20200624 - the option was expressed successfully
    }

    // ---------------------------------------------------------------------------
    // load structure
    xstructure a(input, IOAFLOW_AUTO);
    a.ReScale(1.0);

    // ---------------------------------------------------------------------------
    // get directory info
    if (a.directory.empty()) {
      a.directory = aurostd::getPWD();
    }

    // ---------------------------------------------------------------------------
    // get magnetic moment information
    if (vpflow.flag("DATA::MAGNETIC")) {
      const string magmom_info = vpflow.getattachedscheme("DATA::MAGNETIC");
      ProcessAndAddSpinToXstructure(a, magmom_info); // DX20191108 - condensed into a single function
    }

    // ---------------------------------------------------------------------------
    // get tolerance
    const double tolerance = pflow::getSymmetryTolerance(a, vpflow.getattachedscheme("DATA::TOLERANCE")); // DX20200820 - consolidated setting tolerance into a function

    // ---------------------------------------------------------------------------
    // get space group setting
    const uint setting = pflow::getSpaceGroupSetting(vpflow.getattachedscheme("DATA::SETTING")); // DX20210421 - consolidated space group setting into function

    // ---------------------------------------------------------------------------
    // file type //DX20210226 - string to filetype
    filetype ftype = txt_ft;
    if (XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      ftype = txt_ft;
    } else if (XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
      ftype = json_ft;
    }

    // ---------------------------------------------------------------------------
    // self-consistent tolerance scan
    if (vpflow.flag("DATA::NO_SCAN")) {
      a.sym_eps_no_scan = true;
    } // DX20210406

    // ---------------------------------------------------------------------------
    // SGDATA: do not print all Wyckoff data (useful for web) //DX20210301
    const bool suppress_Wyckoff = (vpflow.flag("DATA::SUPPRESS_WYCKOFF_PRINTING") || XHOST.vflag_control.flag("WWW"));

    // ---------------------------------------------------------------------------
    // add banner // ME20210402 - but not for web
    if (ftype == txt_ft && !XHOST.vflag_control.flag("WWW")) {
      cout << aflow::Banner("BANNER_TINY") << endl;
    }

    // ---------------------------------------------------------------------------
    // perform relevant analysis
    const bool already_calculated = false;
    const bool standalone = true;
    if (smode == "EDATA" || smode == "DATA") {
      oss << PrintData(a, smode, ftype, already_calculated, tolerance, a.sym_eps_no_scan, setting);
    } else if (smode == "REAL_LATTICE") {
      oss << PrintRealLatticeData(a, "EDATA", ftype, standalone, already_calculated, tolerance);
    } else if (smode == "CRYSTAL_POINT_GROUP") {
      oss << PrintCrystalPointGroupData(a, ftype, standalone, already_calculated, tolerance);
    } else if (smode == "SUPERLATTICE") {
      oss << PrintSuperlatticeData(a, ftype, standalone, already_calculated, tolerance);
    } else if (smode == "RECIPROCAL_LATTICE") {
      oss << PrintReciprocalLatticeData(a, ftype, standalone, already_calculated, tolerance);
    } else if (smode == "SGDATA") {
      oss << pflow::PrintSGData(a, ftype, standalone, already_calculated, tolerance, a.sym_eps_no_scan, setting, suppress_Wyckoff);
    }

    return true;
  }
} // namespace pflow

// ***************************************************************************
// pflow::DISP
// ***************************************************************************
namespace pflow {
  void DISP(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 1) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --disp=cutoff < POSCAR");
    }
    double cutoff = 0.0;
    if (!tokens.empty()) {
      cutoff = aurostd::string2utype<double>(tokens[0]);
    }
    const xstructure a(input, IOAFLOW_AUTO);
    cout << aflow::Banner("BANNER_TINY") << endl;
    pflow::PrintDisplacements(a, cutoff, cout);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::DIST
// ***************************************************************************
namespace pflow {
  void DIST(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 1) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --dist=cutoff < POSCAR");
    }
    double cutoff = 0.0;
    if (!tokens.empty()) {
      cutoff = aurostd::string2utype<double>(tokens[0]);
    }
    const xstructure a(input, IOAFLOW_AUTO);
    cout << aflow::Banner("BANNER_TINY") << endl;
    pflow::PrintDistances(a, cutoff, cout);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// // ***************************************************************************
// // pflow::DYNADIEL // CAMILO
// // ***************************************************************************
// namespace pflow {
//   void DYNADIEL(vector<string>& argv) {  // loop pflow::DYNADIEL
//     // modify this for the spectrum
//     string outcar ;
//     xvector<double> real(3), imag(3) ;
//
//     if(argv.size() != 3)
//       {  // user control lines - aflow specific.
//    init::ErrorOption("","pflow::DYNADIEL","aflow --dynadiel OUTCAR*");
//       }
//     outcar = argv.at(2) ;
//     KBIN::GetDynaDiel(outcar, real, imag) ;
//   }  // loop pflow::DYNADIEL
// }

// ***************************************************************************
// pflow::EDOS
// ***************************************************************************
namespace pflow {
  void EDOS(vector<string> argv) {
    // 2008 wahyu setyawan    ws26@duke.edu
    // 2008 roman chepulskyy  rc74@duke.edu
#define _MAXSPECIES_ 10
    //   cout << " udos -h            : help\n"
    //   << " udos 1 s p         : nonspin, includes s and p orbitals\n"
    //   << " udos 2 s p d f     : spin-polarized, includes s,p,d,and f orbitals.\n\n"
    //   << " OUTPUT format:\n"
    //   << " udos 1\n"
    //   << "   energy DOS\n"
    //   << " udos 2\n"
    //   << "   energy DOSup DOSdown\n"
    //   << " udos 2 f\n"
    //   << "   energy fDOSup_elemnt1 fDOSdown_elemnt1 ... fDOSup_elmntN fDOSdown_elmntN\n"
    //   << " udos 2 s d\n"
    //   << "   energy sDOSup_elmnt1 sDOSdown_elmnt1 ... sDOSup_elmntN sDOSdown_elmntN dDOSup_elmnt1 dDOSdown_elmnt1 .. dDOSup_elmntN dDOSdown_elmntN"
    //   << " udos 2 d s\n"
    //   << "   energy sDOSup_elmnt1 sDOSdown_elmnt1 ... sDOSup_elmntN sDOSdown_elmntN dDOSup_elmnt1 dDOSdown_elmnt1 .. dDOSup_elmntN dDOSdown_elmntN"
    //   << " udos 2 s p d f\n"
    //   << "   energy DOSup DOSdown sDOSup_elmnt1 sDOSdown_elmnt1 ...sDOSup_elmntN sDOSdown_elmntN ... fDOSup_elmntN fDOSdown_elmntN\n\n"
    //   << " "
    //   << " note: DOS for spin down is given in (negative) sign.\n"
    //   << "       Splitting of the elements(or species) is according to POSCAR.\n";

    const int argc = argv.size();
    if (argc == 2) {
      return;
    }
    if (!(atoi(&argv.at(2)[0]) == 1 or atoi(&argv.at(2)[0]) == 2)) {
      return;
    }

    bool DO_S;
    bool DO_P;
    bool DO_D;
    bool DO_F;
    int i;
    int j;
    int itmp;
    int ispin;
    int natom;
    int nE;
    int ncol;
    int k;
    int ioffset;
    int iorb;
    int maxOrbital = 0;
    int nspec;
    int species[_MAXSPECIES_ + 1];
    float minE;
    float maxE;
    float Efermi;
    float ftmp;
    string tmpstr;
    ifstream pin;
    ifstream din;

    ispin = atoi(&argv.at(2)[0]);
    if (!(ispin == 1 or ispin == 2)) {
      return;
    }

    DO_S = false;
    DO_P = false;
    DO_D = false;
    DO_F = false;

    if (argc > 3) {
      for (i = 3; i < argc; i++) {
        if (argv.at(i)[0] == 's' or argv.at(i)[0] == 'S') {
          DO_S = true;
        }
        if (argv.at(i)[0] == 'p' or argv.at(i)[0] == 'P') {
          DO_P = true;
        }
        if (argv.at(i)[0] == 'd' or argv.at(i)[0] == 'D') {
          DO_D = true;
        }
        if (argv.at(i)[0] == 'f' or argv.at(i)[0] == 'F') {
          DO_F = true;
        }
      }
    }
    // getting number of each species from POSCAR to accumulate s,p,d,f atomic DOS
    pin.open("POSCAR");
    for (i = 1; i <= 5; i++) {
      getline(pin, tmpstr);
    }
    i = 0;
    j = 0;
    pin >> tmpstr;
    do {
      species[++i] = atoi(&tmpstr[0]);
      pin >> tmpstr;
    } while (atoi(&tmpstr[0]) > 0);
    nspec = i;
    pin.close();
    // processing DOSCAR
    din.open("DOSCAR");
    din >> natom;
    itmp = 0;
    for (i = 1; i <= nspec; i++) {
      itmp += species[i];
    }
    if (itmp != natom) {
      cerr << "DOSCAR is INcompatible with POSCAR\n";
      return;
    }
    getline(din, tmpstr);
    getline(din, tmpstr);
    getline(din, tmpstr);
    getline(din, tmpstr);
    getline(din, tmpstr);
    din >> maxE >> minE >> nE >> Efermi;
    getline(din, tmpstr);
    // energy loop for DOS
    // only DOS will be considered and written out because the Integrated DOS
    // can be calculated from DOS and energy bins
    ncol = 1;
    if (DO_S) {
      ncol = 1 + 1 * nspec;
    }
    if (DO_P) {
      ncol = 1 + 2 * nspec;
    }
    if (DO_D) {
      ncol = 1 + 3 * nspec;
    }
    if (DO_F) {
      ncol = 1 + 4 * nspec;
    }
    if (ispin) {
      ncol = ncol * 2;
    }
    ncol++;
    const xmatrix<float> EDOS(nE, ncol);
    for (i = 1; i <= nE; i++) {
      din >> EDOS[i][1] >> EDOS[i][2];
      if (ispin == 2) {
        din >> EDOS[i][3];
      }
      getline(din, tmpstr);
    }
    // energy loop for DOS for s,p,d,f
    // sum over all atoms for the same SPECIES!
    // note that we read up to the highest between s,p,d,f
    // even though not all of them will be outputed
    // We will output only the orbitals that are requested
    // at the prompt input.
    if (DO_S) {
      maxOrbital = 1;
    }
    if (DO_P) {
      maxOrbital = 2;
    }
    if (DO_D) {
      maxOrbital = 3;
    }
    if (DO_F) {
      maxOrbital = 4;
    }
    if (maxOrbital > 0) {
      ioffset = 0;
      for (i = 1; i <= nspec; i++) {
        ioffset = (i - 1) * maxOrbital * ispin;
        for (j = 1; j <= species[i]; j++) {
          getline(din, tmpstr);
          for (k = 1; k <= nE; k++) {
            din >> ftmp; // discard energy
            if (ispin == 1) {
              for (iorb = 1; iorb <= maxOrbital; iorb++) {
                din >> ftmp;
                EDOS[k][ioffset + iorb + 2] += ftmp; // s
              }
              getline(din, tmpstr);
            }
            if (ispin == 2) {
              for (iorb = 1; iorb <= maxOrbital; iorb++) {
                din >> ftmp;
                EDOS[k][ioffset + (iorb - 1) * 2 + 4] += ftmp; // up
                din >> ftmp;
                EDOS[k][ioffset + (iorb - 1) * 2 + 5] -= ftmp; // down
              }
              getline(din, tmpstr);
            }
          }
        }
      }
    } // if maxOrbital>0
    // output
    for (i = 1; i <= nE; i++) {
      cout << "  " << EDOS[i][1]; // energy
      for (j = 1; j <= ispin; j++) {
        cout << "  " << EDOS[i][j + 1]; // DOS
      }
      if (DO_S) {
        //      cerr << "do S\n";
        for (j = 1; j <= nspec; j++) {
          ioffset = (j - 1) * maxOrbital * ispin;
          for (k = 1; k <= ispin; k++) {
            cout << "  " << EDOS[i][ioffset + k + 1 + ispin];
          }
        }
      }
      if (DO_P) {
        // cerr << "do P\n";
        for (j = 1; j <= nspec; j++) {
          ioffset = (j - 1) * maxOrbital * ispin + ispin;
          for (k = 1; k <= ispin; k++) {
            cout << "  " << EDOS[i][ioffset + k + 1 + ispin];
          }
        }
      }
      if (DO_D) {
        // cerr << "do D\n";
        for (j = 1; j <= nspec; j++) {
          ioffset = (j - 1) * maxOrbital * ispin + 2 * ispin;
          for (k = 1; k <= ispin; k++) {
            cout << "  " << EDOS[i][ioffset + k + 1 + ispin];
          }
        }
      }
      if (DO_F) {
        // cerr << "do F\n";
        for (j = 1; j <= nspec; j++) {
          ioffset = (j - 1) * maxOrbital * ispin + 3 * ispin;
          for (k = 1; k <= ispin; k++) {
            cout << "  " << EDOS[i][ioffset + k + 1 + ispin];
          }
        }
      }
      cout << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::EFFMASS // CAMILO
// ***************************************************************************
namespace pflow {
  void EFFMASS(vector<string>& argv, ostream& oss) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    // aflow --effective_mass directory_name
    // aflow --em             directory_name
    string WorkDir = argv.at(2);
    PrintEffectiveMass(WorkDir, oss);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::EQUIVALENT
// ***************************************************************************
namespace pflow {
  string EQUIVALENT(_aflags& aflags, istream& input, aurostd::xoption& vpflow) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    stringstream message;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    const xstructure _a(input, IOAFLOW_AUTO);
    const bool PRINT_SCREEN = false;
    aflags.QUIET = true;
    //  DEBUG=true;
    ofstream FileMESSAGE("/dev/null");
    // DX+CO START
    _kflags kflags;
    // DX20170815 - Add in consistency checks pflow::PerformFullSymmetry(a,File,aflags,kflags,PRINT_SCREEN,cout);
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_FGROUP = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL = false; // DX20170815 - Add in consistency checks //DX20171205 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON = false; // DX20200129
    kflags.KBIN_SYMMETRY_CALCULATE_SGROUP = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_IATOMS = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_AGROUP = false; // DX20170815 - Add in consistency checks
    const string options = vpflow.getattachedscheme("EQUIVALENT");
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() == 1) {
      if (tokens[0] == "usage" || tokens[0] == "USAGE") {
        init::MessageOption(options, __AFLOW_FUNC__,
                            aurostd::liststring2string("aflow --equivalent|--equiv|--inequivalent|--inequiv|--iatoms|--eatoms[=<tolerance_value>|=tight|=loose] [--no_scan] [--print=txt|json] "
                                                       "[--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR  default: tolerance=(minimum_interatomic_distance)/100.0"));
        return "";
      }
    }
    // DX20170804 - need to rescale, so we make a fast copy and calculate
    xstructure a(_a);
    // DX20180221 - use pwd - START
    if (a.directory.empty()) {
      if (aflags.Directory != "./") {
        a.directory = aflags.Directory;
      } else {
        a.directory = aurostd::getPWD();
        aflags.Directory = a.directory;
      }
    }
    const string print_directory = " [dir=" + a.directory + "]";
    // DX20180221 - use pwd - END
    a.ReScale(1.0);
    // DX20170921 - MAGNETIC SYMMETRY - START
    if (vpflow.flag("SYMMETRY::MAGNETIC")) {
      const string magmom_info = vpflow.getattachedscheme("SYMMETRY::MAGNETIC");
      ProcessAndAddSpinToXstructure(a, magmom_info); // DX20191108 - condensed into a single function
    }
    // DX20170921 - MAGNETIC SYMMETRY - END
    //  get tolerance
    double tolerance = pflow::getSymmetryTolerance(a, vpflow.getattachedscheme("SYMMETRY::TOLERANCE")); // DX20200820 - consolidated setting tolerance into a function
    // DX NOT REALLY NEEDED FOR THIS FUNCTION
    // DX20170803 - Add format flag - START
    string format = "txt";
    if (XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      format = "txt";
    } else if (XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
      format = "json";
    } else { // default is txt
      format = "txt";
    }
    // Perform full scan
    const bool force_perform = true; // if no_scan fails, still return true at default tolerance (even though it cannot be validated)
    if (vpflow.flag("SYMMETRY::NO_SCAN")) {
      a.sym_eps_no_scan = true; // DX20210406
    }

    if (!pflow::PerformFullSymmetry(a, tolerance, a.sym_eps_no_scan, force_perform, FileMESSAGE, aflags, kflags, PRINT_SCREEN, cout)) {
      message << "Could not find commensurate symmetry at tolerance = " << tolerance << " " << print_directory;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_); // CO20200624
    }
    // pflow::PerformFullSymmetry(a,File,aflags,kflags,PRINT_SCREEN,cout); //DX20170815 - Add in consistency checks
    // SYM::CalculatePointGroup(File,a,aflags,true,PRINT_SCREEN,cout);
    // SYM::CalculateFactorGroup(File,a,aflags,true,PRINT_SCREEN,cout);
    // SYM::CalculateSpaceGroup(File,a,aflags,false,PRINT_SCREEN,cout);
    // SYM::CalculateInequivalentAtoms(File,a,aflags,true,PRINT_SCREEN,cout);
    // DX+CO END
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    stringstream oss;
    if (format == "txt") {
      a.write_inequivalent_flag = true;
      stringstream oss;
      oss << a << endl;
      return oss.str();
    } else if (format == "json") {
      const char mode = _IATOMS_;
      KBIN_SymmetryToScreen(a, format, oss, mode);
      return oss.str();
    }
    return oss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::EWALD
// ***************************************************************************
namespace pflow {
  void EWALD(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() >= 2) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --ewald[=eta] < POSCAR");
    }
    // move on
    double eta = -1.0;
    if (!tokens.empty()) {
      eta = aurostd::string2utype<double>(tokens[0]);
    }
    // cout << aflow::Banner("BANNER_TINY") << endl;
    const double SUMTOL = 1.0e-16;
    xstructure str(input, IOAFLOW_AUTO);
    str = GetNiggliStr(str);
    double epoint;
    double ereal;
    double erecip;
    double eewald;
    pflow::Ewald(str, epoint, ereal, erecip, eewald, eta, SUMTOL);
    pflow::PrintEwald(str, epoint, ereal, erecip, eewald, eta, SUMTOL, cout);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::EXTRACT_xcar
// ***************************************************************************
namespace pflow {
  string EXTRACT_xcar(_aflags& aflags, vector<string> argv, string mode, string file) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " mode=" << mode << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " file=" << file << endl;
    }
    if (!argv.empty()) {
      ;
    } // phony just to keep argv busy no complaining about unused
    ofstream FileMESSAGE("/dev/null");
    _kflags kflags;
    kflags.AFLOW_MODE_VASP = true;
    _vflags vflags;
    _xvasp xvasp;
    xvasp.clear();
    if (!aurostd::FileExist(file)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "mode=" + mode + "  file=" + file + " not found", _FILE_CORRUPT_); // CO20200624
    }
    string AflowIn;
    aurostd::file2string(file, AflowIn);
    AflowIn = aurostd::RemoveComments(AflowIn); // NOW Clean AFLOWIN //CO20180502
    aflags.QUIET = true;
    XHOST.QUIET = true;
    if (mode == "POSCAR") {
      vflags = KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn, aflags, kflags);
      KBIN::VASP_Produce_POSCAR(xvasp, AflowIn, FileMESSAGE, aflags, kflags, vflags);
      return xvasp.POSCAR.str();
    }
    if (mode == "INCAR") {
      vflags = KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn, aflags, kflags);
      KBIN::VASP_Produce_INCAR(xvasp, AflowIn, FileMESSAGE, aflags, kflags, vflags);
      KBIN::VASP_Modify_INCAR(xvasp, FileMESSAGE, aflags, kflags, vflags);
      return xvasp.INCAR.str();
    }
    if (mode == "KPOINTS") {
      vflags = KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn, aflags, kflags);
      KBIN::VASP_Produce_POSCAR(xvasp, AflowIn, FileMESSAGE, aflags, kflags, vflags);
      KBIN::VASP_Produce_KPOINTS(xvasp, AflowIn, FileMESSAGE, aflags, kflags, vflags);
      KBIN::VASP_Modify_KPOINTS(xvasp, FileMESSAGE, aflags, vflags);
      return xvasp.KPOINTS.str();
    }
    if (mode == "POTCAR") {
      vflags = KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn, aflags, kflags);
      KBIN::VASP_Produce_POTCAR(xvasp, AflowIn, FileMESSAGE, aflags, kflags, vflags);
      return xvasp.POTCAR.str();
    }
    if (mode == "PARTCAR") { // CO20181226
      const xstructure xstr = pocc::extractPARTCAR(AflowIn);
      stringstream output;
      output << xstr;
      return output.str();
    }
    return mode; // something must go out
  }
} // namespace pflow

// ***************************************************************************
// pflow::EIGCURV // CAMILO
// ***************************************************************************
namespace pflow {
  void EIGCURV(const string& options, ostream& oss) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 1) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --eigcurv=DIRECTORY(with bands)");
    }
    string filename;
    if (!tokens.empty()) {
      filename = (tokens[0]);
    }
    string WorkDir = filename;
    PrintEigCurv(WorkDir, oss);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// DX+CO START
//  ***************************************************************************
//  pflow::PerformFullSymmetry
//  ***************************************************************************
namespace pflow {
  bool PerformFullSymmetry(xstructure& a) {
    ofstream FileMESSAGE("/dev/null");
    _aflags aflags;
    aflags.Directory = ".";
    aflags.QUIET = true;
    _kflags kflags;
    defaultKFlags4SymWrite(kflags, false);
    defaultKFlags4SymCalc(kflags, true);
    const bool osswrite = false;
    const bool see_output = false; // true;
    ostream& oss_empty = cout;
    if (!see_output) {
      oss_empty.setstate(std::ios_base::badbit);
    } // like NULL
    const bool sym_done = PerformFullSymmetry(a, FileMESSAGE, aflags, kflags, osswrite, oss_empty);
    if (!see_output) {
      oss_empty.clear();
    } // clear NULL
    return sym_done;
  }

  // ME20200224 - added directory option
  bool PerformFullSymmetry(xstructure& a, ofstream& FileMESSAGE, const string& directory, _kflags& kflags, bool osswrite, ostream& oss, string format) {
    _aflags aflags;
    aflags.Directory = directory;
    return PerformFullSymmetry(a, FileMESSAGE, aflags, kflags, osswrite, oss, format);
  }

  bool PerformFullSymmetry(xstructure& a, ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, bool osswrite, ostream& oss, string format) {
    double tolerance = a.sym_eps;
    const bool no_scan = false;
    const bool force_perform = true; // if no_scan fails, still return true at default tolerance (even though it cannot be validated)
    return PerformFullSymmetry(a, tolerance, no_scan, force_perform, FileMESSAGE, aflags, kflags, osswrite, oss, format);
  }

  bool PerformFullSymmetry(xstructure& a, double& tolerance, bool no_scan, bool force_perform, ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, bool osswrite, ostream& oss, string format) {
    a.ClearSymmetry(); // CO20190204
    const xstructure b(a); // save for later
    a.ReScale(1.0); // the nuclear option, only way to fix all of the issues with f2c/c2f/ctau/ctrasl/etc.
    const bool LDEBUG = (false || XHOST.DEBUG);
    bool symmetry_commensurate = false;
    // ME20210429 - use a function-wide QUIET flag
    const bool QUIET = (aflags.QUIET || XHOST.QUIET || XHOST.vflag_control.flag("WWW"));
    // ME20210429 - Do not write symmetry information for web; too much output
    osswrite = (osswrite && !XHOST.vflag_control.flag("WWW"));
    ostringstream aus;

    // a.directory = aflags.Directory;
    // DX20180221 - use pwd - START
    if (a.directory.empty()) {
      if (aflags.Directory != "./") {
        a.directory = aflags.Directory;
      } else {
        a.directory = aurostd::getPWD();
        aflags.Directory = a.directory;
      }
    }
    const string print_directory = " [dir=" + a.directory + "]";
    // DX20180221 - use pwd - END

    if (LDEBUG) {
      cerr << XPID << "pflow::PerformFullSymmetry: STRUCTURE" << endl;
      cerr << a << endl;
    }

    if (a.atoms.empty()) {
      cerr << __AFLOW_FUNC__ << " ERROR! No atoms found in the structure" << print_directory << endl;
      return false;
    }

    // MOVED DOWN A BIT if(!aflags.QUIET) aus << XPID << (aflags.QUIET?"":"00000  MESSAGE ") << "Symmetry: starting tolerance " << _EPS_sym_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
    if (true || a.dist_nn_min == AUROSTD_NAN) { // CO20171024 - always recalculate min_dist (SAFE)
      if (LDEBUG) {
        cerr << XPID << "pflow::PerformFullSymmetry: CALCULATING MIN DISTANCE" << print_directory << endl;
      }
      a.MinDist();
      if (LDEBUG) {
        cerr << XPID << "pflow::PerformFullSymmetry: MIN DISTANCE DONE" << print_directory << endl;
      }
    }
    double min_dist = a.dist_nn_min;
    // CO20180420 - we need something tighter here, but this is good to kill POCC
    if (min_dist < _ZERO_TOL_) { //_XPROTO_TOO_CLOSE_ERROR_ perhaps?
      cerr << XPID << "pflow::PerformFullSymmetry: ERROR! Atoms too close (min_dist=" << min_dist << print_directory << endl;
      return false;
    }
    // if(a.sym_eps!=AUROSTD_NAN){ //Tolerance came from user or was calculated
    //   a.sym_eps;
    // }
    // CO, I changed a bit, if tolerance is specified, it should override
    if (tolerance != AUROSTD_NAN) {
      a.sym_eps = tolerance;
    } else if (!a.sym_eps_calculated || a.sym_eps == AUROSTD_NAN) {
      a.sym_eps = SYM::defaultTolerance(a);
    }
    // if(a.sym_eps == AUROSTD_NAN && tolerance == AUROSTD_NAN){
    //   a.sym_eps = SYM::defaultTolerance(a);
    // }
    // else if(!a.sym_eps_calculated && tolerance != AUROSTD_NAN){
    //   a.sym_eps = tolerance;
    // }
    // a.sym_eps = SYM::defaultTolerance(a);
    double orig_tolerance = a.sym_eps;
    if (!aflags.QUIET) {
      aus << XPID << (aflags.QUIET ? "" : "00000  MESSAGE ") << "Symmetry: starting tolerance " << a.sym_eps << " " << Message(__AFLOW_FILE__, aflags) << endl;
    }
    aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, osswrite, oss);

    a.sym_eps_no_scan = no_scan; // DX20210406
    if (LDEBUG) {
      cerr << XPID << "pflow::PerformFullSymmetry: no_scan=" << a.sym_eps_no_scan << endl;
    }

    while (symmetry_commensurate == false) {
      a.sym_eps_calculated = true;
      // Calculate Lattice Point Group
      if (kflags.KBIN_SYMMETRY_CALCULATE_PGROUP) { // DX20170814
        if (!SYM::CalculatePointGroup(FileMESSAGE, a, aflags, kflags.KBIN_SYMMETRY_PGROUP_WRITE, osswrite, oss, format)) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " WARNING: PGROUP calculation is inconsisent." << endl;
          }
          if (!a.sym_eps_no_scan) {
            a.ClearSymmetry();
            if (!SYM::change_tolerance(a, a.sym_eps, min_dist, a.sym_eps_no_scan)) { // CO20200106 - patching for auto-indenting
              a = b; // pretty printing, unmodified structure
              if (force_perform) {
                cerr << __AFLOW_FUNC__ << " Scan failed [0]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
                PerformFullSymmetry(a, orig_tolerance, true, false, FileMESSAGE, aflags, kflags, osswrite, oss, format);
              } else {
                return false;
              }
            }
            if (!QUIET) {
              aus << XPID << "00000  MESSAGE PGROUP Symmetry: changing tolerance to " << a.sym_eps << " " << Message(__AFLOW_FILE__, aflags) << endl;
            }
            aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
            continue;
          }
        }
      } // DX20170814
      if (kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK) { // DX20170814
        if (!SYM::CalculatePointGroupKLattice(FileMESSAGE, a, aflags, kflags.KBIN_SYMMETRY_PGROUPK_WRITE, osswrite, oss, format)) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " WARNING: PGROUPK calculation is inconsisent." << endl;
          }
          if (!a.sym_eps_no_scan) {
            a.ClearSymmetry();
            if (!SYM::change_tolerance(a, a.sym_eps, min_dist, a.sym_eps_no_scan)) { // CO20200106 - patching for auto-indenting
              a = b; // pretty printing, unmodified structure
              if (force_perform) {
                cerr << __AFLOW_FUNC__ << " Scan failed [1]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
                PerformFullSymmetry(a, orig_tolerance, true, false, FileMESSAGE, aflags, kflags, osswrite, oss, format);
              } else {
                return false;
              }
            }
            if (!QUIET) {
              aus << XPID << "00000  MESSAGE PGROUPK Symmetry: changing tolerance to " << a.sym_eps << " " << Message(__AFLOW_FILE__, aflags) << endl;
            }
            aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
            continue;
          }
        }
      } // DX20170814
      // Check for identity element
      if (kflags.KBIN_SYMMETRY_CALCULATE_PGROUP && SYM::CheckForIdentity(a) == false) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " WARNING: Point group does not contain the identity element (impossible for a crystal)." << print_directory << endl;
        }
        if (!a.sym_eps_no_scan) {
          a.ClearSymmetry();
          if (!SYM::change_tolerance(a, a.sym_eps, min_dist, a.sym_eps_no_scan)) { // CO20200106 - patching for auto-indenting
            a = b; // pretty printing, unmodified structure
            if (force_perform) {
              cerr << __AFLOW_FUNC__ << " Scan failed [2]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
              PerformFullSymmetry(a, orig_tolerance, true, false, FileMESSAGE, aflags, kflags, osswrite, oss, format);
            } else {
              return false;
            }
          }
          if (!QUIET) {
            aus << XPID << "00000  MESSAGE PGROUP Symmetry: changing tolerance to " << a.sym_eps << " " << Message(__AFLOW_FILE__, aflags) << endl;
          }
          aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
          continue;
        }
      }
      // Calculate Factor Group
      if (kflags.KBIN_SYMMETRY_CALCULATE_FGROUP) { // DX20170814
        if (!SYM::CalculateFactorGroup(FileMESSAGE, a, aflags, kflags.KBIN_SYMMETRY_FGROUP_WRITE, osswrite, oss, format)) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " WARNING: FGROUP calculation is inconsisent." << endl;
          }
          if (!a.sym_eps_no_scan) {
            a.ClearSymmetry();
            if (!SYM::change_tolerance(a, a.sym_eps, min_dist, a.sym_eps_no_scan)) { // CO20200106 - patching for auto-indenting
              a = b; // pretty printing, unmodified structure
              if (force_perform) {
                cerr << __AFLOW_FUNC__ << " Scan failed [3]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
                PerformFullSymmetry(a, orig_tolerance, true, false, FileMESSAGE, aflags, kflags, osswrite, oss, format);
              } else {
                return false;
              }
            }
            if (!QUIET) {
              aus << XPID << "00000  MESSAGE FGROUP Symmetry: changing tolerance to " << a.sym_eps << " " << Message(__AFLOW_FILE__, aflags) << endl;
            }
            aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
            continue;
          }
        }
      } // DX20170814
      // Calculate Crystallographic Point Group
      if (kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL) { // DX20170814
        if (!SYM::CalculatePointGroupCrystal(FileMESSAGE, a, aflags, kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE, osswrite, oss, format)) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " WARNING: PGROUP_XTAL calculation is inconsisent." << endl;
          }
          if (!a.sym_eps_no_scan) {
            a.ClearSymmetry();
            if (!SYM::change_tolerance(a, a.sym_eps, min_dist, a.sym_eps_no_scan)) { // CO20200106 - patching for auto-indenting
              a = b; // pretty printing, unmodified structure
              if (force_perform) {
                cerr << XPID << "pflow::PerformFullSymmetry: Scan failed [4]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
                PerformFullSymmetry(a, orig_tolerance, true, false, FileMESSAGE, aflags, kflags, osswrite, oss, format);
              } else {
                return false;
              }
            }
            if (!QUIET) {
              aus << XPID << "00000  MESSAGE PGROUP_XTAL Symmetry: changing tolerance to " << a.sym_eps << " " << Message(__AFLOW_FILE__, aflags) << endl;
            }
            aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
            continue;
          }
        }
      } // DX20170814
      // Check if a point group map was found; if not, change tolerance
      if (kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL && a.point_group_Hermann_Mauguin.empty() == true) { // DX20170814
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " WARNING: Point group crystal operations did not match with any Hermann-Mauguin symbols. (i.e. The set of symmetry elements found are not allowed possible for a crystal.) "
               << print_directory << endl;
          ;
        }
        if (!a.sym_eps_no_scan) {
          a.ClearSymmetry();
          if (!SYM::change_tolerance(a, a.sym_eps, min_dist, a.sym_eps_no_scan)) { // CO20200106 - patching for auto-indenting
            a = b; // pretty printing, unmodified structure
            if (force_perform) {
              cerr << __AFLOW_FUNC__ << " Scan failed [5]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
              PerformFullSymmetry(a, orig_tolerance, true, false, FileMESSAGE, aflags, kflags, osswrite, oss, format);
            } else {
              return false;
            }
          }
          if (!QUIET) {
            aus << XPID << "00000  MESSAGE PGROUP_XTAL Symmetry: changing tolerance to " << a.sym_eps << " " << Message(__AFLOW_FILE__, aflags) << endl;
          }
          aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
          continue;
        }
      }
      // Check if factor group is integer multiple of the crystallographic point group; if not, change tolerance
      int multiplicity_of_primitive = 0; // DX20170815 - Added consistency checks, need to initialize and calculate if we have those groups
      if (kflags.KBIN_SYMMETRY_CALCULATE_FGROUP && kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL) { // DX20170814
        multiplicity_of_primitive = a.fgroup.size() / a.pgroup_xtal.size();
        if (a.fgroup.size() % a.pgroup_xtal.size() != 0) { // DX20170814
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " WARNING: Number of factor groups is not an integer multiple of the point group crystal (fgroup: " << a.fgroup.size() << " vs pgroup_xtal: " << a.pgroup_xtal.size()
                 << ")." << print_directory << endl;
          }
          if (!a.sym_eps_no_scan) {
            a.ClearSymmetry();
            if (!SYM::change_tolerance(a, a.sym_eps, min_dist, a.sym_eps_no_scan)) { // CO20200106 - patching for auto-indenting
              a = b; // pretty printing, unmodified structure
              if (force_perform) {
                cerr << __AFLOW_FUNC__ << " Scan failed [6]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
                PerformFullSymmetry(a, orig_tolerance, true, false, FileMESSAGE, aflags, kflags, osswrite, oss, format);
              } else {
                return false;
              }
            }
            if (!QUIET) {
              aus << XPID << "00000  MESSAGE PGROUP_XTAL Symmetry: changing tolerance to " << a.sym_eps << " " << Message(__AFLOW_FILE__, aflags) << endl;
            }
            aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
            continue;
          }
        }
      } // DX20170814
      if (kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL) { // DX20171205 - Added pgroupk_xtal
        if (!SYM::CalculatePointGroupKCrystal(FileMESSAGE, a, aflags, kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE, osswrite, oss, format)) { // DX20180118 - PGROUPK_XTAL not PGROUPK
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " WARNING: PGROUPK_XTAL calculation is inconsisent." << endl;
          }
          if (!a.sym_eps_no_scan) {
            a.ClearSymmetry();
            if (!SYM::change_tolerance(a, a.sym_eps, min_dist, a.sym_eps_no_scan)) { // CO20200106 - patching for auto-indenting
              a = b; // pretty printing, unmodified structure
              if (force_perform) {
                cerr << __AFLOW_FUNC__ << " Scan failed [7]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
                PerformFullSymmetry(a, orig_tolerance, true, false, FileMESSAGE, aflags, kflags, osswrite, oss, format);
              } else {
                return false;
              }
            }
            if (!QUIET) {
              aus << XPID << "00000  MESSAGE PGROUPK_XTAL Symmetry: changing tolerance to " << a.sym_eps << " " << Message(__AFLOW_FILE__, aflags) << endl;
            }
            aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
            continue;
          }
        }
      } // DX20170814
      // Calculate Patterson Point Group
      if (kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON) { // DX20200129
        if (!SYM::CalculatePointGroupKPatterson(FileMESSAGE, a, aflags, kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE, osswrite, oss, format)) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " WARNING: PGROUPK_PATTERSON calculation is inconsisent." << endl;
          }
          if (!a.sym_eps_no_scan) {
            a.ClearSymmetry();
            if (!SYM::change_tolerance(a, a.sym_eps, min_dist, a.sym_eps_no_scan)) { // CO20200106 - patching for auto-indenting
              a = b; // pretty printing, unmodified structure
              if (force_perform) {
                cerr << __AFLOW_FUNC__ << " Scan failed [0]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
                PerformFullSymmetry(a, orig_tolerance, true, false, FileMESSAGE, aflags, kflags, osswrite, oss, format);
              } else {
                return false;
              }
            }
            if (!QUIET) {
              aus << XPID << "00000  MESSAGE PGROUPK_PATTERSON Symmetry: changing tolerance to " << a.sym_eps << " " << Message(__AFLOW_FILE__, aflags) << endl;
            }
            aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
            continue;
          }
        }
      } // DX20200129
      if (kflags.KBIN_SYMMETRY_CALCULATE_SGROUP) { // DX20170814
        if (kflags.KBIN_SYMMETRY_SGROUP_RADIUS > 0.0) {
          if (!QUIET) {
            aus << XPID << "00000  MESSAGE POSCAR SGROUP: found RADIUS=" << kflags.KBIN_SYMMETRY_SGROUP_RADIUS << " " << Message(__AFLOW_FILE__, aflags) << endl;
          }
          aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
        } else {
          kflags.KBIN_SYMMETRY_SGROUP_RADIUS = KBIN_SYMMETRY_SGROUP_RADIUS_DEFAULT;
          if (!QUIET) {
            aus << XPID << "00000  MESSAGE POSCAR SGROUP: Default RADIUS=" << kflags.KBIN_SYMMETRY_SGROUP_RADIUS << " " << Message(__AFLOW_FILE__, aflags) << endl;
          }
          aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
        }
        a.sgroup_radius = kflags.KBIN_SYMMETRY_SGROUP_RADIUS;
        // Calculate Space Group
        if (!SYM::CalculateSpaceGroup(FileMESSAGE, a, aflags, kflags.KBIN_SYMMETRY_SGROUP_WRITE, osswrite, oss, format)) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " WARNING: SGROUP calculation is inconsisent." << endl;
          }
          if (!a.sym_eps_no_scan) {
            a.ClearSymmetry();
            if (!SYM::change_tolerance(a, a.sym_eps, min_dist, a.sym_eps_no_scan)) { // CO20200106 - patching for auto-indenting
              a = b; // pretty printing, unmodified structure
              if (force_perform) {
                cerr << XPID << "pflow::PerformFullSymmetry: Scan failed [8]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
                PerformFullSymmetry(a, orig_tolerance, true, false, FileMESSAGE, aflags, kflags, osswrite, oss, format);
              } else {
                return false;
              }
            }
            if (!QUIET) {
              aus << XPID << "00000  MESSAGE SGROUP Symmetry: changing tolerance to " << a.sym_eps << " " << Message(__AFLOW_FILE__, aflags) << endl;
            }
            aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
            continue;
          }
        }
      } // DX20170814
      // Calculate inequivalent atoms
      if (kflags.KBIN_SYMMETRY_CALCULATE_IATOMS) { // DX20170814
        if (!SYM::CalculateInequivalentAtoms(FileMESSAGE, a, aflags, kflags.KBIN_SYMMETRY_IATOMS_WRITE, osswrite, oss, format)) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " WARNING: IATOMS calculation is inconsisent." << endl;
          }
          if (!a.sym_eps_no_scan) {
            a.ClearSymmetry();
            if (!SYM::change_tolerance(a, a.sym_eps, min_dist, a.sym_eps_no_scan)) { // CO20200106 - patching for auto-indenting
              a = b; // pretty printing, unmodified structure
              if (force_perform) {
                cerr << __AFLOW_FUNC__ << " Scan failed [9]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
                PerformFullSymmetry(a, orig_tolerance, true, false, FileMESSAGE, aflags, kflags, osswrite, oss, format);
              } else {
                return false;
              }
            }
            if (!QUIET) {
              aus << XPID << "00000  MESSAGE IATOMS ATOMS: changing tolerance to " << a.sym_eps << " " << Message(__AFLOW_FILE__, aflags) << endl;
            }
            aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
            continue;
          }
        }
        // Check if number of equivalent atoms is consistent with cell ; if not, change tolerance
        bool iatoms_commensurate = true;
        for (size_t i = 0; i < a.iatoms.size(); i++) {
          if (a.iatoms[i].size() % multiplicity_of_primitive != 0) {
            iatoms_commensurate = false;
            break;
          }
        }
        if (iatoms_commensurate == false) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " WARNING: Number of equivalent atoms is not an integer multiple of the number factor groups." << print_directory << endl;
          }
          if (!a.sym_eps_no_scan) {
            a.ClearSymmetry();
            if (!SYM::change_tolerance(a, a.sym_eps, min_dist, a.sym_eps_no_scan)) { // CO20200106 - patching for auto-indenting
              a = b; // pretty printing, unmodified structure
              if (force_perform) {
                cerr << XPID << "pflow::PerformFullSymmetry: Scan failed [10]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
                PerformFullSymmetry(a, orig_tolerance, true, false, FileMESSAGE, aflags, kflags, osswrite, oss, format);
              } else {
                return false;
              }
            }
            if (!QUIET) {
              aus << XPID << "00000  MESSAGE IATOMS ATOMS: changing tolerance to " << a.sym_eps << " " << Message(__AFLOW_FILE__, aflags) << endl;
            }
            aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
            continue;
          }
        }
      } // DX20170814
      // Calculate site point group
      if (kflags.KBIN_SYMMETRY_CALCULATE_AGROUP) { // DX20170814
        if (!SYM::CalculateSitePointGroup(FileMESSAGE, a, aflags, kflags.KBIN_SYMMETRY_AGROUP_WRITE, osswrite, oss, format)) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " WARNING: AGROUP calculation is inconsisent." << endl;
          }
          if (!a.sym_eps_no_scan) {
            a.ClearSymmetry();
            if (!SYM::change_tolerance(a, a.sym_eps, min_dist, a.sym_eps_no_scan)) { // CO20200106 - patching for auto-indenting
              a = b; // pretty printing, unmodified structure
              if (force_perform) {
                cerr << __AFLOW_FUNC__ << " Scan failed [11]. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << print_directory << endl;
                PerformFullSymmetry(a, orig_tolerance, true, false, FileMESSAGE, aflags, kflags, osswrite, oss, format);
              } else {
                return false;
              }
            }
            if (!QUIET) {
              aus << XPID << "00000  MESSAGE AGROUP Symmetry: changing tolerance to " << a.sym_eps << " " << Message(__AFLOW_FILE__, aflags) << endl;
            }
            aurostd::PrintMessageStream(FileMESSAGE, aus, QUIET, osswrite, oss);
            continue;
          }
        }
      } // DX20170814
      symmetry_commensurate = true; // NOTE: This may not be entirely true if no_scan=true
    }
    a.ReScale(b.scale); // the nuclear option, only way to fix all of the issues with f2c/c2f/ctau/ctrasl/etc.
    a.lattice = b.lattice;
    a.scale = b.scale; // perfect printing, we should also fix cpos of all atoms, but not really worried since we usually print fpos

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " SYMMETRY CALCULATION IS COMPLETE." << endl;
      cerr << __AFLOW_FUNC__ << " # pgroup: " << a.pgroup.size() << endl;
      cerr << __AFLOW_FUNC__ << " # pgroupk: " << a.pgroupk.size() << endl;
      cerr << __AFLOW_FUNC__ << " # fgroup: " << a.fgroup.size() << endl;
      cerr << __AFLOW_FUNC__ << " # pgroup_xtal: " << a.pgroup_xtal.size() << endl;
      cerr << __AFLOW_FUNC__ << " # pgroupk_xtal: " << a.pgroupk_xtal.size() << endl;
      cerr << __AFLOW_FUNC__ << " # pgroupk_Patterson: " << a.pgroupk_Patterson.size() << endl;
      cerr << __AFLOW_FUNC__ << " # sgroup: " << a.sgroup.size() << endl;
      cerr << __AFLOW_FUNC__ << " # agroup: " << a.agroup.size() << endl;
    }
    return symmetry_commensurate;
  }
} // namespace pflow

namespace pflow {
  void defaultKFlags4SymWrite(_kflags& kflags, bool write) {
    kflags.KBIN_SYMMETRY_PGROUP_WRITE = write;
    kflags.KBIN_SYMMETRY_PGROUPK_WRITE = write;
    kflags.KBIN_SYMMETRY_FGROUP_WRITE = write;
    kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE = write;
    kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE = write; // DX20171205 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE = write; // DX20200129
    kflags.KBIN_SYMMETRY_IATOMS_WRITE = write;
    kflags.KBIN_SYMMETRY_AGROUP_WRITE = write;
    kflags.KBIN_SYMMETRY_SGROUP_WRITE = write;
  }
  void defaultKFlags4SymCalc(_kflags& kflags, bool calc) {
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP = calc;
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK = calc;
    kflags.KBIN_SYMMETRY_CALCULATE_FGROUP = calc;
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL = calc;
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL = calc; // DX20171205 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON = calc; // DX20200129
    kflags.KBIN_SYMMETRY_CALCULATE_IATOMS = calc;
    kflags.KBIN_SYMMETRY_CALCULATE_AGROUP = calc;
    kflags.KBIN_SYMMETRY_CALCULATE_SGROUP = calc;
  }
} // namespace pflow

namespace pflow {
  // COMMAND LINE SYMMETRY CALCULATION, calls main function PerformFullSymmetry()!!!!!!!!!!!
  bool CalculateFullSymmetry(istream& input, aurostd::xoption& vpflow, ostream& oss) { // overload

    // DX20201228 - print Python script
    if (XHOST.vflag_control.flag("PRINT_MODE::PYTHON")) {
      SYM::writePythonScript(oss);
      return true;
    }

    xstructure a(input, IOAFLOW_AUTO);

    // default aflags for command line
    _aflags aflags;
    aflags.Directory = ".";
    aflags.QUIET = false;

    // default kflags
    _kflags kflags;
    defaultKFlags4SymWrite(kflags, true);
    defaultKFlags4SymCalc(kflags, true);

    // write output to screen
    const bool osswrite = true;

    return CalculateFullSymmetry(aflags, kflags, a, vpflow, osswrite, oss);
  }

  // COMMAND LINE SYMMETRY CALCULATION, calls main function PerformFullSymmetry()!!!!!!!!!!!
  bool CalculateFullSymmetry(_aflags& aflags, _kflags& kflags, xstructure& _a, aurostd::xoption& vpflow, bool osswrite, ostream& oss) { // main function
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    const string options = vpflow.getattachedscheme("FULLSYMMETRY");
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() == 1) {
      if (tokens[0] == "usage" || tokens[0] == "USAGE") {
        init::MessageOption(options, __AFLOW_FUNC__,
                            aurostd::liststring2string("aflow --aflow-sym|--AFLOW-SYM|--AFLOWSYM|--aflowSYM|--aflowsym|--full_symmetry|--full_sym|--fullsym[=<tolerance_value>|=tight|=loose] [--no_scan] "
                                                       "[--print=txt|json] [--screen_only] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR  default: "
                                                       "tolerance=(minimum_interatomic_distance)/100.0, print=txt")); // DX20200724 - removed return
        return false; // CO20200624 - there has been NO calculation of symmetry
      }
    }

    // DX20201228 - print Python script
    if (XHOST.vflag_control.flag("PRINT_MODE::PYTHON")) {
      SYM::writePythonScript(oss);
      return true;
    }

    // DX20170804 - need to rescale, so we make a fast copy and calculate
    xstructure a(_a);
    a.ReScale(1.0);

    // DX20170921 - MAGNETIC SYMMETRY - START
    if (vpflow.flag("FULLSYMMETRY::MAGNETIC")) {
      const string magmom_info = vpflow.getattachedscheme("FULLSYMMETRY::MAGNETIC");
      ProcessAndAddSpinToXstructure(a, magmom_info); // DX20191108 - condensed into a single function
    }
    // DX20170921 - MAGNETIC SYMMETRY - END

    // get tolerance
    double tolerance = pflow::getSymmetryTolerance(a, vpflow.getattachedscheme("FULLSYMMETRY::TOLERANCE")); // DX20200820 - consolidated setting tolerance into a function

    // DX20170803 - Add format flag - START
    string format = "txt";
    if (XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      format = "txt";
    } else if (XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
      format = "json";
    } else { // default is txt
      format = "txt";
    }
    // DX20170803 - Add format flag - END
    // DX20170803 - Print ouptut to screen - START
    bool print = false;
    if (vpflow.flag("FULLSYMMETRY::SCREEN_ONLY")) {
      print = true;
      defaultKFlags4SymWrite(kflags, false);
      // kflags.KBIN_SYMMETRY_PGROUP_WRITE=false;
      // kflags.KBIN_SYMMETRY_PGROUPK_WRITE=false;
      // kflags.KBIN_SYMMETRY_FGROUP_WRITE=false;
      // kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE=false;
      // kflags.KBIN_SYMMETRY_SGROUP_WRITE=false;
      // kflags.KBIN_SYMMETRY_IATOMS_WRITE=false;
      // kflags.KBIN_SYMMETRY_AGROUP_WRITE=false;
      osswrite = false;
    }
    // DX20170803 - Print ouptut to screen - END
    // DX
    //     if(tokens.size()==0) {tolerance=default_tolerance;}
    //     if(tokens.size()>=1 && tokens.at(0) != "--debug" && tokens.at(0) != "--debug" ) {
    //       if(tokens.at(0).at(0) == 't' || tokens.at(0).at(0) == 'T'){ //Tight
    //         tolerance=default_tolerance;
    //       }
    //       else if(tokens.at(0).at(0) == 'l' || tokens.at(0).at(0) == 'L'){ //Loose
    //         tolerance=default_tolerance*10.0;
    //       }
    //       else {
    //         tolerance=aurostd::string2utype<double>(tokens.at(0));
    //       }
    //     }
    // DX
    //  Perform full scan
    if (vpflow.flag("FULLSYMMETRY::NO_SCAN")) {
      a.sym_eps_no_scan = true; // DX20210406
    }

    bool tocompress = true;
    ofstream FileMESSAGE("/dev/null");
    string directory = aurostd::getPWD();
    if (XHOST.vflag_control.flag("DIRECTORY")) {
      directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");
    }
    aflags.Directory = directory;
    if (!aurostd::FileExist(directory)) {
      oss << "pflow::CalculateFullSymmetry: ERROR: Unable to locate " << directory << ". Exiting." << endl;
      // return oss.str();
      return false;
    }

    // while(symmetry_commensurate==false)
    if (print == false) // DX20170803 - PRINT
    { // CO20200106 - patching for auto-indenting
      if (format == "txt") { // DX20200623 - only overwrite the specified file format
        for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_OUT + XHOST.vext[iext]) ||
              aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT + XHOST.vext[iext]) || // DX20200129
              aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_OUT + XHOST.vext[iext]) || // DX20171205 - Added pgroupk_xtal
              aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_OUT + XHOST.vext[iext]) ||
              aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_OUT + XHOST.vext[iext])) {
            tocompress = true;
          }
        }
        for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_OUT)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUP_OUT + ".*"));
          }
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_OUT)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUP_XTAL_OUT + ".*"));
          }
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT + ".*"));
          } // DX20200129
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_OUT)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_OUT + ".*"));
          }
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_OUT)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_XTAL_OUT + ".*"));
          } // DX20171205 - Added pgroupk_xtal
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_OUT)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_FGROUP_OUT + ".*"));
          }
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_OUT)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_AGROUP_OUT + ".*"));
          }
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_OUT)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_SGROUP_OUT + ".*"));
          }
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_OUT)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_IATOMS_OUT + ".*"));
          }
        }
      } else if (format == "json") { // DX20200623 - only overwrite the specified file format
        for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_JSON + XHOST.vext[iext]) ||
              aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON + XHOST.vext[iext]) || // DX20200129
              aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_JSON + XHOST.vext[iext]) || // DX20171205 - Added pgroupk_xtal
              aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_JSON + XHOST.vext[iext]) ||
              aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_JSON + XHOST.vext[iext])) {
            tocompress = true;
          }
        }
        for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_JSON)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUP_JSON + ".*"));
          }
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_JSON)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUP_XTAL_JSON + ".*"));
          }
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON + ".*"));
          } // DX20200129
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_JSON)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_JSON + ".*"));
          }
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_JSON)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_XTAL_JSON + ".*"));
          } // DX20180516 - added pgroupk_xtal
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_JSON)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_FGROUP_JSON + ".*"));
          }
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_JSON)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_AGROUP_JSON + ".*"));
          }
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_JSON)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_SGROUP_JSON + ".*"));
          }
          if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_JSON)) {
            aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_IATOMS_JSON + ".*"));
          }
        }
      }
    }

    if (!pflow::PerformFullSymmetry(a, tolerance, a.sym_eps_no_scan, true, FileMESSAGE, aflags, kflags, osswrite, oss, format)) {
      return false;
    }

    // DX20170803 - Print to symmetry operators to screen - START
    if (print == true) {
      KBIN_SymmetryToScreen(a, format, oss);
    }
    // DX20170803 - Print to symmetry operators to screen - END

    // BZIP if necessary
    if (tocompress) {
      if (format == "txt") { // DX20170803 - FORMAT
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUP_OUT);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_OUT);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT); // DX20200129
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_OUT);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_OUT); // DX20171205 - Added pgroupk_xtal
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_FGROUP_OUT);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_SGROUP_OUT);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_AGROUP_OUT);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_OUT)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_IATOMS_OUT);
        }
      } else if (format == "json") { // DX20170803 - FORMAT
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUP_JSON);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_JSON);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON); // DX20200129
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_JSON);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_JSON); // DX20171205 - Added pgroupk_xtal
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_FGROUP_JSON);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_SGROUP_JSON);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_AGROUP_JSON);
        }
        if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_JSON)) {
          aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_IATOMS_JSON);
        }
      }
    }
    // return oss.str(); //string("AFLOW "+string(AFLOW_VERSION)+" Symmetry Fixed in "+directory+"  (need aflow>=2948)\n");
    return true; // string("AFLOW "+string(AFLOW_VERSION)+" Symmetry Fixed in "+directory+"  (need aflow>=2948)\n");
  }
} // namespace pflow

// ***************************************************************************
// pflow::fixEmptyAtomNames
// ***************************************************************************
// fixes POSCARs lacking column of atom names
// grabs from POTCAR information of aflow.in
// force_fix: override if a column is present, this may be crucial to fix issues with write_inequivalent_flag
namespace pflow {
  bool fixEmptyAtomNames(xstructure& xstr, bool force_fix) { // CO20190219
    xstr.fixEmptyAtomNames(force_fix);
    return true;
  }
} // namespace pflow
// DX+CO END

// ***************************************************************************
// pflow::EXTRACT_Symmetry
// ***************************************************************************
namespace pflow {
  string EXTRACT_Symmetry(_aflags& aflags, vector<string> argv) {
    string directory = argv.at(2);
    bool tocompress = false;
    directory = aurostd::RemoveSubStringFirst(directory, _AFLOWIN_);
    aflags.Directory = directory;
    if (aurostd::FileExist(directory + "/" + _AFLOWIN_) == false) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "File/Directory not found, nothing to do", _FILE_CORRUPT_); // CO20200624
    }
    for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_OUT + XHOST.vext[iext]) ||
          aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_OUT + XHOST.vext[iext]) ||
          aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_OUT + XHOST.vext[iext]) ||
          aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_OUT + XHOST.vext[iext])) {
        tocompress = true;
      }
    }
    for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
      // txt
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_OUT)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUP_OUT + ".*"));
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_OUT)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_OUT + ".*"));
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_OUT)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUP_XTAL_OUT + ".*"));
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_OUT)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_XTAL_OUT + ".*"));
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_OUT)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_FGROUP_OUT + ".*"));
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_OUT)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_AGROUP_OUT + ".*"));
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_OUT)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_SGROUP_OUT + ".*"));
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_OUT + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_OUT)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_IATOMS_OUT + ".*"));
      }
      // json
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_JSON)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUP_JSON + ".*"));
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_JSON)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_JSON + ".*"));
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_JSON)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUP_XTAL_JSON + ".*"));
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_JSON)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_PGROUPK_XTAL_JSON + ".*"));
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_JSON)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_FGROUP_JSON + ".*"));
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_JSON)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_AGROUP_JSON + ".*"));
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_JSON)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_SGROUP_JSON + ".*"));
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_JSON + XHOST.vext[iext]) || aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_JSON)) {
        aurostd::RemoveFile(directory, std::regex(DEFAULT_AFLOW_IATOMS_JSON + ".*"));
      }
    }
    // LOAD _AFLOWIN_
    ofstream FileMESSAGE("/dev/null");
    _kflags kflags;
    kflags.AFLOW_MODE_VASP = true;
    _vflags vflags;
    _xvasp xvasp;
    xvasp.clear();
    ifstream FileAFLOWIN(string(directory + string("/" + _AFLOWIN_)).c_str());
    aurostd::InFileExistCheck("pflow::EXTRACT_Symmetry", string(directory + string("/" + _AFLOWIN_)).c_str(), FileAFLOWIN);
    string AflowIn;
    AflowIn = "";
    char c;
    while (FileAFLOWIN.get(c)) {
      AflowIn += c; // READ _AFLOWIN_ and put into AflowIn
    }
    FileAFLOWIN.close();
    AflowIn = aurostd::RemoveComments(AflowIn); // NOW Clean AFLOWIN //CO20180502
    aflags.QUIET = true;
    XHOST.QUIET = true;
    vflags = KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn, aflags, kflags);
    KBIN::VASP_Produce_POSCAR(xvasp, AflowIn, FileMESSAGE, aflags, kflags, vflags);
    // CREATE PGROUP/FGROUP/SGROUP
    bool PGROUPWRITE = true;
    bool PGROUPKWRITE = true;
    bool FGROUPWRITE = true;
    bool SGROUPWRITE = false;
    bool IATOMSWRITE = true;
    bool AGROUPWRITE = true;
    const bool OSSWRITE = true; // to FileMESSAGE, does not matter as it is /dev/null
    SYM::CalculatePointGroup(FileMESSAGE, xvasp.str, aflags, PGROUPWRITE, OSSWRITE, cout);
    SYM::CalculatePointGroupKLattice(FileMESSAGE, xvasp.str, aflags, PGROUPKWRITE, OSSWRITE, cout);
    SYM::CalculateFactorGroup(FileMESSAGE, xvasp.str, aflags, FGROUPWRITE, OSSWRITE, cout);
    SYM::CalculateSpaceGroup(FileMESSAGE, xvasp.str, aflags, SGROUPWRITE, OSSWRITE, cout);
    SYM::CalculateInequivalentAtoms(FileMESSAGE, xvasp.str, aflags, IATOMSWRITE, OSSWRITE, cout);
    SYM::CalculateSitePointGroup(FileMESSAGE, xvasp.str, aflags, AGROUPWRITE, OSSWRITE, cout);
    // COMPRESS if necessary
    if (tocompress) {
      // txt
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_OUT)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUP_OUT);
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_OUT)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_OUT);
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_OUT)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_OUT);
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_OUT)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_OUT);
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_OUT)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_FGROUP_OUT);
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_OUT)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_SGROUP_OUT);
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_OUT)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_AGROUP_OUT);
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_OUT)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_IATOMS_OUT);
      }
      // json
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_JSON)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUP_JSON);
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_JSON)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_JSON);
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_JSON)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUP_XTAL_JSON);
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_JSON)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_PGROUPK_XTAL_JSON);
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_FGROUP_JSON)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_FGROUP_JSON);
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_SGROUP_JSON)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_SGROUP_JSON);
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_AGROUP_JSON)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_AGROUP_JSON);
      }
      if (aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_IATOMS_JSON)) {
        aurostd::CompressFile(directory + "/" + DEFAULT_AFLOW_IATOMS_JSON);
      }
    }
    return string("AFLOW " + string(AFLOW_VERSION) + " Symmetry Fixed in " + directory + "  (need aflow>=2948)\n");
  }
} // namespace pflow

// ***************************************************************************
// pflow::FGROUP
// ***************************************************************************
namespace pflow {
  void FGROUP(_aflags& aflags, istream& input) {
    cout << aflow::Banner("BANNER_TINY") << endl;
    aflags.QUIET = true;
    xstructure a(input, IOAFLOW_AUTO);
    const bool WRITE = true;
    ofstream FileMESSAGE("/dev/null");
    // DX20170815 - Add in consistency checks SYM::CalculatePointGroup(FileMESSAGE,a,aflags,WRITE,true,cout);
    // DX20170815 - Add in consistency checks SYM::CalculateFactorGroup(FileMESSAGE,a,aflags,WRITE,true,cout);
    _kflags kflags; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_FGROUP = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL = false; // DX20170815 - Add in consistency checks //DX20171205 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_CALCULATE_SGROUP = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_IATOMS = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_AGROUP = false; // DX20170815 - Add in consistency checks
    pflow::PerformFullSymmetry(a, FileMESSAGE, aflags, kflags, WRITE, cout); // DX20170815 - Add in consistency checks
  }
} // namespace pflow

// ***************************************************************************
// pflow::FINDSYM
// ***************************************************************************
namespace pflow {
  void FINDSYM(aurostd::xoption& vpflow, uint mode, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    string flag_name;
    if (mode == 0) {
      flag_name = "FINDSYM_PRINT";
    }
    if (mode == 1) {
      flag_name = "FINDSYM_EXEC";
    }
    const string options = vpflow.getattachedscheme(flag_name);
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() >= 2) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --findsym[_print][=tolerance_relative] < POSCAR");
    }
    // move on

    xstructure a(input, IOAFLOW_AUTO);
    double tolerance = DEFAULT_FINDSYM_TOL;
    if (vpflow.flag("SG::TOLERANCE")) {
      const string tolerance_string = vpflow.getattachedscheme("SG::TOLERANCE");
      vector<string> tol_tokens;
      aurostd::string2tokens(tolerance_string, tol_tokens, ",");
      if (tol_tokens.empty()) {
        tolerance = DEFAULT_FINDSYM_TOL;
      }
      if (tol_tokens.size() == 1) {
        tolerance = aurostd::string2utype<double>(tol_tokens[0]);
      }
    }
    // Read in input file.
    if (mode == 0) {
      cout << a.findsym2print(tolerance) << endl;
    }
    if (mode == 1) {
      cout << a.findsym2execute(tolerance) << endl;
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::FIXBANDS
// ***************************************************************************
namespace pflow {
  bool FIXBANDS(_aflags& aflags, string opts) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    stringstream aus;
    vector<string> kpoints_old;
    vector<string> eigenval_old;
    vector<string> tokens;
    vector<string> kpoints_new_string;
    stringstream kpoints_new;
    stringstream eigenval_new;
    const string directory = aflags.Directory;
    cout << XPID << "pflow::FIXBANDS: aflow --fix_bands=POSCAR,KPOINTS.bands.old,EIGENVAL.bands.old,KPOINTS.bands.new,EIGENVAL.bands.new" << endl;
    double grid = 0;
    uint points_bz = 0;
    uint eigenval_size = 0;
    uint nbands = 0;
    bool foundBZ;
    double pos4 = 0.0;

    vector<string> vopt;
    aurostd::string2tokens(opts, vopt, ",");
    uint vopt_counter = 0;

    if (vopt.size() < 5) {
      return false;
    }
    // FILE_POSCAR_IN ---------------------------------------
    const string FILE_POSCAR_IN = directory + "/" + vopt.at(vopt_counter++);
    if (aurostd::FileExist(FILE_POSCAR_IN) == false) {
      cout << XPID << "pflow::FIXBANDS: File POSCAR: " << FILE_POSCAR_IN << " not found, nothing to do" << endl;
      return false;
    } else {
      cout << XPID << "pflow::FIXBANDS: Load File POSCAR: " << FILE_POSCAR_IN << endl;
    }
    ifstream File_POSCAR_IN(FILE_POSCAR_IN.c_str());
    const xstructure str(File_POSCAR_IN, IOAFLOW_AUTO);
    File_POSCAR_IN.close();
    // FILE_KPOINTS_BANDS_OLD_IN ---------------------------------------
    const string FILE_KPOINTS_BANDS_OLD_IN = directory + "/" + vopt.at(vopt_counter++);
    if (aurostd::FileExist(FILE_KPOINTS_BANDS_OLD_IN) == false) {
      cout << XPID << "pflow::FIXBANDS: File KPOINTS_BANDS_OLD: " << FILE_KPOINTS_BANDS_OLD_IN << " not found, nothing to do" << endl;
      return false;
    } else {
      cout << XPID << "pflow::FIXBANDS: Load File KPOINTS_BANDS_OLD: " << FILE_KPOINTS_BANDS_OLD_IN << endl;
    }
    aurostd::string2tokens(aurostd::file2string(FILE_KPOINTS_BANDS_OLD_IN), kpoints_old, "\n");
    aus << kpoints_old.at(1);
    aus >> grid;
    // LATTICE ---------------------------------------
    aurostd::string2tokens(kpoints_old.at(0), tokens, " ");
    string LATTICE = tokens.at(0);
    if (LATTICE == "KPOINTS:" && tokens.size() > 1) {
      LATTICE = tokens[1];
    }
    // FILE_EIGENVAL_BANDS_OLD_IN ---------------------------------------
    const string FILE_EIGENVAL_BANDS_OLD_IN = directory + "/" + vopt.at(vopt_counter++);
    if (aurostd::FileExist(FILE_EIGENVAL_BANDS_OLD_IN) == false) {
      cout << XPID << "pflow::FIXBANDS: File EIGENVAL_BANDS_OLD: " << FILE_EIGENVAL_BANDS_OLD_IN << " not found, nothing to do" << endl;
      return false;
    } else {
      cout << XPID << "pflow::FIXBANDS: Load File EIGENVAL_BANDS_OLD: " << FILE_EIGENVAL_BANDS_OLD_IN << endl;
    }
    string tmp_eigenval = aurostd::file2string(FILE_EIGENVAL_BANDS_OLD_IN); // to avoid empty lines
    aurostd::StringSubstInPlace(tmp_eigenval, "\n", " \n"); // to avoid empty lines
    aurostd::string2tokens(tmp_eigenval, eigenval_old, "\n"); // to avoid empty lines
    // loaded eigenval now operate
    aurostd::string2tokens(eigenval_old.at(8), tokens, " ");
    eigenval_size = tokens.size() - 1;
    eigenval_size++;
    aurostd::string2tokens(eigenval_old.at(5), tokens, " ");
    points_bz = aurostd::string2utype<uint>(tokens.at(1));
    aurostd::string2tokens(eigenval_old.at(5), tokens, " ");
    nbands = aurostd::string2utype<uint>(tokens.at(2));
    aurostd::string2tokens(eigenval_old.at(7), tokens, " ");
    pos4 = aurostd::string2utype<double>(tokens.at(3));

    uint step = 7;
    vector<xvector<double>> vbzpoint;
    vector<xmatrix<double>> vbzband;

    for (uint ivbz = 0; ivbz < points_bz; ivbz++) {
      const xvector<double> bzpoint(3);
      aus.clear();
      aus.str(eigenval_old.at(step++));
      aus.precision(20);
      aus >> bzpoint[1] >> bzpoint[2] >> bzpoint[3] >> pos4;
      vbzpoint.push_back(bzpoint);
      const xmatrix<double> bzband(nbands, eigenval_size);
      for (uint ivbands = 0; ivbands < nbands; ivbands++) {
        aus.clear();
        aus.str(eigenval_old.at(step++));
        for (uint ieig = 0; ieig < eigenval_size; ieig++) {
          aus >> bzband[ivbands + 1][ieig + 1];
        }
      }
      vbzband.push_back(bzband);
      step++;
      //    cerr << step << endl;
    }
    for (uint ivbz = 0; ivbz < points_bz; ivbz++) {
      //     cerr << vbzpoint.at(ivbz) << endl;
    }

    ; //   cerr << vbzband.at(0) << endl;
    if (LDEBUG) {
      cerr << "DEBUG eigenval_size=" << eigenval_size << endl;
    }
    if (LDEBUG) {
      cerr << "DEBUG points_bz=" << points_bz << endl;
    }
    if (LDEBUG) {
      cerr << "DEBUG nbands=" << nbands << endl;
    }
    if (LDEBUG) {
      cerr << "DEBUG pos4=" << pos4 << endl;
    }

    cout << XPID << "pflow::FIXBANDS: LATTICE=" << LATTICE << endl;
    string tmp_kpoints = LATTICE::KPOINTS_Directions(LATTICE, str.lattice, grid, str.iomode, foundBZ); // to avoid empty lines
    aurostd::StringSubstInPlace(tmp_kpoints, "\n", " \n"); // to avoid empty lines
    aurostd::string2tokens(tmp_kpoints, kpoints_new_string, "\n"); // to avoid empty lines
    // loaded kpoints now operate
    cout << XPID << "pflow::FIXBANDS: PATH= " << kpoints_new_string.at(0) << endl;
    kpoints_new.clear();
    kpoints_new.str(std::string());
    for (uint i = 0; i <= 3; i++) {
      kpoints_new << kpoints_new_string.at(i) << endl;
    }
    eigenval_new.clear();
    eigenval_new.str(std::string());
    for (uint i = 0; i <= 4; i++) {
      eigenval_new << eigenval_old.at(i) << endl;
    }

    // Patch for EIGENVAL.new (KY adds it)
    int kpoints_NEW;
    kpoints_NEW = (kpoints_new_string.size() - 3) * grid / 3;
    eigenval_new << "   " << int(grid) << "  " << kpoints_NEW << "  " << nbands << endl;
    // eigenval_new << " " << endl; //Comment this, move the blank line before the bands data,
    // Make the total lines of EIGENVAL.aflow exactly same with EIGENVAL.vasp, KY

    for (size_t ikpz = 4; ikpz < kpoints_new_string.size();) {
      aurostd::string2tokens(kpoints_new_string[ikpz], tokens);
      if (tokens.size() < 2) {
        kpoints_new << kpoints_new_string.at(ikpz++) << endl;
      } else {
        // FROM
        xvector<double> kpoint_from(3);
        string string_from;
        aus.clear();
        aus.str(kpoints_new_string[ikpz]);
        string_from = kpoints_new_string.at(ikpz++);
        aus >> kpoint_from[1] >> kpoint_from[2] >> kpoint_from[3];
        // TO
        const xvector<double> kpoint_to(3);
        string string_to;
        aus.clear();
        aus.str(kpoints_new_string[ikpz]);
        string_to = kpoints_new_string.at(ikpz++);
        aus >> kpoint_to[1] >> kpoint_to[2] >> kpoint_to[3];
        // DELTA
        xvector<double> kpoint_delta(3);
        kpoint_delta = (kpoint_to - kpoint_from) / (grid - 1.0);
        const double kpoins_delta_prec = 1.0 * modulus(kpoint_delta);
        for (uint jkpz = 0; jkpz < (uint) grid; jkpz++) {
          // generated and scan
          uint index = 0;
          double err = 1.0e6;
          for (size_t ivbz = 0; ivbz < vbzpoint.size(); ivbz++) {
            if (modulus(vbzpoint[ivbz] - kpoint_from) < err) {
              index = ivbz;
              err = modulus(vbzpoint[ivbz] - kpoint_from);
            }
          }
          if (err > kpoins_delta_prec) {
            cerr << "************************************************" << endl;
            cerr << "kpoint_from: not found (" << kpoint_from << ")" << endl;
            cerr << "closest is vbzpoint.at(index)=" << vbzpoint.at(index) << "   with err=" << err << endl;
            cerr << "************************************************" << endl;
            return false;
          }
          // in index there is the good point

          eigenval_new << endl; // Put the blank line before the bands data, KY
          eigenval_new.precision(7);
          eigenval_new.setf(std::ios::scientific, std::ios::floatfield);
          eigenval_new << "  " << vbzpoint.at(index)[1] << "  " << vbzpoint.at(index)[2] << "  " << vbzpoint.at(index)[3] << "  " << pos4 << endl;
          for (uint ivbands = 1; ivbands <= nbands; ivbands++) {
            eigenval_new << "   " << (int) vbzband.at(index)[ivbands][1];
            eigenval_new.precision(5);
            eigenval_new.setf(std::ios::fixed, std::ios::floatfield);
            for (uint ieig = 2; ieig <= eigenval_size; ieig++) {
              eigenval_new << "   " << vbzband.at(index)[ivbands][ieig];
            }
            eigenval_new << endl;
          }
          // eigenval_new << endl;  //Remove the additional blank line, Put the blank line before the bands data, KY
          //  new one
          kpoint_from = kpoint_from + kpoint_delta;
        }
        // if OK
        kpoints_new << string_from << endl;
        kpoints_new << string_to << endl;
      }
    }
    kpoints_new << endl;
    // FILE_KPOINTS_BANDS_NEW_OUT ---------------------------------------
    const string FILE_KPOINTS_BANDS_NEW_OUT = directory + "/" + vopt.at(vopt_counter++);
    aurostd::stringstream2file(kpoints_new, FILE_KPOINTS_BANDS_NEW_OUT);
    cout << XPID << "pflow::FIXBANDS: Save File KPOINTS_BANDS_NEW: " << FILE_KPOINTS_BANDS_NEW_OUT << endl;
    // FILE_EIGENVAL_BANDS_NEW_OUT ---------------------------------------
    const string FILE_EIGENVAL_BANDS_NEW_OUT = directory + "/" + vopt.at(vopt_counter++);
    aurostd::stringstream2file(eigenval_new, FILE_EIGENVAL_BANDS_NEW_OUT);
    cout << XPID << "pflow::FIXBANDS: Save File EIGENVAL_BANDS_NEW: " << FILE_EIGENVAL_BANDS_NEW_OUT << endl;

    return true;
  }
} // namespace pflow

// ***************************************************************************
// pflow::FRAC
// ***************************************************************************
namespace pflow {
  xstructure FRAC(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    xstructure b(a);
    b.SetCoordinates(_COORDS_FRACTIONAL_);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::FROZSL_ANALYZE
// ***************************************************************************
namespace pflow {
  string FROZSL_ANALYZE(istream& input) {
    ostringstream oss;
    oss.clear();
    oss.setf(std::ios::fixed, std::ios::floatfield);
    const uint _precision_ = 10; // was 16 SC 10 DM
    oss.precision(_precision_);

    vector<string> vinput;
    vector<string> tokens;
    vector<string> Mrefs;
    vector<double> Erefs;
    aurostd::stream2vectorstring(input, tokens);
    // take the good ones
    for (size_t i = 0; i < tokens.size(); i++) {
      if (aurostd::substring2bool(tokens[i], ":")) {
        vinput.push_back(tokens[i]);
      }
    }
    // find m0a0
    for (size_t i = 0; i < vinput.size(); i++) {
      if (aurostd::substring2bool(vinput[i], "m0a0")) {
        aurostd::string2tokens(vinput[i], tokens, "m0a0");
        Mrefs.push_back(tokens.at(0) + "m");
        aurostd::string2tokens(vinput[i], tokens, "=");
        Erefs.push_back(aurostd::string2utype<double>(tokens.at(tokens.size() - 1)));
      }
    }
    //  for(size_t i=0;i<Mrefs.size();i++) cout << Mrefs.at(i) << " " << Erefs.at(i) << endl;
    // now print
    const double frozsl_eps = 1.0e-6;
    const double ev2hartree = 1.0 / 27.211383;
    for (size_t i = 0; i < vinput.size(); i++) {
      for (size_t j = 0; j < Mrefs.size(); j++) {
        if (aurostd::substring2bool(vinput[i], Mrefs[j])) {
          // cout << Mrefs.at(j) << " " << Erefs.at(j) << endl;
          aurostd::string2tokens(vinput[i], tokens, "=");
          if (std::abs(aurostd::string2utype<double>(tokens.at(tokens.size() - 1)) - Erefs.at(j)) > frozsl_eps) {
            oss << aurostd::PaddedPOST(aurostd::utype2string(1000 * ev2hartree * (aurostd::string2utype<double>(tokens.at(tokens.size() - 1)) - Erefs.at(j))), 23, " ") << " ! ";
            aurostd::string2tokens(vinput[i], tokens);
            oss << tokens.at(0) << " - " << Mrefs[j] << "0a0" << endl;
          }
        }
      }
    }
    return oss.str();
  }
} // namespace pflow

namespace pflow {
  string FROZSL_ANALYZE_old(istream& input) {
    // this one was looking for the minimum
    ostringstream oss;
    oss.clear();
    oss.setf(std::ios::fixed, std::ios::floatfield);
    const uint _precision_ = 10; // was 16 SC 10 DM
    oss.precision(_precision_);

    vector<string> vinput;
    vector<string> tokens;
    aurostd::stream2vectorstring(input, tokens);
    // take the good ones
    for (size_t i = 0; i < tokens.size(); i++) {
      if (aurostd::substring2bool(tokens[i], ":")) {
        vinput.push_back(tokens[i]);
      }
    }
    // find minimum
    double minE = 1.0e6;
    for (size_t i = 0; i < vinput.size(); i++) {
      aurostd::string2tokens(vinput[i], tokens, "=");
      minE = min(minE, aurostd::string2utype<double>(tokens.at(tokens.size() - 1)));
    }
    const double frozsl_eps = 1.0e-6;
    const double ev2hartree = 1.0 / 27.211383;
    for (size_t i = 0; i < vinput.size(); i++) {
      aurostd::string2tokens(vinput[i], tokens, "=");
      if (std::abs(aurostd::string2utype<double>(tokens.at(tokens.size() - 1)) - minE) > frozsl_eps) {
        oss << 1000 * ev2hartree * (aurostd::string2utype<double>(tokens.at(tokens.size() - 1)) - minE) << "   ! ";
        aurostd::string2tokens(vinput[i], tokens);
        oss << tokens.at(0) << endl;
      }
    }
    return oss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::FROZSL_INPUT
// ***************************************************************************
namespace pflow {
  string FROZSL_INPUT() {
    _aflags aflags;
    aflags.Directory = "./";
    _kflags kflags;
    ofstream oss;
    string aflowin;
    string MESSAGE = "pflow::FROZSL_INPUT ERROR";
    aflowin = string(aflags.Directory + _AFLOWIN_);
    if (!aurostd::FileExist(aflowin)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "file not found: " + aflowin, _FILE_CORRUPT_); // CO20200624
    }

    string AflowIn;
    aurostd::file2string(aflowin, AflowIn);
    FROZSL::WGET_INPUT(oss, AflowIn, aflags, kflags);
    // cout << oss << endl;
    return "";
  }
} // namespace pflow

// ***************************************************************************
// pflow::FROZSL_OUTPUT
// ***************************************************************************
namespace pflow {
  string FROZSL_OUTPUT() {
    _aflags aflags;
    aflags.Directory = "./";
    _kflags kflags;
    ofstream oss;
    string aflowin;
    string MESSAGE = "pflow::FROZSL_OUTPUT ERROR";
    aflowin = string(aflags.Directory + _AFLOWIN_);
    if (!aurostd::FileExist(aflowin)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "file not found: " + aflowin, _FILE_CORRUPT_); // CO20200624
    }
    FROZSL::WGET_OUTPUT(oss, aflags, kflags);
    stringstream sss;
    aurostd::file2stringstream(aflags.Directory + "/" + DEFAULT_AFLOW_FROZSL_EIGEN_OUT, sss);
    //  cerr << sss.str() << endl;
    return "";
  }
} // namespace pflow

// ***************************************************************************
// pflow::FROZSL_VASPSETUP
// ***************************************************************************
namespace pflow {
  string FROZSL_VASPSETUP(vector<string> argv, int mode) {
    string strout;
    const string strfile = _FROZSL_VASPSETUP_FILE_;
    // cerr << "HERE argv.size()=" << argv.size() << endl;

    if (mode == 0) {
      strout = aurostd::EmbData::get_content("phvaspsetup_AFLOW", "PHVASP");
    }
    if (mode == 1) {
      strout = aurostd::EmbData::get_content("phvaspsetup_POSCAR", "PHVASP");
    }
    if (argv.size() == 2) {
      aurostd::string2file(strout, strfile);
      // cerr << "HERE argv.size()=" << argv.size() << endl;
      aurostd::Chmod(0755, strfile);
      strout = aurostd::execute2string(strfile);
      aurostd::RemoveFile(strfile);
      return strout;
    }
    return strout;
  }
} // namespace pflow

// ***************************************************************************
// pflow::GEOMETRY
// ***************************************************************************
namespace pflow {
  string GEOMETRY(istream& input) { // CO20191110
    xstructure a(input, IOAFLOW_AUTO);
    a.ReScale(1.0); // DX+ME20201229 - included scaling in a,b,c,alpha,beta,gamma
    return aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(Getabc_angles(a.lattice, DEGREES), 5, true), ","); // roundoff
  }
} // namespace pflow

// DX20170927 - get collinear magnetic info - START
//  ***************************************************************************
//  pflow::GetCollinearMagneticInfo
//  ***************************************************************************
namespace pflow {
  bool GetCollinearMagneticInfo(uint num_atoms, const string& magmom_info, vector<double>& vmag) { // DX20191107 - int to uint
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (aurostd::substring2bool(magmom_info, "OUTCAR")) {
      if (aurostd::FileExist(magmom_info)) {
        xOUTCAR outcar;
        outcar.GetPropertiesFile(magmom_info);
        for (size_t i = 0; i < outcar.vmag.size(); i++) {
          vmag.push_back(outcar.vmag[i]);
        }
      } else {
        cerr << XPID << "pflow::GetCollinearMagneticInfo: ERROR: OUTCAR file does not exist." << endl;
        return false;
      }
    } else if (aurostd::substring2bool(magmom_info, "INCAR")) {
      if (aurostd::FileExist(magmom_info)) {
        stringstream sss;
        aurostd::compressfile2stringstream(magmom_info, sss);
        vector<string> vcontent;
        aurostd::string2vectorstring(sss.str(), vcontent);
        bool magmom_found = false;
        for (size_t i = 0; i < vcontent.size(); i++) {
          if (vcontent[i].find("MAGMOM=") != std::string::npos) {
            vcontent[i] = aurostd::RemoveComments(vcontent[i]); // DX20210517 - better way to check/remove comments
            if (vcontent[i].find("MAGMOM=") == std::string::npos) {
              cerr << XPID << "pflow::GetCollinearMagneticInfo: ERROR: MAGMOM tag is commented out." << endl;
              return false;
            } else {
              const string magmom_values = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(aurostd::RemoveCharacter(aurostd::RemoveSubString(vcontent[i], "MAGMOM="), '\n')); // DX20210517 - need to pre-condition string
              vector<string> mag_tokens;
              aurostd::string2tokens(magmom_values, mag_tokens, " "); // DX20210517 - need to split on spaces
              for (size_t m = 0; m < mag_tokens.size(); m++) {
                // INCAR allows multiplication of elements to describe magnetic moment (i.e. 2*2.0)
                if (mag_tokens[m].find("*") != std::string::npos) {
                  vector<string> tokens;
                  aurostd::string2tokens(mag_tokens[m], tokens, "*");
                  const uint num_magmoms = aurostd::string2utype<uint>(tokens[0]); // DX20191107 - int to uint
                  for (uint n = 0; n < num_magmoms; n++) {
                    vmag.push_back(aurostd::string2utype<double>(tokens[1]));
                  }
                } else if (!mag_tokens[m].empty()) { // DX20210517 - ensure the token is not empty
                  vmag.push_back(aurostd::string2utype<double>(mag_tokens[m]));
                }
              }
              magmom_found = true;
            }
          }
        }
        if (!magmom_found) {
          cerr << XPID << "pflow::GetCollinearMagneticInfo: ERROR: MAGMOM tag was not found in the INCAR." << endl;
          return false;
        }
      } else {
        cerr << XPID << "pflow::GetCollinearMagneticInfo: ERROR: INCAR file does not exist." << endl;
        return false;
      }
    } else {
      vector<string> mag_tokens;
      aurostd::string2tokens(magmom_info, mag_tokens, ",");
      for (size_t m = 0; m < mag_tokens.size(); m++) {
        // Allows multiplication of elements to describe magnetic moment (i.e. 2*2.0)
        if (mag_tokens[m].find("*") != std::string::npos) {
          vector<string> tokens;
          aurostd::string2tokens(mag_tokens[m], tokens, "*");
          const uint num_magmoms = aurostd::string2utype<uint>(tokens[0]); // DX20191107 - int to uint
          for (uint n = 0; n < num_magmoms; n++) {
            vmag.push_back(aurostd::string2utype<double>(tokens[1]));
          }
        } else {
          vmag.push_back(aurostd::string2utype<double>(mag_tokens[m]));
        }
      }
    }
    if (vmag.size() != num_atoms) { // DX20191107 - remove (int) type casting
      if (LDEBUG) {
        cerr << XPID << "pflow::GetCollinearMagneticInfo: WARNING: Number of magnetic moments is not equivalent to the number of atoms." << endl;
      }
      return false;
    }
    return true;
  }
} // namespace pflow
// DX20170927 - get collinear magnetic info - END

// DX20171205 - get non-collinear magnetic info - START
//  ***************************************************************************
//  pflow::GetNonCollinearMagneticInfo
//  ***************************************************************************
namespace pflow {
  bool GetNonCollinearMagneticInfo(uint num_atoms, const string& magmom_info, vector<xvector<double>>& vmag_noncoll) { // DX20191107 - int to uint
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (aurostd::substring2bool(magmom_info, "OUTCAR")) {
      if (aurostd::FileExist(magmom_info)) {
        xOUTCAR outcar;
        outcar.GetPropertiesFile(magmom_info);
        for (size_t i = 0; i < outcar.vmag_noncoll.size(); i++) {
          vmag_noncoll.push_back(outcar.vmag_noncoll[i]);
        }
      } else {
        cerr << XPID << "pflow::GetNonCollinearMagneticInfo: ERROR: OUTCAR file does not exist." << endl;
        return false;
      }
    } else if (aurostd::substring2bool(magmom_info, "INCAR")) {
      if (aurostd::FileExist(magmom_info)) {
        stringstream sss;
        aurostd::compressfile2stringstream(magmom_info, sss);
        vector<string> vcontent;
        aurostd::string2vectorstring(sss.str(), vcontent);
        bool magmom_found = false;
        for (size_t i = 0; i < vcontent.size(); i++) {
          if (vcontent[i].find("MAGMOM=") != std::string::npos) {
            vcontent[i] = aurostd::RemoveComments(vcontent[i]); // DX20210517 - better way to check/remove comments
            if (vcontent[i].find("MAGMOM=") == std::string::npos) {
              cerr << XPID << "pflow::GetNonCollinearMagneticInfo: ERROR: MAGMOM tag is commented out." << endl;
              return false;
            } else {
              const string magmom_values = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(aurostd::RemoveCharacter(aurostd::RemoveSubString(vcontent[i], "MAGMOM="), '\n')); // DX20210517 - need to pre-condition string
              vector<string> mag_tokens;
              aurostd::string2tokens(magmom_values, mag_tokens, " "); // DX20210517 - need to split on spaces
              // INCAR allows multiplication of elements to describe magnetic moment (i.e. 2*2.0)
              vector<double> all_magmom_tokens;
              for (size_t m = 0; m < mag_tokens.size(); m++) {
                if (mag_tokens[m].find("*") != std::string::npos) {
                  vector<string> tokens;
                  aurostd::string2tokens(mag_tokens[m], tokens, "*");
                  const uint num_magmoms = aurostd::string2utype<uint>(tokens[0]); // DX20191107 - int to uint
                  for (uint n = 0; n < num_magmoms; n++) {
                    all_magmom_tokens.push_back(aurostd::string2utype<double>(tokens[1]));
                  }
                } else if (!mag_tokens[m].empty()) { // DX20210517 - ensure the token is not empty
                  all_magmom_tokens.push_back(aurostd::string2utype<double>(mag_tokens[m]));
                }
              }
              // non-collinear check (should be divisible by 3)
              if (all_magmom_tokens.size() != 3 * num_atoms) { // DX20191107 - removed (int) type casting
                if (LDEBUG) {
                  cerr << XPID << "pflow::GetNonCollinearMagneticInfo: WARNING: From INCAR. Number of magnetic moments not divisible by 3; not non-collinear system." << endl;
                }
                return false;
              }
              const xvector<double> mag_xyz;
              uint index = 1;
              for (size_t m = 0; m < all_magmom_tokens.size(); m++) {
                mag_xyz(index) = all_magmom_tokens[m];
                index++;
                if (index == 4) {
                  vmag_noncoll.push_back(mag_xyz);
                  index = 1;
                }
              }
              magmom_found = true;
            }
          }
        }
        if (!magmom_found) {
          cerr << XPID << "pflow::GetNonCollinearMagneticInfo: ERROR: MAGMOM tag was not found in the INCAR." << endl;
          return false;
        }
      } else {
        cerr << XPID << "pflow::GetNonCollinearMagneticInfo: ERROR: INCAR file does not exist." << endl;
        return false;
      }
    } else {
      vector<string> mag_tokens;
      aurostd::string2tokens(magmom_info, mag_tokens, ",");
      // Allows multiplication of elements to describe magnetic moment (i.e. 2*2.0)
      vector<double> all_magmom_tokens;
      for (size_t m = 0; m < mag_tokens.size(); m++) {
        if (mag_tokens[m].find("*") != std::string::npos) {
          vector<string> tokens;
          aurostd::string2tokens(mag_tokens[m], tokens, "*");
          const uint num_magmoms = aurostd::string2utype<uint>(tokens[0]); // DX20191107 - int to uint
          for (uint n = 0; n < num_magmoms; n++) {
            all_magmom_tokens.push_back(aurostd::string2utype<double>(tokens[1]));
          }
        } else {
          all_magmom_tokens.push_back(aurostd::string2utype<double>(mag_tokens[m]));
        }
      }
      // non-collinear check (should be divisible by 3)
      if (all_magmom_tokens.size() != 3 * num_atoms) { // DX20191107 - removed (int) type casting
        if (LDEBUG) {
          cerr << XPID << "pflow::GetNonCollinearMagneticInfo: WARNING: From manual input. Number of magnetic moments is not three times the number of atoms; not non-collinear system." << endl;
        }
        return false;
      }
      const xvector<double> mag_xyz;
      uint index = 1;
      for (size_t m = 0; m < all_magmom_tokens.size(); m++) {
        mag_xyz(index) = all_magmom_tokens[m];
        index++;
        if (index == 4) {
          vmag_noncoll.push_back(mag_xyz);
          index = 1;
        }
      }
    }
    return true;
  }
} // namespace pflow
// DX20170927 - get collinear magnetic info - END

// ***************************************************************************
// pflow::ATOMIC_ENVIRONMENT
// ***************************************************************************
namespace pflow { // HE20210331
  void ATOMIC_ENVIRONMENT(const aurostd::xoption& vpflow) {
    string auid;
    uint ae_mode = ATOM_ENVIRONMENT_MODE_1;
    double radius = 4.0;
    if (vpflow.flag("ATOMIC_ENVIRONMENT::AUID")) {
      auid = vpflow.getattachedscheme("ATOMIC_ENVIRONMENT::AUID");
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "missing auid - could not load structure", _INPUT_ERROR_);
    }
    if (vpflow.flag("ATOMIC_ENVIRONMENT::MODE")) {
      ae_mode = vpflow.getattachedutype<uint>("ATOMIC_ENVIRONMENT::MODE");
    }
    if (vpflow.flag("ATOMIC_ENVIRONMENT::RADIUS")) {
      radius = vpflow.getattachedutype<double>("ATOMIC_ENVIRONMENT::RADIUS");
    }
    pflow::outputAtomicEnvironment(auid, ae_mode, radius);
  }
} // namespace pflow

// ***************************************************************************
// pflow::GLASS_FORMING_ABILITY
// ***************************************************************************
namespace pflow { // DF20190329
  void GLASS_FORMING_ABILITY(aurostd::xoption& vpflow) { // DF20190329

    const string alloy = vpflow.getattachedscheme("PFLOW::ALLOY");
    string AE_file_read = vpflow.getattachedscheme("GFA::AE_FILE");
    //[CO20190628 - AFLOWRC default]double fe_cut = aurostd::string2utype<double>(vpflow.getattachedscheme("GFA::FORMATION_ENTHALPY_CUTOFF")); //DF20190619
    double fe_cut = DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF;

    if (!vpflow.flag("PFLOW::ALLOY")) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Must specify an alloy system; e.g.: --alloy=CaMg", _INPUT_MISSING_); // CO20200624
    }

    cout << endl << "Calculating glass-forming ability for " << alloy << "." << endl;

    if (!vpflow.flag("GFA::AE_FILE")) {
      AE_file_read = "none";
      cout << endl
           << "No AE input file? That's ok, will calculate the atomic environments." << endl
           << "alternatively, set --ae_file= to specify a file of atomic environments to read; e.g.: --ae_file=All_atomic_environments_read.dat" << endl;
    } else if (!aurostd::FileExist(AE_file_read)) {
      cout << endl << "Can't find " << AE_file_read << ". Will calculate all atomic environments." << endl;
      AE_file_read = "none";
    } else {
      cout << endl << "Reading atomic environments from file: " << AE_file_read << "." << endl;
    }

    //[CO20190628 - AFLOWRC default]if(!vpflow.flag("GFA::FORMATION_ENTHALPY_CUTOFF")) {  //DF20190619
    //[CO20190628 - AFLOWRC default]  fe_cut = 0.05;
    //[CO20190628 - AFLOWRC default]}
    if (vpflow.flag("GFA::FORMATION_ENTHALPY_CUTOFF")) {
      fe_cut = aurostd::string2utype<double>(vpflow.getattachedscheme("GFA::FORMATION_ENTHALPY_CUTOFF"));
    } // CO20190628
    cout << endl << "Using " << fe_cut << " eV formation enthalpy cutoff for structures to include in the analysis." << endl;

    pflow::CalculateGFA(vpflow, alloy, AE_file_read, fe_cut); // DF20190619
  }
} // namespace pflow

// ***************************************************************************
// pflow::GULP
// ***************************************************************************
namespace pflow {
  void GULP(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    pflow::PrintGulp(a, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::HKL
// ***************************************************************************
namespace pflow {
  void HKL(const string& options, _aflags& aflags, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 3 && tokens.size() != 4) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --hkl=h,k,l[,bond] < POSCAR");
    }
    const xstructure a(input, IOAFLOW_AUTO);
    if (LDEBUG) {
      cerr << XPID << "pflow::HKL: a=" << a << endl;
    }
    vector<vector<double>> planesreducible;
    vector<vector<double>> planesirreducible;

    if (LDEBUG) {
      cerr << XPID << "pflow::HKL: mode=" << tokens.size() << endl;
    }
    const xvector<double> iparams(tokens.size(), 1);
    if (!tokens.empty()) {
      iparams(1) = aurostd::string2utype<double>(tokens[0]);
    }
    if (tokens.size() > 1) {
      iparams(2) = aurostd::string2utype<double>(tokens[1]);
    }
    if (tokens.size() > 2) {
      iparams(3) = aurostd::string2utype<double>(tokens[2]);
    }
    if (tokens.size() > 3) {
      iparams(4) = aurostd::string2utype<double>(tokens[3]);
    }
    if (LDEBUG) {
      cerr << XPID << "pflow::HKL: iparams=" << iparams << endl;
    }

    surface::GetSurfaceHKL(a, aflags, iparams, planesreducible, planesirreducible, cout);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }

} // namespace pflow

// ***************************************************************************
// pflow::HKLSearchSearch Trivial/Simple/Complete
// ***************************************************************************
namespace pflow {
  void HKLSearch(const string& options, _aflags& aflags, istream& input, const string& smode) {
    const bool LDEBUG = true; //(false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");

    if (LDEBUG) {
      cerr << XPID << "pflow::HKLSearch: smode=" << smode << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "pflow::HKLSearch: tokens.size()=" << tokens.size() << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "pflow::HKLSearch: options=" << options << endl;
    }

    if (smode == "HKL_SEARCH_TRIVIAL" && tokens.size() > 3) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --hkl_search[=khlmax[,bond[,step]]] < POSCAR");
    }
    if (smode == "HKL_SEARCH_SIMPLE" && tokens.size() > 4) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --hkl_search_simple[=cutoff[,bond[,khlmax[,step]]]] < POSCAR");
    }
    if (smode == "HKL_SEARCH_COMPLETE" && tokens.size() > 4) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --hkl_search_complete[=cutoff[,bond[,khlmax[,step]]]] < POSCAR");
    }
    const xstructure a(input, IOAFLOW_AUTO);
    //  if(LDEBUG) cerr << XPID << "pflow::HKL: a=" << a << endl;
    vector<vector<double>> planesreducible;
    vector<vector<double>> planesirreducible;
    vector<vector<uint>> planesirreducible_images;

    if (LDEBUG) {
      cerr << XPID << "pflow::HKLSearch: tokens.size()=" << tokens.size() << endl;
    }
    const xvector<double> iparams(tokens.size(), (!tokens.empty()));
    //    cerr << iparams.lrows << endl;
    // cerr << iparams.urows << endl;
    if (LDEBUG) {
      cerr << XPID << "pflow::HKLSearch: iparams.rows=" << iparams.rows << endl;
    }

    if (!tokens.empty()) {
      iparams(1) = aurostd::string2utype<double>(tokens[0]);
    }
    if (tokens.size() > 1) {
      iparams(2) = aurostd::string2utype<double>(tokens[1]);
    }
    if (tokens.size() > 2) {
      iparams(3) = aurostd::string2utype<double>(tokens[2]);
    }
    if (tokens.size() > 3) {
      iparams(4) = aurostd::string2utype<double>(tokens[3]);
    }
    if (LDEBUG) {
      cerr << XPID << "pflow::HKLSearch: iparams=" << iparams << endl;
    }

    surface::GetSurfaceHKLSearch(a, aflags, iparams, planesreducible, planesirreducible, planesirreducible_images, cout, smode);
    cout << "REDUCIBLE_SIZE " << planesreducible.size() << endl;
    cout << "IRREDUCIBLE_SIZE " << planesirreducible.size() << endl;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// Function pflow::HNF
// ***************************************************************************
namespace pflow {
  bool setPOCCTOL(xstructure& xstr, const string& pocc_tol_string) { // CO20181226
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
      cerr << __AFLOW_FUNC__ << " pocc_tol_string=" << pocc_tol_string << endl;
    }
    if (pocc_tol_string.empty()) {
      return false;
    }
    stringstream message;
    if (aurostd::substring2bool(pocc_tol_string, ",")) {
      message << "Cannot handle more than one pocc_tol specification";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    vector<string> tokens;
    aurostd::string2tokens(pocc_tol_string, tokens, ":");
    double tol;
    if (!tokens.empty()) {
      if (!aurostd::isfloat(tokens[0])) {
        message << "Input is not a float [input=" << tokens[0] << "]";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      tol = aurostd::string2utype<double>(tokens[0]);
      if (aurostd::isequal(tol, 0.0, _AFLOW_POCC_ZERO_TOL_)) {
        message << "Tolerance is too small (0) [input=" << tokens[0] << "]";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      if (std::signbit(tol)) {
        xstr.neg_scale_second = true;
        xstr.partial_occupation_HNF = -1 * (int) tol;
      } else {
        // defaults both tols, stoich_tol to be overwritten if another tolerance is provided
        xstr.partial_occupation_site_tol = xstr.partial_occupation_stoich_tol = tol;
      }
    }
    if (tokens.size() > 1) {
      if (!aurostd::isfloat(tokens[1])) {
        message << "Input is not a float [input=" << tokens[1] << "]";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      tol = aurostd::string2utype<double>(tokens[1]);
      if (aurostd::isequal(tol, 0.0, _AFLOW_POCC_ZERO_TOL_)) {
        message << "Tolerance is too small (0) [input=" << tokens[1] << "]";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      if (std::signbit(tol)) {
        message << "Tolerance cannot be negative [input=" << tokens[1] << "]";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      xstr.scale_third.isentry = true;
      xstr.partial_occupation_stoich_tol = tol;
    }
    return true;
  }
} // namespace pflow

// ***************************************************************************
// Function pflow::POCC_COMMAND_LINE
// ***************************************************************************
namespace pflow {
  bool POCC_COMMAND_LINE(aurostd::xoption& vpflow, istream& input, ostream& oss) { // CO20181226
    // ME20210429 - Do not print banner for web output
    if (!XHOST.vflag_control.flag("WWW")) {
      oss << aflow::Banner("BANNER_NORMAL");
    }
    xstructure xstr(input, IOAFLOW_AUTO);
    if (vpflow.flag("POCC_TOL")) {
      setPOCCTOL(xstr, vpflow.getattachedscheme("POCC_TOL"));
    }
    pocc::POccCalculator pcalc(xstr, vpflow, oss); // CO20190319
    pcalc.m_count_unique_fast = vpflow.flag("POCC_COUNT_UNIQUE_FAST");
    if (!pcalc.m_initialized) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "POccCalculator failed to initialized");
    }
    pcalc.calculateHNF();
    if (vpflow.flag("HNF")) {
      return pcalc.m_initialized;
    }
    pcalc.getTotalPermutationsCount();
    if (vpflow.flag("POCC_COUNT_TOTAL")) {
      return pcalc.m_initialized;
    }
    pcalc.calculate();
    if (vpflow.flag("POCC_COUNT_UNIQUE") || vpflow.flag("POCC_COUNT_UNIQUE_FAST")) {
      return pcalc.m_initialized;
    }
    oss << AFLOWIN_SEPARATION_LINE << endl;
    oss << "Creating list of unique derivative supercells." << endl; // CO20190116
    for (unsigned long long int i = 0; i < pcalc.getUniqueSuperCellsCount(); i++) {
      // populate POCC_UNIQUE_DERIVATIVE_STRUCTURES_FILE
      oss << AFLOWIN_SEPARATION_LINE << endl;
      oss << _VASP_POSCAR_MODE_EXPLICIT_START_ << endl; // ." << ss_pocc_count.str() << endl;
      oss << pcalc.getUniqueSuperCell(i);
      oss << _VASP_POSCAR_MODE_EXPLICIT_STOP_ << endl; // ." << ss_pocc_count.str() << endl;
      oss << AFLOWIN_SEPARATION_LINE << endl;
    }
    return pcalc.m_initialized;
  }
} // namespace pflow

namespace pflow {
  string HNF(vector<string> argv, istream& input, ostream& oss) {
    oss << AFLOWIN_SEPARATION_LINE << endl; // --------------------------------
    oss << "[AFLOW] " << aflow::Banner("BANNER_TINY") << endl;
    //  oss << "argv.size()=" << argv.size() << endl;
    xstructure str(input, IOAFLOW_AUTO);
    str.BringInCell();
    //  str.partial_occupation_sum=0;
    // GET SYMMETRY
    str.CalculateSymmetryPointGroup();
    //  str.CalculateSymmetryFactorGroup();
    oss << "[AFLOW] CALCULATION OF HNF xstructures" << endl;
    oss << "[AFLOW] pgroup operations = " << str.pgroup.size() << endl;
    oss.flush();
    // oss << "fgroup operations  = " << str.fgroup.size() << endl;
    // oss << "pgroupk operations = " << str.pgroupk.size() << endl;

    int choice = 2;
    if (argv.size() >= 3) {
      choice = aurostd::args2utype(argv, "--hnf", 2);
    }
    vector<int> vn;
    if (choice > 0) {
      vn.push_back(choice);
    }
    if (choice < 0) {
      for (int i = 2; i <= -choice; i++) {
        vn.push_back(i);
      }
    }
    if (choice == 0) {
      oss << "[AFLOW] pflow::HNF n=0" << endl << AFLOWIN_SEPARATION_LINE << endl;
    }

    vector<xmatrix<double>> sHNF; // contains only the supercells HNF that are uniques
    for (size_t in = 0; in < vn.size(); in++) {
      const int n = vn[in];
      oss << AFLOWIN_SEPARATION_LINE << endl; // --------------------------------
      oss << "[AFLOW] n = " << n << endl;
      sHNF = CalculateHNF(str, n);
      // PHYSICAL REVIEW B 77, 224115 2008
      // Algorithm for generating derivative structures
      // Gus L. W. Hart and Rodney W. Forcade
      // xmatrix<double> HNF(3,3);
      // vector<xmatrix<double> > vHNF;
      // xmatrix<double> A(3,3),B(3,3),Bi(3,3),Bj(3,3),R(3,3),H(3,3);A=trasp(str.lattice);
      // vector<xmatrix<double> > vB;   // contains lattices of the potential xstructures
      // vector<xmatrix<double> > sHNF; // contains only the supercells HNF that are uniques
      // long long int i=0;
      // bool found=false;
      // for(int a=1;a<=n;a++)
      //  for(int c=1;c<=n/a;c++)
      //    for(int f=1;f<=n/a/c;f++)
      //      if(a*c*f==n)
      //        for(int b=0;b<c;b++)
      //          for(int d=0;d<f;d++)
      //            for(int e=0;e<f;e++) {
      //              i++;
      //              HNF(1,1)=a;HNF(1,2)=0;HNF(1,3)=0;
      //              HNF(2,1)=b;HNF(2,2)=c;HNF(2,3)=0;
      //              HNF(3,1)=d;HNF(3,2)=e;HNF(3,3)=f;
      //              vHNF.push_back(HNF);
      //            }
      // oss << "[AFLOW] HNF = " << vHNF.size() << endl;
      // oss.flush();
      // for(i=0;i<(int) vHNF.size();i++) {
      //  Bi=A*vHNF.at(i);
      //  found=false;
      //  for(size_t istr=0;istr<vB.size()&&found==false;istr++)  {                     // cycle through the unique vBs
      //    Bj=vB.at(istr);
      //    for(size_t pgroup=0;pgroup<str.pgroup.size()&&found==false;pgroup++) {      // cycle through the pgroup of str
      //      R=trasp(str.pgroup.at(pgroup).Uc);
      //      H=aurostd::inverse(Bj)*aurostd::inverse(R)*Bi;
      //      if(aurostd::isinteger(H)) found=true;
      //    }
      //  }
      //  if(found==false) { // not found, then plug
      //    vB.push_back(Bi);
      //    sHNF.push_back(vHNF.at(i));
      //    //      cerr << "+";
      //  }
      //}
      //  oss << endl;
      // oss << "[AFLOW] supercells = " << vB.size() << endl;
      oss << "[AFLOW] supercells = " << sHNF.size() << endl;
      oss << AFLOWIN_SEPARATION_LINE << endl; // --------------------------------
      oss.flush();
      for (size_t j = 0; j < sHNF.size(); j++) {
        xstructure strj = GetSuperCell(str, trasp(sHNF[j]));
        // strj.Standard_Conventional_UnitCellForm();
        // strj.iomode=IOVASP_ABCCAR;
        strj.title += "  [HNF(" + aurostd::utype2string(n) + ")=" + aurostd::utype2string(j + 1) + "/" + aurostd::utype2string(sHNF.size()) + "=";
        for (uint ii = 1; ii <= 3; ii++) {
          for (uint jj = 1; jj <= 3; jj++) {
            strj.title += aurostd::utype2string(sHNF[j](ii, jj)) + (ii * jj < 9 ? "," : "]");
          }
        }
        oss << "[VASP_POSCAR_MODE_EXPLICIT]START.HNF_" << n << "_" << j + 1 << "_" << sHNF.size() << endl;
        oss << strj; // << endl;
        // oss << strj.title << endl;
        oss << "[VASP_POSCAR_MODE_EXPLICIT]STOP.HNF_" << n << "_" << j + 1 << "_" << sHNF.size() << endl;
        oss << AFLOWIN_SEPARATION_LINE << endl; // --------------------------------
        oss.flush();
      }
    }
    cerr << aurostd::ostream2string(oss) << endl;
    return aurostd::ostream2string(oss);
  }
} // namespace pflow

// ***************************************************************************
// Function pflow::HNFTOL
// ***************************************************************************
namespace pflow {
  string HNFTOL(vector<string> argv, istream& input, ostream& oss) {
    oss << AFLOWIN_SEPARATION_LINE << endl; // --------------------------------
    oss << "[AFLOW] " << aflow::Banner("BANNER_TINY") << endl;
    //  oss << "argv.size()=" << argv.size() << endl;
    xstructure str(input, IOAFLOW_AUTO);
    str.BringInCell();
    //  str.partial_occupation_sum=0;
    // GET SYMMETRY
    // str.CalculateSymmetryPointGroup(); //CO, cannot calculate symmetry of a structure with overlapping atoms
    //  str.CalculateSymmetryFactorGroup();
    oss << "[AFLOW] CALCULATION OF HNFTOL xstructures" << endl;
    // oss << "[AFLOW] pgroup operations = " << str.pgroup.size() << endl;
    oss.flush();
    // oss << "fgroup operations  = " << str.fgroup.size() << endl;
    // oss << "pgroupk operations = " << str.pgroupk.size() << endl;
    uint digits = 6;
    uint digits1 = 4;
    uint digits2 = 2 * digits + 11;

    double tolerance = DEFAULT_POCC_SITE_TOL; // DEFAULT_PARTIAL_OCCUPATION_TOLERANCE; //CO20181226
    if (argv.size() >= 3) {
      tolerance = aurostd::args2utype(argv, "--hnftol", (double) DEFAULT_POCC_SITE_TOL); // DEFAULT_PARTIAL_OCCUPATION_TOLERANCE); //CO20181226
    }
    if (argv.size() == 2) {
      tolerance = str.partial_occupation_site_tol;
    } // CO20180409
    oss << "[AFLOW] hnf_tolerance=" << tolerance << endl;
    oss << AFLOWIN_SEPARATION_LINE << endl; // --------------------------------
    oss << _VASP_POSCAR_MODE_EXPLICIT_START_ << endl;
    oss << str; // << endl;
    oss << _VASP_POSCAR_MODE_EXPLICIT_STOP_ << endl;
    oss << AFLOWIN_SEPARATION_LINE << endl; // --------------------------------
    oss.precision(digits);
    double error = 1.0;
    double eps;
    int HNFi = 0;
    vector<double> error_iatom;
    vector<double> effective_pocc_iatom;
    vector<uint> M_iatom;
    for (size_t iatom = 0; iatom < str.atoms.size(); iatom++) {
      error_iatom.push_back(1.0);
    }
    for (size_t iatom = 0; iatom < str.atoms.size(); iatom++) {
      effective_pocc_iatom.push_back(0.0);
    }
    for (size_t iatom = 0; iatom < str.atoms.size(); iatom++) {
      M_iatom.push_back(0);
    }

    oss << aurostd::PaddedPRE(aurostd::utype2string(0), digits1) << "  " << "--------" << " ";
    for (size_t iatom = 0; iatom < str.atoms.size(); iatom++) {
      if (str.atoms[iatom].partial_occupation_value < 1.0) {
        oss << "| " + aurostd::PaddedPOST("iatom=" + aurostd::utype2string(iatom) + "/" + aurostd::utype2string(str.atoms.size()), digits2) << " ";
      }
    }
    oss << " | " << "error" << endl;
    while (error > tolerance) {
      HNFi++;
      error = 1.0 / HNFi;
      oss << aurostd::PaddedPRE(aurostd::utype2string(HNFi), digits1) << "  " << error << " "; // endl;
      for (size_t iatom = 0; iatom < str.atoms.size(); iatom++) {
        error_iatom.at(iatom) = 1.0;
        M_iatom.at(iatom) = 0;
        for (int j = 0; j <= HNFi; j++) {
          eps = std::abs((double) j / HNFi - str.atoms[iatom].partial_occupation_value);
          if (eps < error_iatom.at(iatom)) {
            error_iatom.at(iatom) = eps;
            M_iatom.at(iatom) = (uint) j;
            effective_pocc_iatom.at(iatom) = (double) M_iatom.at(iatom) / HNFi;
          }
        }
      }
      for (size_t iatom = 0; iatom < str.atoms.size(); iatom++) {
        if (str.atoms[iatom].partial_occupation_value < 1.0) {
          stringstream aus;
          aus.precision(digits);
          aus.setf(std::ios::fixed, std::ios::floatfield);
          aus << effective_pocc_iatom.at(iatom) << "," << error_iatom.at(iatom) << "," << M_iatom.at(iatom) << "/" << HNFi;
          oss << "| " << aurostd::PaddedPOST(aus.str(), digits2) << " ";
        }
      }
      error = aurostd::max(error_iatom);
      oss << " | " << error << endl;
    }
    oss << AFLOWIN_SEPARATION_LINE << endl; // --------------------------------
    // DONE
    return aurostd::ostream2string(oss);
  }
} // namespace pflow

// ***************************************************************************
// pflow::IDENTICAL
// ***************************************************************************
namespace pflow {
  xstructure IDENTICAL(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    xstructure b(a);
    b = IdenticalAtoms(b);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::INCELL
// ***************************************************************************
namespace pflow {
  xstructure INCELL(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    xstructure b(a);
    b = BringInCell(b);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::INCOMPACR
// ***************************************************************************
namespace pflow {
  xstructure INCOMPACT(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    xstructure b;
    b = PutInCompact(a);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::INTPOL
// ***************************************************************************
namespace pflow {
  void INTPOL(const string& options) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 4) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --intpol=file1,file2,nimages,nearest_image_flag");
    }
    if (!aurostd::FileExist(tokens.at(0))) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "file not found: " + tokens[0], _FILE_CORRUPT_); // CO20200624
    }
    if (!aurostd::FileExist(tokens.at(1))) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "file not found: " + tokens[1], _FILE_CORRUPT_); // CO20200624
    }
    const xstructure strA(tokens.at(0), IOAFLOW_AUTO);
    const xstructure strB(tokens.at(1), IOAFLOW_AUTO);

    cout << aflow::Banner("BANNER_TINY") << endl;
    const int nimages = aurostd::string2utype<int>(tokens.at(2));
    const string nearest_image_flag = tokens.at(3);
    PrintImages(strA, strB, nimages, nearest_image_flag);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::INWS
// ***************************************************************************
namespace pflow {
  xstructure INWS(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    xstructure b;
    b = BringInWignerSeitz(a);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::JMOL
// ***************************************************************************
// Stefano Curtarolo (Dec-2009) //modification Richard Taylor // modification SC Oct2014
namespace pflow {
  void JMOL(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    string COLOR = "white"; // default
    bool WRITE = false; // default
    const xvector<int> ijk(3); // default
    ijk[1] = 1;
    ijk[2] = 1;
    ijk[3] = 1; // default

    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");

    if (tokens.size() > 5) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --jmol[=n1[,n2[,n3[,color[,true/false]]]]] < POSCAR");
    }

    if (!tokens.empty()) {
      ijk[1] = aurostd::string2utype<int>(tokens[0]);
    }
    if (tokens.size() >= 2) {
      ijk[2] = aurostd::string2utype<int>(tokens[1]);
    }
    if (tokens.size() >= 3) {
      ijk[3] = aurostd::string2utype<int>(tokens[2]);
    }
    if (tokens.size() >= 4) {
      COLOR = tokens[3];
    }
    if (tokens.size() >= 5) {
      if (tokens[4] == "true" || tokens[4] == "1") {
        WRITE = true;
        cerr << XPID << "pflow::JMOL: saving output..." << endl;
      }
    }
    // ofstream script; //file to determine color of background, bond width etc
    ofstream ras;
    const string rasFILE = "ras." + XHOST.ostrPID.str() + "." + XHOST.ostrTID.str(); // CO20200502 - threadID
    ras.open(rasFILE.c_str(), std::ios::out);
    ras << "background " << COLOR << endl << "wireframe " << 0.05 << endl; // radius of bond in angstroms
    // anymore? (http://jmol.sourceforge.net/docs/JmolUserGuide/ch04.html)
    ras.close();

    ofstream FileOUTPUT;
    const string FileOUTPUTName = "aflow.jmol.xyz." + XHOST.ostrPID.str() + "." + XHOST.ostrTID.str(); // CO20200502 - threadID
    FileOUTPUT.open(FileOUTPUTName.c_str(), std::ios::out);
    if (FileOUTPUT) {
      xvector<int> _ijk(vabs(ijk));
      if (max(_ijk) == 0) {
        _ijk.set(1);
      }
      ostringstream aus_exec;
      xstructure a(input, IOAFLOW_AUTO);
      if (_ijk(1) != 1 || _ijk(2) != 1 || _ijk(3) != 1) {
        a = GetSuperCell(a, _ijk);
      }
      pflow::PrintCIF(FileOUTPUT, a);
      FileOUTPUT.close();
      // cerr << _ijk << endl;
      aus_exec << XHOST.command("jmol") << " -s " << rasFILE << " " << FileOUTPUTName << endl; // " && rm -f " << FileOUTPUTName << endl;
      aurostd::execute(aus_exec);
      // aurostd::StringstreamClean(aus_exec);
      if (WRITE == false) {
        aurostd::RemoveFile(FileOUTPUTName);
      }
      aurostd::RemoveFile(rasFILE);
    } else {
      FileOUTPUT.close();
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::KBAND
// ***************************************************************************
// Wahyu Setyawan (Nov-2009)
namespace pflow {
  void KBAND(vector<string> argv) {
    // 2008 wahyu setyawan    ws26@duke.edu
    // 2008 roman chepulskyy  rc74@duke.edu
    //   cout << " uband -h     : help\n"
    //   << " uband 1      : nonspin\n"
    //   << " uband 2      : spin-polarized.\n";
#define _NSEGMAX_ 50

    int const argc = argv.size();
    ifstream kin;
    ifstream ein;
    ifstream din;

    if (argc != 3) {
      return;
    }
    if (!(atoi(&argv.at(2)[0]) == 1 or atoi(&argv.at(2)[0]) == 2)) {
      return;
    }

    int i;
    int ind;
    int ispin;
    int iseg;
    int itmp;
    int j;
    int ngrid;
    int nseg;
    int nband;
    float ftmp;
    float kseg[_NSEGMAX_][4];
    float kx1;
    float ky1;
    float kz1;
    float kx2;
    float ky2;
    float kz2;
    float dkx;
    float dky;
    float dkz;
    float dk;
    float koffset;
    string tmpstr;

    ispin = atoi(&argv.at(2)[0]);
    kin.open("KPOINTS");
    getline(kin, tmpstr);
    kin >> ngrid;
    getline(kin, tmpstr);
    getline(kin, tmpstr);
    getline(kin, tmpstr);
    i = 0;
    while (!kin.eof()) {
      i++;
      kin >> kseg[i][1] >> kseg[i][2] >> kseg[i][3];
      getline(kin, tmpstr);
    }
    nseg = i / 2;
    kin.close();
    // constructing kline from the norm of kpoints cascaded for continuous plot
    const xvector<float> kline(nseg * ngrid);
    koffset = 0;
    ind = 1;
    for (iseg = 1; iseg <= nseg; iseg++) {
      i = 2 * iseg - 1;
      kx1 = kseg[i][1];
      ky1 = kseg[i][2];
      kz1 = kseg[i][3];
      kx2 = kseg[i + 1][1];
      ky2 = kseg[i + 1][2];
      kz2 = kseg[i + 1][3];
      //    knorm1=sqrt(kx1*kx1 + ky1*ky1 + kz1*kz1); //kstart
      // knorm2=koffset+sqrt(kx2*kx2 + ky2*ky2 + kz2*kz2); //kstop
      dkx = kx2 - kx1;
      dky = ky2 - ky1;
      dkz = kz2 - kz1;
      dk = (sqrt(dkx * dkx + dky * dky + dkz * dkz)) / (ngrid - 1);
      kline[ind++] = koffset;
      for (j = 2; j <= ngrid; j++) {
        kline[ind++] = koffset + (j - 1) * dk;
      }
      koffset = koffset + dk * (ngrid - 1);
    }
    // creating overall kband data
    // format:
    //              spinup                 spindown
    // kline[1] band1 band2 band3 ...  band1 band2 band3 ...
    // kline[2] band1 band2 band3 ...  band1 band2 band3 ...
    //...
    // this format is nice for plotting with kpoints as the x-axis
    ein.open("EIGENVAL");
    getline(ein, tmpstr);
    getline(ein, tmpstr);
    getline(ein, tmpstr);
    getline(ein, tmpstr);
    getline(ein, tmpstr);
    ein >> itmp >> itmp;
    if (nseg != itmp / ngrid) {
      cerr << "error nseg inconsistent between KPOINTS and EIGENVAL\n";
      return;
    }
    ein >> nband;
    getline(ein, tmpstr);
    const xmatrix<float> kband(nseg * ngrid, ispin * nband + 1);
    for (i = 1; i <= nseg * ngrid; i++) {
      ein >> ftmp >> ftmp >> ftmp >> ftmp;
      kband[i][1] = kline[i];
      for (j = 1; j <= nband; j++) {
        ein >> ftmp >> kband[i][j + 1];
        if (ispin == 2) {
          ein >> kband[i][nband + j + 1];
        }
      }
    }
    ein.close();
    // output
    for (i = 1; i <= nseg * ngrid; i++) {
      for (j = 1; j <= nband * ispin + 1; j++) {
        cout << "  " << kband[i][j];
      }
      cout << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::KPATH
// ***************************************************************************
namespace pflow {
  void KPATH(istream& input, bool WWW) { // CO20200329
    const xstructure str_in(input, IOAFLOW_AUTO);
    double grid = DEFAULT_BANDS_GRID;
    if (str_in.num_each_type.size() == 1) {
      grid = DEFAULT_UNARY_BANDS_GRID;
    }
    return KPATH(input, grid, WWW);
  }
  void KPATH(istream& input, double grid, bool WWW) {
    const xstructure str_in(input, IOAFLOW_AUTO);
    if (std::signbit(grid)) { // CO20200329
      grid = DEFAULT_BANDS_GRID;
      if (str_in.num_each_type.size() == 1) {
        grid = DEFAULT_UNARY_BANDS_GRID;
      }
    }
    xstructure str_sp;
    xstructure str_sc;
    const bool full_sym = false; // DX20170829 - Speed increase
    LATTICE::Standard_Lattice_Structure(str_in, str_sp, str_sc, full_sym); // DX20170829 - Speed increase
    string lattice_type;
    stringstream oss;
    string sss;
    bool foundBZ;
    //    if(grid<=0.0) grid=16.0; // no more default
    lattice_type = str_sp.bravais_lattice_variation_type;
    //  if(WWW) cout << "<pre>" << endl;
    cout << "! STRUCTURE TO RUN ***********************************************************************" << endl; // AZ202402022 switched from // to ! so this works in qe and vasp
    oss.clear();
    oss << str_sp;
    sss = oss.str();
    if (WWW) {
      aurostd::StringSubstInPlace(sss, "<", "&#60;");
      aurostd::StringSubstInPlace(sss, ">", "&#62;");
    }
    cout << sss;
    cout << "! KPOINTS TO RUN *************************************************************************" << endl; // AZ202402022 switched from // to ! so this works in qe and vasp
    sss = LATTICE::KPOINTS_Directions(lattice_type, str_sp.lattice, grid, str_sp.iomode, foundBZ);
    if (WWW) {
      aurostd::StringSubstInPlace(sss, "<", "&#60;");
      aurostd::StringSubstInPlace(sss, ">", "&#62;");
    }
    cout << sss;
    // cout << LATTICE::KPOINTS_Directions(lattice_type,str_sp.lattice,grid,foundBZ);
    //  if(WWW) cout << "</pre>" << endl;
    if (WWW) {
      cout << "// PATH ***********************************************************************************" << endl;
      //   cout << lattice_type << endl;
      //   cout << str_sp.lattice << endl;
      cout << "</pre>" << endl;
      cout << "<img src=https://" << XHOST.AFLOW_MATERIALS_SERVER << "/SCIENCE/images/brillouin/" << lattice_type << ".PNG><br>" << endl; // ME+DX20210521 - remove height, let the CSS handle it
      // cout << "<br> [ ";
      // cout << "<a href=http://" << XHOST.AFLOW_MATERIALS_SERVER << "/SCIENCE/images/brillouin/" << lattice_type << ".PNG>png</a>";
      // cout << " | ";
      // cout << "<a href=http://" << XHOST.AFLOW_MATERIALS_SERVER << "/SCIENCE/images/brillouin/" << lattice_type << ".EPS>eps</a>";
      // cout << " ] ";
      cout << "<pre>" << endl;
    }
    cout << "! END ************************************************************************************" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::KPOINTS
// ***************************************************************************
namespace pflow {
  xstructure KPOINTS(const string& options, istream& input, ostream& oss) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl; // CO20190520
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 1) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --kpoints=KDENS [or --kppra=KDENS,-k=KDENS] < POSCAR"); // CO20190520
    }
    int NK = aurostd::string2utype<int>(tokens.at(0));

    oss << aflow::Banner("BANNER_TINY") << endl;
    xstructure str(input, IOAFLOW_AUTO);
    if (NK == 0) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "--kpoints=NUMBER: NUMBER must be bigger than zero", _INPUT_ILLEGAL_); // CO20200624
    } else {
      const ofstream FileDUMMY;
      if (NK > 10000000) {
        NK = 10000000;
      }
      oss << kintro << "KPPRA     = " << NK << " (requested) " << endl;
      const double NK_tmp = (int) ((double) NK / str.atoms.size() + 0.5);
      if (NK < 1) {
        NK = 1; // CO20180226 - this is done inside KPPRA(), so don't overwrite NK
      }
      oss << kintro << "KPPRA/#AT = " << NK_tmp << endl; // CO20180226
      KPPRA(str, NK);
      if (LDEBUG) { // CO20190401
        cerr << __AFLOW_FUNC__ << " str.kpoints_k1=" << str.kpoints_k1 << endl;
        cerr << __AFLOW_FUNC__ << " str.kpoints_k2=" << str.kpoints_k2 << endl;
        cerr << __AFLOW_FUNC__ << " str.kpoints_k3=" << str.kpoints_k3 << endl;
      }
      //	oss << kintro << "KPOINTS   = [" << str.kpoints_k1 << "," << str.kpoints_k2 << "," << str.kpoints_k3 << "]" << endl;
      oss << kintro << "DKGRID    = [" << modulus(str.klattice(1)) / str.kpoints_k1 << "," << modulus(str.klattice(2)) / str.kpoints_k2 << "," << modulus(str.klattice(3)) / str.kpoints_k3 << "]" << endl;
      oss << kintro << "KPPRA     = " << str.kpoints_kppra << " (found) " << endl;
      // oss << kintro << "next line for automatic scripting (with cat POSCAR | aflow --kpoints | grep -i AUTO | sed \"s/AUTO//g\")" << endl;
      oss << kintro << "KPOINTS   = " << str.kpoints_k1 << " " << str.kpoints_k2 << " " << str.kpoints_k3 << endl;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " END" << endl; // CO20190520
      }
      return str;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::KPOINTS_DELTA
// ***************************************************************************
namespace pflow {
  xstructure KPOINTS_DELTA(aurostd::xoption& vpflow, istream& input, ostream& oss) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    const double DK = aurostd::string2utype<double>(vpflow.getattachedscheme("FLAG::XVASP_KPOINTS_DELTA")); // tokens.at(0));  //CO20171010

    oss << aflow::Banner("BANNER_TINY") << endl;
    xstructure str(input, IOAFLOW_AUTO);
    if (DK <= 1.0E-6) {
      stringstream message;
      message << "--delta_kpoints=NUMBER: NUMBER must be bigger than zero (dk=" << DK << ")";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_); // CO20200624
    } else {
      const ofstream FileDUMMY;
      oss << kintro << "DK     = " << DK << " (requested) " << endl;
      KPPRA_DELTA(str, DK);
      //	oss << kintro << "KPOINTS   = [" << str.kpoints_k1 << "," << str.kpoints_k2 << "," << str.kpoints_k3 << "]" << endl;
      oss << kintro << "DKGRID    = [" << modulus(str.klattice(1)) / str.kpoints_k1 << "," << modulus(str.klattice(2)) / str.kpoints_k2 << "," << modulus(str.klattice(3)) / str.kpoints_k3 << "]" << endl;
      oss << kintro << "KPPRA     = " << str.kpoints_kppra << " (found) " << endl;
      // oss << kintro << "next line for automatic scripting (with cat POSCAR | aflow --kpoints | grep -i AUTO | sed \"s/AUTO//g\")" << endl;
      oss << kintro << "KPOINTS   = " << str.kpoints_k1 << " " << str.kpoints_k2 << " " << str.kpoints_k3 << endl;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " END" << endl;
      }
      return str;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::JOINSTRLIST
// ***************************************************************************
namespace pflow {
  void JOINSTRLIST(vector<string> argv) {
    ifstream list1_inf(argv.at(2).c_str());
    aurostd::InFileExistCheck("aflow", argv.at(2).c_str(), list1_inf);
    ifstream list2_inf(argv.at(3).c_str());
    aurostd::InFileExistCheck("aflow", argv.at(3).c_str(), list2_inf);
    std::vector<xstructure> str_vec_1(0);
    std::vector<xstructure> str_vec_2(0);
    std::vector<xstructure> str_vec_tot(0);
    pflow::ReadInStrVec(str_vec_1, list1_inf);
    pflow::ReadInStrVec(str_vec_2, list2_inf);
    pflow::JoinStrVec(str_vec_1, str_vec_2, str_vec_tot);
    pflow::PrintStrVec(str_vec_tot, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::LATTICEREDUCTION
// ***************************************************************************
namespace pflow {
  xstructure LATTICEREDUCTION(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    return LatticeReduction(a);
  }
} // namespace pflow

// ***************************************************************************
// pflow::LATTICE_TYPE
// ***************************************************************************
namespace pflow {
  string LATTICE_TYPE(istream& input, aurostd::xoption& vpflow) { // DX20200820
    stringstream sss;
    xstructure a(input, IOAFLOW_AUTO);
    const double tolerance = pflow::getSymmetryTolerance(a, vpflow.getattachedscheme("LATTICE::TOLERANCE")); // DX20200820 - consolidated setting tolerance into a function
    a.sym_eps = tolerance; // DX20200820 - add tolerance
    xstructure str_sp;
    xstructure str_sc; // DX20170824 - Speed increase
    const bool full_sym = false; // DX20170829 - Speed increase
    LATTICE::Standard_Lattice_Structure(a, str_sp, str_sc, full_sym); // DX20170829 - Speed increase
    // sss << " Real space lattice primitive           = " << a.bravais_lattice_type << endl;
    // sss << " Real space lattice variation           = " << a.bravais_lattice_variation_type << endl;//WSETYAWAN mod
    // sss << " Real space conventional lattice        = " << a.bravais_conventional_lattice_type << endl;
    // sss << " Real space Pearson symbol              = " << a.pearson_symbol << endl;
    sss << str_sp.bravais_lattice_type << "," << str_sp.bravais_lattice_variation_type << endl; // WSETYAWAN mod //DX20170824 - a to str_sp
    //  sss << a.bravais_lattice_type << "," << a.bravais_conventional_lattice_type << endl;
    return sss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::LATTICE_LATTICE_TYPE
// ***************************************************************************
namespace pflow {
  string LATTICE_LATTICE_TYPE(istream& input, aurostd::xoption& vpflow) { // DX20200820
    stringstream sss;
    xstructure a(input, IOAFLOW_AUTO);
    const double tolerance = pflow::getSymmetryTolerance(a, vpflow.getattachedscheme("LATTICE_LATTICE::TOLERANCE")); // DX20200820 - consolidated setting tolerance into a function
    a.sym_eps = tolerance; // DX20200820 - add tolerance
    xstructure str_sp;
    xstructure str_sc; // DX20170824 - Speed increase
    const bool full_sym = false; // DX20170829 - Speed increase
    LATTICE::Bravais_Lattice_StructureDefault(a, str_sp, str_sc, full_sym); // DX20170829 - Speed increase
    // sss << " Real space lattice primitive           = " << a.bravais_lattice_lattice_type << endl;
    // sss << " Real space lattice variation           = " << a.bravais_lattice_lattice_variation_type << endl;//WSETYAWAN mod
    // sss << " Real space conventional lattice        = " << a.bravais_conventional_lattice_lattice_type << endl;
    // sss << " Real space Pearson symbol              = " << a.pearson_symbol << endl;
    sss << str_sp.bravais_lattice_lattice_type << "," << str_sp.bravais_lattice_lattice_variation_type << endl; // WSETYAWAN mod //DX20170824 - a to str_sp
    //  sss << a.bravais_lattice_lattice_type << "," << a.bravais_conventional_lattice_lattice_type << endl;
    return sss.str();
  }
} // namespace pflow

// DX20190109 - list prototype labels - START
//  ***************************************************************************
//  pflow::listPrototypeLabels()
//  ***************************************************************************
namespace pflow {
  string listPrototypeLabels(aurostd::xoption& vpflow) {
    const anrl::ProtoData pd = anrl::ProtoData::get();

    // ----- mode -----//
    enum Mode {
      all, //< get all prototypes
      n_nary, //< get n-nary prototypes
      stoichiometry, //< get prototypes with given stoichiometry
      symmetry //< get prototypes with given symmetry
    };
    Mode mode = all;

    // ----- library name ----- //
    string library = "all"; // default: all; legacy; encyclopedia
    if (vpflow.flag("LIST_PROTOTYPE_LABELS::LIBRARY")) {
      library = vpflow.getattachedscheme("LIST_PROTOTYPE_LABELS::LIBRARY");
    }

    // ----- arity/number of species ----- //
    uint arity = 0; // default: 0=all; other options: 1=unary, 2=binary, etc.
    if (vpflow.flag("LIST_PROTOTYPE_LABELS::ARITY")) {
      mode = Mode::n_nary;
      arity = aurostd::string2utype<uint>(vpflow.getattachedscheme("LIST_PROTOTYPE_LABELS::ARITY"));
    }

    // ----- stoichiometry ----- //
    std::string stoichiometry_key;
    // default: <empty>=any; other options: x,x,...
    if (vpflow.flag("LIST_PROTOTYPE_LABELS::STOICHIOMETRY")) {
      vector<string> tokens;
      stringstream message;
      vector<uint> stoichiometry;
      mode = Mode::stoichiometry;
      aurostd::string2tokens(vpflow.getattachedscheme("LIST_PROTOTYPE_LABELS::STOICHIOMETRY"), tokens, ",");
      for (size_t i = 0; i < tokens.size(); i++) {
        stoichiometry.push_back(aurostd::string2utype<uint>(tokens[i]));
      }
      if (arity != 0 && arity != stoichiometry.size()) {
        message << "arity=" << arity << " and stoichiometry=" << aurostd::joinWDelimiter(stoichiometry, ",") << " do not match.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); // DX20191107 - exit -> throw
      }
      std::sort(stoichiometry.begin(), stoichiometry.end(), std::greater<>()); // must sort to properly filter
      stoichiometry_key = aurostd::joinWDelimiter(stoichiometry, ':');
    }

    // ----- space group number ----- //
    uint space_group_number = 0; // default: 0=any; other options: 1-230.
    if (vpflow.flag("LIST_PROTOTYPE_LABELS::SPACE_GROUP")) {
      mode = Mode::symmetry;
      space_group_number = aurostd::string2utype<uint>(vpflow.getattachedscheme("LIST_PROTOTYPE_LABELS::SPACE_GROUP"));
    }

    //-------------------------------------------------
    // get prototype uids
    vector<string> prototype_uid;

    // get all prototype labels
    switch (mode) {
      case all: prototype_uid = pd.lookup[library]; break;
      // get n-nary prototype labels
      case n_nary: prototype_uid = pd.lookup["number_of_species"][arity]; break;
      // get prototype labels via stoichiometry
      case stoichiometry: prototype_uid = pd.lookup["stoichiometry"][stoichiometry_key]; break;
      // get prototype labels via symmetry
      case symmetry: prototype_uid = pd.lookup["space_group_number"][space_group_number]; break;
    }
    std::string output;
    if (vpflow.flag("LIST_PROTOTYPE_LABELS::JSON")) {
      const aurostd::JSON::object output_json(aurostd::JSON::object_types::DICTIONARY);
      for (const auto& uid : prototype_uid) {
        output_json[uid] = pd.content[uid]["label"];
      }
      output = output_json.toString();
    } else {
      for (const auto& uid : prototype_uid) {
        output += static_cast<string>(pd.content[uid]["label"]) + "(" + uid + ")\n";
      }
    }
    return output;
  }
} // namespace pflow
// DX20190109 - list prototype labels - END

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// START - all relevent functions for loading entries here
// Added by Corey Oses - May 2017
// load entries is heavily overloaded, mostly to accommodate entries separated as
// vector<vector<vector<> > > entries (unaries vs. binaries, then species-specific, good for convex hull),
// vector<vector<> > entries (unaries vs. binaries), OR
// vector<> entries (all together)
//  ***************************************************************************
//  pflow::loadEntries(aurostd::xoption& vpflow,vector<string>& velements,string server,
//  ofstream& FileMESSAGE, ostream& oss)
//  ***************************************************************************
namespace pflow {
  string arity_string(uint arity, bool capital, bool plural) { // for printing, so we don't need to make zillions of if statements elsewhere
    if (arity == 0) {
      return (capital ? string("N") : string("n")) + string("ullar") + (plural ? string("ies") : string("y"));
    }
    if (arity == 1) {
      return (capital ? string("U") : string("u")) + string("nar") + (plural ? string("ies") : string("y"));
    }
    if (arity == 2) {
      return (capital ? string("B") : string("b")) + string("inar") + (plural ? string("ies") : string("y"));
    }
    if (arity == 3) {
      return (capital ? string("T") : string("t")) + string("ernar") + (plural ? string("ies") : string("y"));
    }
    if (arity == 4) {
      return (capital ? string("Q") : string("q")) + string("uaternar") + (plural ? string("ies") : string("y"));
    }
    if (arity == 5) {
      return (capital ? string("Q") : string("q")) + string("uinar") + (plural ? string("ies") : string("y"));
    }
    if (arity == 6) {
      return (capital ? string("S") : string("s")) + string("enar") + (plural ? string("ies") : string("y"));
    }
    if (arity == 7) {
      return (capital ? string("S") : string("s")) + string("eptenar") + (plural ? string("ies") : string("y"));
    }
    if (arity == 8) {
      return (capital ? string("O") : string("o")) + string("ctonar") + (plural ? string("ies") : string("y"));
    }
    if (arity == 9) {
      return (capital ? string("N") : string("n")) + string("ovenar") + (plural ? string("ies") : string("y"));
    }
    if (arity == 10) {
      return (capital ? string("D") : string("d")) + string("enar") + (plural ? string("ies") : string("y"));
    }
    return aurostd::utype2string(arity) + string("-ar") + (plural ? string("ies") : string("y"));
  }
} // namespace pflow
namespace pflow {
  bool loadEntries(vector<string>& velements, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries, ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadEntries(velements, entries, FileMESSAGE, oss);
  }
  bool loadEntries(vector<string>& velements, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    string server;
    if (XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
      server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
    } else {
      server = AFLOWLIB_SERVER_DEFAULT;
    }
    return loadEntries(velements, server, entries, FileMESSAGE, oss);
  }
  bool loadEntries(vector<string>& velements, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries, ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadEntries(velements, server, entries, FileMESSAGE, oss);
  }
  bool loadEntries(vector<string>& velements, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    aurostd::xoption vpflow;
    pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, std::string("A"));
    return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries,
                   ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadEntries(vpflow, velements, entries, FileMESSAGE, oss);
  }
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    string server;
    if (XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
      server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
    } else {
      server = AFLOWLIB_SERVER_DEFAULT;
    }
    return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries,
                   ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries, ofstream& FileMESSAGE, ostream& oss) { // main function

    stringstream message;

    // get directory info
    _aflags aflags;
    aflags.Directory = chull::getPath(); // default
    if (vpflow.flag("CHULL::PATH")) {
      aflags.Directory = chull::getPath(vpflow, FileMESSAGE, oss);
    }

    vpflow.flag("PFLOW::LOAD_ENTRIES_COMING_FROM_LOADENTRIESX", true); // silence some output

    message << "Loading entries for: " << aurostd::joinWDelimiter(velements, ",");
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);

    if (vpflow.flag("PFLOW::LOAD_ENTRIES_NON_ALPHABETICAL")) {
      message << "Loading NON-alphabetized entries as well";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_OPTION_);
    }

    vector<vector<string>> combinations;
    vector<vector<vector<aflowlib::_aflowlib_entry>>> _entries; // unsorted

    if (false) {
      // anticipating AFLUX integration

    } else {
      const bool load_from_common = pflow::loadFromCommon(vpflow);

      if (load_from_common) {
        message << "Loading entries from COMMON";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_OPTION_);
      } else {
        if ((server == AFLOW_MATERIALS_SERVER_DEFAULT) || (server == AFLOWLIB_SERVER_DEFAULT)) {
          message << "Using " << server << " as server";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_OPTION_);
        } else {
          message << "Server must be either " << AFLOW_MATERIALS_SERVER_DEFAULT << " or " << AFLOWLIB_SERVER_DEFAULT;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
          return false; // entries;
        }
      }

      string lib_name;
      string lib_count_string;
      string load_lib_flag_name;

      //////////////////////////////////////////////////////////////////////////////
      // START loadLIBX
      //////////////////////////////////////////////////////////////////////////////

      for (size_t lib = 1; lib <= velements.size() && lib <= _AFLOW_LIB_MAX_; lib++) {
        lib_count_string = aurostd::utype2string(lib);
        lib_name = std::string("LIB") + lib_count_string;
        load_lib_flag_name = "PFLOW::LOAD_ENTRIES_LOAD_" + lib_name;
        if (vpflow.flag(load_lib_flag_name)) {
          message << "Loading " + lib_name;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
          combinations = pflow::elementalCombinations(velements, lib);
          for (size_t i = 0; i < combinations.size(); i++) {
            if (!loadAndMergeLIBX(vpflow, combinations[i], lib_count_string, server, _entries, FileMESSAGE, oss)) {
              message << "Merging entries for " + lib_name + " failed";
              pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
              return false;
            }
          }
        }
      }

      //////////////////////////////////////////////////////////////////////////////
      // END loadLIBX
      //////////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////////
      // START loadICSD
      //////////////////////////////////////////////////////////////////////////////

      lib_name = std::string("ICSD");
      load_lib_flag_name = "PFLOW::LOAD_ENTRIES_LOAD_" + lib_name;
      if (vpflow.flag(load_lib_flag_name)) {
        message << "Loading " + lib_name;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        if (!loadAndMergeLIBX(vpflow, velements, lib_name, server, _entries, FileMESSAGE, oss)) {
          message << "Merging entries for " + lib_name + " failed";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
          return false;
        }
      }

      //////////////////////////////////////////////////////////////////////////////
      // END loadICSD
      //////////////////////////////////////////////////////////////////////////////
    }

    //////////////////////////////////////////////////////////////////////////////
    // START Sort
    //////////////////////////////////////////////////////////////////////////////

    bool found;
    vector<string> vspecies;
    // make space and sort correctly
    for (size_t i = 0; i < velements.size(); i++) {
      entries.emplace_back(0);
      combinations = pflow::elementalCombinations(velements, i + 1);
      for (size_t j = 0; j < combinations.size(); j++) {
        entries[i].emplace_back(0);
      }
      for (size_t j = 0; j < combinations.size(); j++) {
        found = false;
        if (_entries.size() > i) {
          for (size_t k = 0; k < _entries[i].size() && !found; k++) {
            if (!_entries[i][k].empty()) {
              vspecies = _entries[i][k][0].vspecies;
              if (vpflow.flag("PFLOW::LOAD_ENTRIES_NON_ALPHABETICAL")) {
                sort(vspecies.begin(), vspecies.end());
              }
              if (vspecies == combinations[j]) {
                entries[i][j] = _entries[i][k];
                found = true;
              }
            }
          }
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////////
    // END Sort
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // START Print loaded entries summary
    //////////////////////////////////////////////////////////////////////////////

    vector<uint> sizes;
    uint totalNum = 0;

    for (size_t i = 0; i < entries.size(); i++) {
      sizes.push_back(0);
      for (size_t j = 0; j < entries[i].size(); j++) {
        for (size_t k = 0; k < entries[i][j].size(); k++) {
          sizes.at(i)++;
        }
      }
      totalNum += sizes.at(i);
    }

    for (size_t i = 0; i < entries.size(); i++) {
      message << "Loaded ";
      message << entries[i].size() << " " << pflow::arity_string(i + 1, false, true) << ": " << sizes.at(i) << " entries";
      // if(i == 0) {
      //   message << "unaries: ";
      // } else if(i == 1) {
      //   message << "binaries: ";
      // } else if(i == 2) {
      //   message << "ternaries: ";
      // } else {
      //   message << (i + 1) << "-naries: ";
      // }
      // message << sizes.at(i) << " entries";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
    }
    message << "Loaded " << totalNum << " entries total";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);

    //////////////////////////////////////////////////////////////////////////////
    // END Print loaded entries summary
    //////////////////////////////////////////////////////////////////////////////

    return true; // entries;
  }
} // namespace pflow

namespace pflow {
  // ***************************************************************************
  // pflow::loadEntries(aurostd::xoption& vpflow,vector<string>& velements,string
  // server,ostream& oss,ofstream& FileMESSAGE)
  // ***************************************************************************
  bool loadEntries(vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadEntries(velements, entries, FileMESSAGE, oss);
  }
  bool loadEntries(vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    string server;
    if (XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
      server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
    } else {
      server = AFLOWLIB_SERVER_DEFAULT;
    }
    return loadEntries(velements, server, entries, FileMESSAGE, oss);
  }
  bool loadEntries(vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries,
                   ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadEntries(velements, server, entries, FileMESSAGE, oss);
  }
  bool loadEntries(vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    aurostd::xoption vpflow;
    pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, std::string("A"));
    return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries,
                   ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadEntries(vpflow, velements, entries, FileMESSAGE, oss);
  }
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    string server;
    if (XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
      server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
    } else {
      server = AFLOWLIB_SERVER_DEFAULT;
    }
    return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries,
                   ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss) { // main function

    vector<vector<vector<aflowlib::_aflowlib_entry>>> naries;
    if (!loadEntries(vpflow, velements, server, naries, FileMESSAGE, oss)) {
      return false;
    }
    if (!aflowlib::mergeEntries(entries, naries, true)) {
      return false;
    }
    return true;
  }
} // namespace pflow

namespace pflow {
  // ***************************************************************************
  // pflow::loadEntries(aurostd::xoption& vpflow,vector<string>
  // velements,string server,ostream& oss,ofstream& FileMESSAGE)
  // ***************************************************************************
  bool loadEntries(vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadEntries(velements, entries, FileMESSAGE, oss);
  }
  bool loadEntries(vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    string server;
    if (XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
      server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
    } else {
      server = AFLOWLIB_SERVER_DEFAULT;
    }
    return loadEntries(velements, server, entries, FileMESSAGE, oss);
  }
  bool loadEntries(vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadEntries(velements, server, entries, FileMESSAGE, oss);
  }
  bool loadEntries(vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    aurostd::xoption vpflow;
    pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, std::string("A"));
    return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries,
                   ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadEntries(vpflow, velements, entries, FileMESSAGE, oss);
  }
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    string server;
    if (XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
      server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
    } else {
      server = AFLOWLIB_SERVER_DEFAULT;
    }
    return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries,
                   ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadEntries(vpflow, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss) { // main function

    vector<vector<vector<aflowlib::_aflowlib_entry>>> naries;
    if (!loadEntries(vpflow, velements, server, naries, FileMESSAGE, oss)) {
      return false;
    }
    if (!aflowlib::mergeEntries(entries, naries)) {
      return false;
    }
    return true;
  }
} // namespace pflow

namespace pflow {
  // ***************************************************************************
  // pflow::loadFromCommon(aurostd::xoption& vpflow)
  // simple function for determining if we can load from common
  // it's a function because we call it a few times
  // ***************************************************************************
  bool loadFromCommon(aurostd::xoption& vpflow) {
    const bool load_from_common =
        (!vpflow.flag("PFLOW::LOAD_API") && (aurostd::substring2bool(XHOST.hostname, "nietzsche") || aurostd::substring2bool(XHOST.hostname, "aflowlib") || aurostd::substring2bool(XHOST.hostname, "habana")));
    return load_from_common;
  }
} // namespace pflow

namespace pflow {
  // ***************************************************************************
  // pflow::loadAndMergeLIBX(aurostd::xoption& vpflow, vector<string>& combination,
  // string LIB, string server,vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries,
  // ofstream& FileMESSAGE, ostream& oss
  // ***************************************************************************
  // helper function for loadEntries, loads select library, LIB can be "LIB2" or "2" or "ICSD"
  bool loadAndMergeLIBX(string combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ostream& oss) {
    ofstream FileMESSAGE;
    return loadAndMergeLIBX(combination, LIB, server, naries, FileMESSAGE, oss);
  }
  bool loadAndMergeLIBX(string _combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ofstream& FileMESSAGE, ostream& oss) {
    vector<string> combination;
    combination.push_back(_combination);
    return loadAndMergeLIBX(combination, LIB, server, naries, FileMESSAGE, oss);
  }
  bool loadAndMergeLIBX(vector<string>& combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ostream& oss) {
    ofstream FileMESSAGE;
    return loadAndMergeLIBX(combination, LIB, server, naries, FileMESSAGE, oss);
  }
  bool loadAndMergeLIBX(vector<string>& combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ofstream& FileMESSAGE, ostream& oss) {
    aurostd::xoption vpflow;
    pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, LIB);
    return loadAndMergeLIBX(vpflow, combination, LIB, server, naries, oss);
  }
  bool loadAndMergeLIBX(aurostd::xoption& vpflow, string combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ostream& oss) {
    ofstream FileMESSAGE;
    return loadAndMergeLIBX(vpflow, combination, LIB, server, naries, FileMESSAGE, oss);
  }
  bool loadAndMergeLIBX(aurostd::xoption& vpflow, string _combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ofstream& FileMESSAGE, ostream& oss) {
    vector<string> combination;
    combination.push_back(_combination);
    return loadAndMergeLIBX(vpflow, combination, LIB, server, naries, FileMESSAGE, oss);
  }
  bool loadAndMergeLIBX(aurostd::xoption& vpflow, vector<string>& combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ostream& oss) {
    ofstream FileMESSAGE;
    return loadAndMergeLIBX(vpflow, combination, LIB, server, naries, FileMESSAGE, oss);
  }
  bool loadAndMergeLIBX(aurostd::xoption& vpflow, vector<string>& combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ofstream& FileMESSAGE, ostream& oss) {
    stringstream message;

    // get directory info
    _aflags aflags;
    aflags.Directory = chull::getPath(); // default
    if (vpflow.flag("CHULL::PATH")) {
      aflags.Directory = chull::getPath(vpflow, FileMESSAGE, oss);
    }

    LIB = aurostd::toupper(LIB); // removes ALL ambiguity with case
    vector<vector<aflowlib::_aflowlib_entry>> v_temp;
    if (!loadLIBX(vpflow, LIB, combination, server, v_temp, FileMESSAGE, oss)) {
      return false;
    }
    if (vpflow.flag("PFLOW::LOAD_ENTRIES_NARIES_MINUS_ONE")) {
      for (size_t i = 0; i < v_temp.size() - 1; i++) {
        if (!aflowlib::mergeEntries(naries, v_temp[i], false)) { // e.g. for ternary MnPdPt, there are 3 binary combinations
          if (LIB == "ICSD") // or LIB == "icsd")
          { // CO20200106 - patching for auto-indenting
            message << "Merging entries for ICSD (" + aurostd::utype2string(i + 1) + "-naries) failed";
          } else {
            message << "Merging entries for LIB" + LIB + " (" + aurostd::utype2string(i + 1) + "-naries) failed";
          }
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
          return false;
        }
      }
    }
    if (!aflowlib::mergeEntries(naries, v_temp[v_temp.size() - 1], true)) { // e.g. for ternary MnPdPt, there is only ONE ternary combination
      if (LIB == "ICSD") // or LIB == "icsd")
      { // CO20200106 - patching for auto-indenting
        message << "Merging entries for ICSD (" + aurostd::utype2string(v_temp.size()) + "-naries) failed";
      } else {
        message << "Merging entries for LIB" + LIB + " (" + aurostd::utype2string(v_temp.size()) + "-naries) failed";
      }
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;
    }
    return true;
  }
} // namespace pflow

namespace pflow {
  uint SubLayersRestAPILS(const string& url, vector<string>& vsuburl) { // CO20200731 - mimicking SubDirectoryLS() in aurostd for REST-API
    vector<string> vtokens;
    vector<string> _vtokens;
    aurostd::httpGetTokens(url + "/?aflowlib_entries", vtokens, ",");
    for (size_t i = 0; i < vtokens.size(); i++) {
      aurostd::httpGetTokens(url + "/" + vtokens[i] + "/?aflowlib_entries", _vtokens, ",");
      if (!_vtokens.empty()) {
        vsuburl.push_back(url + "/" + vtokens[i]);
        SubLayersRestAPILS(url + "/" + vtokens[i], vsuburl);
      }
    }
    return vsuburl.size();
  }
} // namespace pflow

namespace pflow {
  // ***************************************************************************
  // pflow::loadLIBX(aurostd::xoption& vpflow,string elements,string server,
  // ofstream& FileMESSAGE, ostream& oss)
  // ***************************************************************************
  // there are MANY loadLIBX overloads, basically variations of string/vector elements
  // and whether entries are loaded into vector<> or vector<vector<> > (unaries vs. binaries, etc.)
  // loadLIBX string elements
  bool loadLIBX(string LIB, string elements, vector<aflowlib::_aflowlib_entry>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(LIB, elements, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(string LIB, string elements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    string server;
    if (XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
      server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
    } else {
      server = AFLOWLIB_SERVER_DEFAULT;
    }
    return loadLIBX(LIB, elements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(string LIB, string elements, string server, vector<aflowlib::_aflowlib_entry>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(LIB, elements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(string LIB, string elements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    aurostd::xoption vpflow;
    pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, LIB);
    return loadLIBX(vpflow, LIB, elements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, vector<aflowlib::_aflowlib_entry>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(vpflow, LIB, elements, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    string server;
    if (XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
      server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
    } else {
      server = AFLOWLIB_SERVER_DEFAULT;
    }
    return loadLIBX(vpflow, LIB, elements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, string server, vector<aflowlib::_aflowlib_entry>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(vpflow, LIB, elements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload

    vector<string> velements = aurostd::getElements(elements, composition_string, FileMESSAGE, true, true, false, oss); // clean and sort, do not keep_pp  //CO20190712
    return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
  }
  // ***************************************************************************
  // pflow::loadLIBX(aurostd::xoption& vpflow,string LIB,vector<string>& velements,string server,
  // ofstream& FileMESSAGE, ostream& oss)
  // ***************************************************************************
  // loadLIBX vector elements
  bool loadLIBX(string LIB, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(LIB, velements, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(string LIB, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    string server;
    if (XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
      server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
    } else {
      server = AFLOWLIB_SERVER_DEFAULT;
    }
    return loadLIBX(LIB, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(string LIB, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(LIB, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(string LIB, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    aurostd::xoption vpflow;
    pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, LIB);
    return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(vpflow, LIB, velements, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    string server;
    if (XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
      server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
    } else {
      server = AFLOWLIB_SERVER_DEFAULT;
    }
    return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss) { // main function

    vector<vector<aflowlib::_aflowlib_entry>> naries; // most intuitive structure from LIBs construction (unary, binaries, etc.), always start here and merge to get other variants

    if (!loadLIBX(vpflow, LIB, velements, server, naries, FileMESSAGE, oss)) {
      return false;
    }
    if (!aflowlib::mergeEntries(entries, naries)) {
      return false;
    }
    return true;
  }
  // ***************************************************************************
  // pflow::loadLIBX(aurostd::xoption& vpflow,string elements,string server,
  // ofstream& FileMESSAGE, ostream& oss)
  // ***************************************************************************
  // loadLIBX string elements
  bool loadLIBX(string LIB, string elements, vector<vector<aflowlib::_aflowlib_entry>>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(LIB, elements, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(string LIB, string elements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    string server;
    if (XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
      server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
    } else {
      server = AFLOWLIB_SERVER_DEFAULT;
    }
    return loadLIBX(LIB, elements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(string LIB, string elements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(LIB, elements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(string LIB, string elements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    aurostd::xoption vpflow;
    pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, LIB);
    return loadLIBX(vpflow, LIB, elements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, vector<vector<aflowlib::_aflowlib_entry>>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(vpflow, LIB, elements, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    string server;
    if (XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
      server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
    } else {
      server = AFLOWLIB_SERVER_DEFAULT;
    }
    return loadLIBX(vpflow, LIB, elements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(vpflow, LIB, elements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload

    vector<string> velements = aurostd::getElements(elements, composition_string, FileMESSAGE, true, true, false, oss); // clean and sort, do not keep_pp  //CO20190712
    return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
  }
  // ***************************************************************************
  // pflow::loadLIBX(aurostd::xoption& vpflow,string LIB,vector<string>& velements,string server,
  // ofstream& FileMESSAGE, ostream& oss)
  // ***************************************************************************
  // loadLIBX vector elements
  bool loadLIBX(string LIB, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(LIB, velements, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(string LIB, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    string server;
    if (XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
      server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
    } else {
      server = AFLOWLIB_SERVER_DEFAULT;
    }
    return loadLIBX(LIB, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(string LIB, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(LIB, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(string LIB, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    aurostd::xoption vpflow;
    pflow::defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, LIB);
    return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(vpflow, LIB, velements, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss) { // overload
    string server;
    if (XHOST.vflag_control.flag("AFLOWLIB_SERVER")) {
      server = XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");
    } else {
      server = AFLOWLIB_SERVER_DEFAULT;
    }
    return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries,
                ostream& oss) { // overload
    ofstream FileMESSAGE;
    return loadLIBX(vpflow, LIB, velements, server, entries, FileMESSAGE, oss);
  }
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss) { // main function
    const bool LDEBUG = (false || XHOST.DEBUG);
    stringstream message;

    // get directory info
    _aflags aflags;
    aflags.Directory = chull::getPath(); // default
    if (vpflow.flag("CHULL::PATH")) {
      aflags.Directory = chull::getPath(vpflow, FileMESSAGE, oss);
    }

    // vector<vector<aflowlib::_aflowlib_entry> > entries;
    for (size_t i = 0; i < velements.size(); i++) {
      entries.emplace_back(0);
    }

    // categorize LIB
    string lib_name;
    const vector<string> symmetries{"BCC", "BCT", "CUB", "FCC", "HEX", "MCL", "MCLC", "ORC", "ORCC", "ORCF", "ORCI", "RHL", "TET", "TRI"};
    LIB = aurostd::toupper(LIB); // removes ALL ambiguity with case
    const bool isICSD = (LIB == "ICSD"); //(LIB == "ICSD" || LIB == "icsd");
    if (isICSD) {
      lib_name = LIB;
      // add here other special LIBS, i.e., not LIB2, LIB3, etc.
    } else {
      string lib_count_string;
      // make robust so input can be "LIB2" or "2"
      if (aurostd::substring2bool(LIB, "LIB")) {
        lib_name = LIB;
        lib_count_string = aurostd::RemoveSubString(LIB, "LIB");
      } else {
        lib_name = "LIB" + LIB;
        lib_count_string = LIB;
      }
      // check validity of input
      for (size_t i = 0; i < lib_count_string.size(); i++) {
        if (!isdigit(lib_count_string[i])) {
          message << "Unknown LIB specification (" << LIB << R"(), should be "LIB1", "LIB2", or "1", "2", etc.)";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
          return false; // entries;
        }
      }
      const uint lib_count_uint = aurostd::string2utype<int>(lib_count_string);
      if (velements.size() != lib_count_uint) {
        message << "LIB" << lib_count_uint << " loads " << lib_count_uint << "-naries ONLY";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
        return false; // entries;
      }
    }

    message << "Loading " << lib_name << " entries for: " << aurostd::joinWDelimiter(velements, "");
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);

    bool load_from_common = pflow::loadFromCommon(vpflow);
    bool override_load_from_common = false;

    string LIB_path; // will be actual directory path or URL, depending on machine
    if (load_from_common) {
      if (!vpflow.flag("PFLOW::LOAD_ENTRIES_COMING_FROM_LOADENTRIESX")) {
        message << "Loading entries from COMMON";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_OPTION_);
      }
      LIB_path = "/common/" + lib_name + "/RAW";
      if (!aurostd::IsDirectory(LIB_path)) {
        load_from_common = false;
        message << LIB_path << " does not exist! Cannot load from COMMON, switching to API";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
        override_load_from_common = true;
      }
    }
    if (!load_from_common) {
      if ((server == AFLOW_MATERIALS_SERVER_DEFAULT) || (server == AFLOWLIB_SERVER_DEFAULT)) {
        if (override_load_from_common || (!vpflow.flag("PFLOW::LOAD_ENTRIES_COMING_FROM_LOADENTRIESX"))) {
          message << "Using " << server << " as server";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_OPTION_);
        }
      } else {
        message << "Server must be either " << AFLOW_MATERIALS_SERVER_DEFAULT << " or " << AFLOWLIB_SERVER_DEFAULT;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
        return false; // entries;
      }
      if (isICSD) {
        LIB_path = server + "/AFLOWDATA/ICSD_WEB";
      } else {
        LIB_path = server + "/AFLOWDATA/" + lib_name + "_RAW";
      }
    }

    //////////////////////////////////////////////////////////////////////////
    // START Finding matches
    //////////////////////////////////////////////////////////////////////////
    aflowlib::_aflowlib_entry _aflowlib_tmp;
    uint nary = 0;
    uint total_count = 0;
    string::size_type loc = 0;
    vector<string> input_velements;
    vector<double> input_vcomposition;
    vector<string> vloadpaths;
    if (isICSD) {
      const bool double_check_icsd = false; // NOT NECESSARY for ICSD since we load from the calculation layer
      string symmetry_path;
      string clean_icsd;
      vector<string> icsds;
      vector<string> tokens;
      vector<string> elements;
      for (size_t i = 0; i < symmetries.size(); i++) {
        symmetry_path = LIB_path + "/" + symmetries[i];
        if (load_from_common && !aurostd::IsDirectory(symmetry_path)) {
          continue;
        }
        // not many symmetries, so this boolean placement won't impact much
        if (load_from_common) {
          aurostd::DirectoryLS(symmetry_path, icsds);
        } else {
          aurostd::httpGetTokens(symmetry_path + "/?aflowlib_entries", icsds, ",");
        }
        if (!icsds.empty() && icsds[0] == "<!DOCTYPE") { // CO20180627
          message << "REST-API appears to be down";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
          return false;
        }
        for (size_t j = 0; j < icsds.size(); j++) {
          loc = icsds[j].find("_ICSD_");
          if (loc != string::npos) { // found a good ICSD
            clean_icsd = icsds[j].substr(0, loc); // get just compound
            if (compoundsBelong(velements, clean_icsd, input_velements, input_vcomposition, FileMESSAGE, oss, false, vpflow.flag("PFLOW::LOAD_ENTRIES_NON_ALPHABETICAL"),
                                composition_string)) { // NOT recommended, these are BS entries // no need to clean, ICSD compound format is already clean of pps
              // url2aflowlib has bad printing to screen options, so we mimic here
              vloadpaths.clear();
              if (load_from_common) {
                aurostd::SubDirectoryLS(symmetry_path + "/" + icsds[j], vloadpaths);
              } else {
                pflow::SubLayersRestAPILS(symmetry_path + "/" + icsds[j], vloadpaths);
              }
              if (LDEBUG) {
                cerr << __AFLOW_FUNC__ << " vloadpaths=" << aurostd::joinWDelimiter(vloadpaths, ",") << endl;
              }
              if (vloadpaths.size() > 1 && vloadpaths[1] == "<!DOCTYPE") { // CO20180627
                message << "REST-API appears to be down";
                pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
                return false;
              }
              for (size_t k = 0; k < vloadpaths.size(); k++) {
                const string& load_path = vloadpaths[k];
                message << "Loading " << (load_from_common ? "path" : "url") << "=" << load_path;
                pflow::logger(__AFLOW_FILE__, std::string("aurostd::") + std::string((load_from_common ? "file" : "url")) + std::string("2string():"), message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_,
                              !vpflow.flag("PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT")); // print to log anyway
                // complicated ternary operator, returns bool, but necessary to avoid double code
                if ((load_from_common) ?
                                       // if load_from_common, check this bool
                                       // this part of the ternary operator will first look for a non-empty aflowlib.out,
                                       // then try to load it
                                       // the last part checks that the loaded entry has a set of species that belong to the original query (protection)
                        (
                            //(aurostd::FileExist(load_path + "/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) &&
                            (aurostd::FileNotEmpty(load_path + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) && (_aflowlib_tmp.file2aflowlib(load_path + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_OUT, message) > 0) &&
                            (double_check_icsd ? compoundsBelong(velements, _aflowlib_tmp.vspecies, FileMESSAGE, oss, vpflow.flag("PFLOW::LOAD_ENTRIES_NON_ALPHABETICAL"))
                                               : true) // sometimes we find odd entries in the wrong LIBS, better to be safe, NOT NECESSARY for ICSD since we load from the calculation layer
                            )
                                       :
                                       // if load_from_api, check this bool
                                       // this part of the ternary operator loads in the url
                                       // the last part checks that the loaded entry has a set of species that belong to the original query (protection)
                        ((_aflowlib_tmp.url2aflowlib(load_path, message, false) > 0) &&
                         (double_check_icsd ? compoundsBelong(velements, _aflowlib_tmp.vspecies, FileMESSAGE, oss, vpflow.flag("PFLOW::LOAD_ENTRIES_NON_ALPHABETICAL"))
                                            : true) // sometimes we find odd entries in the wrong LIBS, better to be safe, NOT NECESSARY for ICSD since we load from the calculation layer
                         )) {
                  nary = _aflowlib_tmp.vspecies.size();
                  if (entries.size() < nary) {
                    continue;
                  } // this entry is garbage (wrong directory)
                  entries[nary - 1].push_back(_aflowlib_tmp);
                  total_count++;
                  if (vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_XSTRUCTURES")) {
                    if (!loadXstructures(entries[nary - 1].back(), FileMESSAGE, oss, vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_XSTRUCTURES_RELAXED_ONLY"))) { // CO20171202 - let it decided whether to load from common or not
                      entries[nary - 1].pop_back();
                      total_count--;
                    }
                  }
                }
                message.str("");
              }
            }
          }
        }
      }
    } else {
      const bool double_check_lib = true; // while NOT NECESSARY for ICSD, LIBs are screwed up sometimes
      string clean_system;
      string calculation_path;
      string clean_calculation;
      vector<string> systems;
      vector<string> calculations;
      if (load_from_common) {
        aurostd::DirectoryLS(LIB_path, systems);
      } else {
        aurostd::httpGetTokens(LIB_path + "/?aflowlib_entries", systems, ",");
      }
      if (!systems.empty() && systems[0] == "<!DOCTYPE") { // CO20180627
        message << "REST-API appears to be down";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
        return false;
      }
      for (size_t i = 0; i < systems.size(); i++) {
        calculation_path = LIB_path + "/" + systems[i];
        if (load_from_common && !aurostd::IsDirectory(calculation_path)) {
          continue;
        }
        loc = systems[i].find(":");
        clean_system = systems[i].substr(0, loc); // even if we don't find it, simply copy string
        if (compoundsBelong(velements, clean_system, input_velements, input_vcomposition, FileMESSAGE, oss, false, vpflow.flag("PFLOW::LOAD_ENTRIES_NON_ALPHABETICAL"), pp_string,
                            true)) { // NOT recommended, these are BS entries  // no need to clean, already done  //use short_pp_string_AFLOW_database
          if (load_from_common) {
            aurostd::DirectoryLS(calculation_path, calculations);
          } else {
            aurostd::httpGetTokens(calculation_path + "/?aflowlib_entries", calculations, ",");
          }
          if (!calculations.empty() && calculations[0] == "<!DOCTYPE") { // CO20180627
            message << "REST-API appears to be down";
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
            return false;
          }
          for (size_t j = 0; j < calculations.size(); j++) {
            if (load_from_common && !aurostd::IsDirectory(calculation_path + "/" + calculations[j])) {
              continue;
            }
            // if(vpflow.flag("PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT"))
            vloadpaths.clear();
            if (load_from_common) {
              aurostd::SubDirectoryLS(calculation_path + "/" + calculations[j], vloadpaths);
            } else {
              pflow::SubLayersRestAPILS(calculation_path + "/" + calculations[j], vloadpaths);
            }
            if (LDEBUG) {
              cerr << __AFLOW_FUNC__ << " vloadpaths=" << aurostd::joinWDelimiter(vloadpaths, ",") << endl;
            }
            if (vloadpaths.size() > 1 && vloadpaths[1] == "<!DOCTYPE") { // CO20180627
              message << "REST-API appears to be down";
              pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
              return false;
            }
            for (size_t k = 0; k < vloadpaths.size(); k++) {
              const string& load_path = vloadpaths[k];
              message << "Loading " << (load_from_common ? "path" : "url") << "=" << load_path;
              pflow::logger(__AFLOW_FILE__, std::string("aurostd::") + std::string((load_from_common ? "file" : "url")) + std::string("2string():"), message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_,
                            !vpflow.flag("PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT")); // print to log anyway
              // complicated ternary operator, returns bool, but necessary to avoid double code
              if ((load_from_common) ?
                                     // if load_from_common, check this bool
                                     // this part of the ternary operator will first look for a non-empty aflowlib.out,
                                     // then try to load it
                                     // the last part checks that the loaded entry has a set of species that belong to the original query (protection)
                      (
                          //(aurostd::FileExist(load_path+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) &&
                          (aurostd::FileNotEmpty(load_path + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) && (_aflowlib_tmp.file2aflowlib(load_path + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_OUT, message) > 0) &&
                          (double_check_lib ? compoundsBelong(velements, _aflowlib_tmp.vspecies, FileMESSAGE, oss, vpflow.flag("PFLOW::LOAD_ENTRIES_NON_ALPHABETICAL")) : true) // sometimes we find odd entries in the wrong LIBS, better to be safe
                          )
                                     :
                                     // if load_from_api, check this bool
                                     // this part of the ternary operator loads in the url
                                     // the last part checks that the loaded entry has a set of species that belong to the original query (protection)
                      ((_aflowlib_tmp.url2aflowlib(load_path, message, false) > 0) &&
                       (double_check_lib ? compoundsBelong(velements, _aflowlib_tmp.vspecies, FileMESSAGE, oss, vpflow.flag("PFLOW::LOAD_ENTRIES_NON_ALPHABETICAL")) : true) // sometimes we find odd entries in the wrong LIBS, better to be safe
                       )) {
                nary = _aflowlib_tmp.vspecies.size();
                if (entries.size() < nary) {
                  continue;
                } // this entry is garbage (wrong directory)
                entries[nary - 1].push_back(_aflowlib_tmp);
                total_count++;
                if (vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_XSTRUCTURES")) {
                  if (!loadXstructures(entries[nary - 1].back(), FileMESSAGE, oss, vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_XSTRUCTURES_RELAXED_ONLY"))) { // CO20171202 - let it decided whether to load from common or not
                    entries[nary - 1].pop_back();
                    total_count--;
                  }
                }
              }
              message.str("");
            }
          }
        }
      }
    }
    //////////////////////////////////////////////////////////////////////////
    // END Finding matches
    //////////////////////////////////////////////////////////////////////////

    // Raises PureA and PureB flags
    if (velements.size() == 2) {
      if (!entries.empty() && !entries[0].empty()) {
        for (size_t j = 0; j < entries[0].size(); j++) {
          if (entries[0][j].vstoichiometry.size() == 1) {
            if (entries[0][j].species == velements[0]) {
              entries[0][j].pureA = true;
            } else if (entries[0][j].species == velements[1]) {
              entries[0][j].pureB = true;
            }
          }
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////////
    // START Print loaded entries summary
    //////////////////////////////////////////////////////////////////////////////

    message << "Loaded " << total_count << " entries";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);

    //////////////////////////////////////////////////////////////////////////////
    // END Print loaded entries summary
    //////////////////////////////////////////////////////////////////////////////

    return true; // entries;
  }
} // namespace pflow

namespace pflow {
  // ***************************************************************************
  // pflow::getCombination(vector<string>& velements,uint nary)
  // ***************************************************************************
  // for given set of elements, will return nary combinations
  // binary combinations of MnPdPt: MnPd, MnPt, PdPt

  vector<vector<string>> elementalCombinations(const vector<string>& velements, uint nary) {
    // for given set of elements, will return nary combinations
    // binary combinations of MnPdPt: MnPd, MnPt, PdPt
    const bool LDEBUG = (false || XHOST.DEBUG);
    vector<vector<string>> combos;
    aurostd::xcombos xc(velements.size(), nary, 'C');
    while (xc.increment()) {
      combos.push_back(xc.applyCombo(velements));
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " combos.back()=";
        for (size_t i = 0; i < combos.back().size(); i++) {
          cerr << combos.back()[i];
        }
        cerr << endl;
      }
    }
    return combos;
  }
} // namespace pflow

namespace pflow {
  // ***************************************************************************
  // pflow::compoundsBelong(vector<string>& velements, vector<string>& elements)
  // ***************************************************************************
  // let's you know if input (or elements) belongs on hull of velements
  // if sort_elements==True, will sort(elements) first before matching with velements,
  // sorting is generally NOT preferred as it would match unsorted compounds from database (NOT good)

  // vector
  bool compoundsBelong(const vector<string>& velements2check, const string& input, ostream& oss, bool clean, bool sort_elements, elements_string_type e_str_type, bool shortcut_pp_string_AFLOW_database) {
    ofstream FileMESSAGE;
    return compoundsBelong(velements2check, input, FileMESSAGE, oss, clean, sort_elements, e_str_type, shortcut_pp_string_AFLOW_database);
  }
  bool compoundsBelong(const vector<string>& velements2check, const string& input, vector<string>& input_velements, vector<double>& input_vcomposition, ostream& oss, bool clean, bool sort_elements, elements_string_type e_str_type, bool shortcut_pp_string_AFLOW_database) {
    ofstream FileMESSAGE;
    return compoundsBelong(velements2check, input, input_velements, input_vcomposition, FileMESSAGE, oss, clean, sort_elements, e_str_type, shortcut_pp_string_AFLOW_database);
  }
  bool compoundsBelong(const vector<string>& velements2check, const string& input, ofstream& FileMESSAGE, ostream& oss, bool clean, bool sort_elements, elements_string_type e_str_type, bool shortcut_pp_string_AFLOW_database) {
    vector<string> input_velements;
    vector<double> input_vcomposition;
    return compoundsBelong(velements2check, input, input_velements, input_vcomposition, FileMESSAGE, oss, clean, sort_elements, e_str_type, shortcut_pp_string_AFLOW_database);
  }
  bool compoundsBelong(const vector<string>& velements2check,
                       const string& input,
                       vector<string>& input_velements,
                       vector<double>& input_vcomposition,
                       ofstream& FileMESSAGE,
                       ostream& oss,
                       bool clean,
                       bool sort_elements,
                       elements_string_type e_str_type,
                       bool shortcut_pp_string_AFLOW_database) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (e_str_type == pp_string && shortcut_pp_string_AFLOW_database == true) {
      // pp_string parsing is slow because of LONG list of strings to substitute in VASP_PseudoPotential_CleanName()
      // instead, we are safe with faster composition_string parsing IFF we only remove _GW, which will confuse elementsFromCompositionString()
      // all other characters are eliminated because we only look for A-Z a-z from substrings created from capital letters
      // these must be directory strings from REST-API/AFLUX
      // otherwise, we also need to be careful of "potpaw_PBE/Mn_pv" species (for example) inside aflow.in (LEGACY)
      // compoundsBelong is only used for loadEntries() (as of 20190712)
      string input_new = input;
      // aurostd::RemoveSubStringInPlace(input_new,"_GW");  //CO20190712
      aurostd::VASP_PseudoPotential_CleanName_InPlace(input_new, true); // capital_letters_only==true
      input_velements = aurostd::getElements(input_new, input_vcomposition, composition_string, FileMESSAGE, clean, sort_elements, false, oss); // use composition_string (FASTER) //do not keep_pp
    } else { // default
      input_velements = aurostd::getElements(input, input_vcomposition, e_str_type, FileMESSAGE, clean, sort_elements, false, oss); // do not keep_pp
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " input=\"" << input << "\", elements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(input_velements, "\""), ",") << endl;
    }
    return compoundsBelong(velements2check, input_velements, FileMESSAGE, oss, false); // already sorted
  }
  bool compoundsBelong(const vector<string>& velements2check, const vector<string>& input_velements, ostream& oss, bool sort_elements) {
    ofstream FileMESSAGE;
    return compoundsBelong(velements2check, input_velements, FileMESSAGE, oss, sort_elements);
  }
  bool compoundsBelong(const vector<string>& velements2check, const vector<string>& input_velements, ofstream& FileMESSAGE, ostream& oss, bool sort_elements) {
    stringstream message;
    if (velements2check.empty() || input_velements.empty()) { // null case, simply return false
      message << "Invalid input (velements2check.size()==" << velements2check.size() << ",input_velements.size()==" << input_velements.size() << ")";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;
    }
    if (input_velements.size() > velements2check.size()) {
      return false;
    } // fast return
    if (sort_elements) { // call recursion if sort needed
      vector<string> velements2check_sorted;
      vector<string> input_velements_sorted;
      for (size_t i = 0; i < velements2check.size(); i++) {
        velements2check_sorted.push_back(velements2check[i]);
      }
      for (size_t i = 0; i < input_velements.size(); i++) {
        input_velements_sorted.push_back(input_velements[i]);
      }
      sort(velements2check_sorted.begin(), velements2check_sorted.end());
      sort(input_velements_sorted.begin(), input_velements_sorted.end());
      return compoundsBelong(velements2check_sorted, input_velements_sorted, FileMESSAGE, oss, false);
    }
    // check if all input_velements are in velements2check (in order)
    bool found;
    uint starting_index = 0; // ensures we search in order
    for (size_t i = 0; i < input_velements.size(); i++) {
      found = false;
      for (size_t j = starting_index; j < velements2check.size() && !found; j++) {
        if (input_velements[i] == velements2check[j]) {
          found = true;
          starting_index = j + 1;
        } // start search at the next velement
      }
      if (!found) {
        return false;
      }
    }
    return true;
  }

} // namespace pflow

// functions for loading entries
namespace pflow {
  // ***************************************************************************
  // pflow::loadXstructures(aflowlib::_aflowlib_entry& entry,string path,bool
  // relaxed_only,ostream& oss,ofstream& FileMESSAGE)
  // ***************************************************************************
  // loads xstructures
  bool loadXstructures(aflowlib::_aflowlib_entry& entry, ostream& oss, bool relaxed_only, string path, bool is_url_path) {
    ofstream FileMESSAGE;
    vector<string> structure_files; // DX20200224
    return loadXstructures(entry, structure_files, FileMESSAGE, oss, relaxed_only, path, is_url_path); // DX20200224
  }
  bool loadXstructures(aflowlib::_aflowlib_entry& entry, ofstream& FileMESSAGE, ostream& oss, bool relaxed_only, string path, bool is_url_path) { // DX20200224
    vector<string> structure_files; // DX20200224
    return loadXstructures(entry, structure_files, FileMESSAGE, oss, relaxed_only, path, is_url_path); // DX20200224
  }
  bool loadXstructures(aflowlib::_aflowlib_entry& entry, vector<string>& structure_files, ofstream& FileMESSAGE, ostream& oss, bool relaxed_only, string path, bool is_url_path) { // DX20200224 - added structure_files as input
    const bool LDEBUG = (false || XHOST.DEBUG);
    stringstream message;

    // get path if not provided
    if (path.empty()) {
      path = entry.getPathAURL(FileMESSAGE, oss, true);
      is_url_path = false;
      if (!(!path.empty() && aurostd::IsDirectory(path))) {
        path = "";
      } // reset
    }
    if (path.empty()) {
      path = entry.getPathAURL(FileMESSAGE, oss, false);
      is_url_path = true;
    }
    if (path.empty()) {
      message << "Path cannot be loaded from entry, skipping loadXstructure()";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_WARNING_);
      return false;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " path=" << path << endl;
    }

    vector<string> vspecies = entry.vspecies;
    //[CO20221110 - fails for unaries in LIB2]if(vspecies.empty()){vspecies=entry.getSpeciesAURL(FileMESSAGE,oss);}
    if (vspecies.empty()) {
      vspecies = entry.getSpecies();
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vspecies=" << aurostd::joinWDelimiter(vspecies, ",") << endl;
    }

    xstructure xstrAux;
    stringstream ss;
    vector<string> files;

    if (is_url_path) {
      aurostd::httpGetTokens(path + "/?files", files, ",");
    } else {
      aurostd::DirectoryLS(path, files);
    }

    // CO20200404 - creating vectors for for-loops to reduce code
    const vector<string> poscar_files_orig{"POSCAR.orig", "POSCAR.relax1"}; // flipped order so orig is preferred if available //POSCAR.relax1,POSCAR.orig
    const vector<string> poscar_files_relax_mid{"POSCAR.relax2", "CONTCAR.relax1"};
    // flipped order for speed and likelihood of existence //POSCAR.bands,CONTCAR.bands,POSCAR.static,CONTCAR.static,CONTCAR.relax2,CONTCAR.relax
    const vector<string> poscar_files_relax_final{"CONTCAR.relax", "CONTCAR.relax2", "POSCAR.static", "POSCAR.bands", "CONTCAR.static", "CONTCAR.bands"};

    // CO20200223 - substring2bool - > CompressFileWithinList()
    string efile;
    if ((!aurostd::CompressFileWithinList(files, "POSCAR.relax1", efile) && !aurostd::CompressFileWithinList(files, "POSCAR.orig", efile) && !relaxed_only) ||
        (!aurostd::CompressFileWithinList(files, "POSCAR.relax2", efile) && !aurostd::CompressFileWithinList(files, "CONTCAR.relax1", efile) && !relaxed_only) ||
        (!aurostd::CompressFileWithinList(files, "POSCAR.bands", efile) && !aurostd::CompressFileWithinList(files, "CONTCAR.bands", efile) && !aurostd::CompressFileWithinList(files, "POSCAR.static", efile) &&
         !aurostd::CompressFileWithinList(files, "CONTCAR.static", efile) && !aurostd::CompressFileWithinList(files, "CONTCAR.relax2", efile) && !aurostd::CompressFileWithinList(files, "CONTCAR.relax", efile) && true)) {
      if (entry.prototype.find("POCC") != string::npos) { // POCC entries have no composition
        message << "path=" << path << " is a POCC-structure and has no structure files. Ignoring entry";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        return false;
      } else {
        message << "path=" << path << " missing structure files. Ignoring entry";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_WARNING_);
        return false;
      }
    }
    if (!relaxed_only) {
      //////////////////////////////////////////////////////////////////////////
      // START Loading original structures
      //////////////////////////////////////////////////////////////////////////

      for (size_t i = 0; i < poscar_files_orig.size() && xstrAux.atoms.empty(); i++) {
        if (aurostd::CompressFileWithinList(files, poscar_files_orig[i], efile)) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " looking for " << path + "/" + efile << endl;
          }
          ss.str("");
          if (is_url_path) {
            ss << aurostd::httpGetCompressedFileContent(path + "/" + efile);
          } else {
            aurostd::compressfile2stringstream(path + "/" + efile, ss);
          }
          if (!ss.str().empty()) {
            try {
              xstrAux = xstructure(ss, IOVASP_AUTO);
              if (xstrAux.GetElementsFromAtomNames().empty()) {
                if (vspecies.empty()) {
                  throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "could not extract species from AURL", _INPUT_ERROR_);
                }
                xstrAux.SetSpecies(aurostd::vector2deque(vspecies));
              }
              if (LDEBUG) {
                cerr << __AFLOW_FUNC__ << " xstrAux=" << endl;
                cerr << xstrAux << endl;
              }
              structure_files.push_back(poscar_files_orig[i]);
            } // DX20191210 - added try-catch //DX20200224 - added structure_files tag
            catch (aurostd::xerror& excpt) {
              xstrAux.clear(); // clear it if it is garbage //DX20191220 - uppercase to lowercase clear
              message << poscar_files_orig[i] << ": Path exists, but could not load structure (e.g., URL timeout or bad structure file).";
              pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, _LOGGER_WARNING_);
            } // DX20191210
          }
        }
      }
      if (xstrAux.atoms.empty()) {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Cannot load original structure", FileMESSAGE, oss, _LOGGER_WARNING_);
        return false;
      }
      entry.vstr.push_back(xstrAux);
      xstrAux.clear(); // DX20191220 - uppercase to lowercase clear
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " loaded ORIGINAL structure" << endl;
      }

      //////////////////////////////////////////////////////////////////////////
      // END Loading original structures
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      // START Loading singularly-relaxed structures
      //////////////////////////////////////////////////////////////////////////

      for (size_t i = 0; i < poscar_files_relax_mid.size() && xstrAux.atoms.empty(); i++) {
        if (aurostd::CompressFileWithinList(files, poscar_files_relax_mid[i], efile)) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " looking for " << path + "/" + efile << endl;
          }
          ss.str("");
          if (is_url_path) {
            ss << aurostd::httpGetCompressedFileContent(path + "/" + efile);
          } else {
            aurostd::compressfile2stringstream(path + "/" + efile, ss);
          }
          if (!ss.str().empty()) {
            try {
              xstrAux = xstructure(ss, IOVASP_AUTO);
              if (xstrAux.GetElementsFromAtomNames().empty()) {
                if (vspecies.empty()) {
                  throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "could not extract species from AURL", _INPUT_ERROR_);
                }
                xstrAux.SetSpecies(aurostd::vector2deque(vspecies));
              }
              if (LDEBUG) {
                cerr << __AFLOW_FUNC__ << " xstrAux=" << endl;
                cerr << xstrAux << endl;
              }
              structure_files.push_back(poscar_files_relax_mid[i]);
            } // DX20191210 - added try-catch //DX20200224 - added structure_files tag
            catch (aurostd::xerror& excpt) {
              xstrAux.clear(); // clear it if it is garbage //DX20191220 - uppercase to lowercase clear
              message << poscar_files_relax_mid[i] << ": Path exists, but could not load structure (e.g., URL timeout or bad structure file).";
              pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, _LOGGER_WARNING_);
            } // DX20191210
          }
        }
      }
      if (xstrAux.atoms.empty()) {
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Cannot load mid-relaxed structure", FileMESSAGE, oss, _LOGGER_WARNING_);
        return false;
      }
      entry.vstr.push_back(xstrAux);
      xstrAux.clear(); // DX20191220 - uppercase to lowercase clear
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " loaded RELAX1 structure" << endl;
      }

      //////////////////////////////////////////////////////////////////////////
      // END Loading singularly-relaxed structures
      //////////////////////////////////////////////////////////////////////////
    }

    ////////////////////////////////////////////////////////////////////////////
    // START Loading fully-relaxed structures
    ////////////////////////////////////////////////////////////////////////////

    for (size_t i = 0; i < poscar_files_relax_final.size() && xstrAux.atoms.empty(); i++) {
      if (aurostd::CompressFileWithinList(files, poscar_files_relax_final[i], efile)) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " looking for " << path + "/" + efile << endl;
        }
        ss.str("");
        if (is_url_path) {
          ss << aurostd::httpGetCompressedFileContent(path + "/" + efile);
        } else {
          aurostd::compressfile2stringstream(path + "/" + efile, ss);
        }
        if (!ss.str().empty()) {
          try {
            xstrAux = xstructure(ss, IOVASP_AUTO);
            if (xstrAux.GetElementsFromAtomNames().empty()) {
              if (vspecies.empty()) {
                throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "could not extract species from AURL", _INPUT_ERROR_);
              }
              xstrAux.SetSpecies(aurostd::vector2deque(vspecies));
            }
            if (LDEBUG) {
              cerr << __AFLOW_FUNC__ << " xstrAux=" << endl;
              cerr << xstrAux << endl;
            }
            structure_files.push_back(poscar_files_relax_final[i]);
          } // DX20191210 - added try-catch //DX20200224 - added structure_files tag
          catch (aurostd::xerror& excpt) {
            xstrAux.clear(); // clear it if it is garbage //DX20191220 - uppercase to lowercase clear
            message << poscar_files_relax_final[i] << ": Path exists, but could not load structure (e.g., URL timeout or bad structure file).";
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, _LOGGER_WARNING_);
          } // DX20191210
        }
      }
    }
    if (xstrAux.atoms.empty()) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Cannot load fully-relaxed structure", FileMESSAGE, oss, _LOGGER_WARNING_);
      return false;
    }
    entry.vstr.push_back(xstrAux);
    xstrAux.clear(); // DX20191220 - uppercase to lowercase clear
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " loaded FULLY-RELAXED structure" << endl;
    }

    ////////////////////////////////////////////////////////////////////////////
    // END Loading fully-relaxed structures
    ////////////////////////////////////////////////////////////////////////////

    return true;
  }

  /// @brief load a xstructure directly from a AFLOW lib entry
  /// @param entry AFLOW lib entry
  /// @param new_structure save-to structure
  /// @return xstructure
  /// @note this is always the final structure
  /// @note does not add the structure to entry.vstr
  bool loadXstructureLibEntry(const aflowlib::_aflowlib_entry& entry, xstructure& new_structure) { // HE20220606
    if (entry.positions_cartesian.empty()) {
      return false;
    }
    new_structure.title = "Relaxed AFLUX: " + entry.aurl;
    new_structure.scale = 1.0;
    new_structure.neg_scale = false;
    new_structure.iomode = IOAFLOW_AUTO;
    std::vector<double> geometry;
    aurostd::string2tokens(entry.geometry, geometry, ",");
    new_structure.a = geometry[0];
    new_structure.b = geometry[1];
    new_structure.c = geometry[2];
    new_structure.alpha = geometry[3];
    new_structure.beta = geometry[4];
    new_structure.gamma = geometry[5];
    new_structure.lattice = GetClat(new_structure.a, new_structure.b, new_structure.c, new_structure.alpha, new_structure.beta, new_structure.gamma);
    new_structure.spacegroupnumber = entry.spacegroup_relax;
    new_structure.spacegrouplabel = GetSpaceGroupLabel(new_structure.spacegroupnumber);
    if (new_structure.spacegroupnumber > 0 && new_structure.spacegroupnumber <= 230) {
      new_structure.spacegroup = GetSpaceGroupName(new_structure.spacegroupnumber, new_structure.directory);
    } else {
      new_structure.spacegroup = "";
    }
    new_structure.coord_flag = _COORDS_FRACTIONAL_;
    strcpy(new_structure.coord_type, "D");
    new_structure.FixLattices();

    std::vector<int> composition;
    aurostd::string2tokens(entry.composition, composition, ",");
    uint num_atoms = 0;
    std::vector<std::vector<double>> positions_fractional;
    std::vector<double> positions_fractional_raw;
    std::vector<std::string> positions_fractional_raw_str;
    aurostd::string2tokens(entry.positions_fractional, positions_fractional_raw_str, ";");
    for (std::vector<std::string>::const_iterator pf_str = positions_fractional_raw_str.begin(); pf_str != positions_fractional_raw_str.end(); pf_str++) {
      positions_fractional_raw.clear();
      aurostd::string2tokens(*pf_str, positions_fractional_raw, ",");
      positions_fractional.emplace_back(positions_fractional_raw);
    }

    std::vector<std::string> species;
    aurostd::string2tokens(entry.species, species, ",");
    xvector<double> v(3);
    _atom atom;
    for (size_t type_idx = 0; type_idx < composition.size(); type_idx++) {
      for (int atom_idx = 0; atom_idx < composition[type_idx]; atom_idx++) {
        atom.clear();
        v.clear();
        v(1) = positions_fractional[num_atoms][0];
        v(2) = positions_fractional[num_atoms][1];
        v(3) = positions_fractional[num_atoms][2];
        atom.fpos = v;
        atom.cpos = new_structure.f2c * atom.fpos;
        atom.basis = num_atoms;
        atom.ijk(1) = 0;
        atom.ijk(2) = 0;
        atom.ijk(3) = 0;
        atom.corigin(1) = 0.0;
        atom.corigin(2) = 0.0;
        atom.corigin(3) = 0.0; // inside the zero cell
        atom.coord(1) = 0.0;
        atom.coord(2) = 0.0;
        atom.coord(3) = 0.0; // inside the zero cell
        atom.spin = 0.0;
        atom.noncoll_spin.clear();
        atom.type = type_idx;
        atom.order_parameter_value = 0;
        atom.order_parameter_atom = false;
        atom.partial_occupation_value = 1.0;
        atom.partial_occupation_flag = false;
        atom.name = species[type_idx];
        atom.cleanname = species[type_idx];
        atom.name_is_given = true;
        new_structure.AddAtom(atom, false); // CO20230319 - add by type
        num_atoms++;
      }
    }
    return true;
  }
} // namespace pflow

namespace pflow {
  //***************************************************************************//
  // pflow::defaultLoadEntriesFlags(const string& input)
  //***************************************************************************//
  // sets some default flags
  void defaultLoadEntriesFlags(aurostd::xoption& vpflow, ostream& oss, string input, bool entry_output, bool silent) { // overload
    ofstream FileMESSAGE;
    return defaultLoadEntriesFlags(vpflow, FileMESSAGE, oss, input, entry_output, silent);
  }
  void defaultLoadEntriesFlags(aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& oss, string input, bool entry_output, bool silent) { // main function
    stringstream message;

    // get directory info
    _aflags aflags;
    aflags.Directory = chull::getPath(); // default
    if (vpflow.flag("CHULL::PATH")) {
      aflags.Directory = chull::getPath(vpflow, FileMESSAGE, oss);
    }

    // common
    if (entry_output) {
      vpflow.flag("PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT", true);
      if (!silent) {
        message << "PFLOW::LOAD_ENTRIES_ENTRY_OUTPUT set to true";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_OPTION_);
      }
    }
    //[CO20190715 - LOAD_ENTRIES_ONLY_ALPHABETICAL -> LOAD_ENTRIES_NON_ALPHABETICAL]vpflow.flag("PFLOW::LOAD_ENTRIES_ONLY_ALPHABETICAL", true);  // un-alphabetized entries are crap
    if (!silent) {
      message << "PFLOW::LOAD_ENTRIES_NON_ALPHABETICAL set to true";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_OPTION_);
    }
    vpflow.flag("PFLOW::LOAD_ENTRIES_NARIES_MINUS_ONE", true); // if loading ternary, also load relevant binaries and unaries
    if (!silent) {
      message << "PFLOW::LOAD_ENTRIES_NARIES_MINUS_ONE set to true";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_OPTION_);
    }
    input = aurostd::toupper(input); // removes ALL ambiguity with case
    // must be specific LIB1, LIB2, etc.
    if (!(input == "A" || input == "ICSD")) //|| input == "icsd"))  //put other special libs here
    { // CO20200106 - patching for auto-indenting
      string lib_name;
      string load_lib_flag_name;
      // make robust so input can be "LIB2" or "2"
      if (aurostd::substring2bool(input, "LIB")) {
        lib_name = input;
      } else {
        lib_name = "LIB" + input;
      }
      load_lib_flag_name = "PFLOW::LOAD_ENTRIES_LOAD_" + lib_name;
      vpflow.flag(load_lib_flag_name, true);
      if (!silent) {
        message << load_lib_flag_name << " set to true";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_OPTION_);
      }
    } else {
      if (input == "A") {
        // get appropriate size, slightly inefficient (as we got this before), but it's cheap
        // string system = vpflow.getattachedscheme("PFLOW::ALLOY");  //CO20170908 - don't want to have to set this everytime
        // vector<string> velements = aurostd::getElements(system, oss, FileMESSAGE);  //un-sorted, okay
        string lib_count_string;
        string load_lib_flag_name;
        // for (size_t i = 0; i < velements.size() && i <= _AFLOW_LIB_MAX_; i++) //CO20170908 - simply load all, LoadEntries() limits appropriately by velements.size()
        for (uint lib = 1; lib <= _AFLOW_LIB_MAX_; lib++) { // CO20200106 - patching for auto-indenting
          // if(i == 0) { continue; }  //skip LIB1 by default  //CO20180316 - DON'T skip LIB1, Mn unary can be very skewed
          lib_count_string = aurostd::utype2string(lib);
          load_lib_flag_name = "PFLOW::LOAD_ENTRIES_LOAD_LIB" + lib_count_string;
          vpflow.flag(load_lib_flag_name, true);
          if (!silent) {
            message << load_lib_flag_name << " set to true";
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_OPTION_);
          }
        }
      }
      if (input == "A" || input == "ICSD") //|| input == "icsd")
      { // CO20200106 - patching for auto-indenting
        vpflow.flag("PFLOW::LOAD_ENTRIES_LOAD_ICSD", true);
        if (!silent) {
          message << "PFLOW::LOAD_ENTRIES_LOAD_ICSD set to true";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_OPTION_);
        }
      }
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::prototypeMatch(string p1,string p2)
// ***************************************************************************
// better than an exact match of proto_database and proto_search
// we match for 549 and 549.bis
// proto_database is "549.bis" and proto_search is "549"
// only knows "." and ":" for now
namespace pflow {
  bool prototypeMatch(string proto_database, string proto_search) {
    if (proto_database == proto_search) {
      return true;
    }
    if (proto_database.size() < (proto_search.size() + 1)) {
      return false;
    } // we only match "549" + something
    if (proto_database.substr(0, proto_search.size()) != proto_search) {
      return false;
    } // not something + "549"
    if (proto_database[proto_search.size()] == '.' || proto_database[proto_search.size()] == ':') {
      return true;
    }
    // cerr << "What is this: " << proto_database << endl;
    return false;
  }
} // namespace pflow
// Added by Corey Oses - May 2017
// END - all relevent functions for loading entries here
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace pflow {

  // ME20200512 - Created by CO for POCC.
  const int BAR_WIDTH = 70;
  void updateProgressBar(unsigned long long int current, unsigned long long int end, ostream& oss) {
    // CO20220630 - decide if we should print
    // this is one of the few functions that escapes pflow::logger() and aurostd::PrintXXStream() for printing
    // all other printing should route through pflow::logger() (routing directly through aurostd::PrintXXStream() is obsolete)
    if (XHOST.QUIET || XHOST.QUIET_GLOBAL || XHOST.vflag_control.flag("WWW")) {
      return;
    } // CO20190520 - no progress bar for web stuff  //CO20200404 - new web flag // ME20210428 - do not update when quiet either
    if ((&oss == &cout) && XHOST.QUIET_COUT) {
      return;
    }
    if ((&oss == &cerr) && XHOST.QUIET_CERR) {
      return;
    }

    const double progress = (double) current / (double) end;
    const int pos = BAR_WIDTH * progress;

    oss << "[";
    for (int i = 0; i < BAR_WIDTH; ++i) {
      if (i < pos) {
        oss << "=";
      } else if (i == pos) {
        oss << ">";
      } else {
        oss << " ";
      }
    }
    oss << "] " << int(progress * 100.0) << " %\r";
    oss.flush();
    if (current == end) {
      oss << endl;
    }
  }

  //***************************************************************************//
  // pflow::logger(const char& type,const string& function_name,string
  // _message,bool silent,ostream& oss,ofstream& FileMESSAGE)
  //***************************************************************************//
  // added by Corey Oses - May 2017
  // effectively logs EVERYTHING, deals with colors, oss, and log files (very robust)
  // types include error, warning, options, and more as needed (see below)
  // raw is special, prints raw message
  // message can be string or stringstream (if stringstream, it gets cleared out)
  // oss = cout, cerr (as you prefer)
  // FileMESSAGE - logfile
  void logger(const string& filename, const string& function_name, stringstream& message, const char& type, ostream& oss, bool silent, const string& message_metadata) { // overload
    const string _message = message.str();
    logger(filename, function_name, _message, oss, type, silent, message_metadata);
    aurostd::StringstreamClean(message);
  }

  void logger(const string& filename, const string& function_name, stringstream& message, ostream& oss, const char& type, bool silent, const string& message_metadata) { // overload
    const string _message = message.str();
    logger(filename, function_name, _message, oss, type, silent, message_metadata);
    aurostd::StringstreamClean(message);
  }

  void logger(const string& filename, const string& function_name, stringstream& message, ofstream& FileMESSAGE, ostream& oss, const char& type, bool silent, const string& message_metadata) { // overload
    const string _message = message.str();
    logger(filename, function_name, _message, FileMESSAGE, oss, type, silent, message_metadata);
    aurostd::StringstreamClean(message);
  }

  void logger(const string& filename, const string& function_name, stringstream& message, const string& directory, ostream& oss, const char& type, bool silent, const string& message_metadata) { // overload
    const string _message = message.str();
    logger(filename, function_name, _message, directory, oss, type, silent, message_metadata);
    aurostd::StringstreamClean(message);
  }

  void logger(const string& filename, const string& function_name, stringstream& message, const string& directory, ofstream& FileMESSAGE, ostream& oss, const char& type, bool silent, const string& message_metadata) { // overload
    const string _message = message.str();
    logger(filename, function_name, _message, directory, FileMESSAGE, oss, type, silent, message_metadata);
    aurostd::StringstreamClean(message);
  }

  void logger(const string& filename, const string& function_name, stringstream& message, const _aflags& aflags, ostream& oss, const char& type, bool silent, const string& message_metadata) { // overload
    const string _message = message.str();
    logger(filename, function_name, _message, aflags, oss, type, silent, message_metadata);
    aurostd::StringstreamClean(message);
  }

  void logger(const string& filename, const string& function_name, stringstream& message, const _aflags& aflags, ofstream& FileMESSAGE, ostream& oss, const char& type, bool silent, const string& message_metadata) { // overload
    const string _message = message.str();
    logger(filename, function_name, _message, aflags, FileMESSAGE, oss, type, silent, message_metadata);
    aurostd::StringstreamClean(message);
  }

  void logger(const string& filename, const string& function_name, const string& _message, const char& type, ostream& oss, bool silent, const string& message_metadata) { // overload
    ofstream FileMESSAGE;
    logger(filename, function_name, _message, FileMESSAGE, oss, type, silent, message_metadata);
  }

  void logger(const string& filename, const string& function_name, const string& _message, ostream& oss, const char& type, bool silent, const string& message_metadata) { // overload
    ofstream FileMESSAGE;
    logger(filename, function_name, _message, FileMESSAGE, oss, type, silent, message_metadata);
  }

  void logger(const string& filename, const string& function_name, const string& _message, ofstream& FileMESSAGE, ostream& oss, const char& type, bool silent, const string& message_metadata) { // overload
    _aflags aflags;
    if (XHOST.vflag_control.flag("DIRECTORY_CLEAN")) {
      aflags.Directory = XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN");
    } // CO20190402
    logger(filename, function_name, _message, aflags, FileMESSAGE, oss, type, silent, message_metadata);
  }

  void logger(const string& filename, const string& function_name, const string& _message, const string& directory, ostream& oss, const char& type, bool silent, const string& message_metadata) { // overload
    ofstream FileMESSAGE;
    logger(filename, function_name, _message, directory, FileMESSAGE, oss, type, silent, message_metadata);
  }

  void logger(const string& filename, const string& function_name, const string& _message, const _aflags& aflags, ostream& oss, const char& type, bool silent, const string& message_metadata) { // overload
    ofstream FileMESSAGE;
    logger(filename, function_name, _message, aflags, FileMESSAGE, oss, type, silent, message_metadata);
  }

  void logger(const string& filename, const string& function_name, const string& _message, const string& directory, ofstream& FileMESSAGE, ostream& oss, const char& type, bool silent, const string& message_metadata) { // overload
    _aflags aflags;
    aflags.Directory = directory;
    logger(filename, function_name, _message, aflags, FileMESSAGE, oss, type, silent, message_metadata);
  }

  void logger(const string& filename, const string& function_name, const string& _message, const _aflags& aflags, ofstream& FileMESSAGE, ostream& oss, const char& type, bool silent, const string& message_metadata) { // main function
    // five types:  M - Message, W - Warning, E - Error, O - Option, C - Complete, R - RAW
    // treat E, W, O, C, and R as special, otherwise treat as a message
    // no need for function name for R, you can put ""

    const string message = aurostd::RemoveWhiteSpacesFromTheBack(_message);
    if (message.empty()) {
      return;
    }
    // CO20181226 - split by newlines and print separately
    vector<string> message_parts;
    vector<string> _message_parts;
    aurostd::string2vectorstring(message, _message_parts);
    for (size_t i = 0; i < _message_parts.size(); i++) {
      if (!aurostd::RemoveWhiteSpacesFromTheBack(_message_parts[i]).empty()) {
        message_parts.push_back(_message_parts[i]);
      }
    }
    if (message_parts.empty()) {
      return;
    }

    string tag_code = "00000";
    string tag_message = "MESSAGE";
    if (type == _LOGGER_ERROR_) {
      tag_code = "EEEEE";
      tag_message = "ERROR";
    } else if (type == _LOGGER_WARNING_) {
      tag_code = "WWWWW";
      tag_message = "WARNING";
    } else if (type == _LOGGER_COMPLETE_) {
      tag_code = "CCCCC";
      tag_message = "COMPLETE";
    } else if (type == _LOGGER_OPTION_) {
      tag_code = "-OPT-";
      tag_message = "MESSAGE-OPTION";
    } else if (type == _LOGGER_RAW_) {
      tag_code = "";
      tag_message = "";
    } else if (type == _LOGGER_NOTICE_) { // ME20200514
      tag_code = "00000";
      tag_message = "NOTICE";
    }

    ostringstream stream;
    if (type == _LOGGER_RAW_) {
      for (size_t i = 0; i < message_parts.size(); i++) {
        stream << message_parts[i] << endl;
      } // CO20181226 //message;
    } else if (XHOST.vflag_control.flag("WWW")) { // ME20210429 - remove metadata for web output
      for (size_t i = 0; i < message_parts.size(); i++) {
        stream << tag_code;
        stream << "  ";
        stream << tag_message;
        stream << " ";
        stream << message_parts[i];
        stream << endl;
      }
    } else {
      for (size_t i = 0; i < message_parts.size(); i++) {
        stream << XPID; // CO20200524
        stream << tag_code;
        stream << "  ";
        stream << tag_message;
        stream << " ";
        stream << function_name; // CO20230319
        stream << " ";
        stream << message_parts[i]; // CO20181226 //message;
        stream << Message(filename, aflags, message_metadata);
        stream << endl;
      }
    }
    // HE+ME20220503
    //  Be permissive and search for substrings to allow for white/black listing
    //  of groups of functions or namespaces without implementing regexes
    //  The white and black lists should be treated like a stack: only push and pop
    bool quiet = XHOST.QUIET;
    if (XHOST.QUIET) {
      quiet = !aurostd::substringlist2bool(function_name, XHOST.LOGGER_WHITELIST, false);
    } else {
      quiet = aurostd::substringlist2bool(function_name, XHOST.LOGGER_BLACKLIST, false);
    }

    // CO20220630 - note about osswrite, it is redundant with quiet, so it would be nice to get rid of it in the future
    // but it would require a full overhaul of many critical aflow printing functions
    // better not to touch and leave the overloads
    const bool osswrite = !silent;
    if (type == _LOGGER_ERROR_) {
      aurostd::PrintErrorStream(FileMESSAGE, stream, quiet, osswrite);
    } // oss - DEFAULT TO cerr
    else if (type == _LOGGER_WARNING_) {
      aurostd::PrintWarningStream(FileMESSAGE, stream, quiet, osswrite);
    } // oss - DEFAULT TO cerr
    else {
      aurostd::PrintMessageStream(FileMESSAGE, stream, quiet, osswrite, oss);
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::LTCELL
// ***************************************************************************
namespace pflow {
  xstructure LTCELL(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    const xstructure str(input, IOAFLOW_AUTO);
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    const xmatrix<double> mlt(3, 3);

    if (tokens.size() != 9 && tokens.size() != 3 && tokens.size() != 1 && tokens.size() != 4) {
      init::ErrorOption(options, __AFLOW_FUNC__,
                        aurostd::liststring2string("aflow --ltcell=a11,a12,a13,a21,a22,a23,a31,a32,a33 < POSCAR", "aflow --ltcell=a11,a22,a33 < POSCAR", "aflow --ltcell=file < POSCAR",
                                                   "aflow --ltcellfv=v1,v2,v3,phi < POSCAR"));
    }

    if (tokens.size() == 9) {
      if (LDEBUG) {
        cerr << XPID << "pflow::LTCELL: 9 entries" << endl;
      }
      mlt(1, 1) = aurostd::string2utype<double>(tokens[0]);
      mlt(1, 2) = aurostd::string2utype<double>(tokens[1]);
      mlt(1, 3) = aurostd::string2utype<double>(tokens[2]);
      mlt(2, 1) = aurostd::string2utype<double>(tokens[3]);
      mlt(2, 2) = aurostd::string2utype<double>(tokens[4]);
      mlt(2, 3) = aurostd::string2utype<double>(tokens[5]);
      mlt(3, 1) = aurostd::string2utype<double>(tokens[6]);
      mlt(3, 2) = aurostd::string2utype<double>(tokens[7]);
      mlt(3, 3) = aurostd::string2utype<double>(tokens[8]);
      if (std::abs(det(mlt)) < 0.01) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "singular ltcell matrix", _INPUT_ILLEGAL_); // CO20200624
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " END" << endl;
      }
      return GetLTCell(mlt, str);
    }

    if (tokens.size() == 3) {
      if (LDEBUG) {
        cerr << XPID << "pflow::LTCELL: 3 entries" << endl;
      }
      mlt(1, 1) = aurostd::string2utype<double>(tokens[0]);
      mlt(2, 2) = aurostd::string2utype<double>(tokens[1]);
      mlt(3, 3) = aurostd::string2utype<double>(tokens[2]);
      if (std::abs(det(mlt)) < 0.01) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "singular ltcell matrix", _INPUT_ILLEGAL_); // CO20200624
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " END" << endl;
      }
      return GetLTCell(mlt, str);
    }

    if (tokens.size() == 1) {
      if (LDEBUG) {
        cerr << XPID << "pflow::LTCELL: 1 entries" << endl;
      }
      ifstream infile(tokens[0].c_str());
      aurostd::InFileExistCheck("pflow::LTCELL", tokens[0].c_str(), infile);
      for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
          infile >> mlt(i, j);
        }
      }
      if (std::abs(det(mlt)) < 0.01) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "singular ltcell matrix", _INPUT_ILLEGAL_); // CO20200624
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " END" << endl;
      }
      return GetLTCell(mlt, str);
    }

    if (tokens.size() == 4) {
      if (LDEBUG) {
        cerr << XPID << "pflow::LTCELL: 4 entries => LTCELLFV mode" << endl;
      }
      const xvector<double> nvec(3);
      nvec(1) = aurostd::string2utype<double>(tokens[0]);
      nvec(2) = aurostd::string2utype<double>(tokens[1]);
      nvec(3) = aurostd::string2utype<double>(tokens[2]);
      const double angle = aurostd::string2utype<double>(tokens[3]) / rad2deg;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " END" << endl;
      }
      return GetLTFVCell(nvec, angle, str);
    }
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::MagneticParameters
// ***************************************************************************
namespace pflow {
  void MagneticParameters(string _directory, ostream& oss) {
    string directory(_directory);
    if (directory.empty()) {
      directory = "./";
    }
    directory = aurostd::CleanFileName(directory);

    bool found = false;
    const vector<string> vtype{".static", ".relax1", ".relax2"};

    vector<string> vfiles;
    vfiles.clear();
    vfiles.push_back(aurostd::CleanFileName(directory + "/OUTCAR"));
    for (size_t i = 0; i < vtype.size(); i++) {
      for (size_t iext = 0; iext < XHOST.vext.size(); iext++) {
        vfiles.push_back(aurostd::CleanFileName(directory + "/OUTCAR" + vtype[i] + XHOST.vext[iext]));
      }
    }

    xOUTCAR outcar;
    for (size_t i = 0; i < vfiles.size() && !found; i++) {
      if (aurostd::FileExist(vfiles[i])) {
        found = true;
        outcar.GetPropertiesFile(vfiles[i]);
        oss << "pflow::MagneticParameters  : using file=" << vfiles[i] << endl;

        oss << "MAGNETIC MOMENTUM CELL : " << outcar.mag_cell << endl;
        oss << "MAGNETIC MOMENTUM ATOM : " << outcar.mag_atom << endl;
        oss << "VOLUME CELL            : " << outcar.volume_cell << endl;
        oss << "VOLUME ATOM            : " << outcar.volume_atom << endl;
        if (outcar.vmag.size() > 1) {
          oss << "SPIN DECOMPOSITION     : " << outcar.vmag[0];
          for (size_t i = 1; i < outcar.vmag.size(); i++) {
            oss << "," << outcar.vmag[i];
          }
          oss << endl;
        }
      }
    }

    found = false;
    vfiles.clear();
    for (size_t iext = 0; iext < XHOST.vext.size(); iext++) {
      vfiles.push_back(aurostd::CleanFileName(directory + "/DOSCAR.static" + XHOST.vext[iext]));
    }
    vfiles.push_back(aurostd::CleanFileName(directory + "/DOSCAR"));

    xDOSCAR doscar;
    for (size_t i = 0; i < vfiles.size() && !found; i++) {
      if (aurostd::FileExist(vfiles[i])) {
        found = true;
        doscar.GetPropertiesFile(vfiles[i]);
        oss << "pflow::MagneticParameters  : using file=" << vfiles[i] << endl;
        oss << "POLARIZATION FERMI     : " << doscar.spinF << endl;
      }
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::MAKESTRLIST
// ***************************************************************************
namespace pflow {
  void MAKESTRLIST(vector<string> argv) {
    ifstream outcar_inf(argv.at(2).c_str());
    aurostd::InFileExistCheck("convasp", argv.at(2).c_str(), outcar_inf);
    ifstream xdatcar_inf(argv.at(3).c_str());
    aurostd::InFileExistCheck("convasp", argv.at(3).c_str(), xdatcar_inf);
    const vector<xstructure> str_vec = pflow::GetStrVecFromOUTCAR_XDATCAR(outcar_inf, xdatcar_inf);
    pflow::PrintStrVec(str_vec, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::MINKOWSKIBASISREDUCTION
// ***************************************************************************
namespace pflow {
  xstructure MINKOWSKIBASISREDUCTION(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    xstructure b(a);
    b.MinkowskiBasisReduction();
    // b.BringInCell(); // not necessary
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::MISCIBILITY
// ***************************************************************************
namespace pflow {
  string MISCIBILITY(vector<string> argv) {
    stringstream output;
    for (size_t i = 2; i < argv.size(); i++) {
      string system = string(argv[i]);
      XATOM_AlphabetizationSpecies(system);
      //  XATOM_AlphabetizationCompound(system);
      const int mix = MiscibilityCheck(system);
      if (mix == MISCIBILITY_SYSTEM_NOMIX) {
        output << system << " " << "MISCIBILITY_SYSTEM_NOMIX" << endl;
      }
      if (mix == MISCIBILITY_SYSTEM_MISCIBLE) {
        output << system << " " << "MISCIBILITY_SYSTEM_MISCIBLE" << endl;
      }
      if (mix == MISCIBILITY_SYSTEM_UNKNOWN) {
        output << system << " " << "MISCIBILITY_SYSTEM_UNKNOWN" << endl;
      }
    }
    return output.str();
  };
} // namespace pflow

// ***************************************************************************
// pflow::MOM
// ***************************************************************************
namespace pflow {
  void MOM(istream& input) {
    const xstructure str(input, IOAFLOW_AUTO);
    xvector<double> m(3);
    m = GetMom1(str);
    cout.setf(std::ios::left, std::ios::adjustfield);
    cout.setf(std::ios::fixed, std::ios::floatfield);
    cout << " Moment 1 : " << setw(15) << setprecision(10) << m(1) << setw(15) << setprecision(10) << m(2) << setw(15) << setprecision(10) << m(3) << endl;
    m = m / str.scale;
    cout << " (unscaled) " << setw(15) << setprecision(10) << m(1) << setw(15) << setprecision(10) << m(2) << setw(15) << setprecision(10) << m(3) << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::MSI
// ***************************************************************************
namespace pflow {
  void MSI(istream& input) {
    const xstructure str(input, IOAFLOW_AUTO);
    PrintMSI(str, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::NATOMS
// ***************************************************************************
namespace pflow {
  uint NATOMS(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    return a.atoms.size();
  }
} // namespace pflow

// ***************************************************************************
// pflow::NBONDXX
// ***************************************************************************
namespace pflow {
  // CO20171025
  string NBONDXX(istream& input, bool aflowlib_legacy_format) {
    const xstructure a(input, IOAFLOW_AUTO);
    return NBONDXX(a, aflowlib_legacy_format);
  }

  string NBONDXX(const xstructure& a, bool aflowlib_legacy_format) {
    vector<double> nbondxx = GetNBONDXX(a);

    // return aurostd::joinWDelimiter(aurostd::vecDouble2vecString(nbondxx,9),',');
    // print nicely

    // are there names in the atoms?
    bool atom_names = true;
    for (size_t i = 0; i < a.atoms.size() && atom_names; i++) {
      if (a.atoms[i].name.empty()) {
        atom_names = false;
      }
    }

    // get names
    uint iat = 0;
    vector<int> first_itypes;
    first_itypes.push_back(0);
    bool found;
    for (size_t i = 1; i < a.atoms.size(); i++) {
      found = false;
      for (size_t j = 0; j < first_itypes.size() && !found; j++) {
        if (a.atoms[i].type == a.atoms[first_itypes[j]].type) {
          found = true;
        }
      }
      if (!found) {
        first_itypes.push_back(i);
      }
    }

    if (aflowlib_legacy_format) {
      atom_names = false;
    }

    vector<string> names;
    stringstream ss;
    if (atom_names) {
      for (size_t i = 0; i < first_itypes.size(); i++) {
        names.push_back(a.atoms[first_itypes[i]].name);
      }
    } else {
      for (size_t i = 0; i < first_itypes.size(); i++) {
        ss.str("");
        ss << char('A' + iat++);
        names.push_back(ss.str());
      }
    }

    iat = 0;
    stringstream output;
    output.str("");
    if (aflowlib_legacy_format) {
      for (size_t itype = 0; itype < first_itypes.size(); itype++) {
        for (size_t jtype = itype; jtype < first_itypes.size(); jtype++) {
          output << aurostd::PaddedPOST("BOND_" + names[itype] + names[jtype], 11) << " " << aurostd::utype2string(nbondxx[iat++], 6, FIXED_STREAM) << " [Angst]" << endl; // CO20220627
        }
      }
    } else {
      for (size_t itype = 0; itype < first_itypes.size(); itype++) {
        for (size_t jtype = itype; jtype < first_itypes.size(); jtype++) {
          output << aurostd::PaddedPOST(names[itype] + "-" + names[jtype] + ":", 8) << " " << aurostd::utype2string(nbondxx[iat++], 6, FIXED_STREAM) << " [Angst]" << endl; // CO20220627
        }
      }
    }

    return output.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::NNDIST
// ***************************************************************************
namespace pflow {
  double NNDIST(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    return (NearestNeighbor(a));
  }
} // namespace pflow

// ***************************************************************************
// pflow::NSPECIES
// ***************************************************************************
namespace pflow {
  uint NSPECIES(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    return a.num_each_type.size();
  }
} // namespace pflow

// ***************************************************************************
// pflow::NAMES
// ***************************************************************************
namespace pflow {
  xstructure NAMES(vector<string> argv, istream& input) {
    xstructure a(input, IOAFLOW_AUTO);
    if (argv.size() != a.num_each_type.size() + 2) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "you need to specify as many names as atom types", _INPUT_ILLEGAL_); // CO20200624
    }
    xstructure b = a;
    int iatom = 0;
    for (size_t itype = 0; itype < a.num_each_type.size(); itype++) {
      const string species = string(argv.at(2 + b.atoms.at(iatom).type));
      b.species.at(itype) = species;
      for (int j = 0; j < a.num_each_type[itype]; j++) {
        b.atoms.at(iatom).name = species; // CONVASP_MODE
        b.atoms.at(iatom).CleanName();
        b.atoms.at(iatom).CleanSpin();
        b.atoms.at(iatom).name_is_given = true;
        iatom++;
      }
    }
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::NANOPARTICLE
// ***************************************************************************
namespace pflow {
  xstructure NANOPARTICLE(istream& input, const xvector<double>& iparams) {
    const bool LDEBUG = (false || XHOST.DEBUG); // CO20180226
    //  cout << aflow::Banner("BANNER_TINY") << endl;
    double radius = NANOPARTICLE_RADIUS_DEFAULT;
    double distance = NANOPARTICLE_DISTANCE_DEFAULT;
    if (iparams.rows >= 1) {
      radius = iparams[1];
    }
    if (iparams.rows >= 2) {
      distance = iparams[2];
    }
    if (LDEBUG) { // CO20180226
      cerr << __AFLOW_FUNC__ << " radius=" << radius << endl;
      cerr << __AFLOW_FUNC__ << " distance=" << distance << endl;
    }
    xstructure a(input, IOAFLOW_AUTO);
    xstructure b;
    xvector<int> dims(3);
    xvector<double> fshift(3);
    xvector<double> cshift(3);
    xvector<double> ucell(3);
    a.ReScale(1.0);

    b = a; // COPY
    while (!b.atoms.empty()) {
      b.RemoveAtom(0);
    } // EMPTY
    dims = LatticeDimensionSphere(a.lattice, 2 * radius + distance);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " dims=" << dims << endl;
    } // CO20180226
    for (int i = 1; i <= 3; i++) {
      ucell(i) = ceil((2.0 * radius + distance) / modulus(b.lattice(i)));
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " ucell(i=" << i << ")=" << ucell(i) << endl;
      } // CO20180226
      for (int j = 1; j <= 3; j++) {
        b.lattice(i, j) = b.lattice(i, j) * ucell(i);
      }
    }
    b.FixLattices();
    for (size_t iat = 0; iat < a.atoms.size(); iat++) { // SHIFT
      for (int i = -dims(1); i <= dims(1); i++) {
        for (int j = -dims(2); j <= dims(2); j++) {
          for (int k = -dims(3); k <= dims(3); k++) {
            _atom atom = a.atoms[iat];
            atom.cpos = ((double) i) * a.lattice(1) + ((double) j) * a.lattice(2) + ((double) k) * a.lattice(3) + a.atoms[iat].cpos;
            atom.fpos = C2F(b.lattice, atom.cpos); // put in fractional of new basis
            if (modulus(atom.cpos) <= radius) {
              b.AddAtom(atom, false); // CO20230319 - add by species
            }
          }
        }
      }
    }
    cerr << "atoms=" << b.atoms.size() << endl;
    //  b.SetCoordinates(_COORDS_FRACTIONAL_);
    b.SetCoordinates(_COORDS_CARTESIAN_);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::NDATA
// ***************************************************************************
namespace pflow {
  void NDATA(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    PrintNdata(a, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::NIGGLI
// ***************************************************************************
namespace pflow {
  xstructure NIGGLI(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    return GetNiggliStr(a);
  }
} // namespace pflow

// ***************************************************************************
// pflow::NOORDERPARAMETER
// ***************************************************************************
namespace pflow {
  xstructure NOORDERPARAMETER(istream& input) {
    xstructure a(input, IOAFLOW_AUTO);
    a.order_parameter_structure = false;
    a.order_parameter_atoms.clear();
    for (size_t i = 0; i < a.atoms.size(); i++) {
      a.atoms[i].order_parameter_atom = false;
      a.atoms[i].order_parameter_value = 0;
    }
    return a;
  }
} // namespace pflow

// ***************************************************************************
// pflow::NOSD
// ***************************************************************************
namespace pflow {
  xstructure NOSD(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    xstructure b = a;
    // Read in input file.
    b.isd = false;
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::NUMNAMES
// ***************************************************************************
namespace pflow {
  xstructure NUMNAMES(vector<string> argv, istream& input) {
    xstructure a(input, IOAFLOW_AUTO);
    if (argv.size() != a.num_each_type.size() + 2) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "you need to specify as many names as atom types", _INPUT_ILLEGAL_); // CO20200624
    }
    xstructure b = a;
    int iatom = 0;
    for (size_t itype = 0; itype < a.num_each_type.size(); itype++) {
      const string species = string(argv.at(2 + b.atoms.at(iatom).type));
      b.species.at(itype) = species;
      for (int j = 0; j < a.num_each_type[itype]; j++) {
        ostringstream aus;
        aus << species << j + 1; // CONVASP_MODE
        b.atoms.at(iatom).name = aus.str();
        b.atoms.at(iatom).CleanName();
        b.atoms.at(iatom).CleanSpin();
        b.atoms.at(iatom).name_is_given = true;
        iatom++;
      }
    }
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::PDB
// ***************************************************************************
namespace pflow {
  void PDB(istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    PrintPDB(a, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::PDOS
// ***************************************************************************
namespace pflow {
  void PDOS(vector<string> argv) {
    cerr << "# WARNING: THIS REQUIRES AN ALTERED VERSION OF VASP - SEE aflow --help" << endl;
    const int only_occ = 0; // Use unoccupied states in PDOS.
    pflow::projdata prd;
    pflow::pdosdata pdd;
    prd.PROOUTinfile = argv.at(3);
    pdd.PDOSinfile = argv.at(2);
    pflow::ReadInProj(prd);
    pflow::CalcNeatProj(prd, only_occ);
    pflow::ReadInPDOSData(prd, pdd);
    pflow::CalcPDOS(prd, pdd);
    pdd.PrintPDOS(cout, prd.sp);
    if (pdd.print_params) {
      pdd.PrintParams(cout, prd.LLMnames);
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::PEARSON_SYMBOL
// ***************************************************************************
// DX20210610 - updated function: added overload and added options
namespace pflow {
  string PEARSON_SYMBOL(istream& input) { // DX20210611 - overloaded
    const aurostd::xoption vpflow;
    return PEARSON_SYMBOL(input, vpflow);
  }
} // namespace pflow
namespace pflow {
  string PEARSON_SYMBOL(istream& input, const aurostd::xoption& vpflow) { // DX20210611 - added xoption
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }

    // ---------------------------------------------------------------------------
    // usage
    if (vpflow.flag("PEARSON_SYMBOL::USAGE")) {
      const string usage = "aflow --pearson_symbol|--pearson|--Pearson_symbol|--Pearson[=<tolerance_value>|=tight|=loose]";
      const string options = "options: [--no_scan] [--tolerance_spectrum=<start>,<stop>,<nsteps>]";
      vector<string> voptions;
      aurostd::string2tokens(options, voptions, " ");
      voptions.insert(voptions.begin(), usage);
      init::MessageOption("--usage", "pflow::PEARSON_SYMBOL()", voptions);
      return "";
    }

    stringstream sss;

    // ---------------------------------------------------------------------------
    // load structure
    xstructure a(input, IOAFLOW_AUTO);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " X1" << endl;
    }

    // ---------------------------------------------------------------------------
    // get tolerance
    const double tolerance = pflow::getSymmetryTolerance(a, vpflow.getattachedscheme("PEARSON_SYMBOL::TOLERANCE"));

    // ---------------------------------------------------------------------------
    // self-consistent tolerance scan
    if (vpflow.flag("PEARSON_SYMBOL::NO_SCAN")) {
      a.sym_eps_no_scan = true;
    } // DX20210406

    // ---------------------------------------------------------------------------
    // tolerance spectrum
    bool tolerance_spectrum_analysis = false;
    vector<double> tolerance_spectrum;
    if (vpflow.flag("PEARSON_SYMBOL::TOLERANCE_SPECTRUM")) {
      tolerance_spectrum_analysis = true;
      tolerance_spectrum = pflow::getSymmetryToleranceSpectrum(vpflow.getattachedscheme("PEARSON_SYMBOL::TOLERANCE_SPECTRUM"));
    } else if (vpflow.flag("PEARSON_SYMBOL::TOLERANCE") && vpflow.flag("PEARSON_SYMBOL::TOLERANCE_SPECTRUM")) {
      sss << __AFLOW_FUNC__ << " ERROR: Cannot specify a single tolerance value and perform the tolerance spectrum at the same time. Please choose one or the other.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, sss, _INPUT_ILLEGAL_);
    }

    // ---------------------------------------------------------------------------
    // get magnetic moment information
    if (vpflow.flag("PEARSON_SYMBOL::MAGNETIC")) {
      const string magmom_info = vpflow.getattachedscheme("PEARSON_SYMBOL::MAGNETIC");
      ProcessAndAddSpinToXstructure(a, magmom_info);
    }

    if (!tolerance_spectrum_analysis) {
      a.GetRealLatticeType(tolerance); // DX20210611 - more concise method for calling
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " X2" << endl;
        cerr << " Real space lattice primitive           = " << a.bravais_lattice_type << endl;
        cerr << " Real space lattice variation           = " << a.bravais_lattice_variation_type << endl;
        cerr << " Real space Pearson symbol              = " << a.pearson_symbol << endl;
      }
      sss << a.pearson_symbol << endl;
    } else {
      // determine Pearson symbol through a range of tolerances
      const xstructure a_orig = a; // store clean xstructure (no symmetry), but keeps options from command line
      const uint ntolerances = tolerance_spectrum.size();
      for (uint i = 0; i < ntolerances; i++) {
        a = a_orig;
        a.GetRealLatticeType(tolerance_spectrum[i]);
        sss << "tol=" << tolerance_spectrum[i] << ": " << a.pearson_symbol << endl;
      }
    }
    return sss.str();
  }
} // namespace pflow

// ***************************************************************************
// pflow::PLANEDENS
// ***************************************************************************
namespace pflow {
  void PLANEDENS(vector<string> argv) {
    // Read in charge
    xstructure str;
    ifstream chgfile(argv.at(3).c_str());
    vector<int> ngrid(3);
    vector<double> chg_tot;
    vector<double> chg_diff;
    pflow::ReadChg(str, ngrid, chg_tot, chg_diff, chgfile);
    // Read in planar density parameters.
    pflow::pd_params pdp;
    ifstream pdfile(argv.at(2).c_str());
    pflow::ReadPlaneDensParams(str, pdp, pdfile);
    // Get planar charge density
    vector<double> dens2d_tot;
    vector<double> dens2d_diff;
    pflow::GetPlaneDens(pdp, dens2d_tot, dens2d_diff, str, ngrid, chg_tot, chg_diff);
    pflow::PrintPlaneDens(pdp, dens2d_tot, dens2d_diff, str);
  }
} // namespace pflow

// ***************************************************************************
// pflow::PLATON
// ***************************************************************************
namespace pflow {
  string PLATON(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() == 1) {
      if (tokens[0] == "usage" || tokens[0] == "USAGE") {
        init::MessageOption(options, __AFLOW_FUNC__,
                            aurostd::liststring2string("aflow --platonSG[_label,_number][=EQUAL| EXACT][,ang,d1,d2,d3] < POSCAR  default:" + string("EQUAL=") + aurostd::utype2string<int>(DEFAULT_PLATON_P_EQUAL) +
                                                       "," + string("EXACT=") + aurostd::utype2string<int>(DEFAULT_PLATON_P_EXACT) + "," + aurostd::utype2string(DEFAULT_PLATON_P_ANG, 5) + "," +
                                                       aurostd::utype2string(DEFAULT_PLATON_P_D1, 5) + "," + aurostd::utype2string(DEFAULT_PLATON_P_D2, 5) + "," + aurostd::utype2string(DEFAULT_PLATON_P_D3, 5)));
        return ""; // CO20200624 - the option was expressed successfully
      }
    }
    // move on

    if (LDEBUG) {
      cerr << XPID << "pflow::PLATON: tokens.size()=" << tokens.size() << endl;
    }

    xstructure a(input, IOAFLOW_AUTO);
    bool Platon_EQUAL = DEFAULT_PLATON_P_EQUAL;
    bool Platon_EXACT = DEFAULT_PLATON_P_EXACT;
    double Platon_ang = DEFAULT_PLATON_P_ANG;
    double Platon_d1 = DEFAULT_PLATON_P_D1;
    double Platon_d2 = DEFAULT_PLATON_P_D2;
    double Platon_d3 = DEFAULT_PLATON_P_D3;
    //     if(argv.size()==2) Platon_ang=1.0e-2;
    // if(argv.size()==3) Platon_ang=atof(tokens.at(0));
    // Read in input file.
    if ((!tokens.empty() && tokens[0] == "EQUAL") || (tokens.size() >= 2 && tokens[1] == "EQUAL")) {
      Platon_EQUAL = true;
    }
    if ((!tokens.empty() && tokens[0] == "EXACT") || (tokens.size() >= 2 && tokens[1] == "EXACT")) {
      Platon_EXACT = true;
    }
    if (tokens.size() >= 3) {
      Platon_ang = aurostd::string2utype<double>(tokens.at(tokens.size() - 4));
      Platon_d1 = aurostd::string2utype<double>(tokens.at(tokens.size() - 3));
      Platon_d2 = aurostd::string2utype<double>(tokens.at(tokens.size() - 2));
      Platon_d3 = aurostd::string2utype<double>(tokens.at(tokens.size() - 1));
    }
    a.DecorateWithElements(); // DX20200727 - FakeNames() -> DecorateWithElements();
    return a.platon2print(Platon_EQUAL, Platon_EXACT, Platon_ang, Platon_d1, Platon_d2, Platon_d3);
    // aflow --platon [EQUAL] [EXACT] [ang d1 d2 d3]
    // CALC ADDSYM (EQUAL) (EXACT) (ang d1 d2 d3)
    // EQUAL - Search with all atom type treated as equivalent.
    // EXACT - All atoms should fit for given criteria.
    // ang - Angle criterium in search for metrical symmetry of the
    //       lattice (default 1.0 degree).
    // d1 - Distance criterium for coinciding atoms for non-inversion
    //      (pseudo)symmetry elements (default 0.25 Angstrom).
    // d2 - Distance criterium for coinciding atoms for (pseudo)
    //      inversion symmetry (default 0.45, 0.25 Angstrom).
    // d3 - Distance criterium for coinciding atoms for (pseudo)
    //      translation symmetry (default 0.45, 0.25 Angstrom).
    // SEE: http://www.cryst.chem.uu.nl/platon/pl000401.html
  }
} // namespace pflow

// ***************************************************************************
// pflow::POCC
// ***************************************************************************
namespace pflow {
  void POCC(vector<string> argv) {
    cerr << "# WARNING: THIS REQUIRES AN ALTERED VERSION OF VASP - SEE aflow --help" << endl;
    const int only_occ = 1; // Use occupied states only in projections.
    pflow::projdata proj_dat;
    proj_dat.PROOUTinfile = argv.at(2);
    pflow::ReadInProj(proj_dat);
    pflow::CalcNeatProj(proj_dat, only_occ);
    PrintNeatProj(proj_dat, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::RENDER //HE20250324
// ***************************************************************************
namespace pflow {
  /// @brief allow to render any structure from the command line
  /// @param input structure content
  /// @param output output folder path
  void RENDER(istream& input, const fs::path& output) {
    const xstructure xstr(input, IOAFLOW_AUTO);
    const fs::path save_path = output / "structure.html";
    aurostd::x3DWriter w = xstr.render();
    const std::string html = w.toHTML();
    aurostd::string2file(html, save_path);
    for (int axis_idx = 0; axis_idx < 3; axis_idx++) {
      w.tachyon_lattice_views_idx = axis_idx;
      aurostd::string2file(w.toTachyon(), output / ("structure_view_" + std::to_string(axis_idx + 1) + ".dat"));
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::POSCAR
// ***************************************************************************
namespace pflow {
  xstructure POSCAR(istream& input) {
    xstructure str(input, IOAFLOW_AUTO);
    str.iomode = IOVASP_POSCAR;
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::POSCAR2WYCKOFF
// ***************************************************************************
namespace pflow {
  void POSCAR2WYCKOFF(istream& input) {
    // Call findsym to output the wyckoff position from POSCAR
    // Note that findsym must be installed properly
    const string findsym_in = aurostd::TmpFileCreate("findsym.in");

    stringstream oss;
    xstructure str(input, IOAFLOW_AUTO);
    oss << str.findsym2print();
    aurostd::stringstream2file(oss, findsym_in);
    vector<string> stroutput;
    FINDSYM::Write("data_space.txt", "./");
    FINDSYM::Write("data_wyckoff.txt", "./");
    aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("findsym") + " < " + findsym_in), stroutput);
    FROZSL::Delete("data_space.txt", "./");
    FROZSL::Delete("data_wyckoff.txt", "./");

    bool flag = false;
    for (size_t i = 0; i < stroutput.size(); i++) {
      if (aurostd::substring2bool(stroutput[i], "----")) {
        flag = !flag;
      }
      if (flag) {
        if (!aurostd::substring2bool(stroutput[i], "----")) {
          cout << stroutput[i] << endl;
        }
      }
    }
    aurostd::RemoveFile("findsym.log " + findsym_in);
  }
} // namespace pflow

// ***************************************************************************
// pflow::PGROUP
// ***************************************************************************
namespace pflow {
  void PGROUP(_aflags& aflags, istream& input) {
    cout << aflow::Banner("BANNER_TINY") << endl;
    aflags.QUIET = true;
    xstructure a(input, IOAFLOW_AUTO);
    const bool WRITE = true;
    ofstream File("/dev/null");
    _kflags kflags; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_FGROUP = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL = false; // DX20170815 - Add in consistency checks //DX20171205 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_CALCULATE_SGROUP = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_IATOMS = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_AGROUP = false; // DX20170815 - Add in consistency checks
    pflow::PerformFullSymmetry(a, File, aflags, kflags, WRITE, cout); // DX20170815 - Add in consistency checks
    // DX20170815 - Add in consistency checks - SYM::CalculatePointGroup(File,a,aflags,WRITE,true,cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::PGROUPK
// ***************************************************************************
namespace pflow {
  void PGROUPK(_aflags& aflags, istream& input) {
    cout << aflow::Banner("BANNER_TINY") << endl;
    aflags.QUIET = true;
    xstructure a(input, IOAFLOW_AUTO);
    const bool WRITE = true;
    ofstream File("/dev/null");
    // DX20170815 SYM::CalculatePointGroupKLattice(File,a,aflags,WRITE,true,cout);
    _kflags kflags; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_FGROUP = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL = false; // DX20170815 - Add in consistency checks //DX20171205 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_CALCULATE_SGROUP = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_IATOMS = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_AGROUP = false; // DX20170815 - Add in consistency checks
    pflow::PerformFullSymmetry(a, File, aflags, kflags, WRITE, cout); // DX20170815 - Add in consistency checks
  }
} // namespace pflow

// ***************************************************************************
// pflow::PGROUPXTAL
// ***************************************************************************
namespace pflow {
  void PGROUPXTAL(_aflags& aflags, istream& input) {
    //  cout << aflow::Banner("BANNER_TINY") << endl;
    aflags.QUIET = true;
    xstructure a(input, IOAFLOW_AUTO);
    const bool WRITE = true;
    ofstream File("/dev/null");
    // DX20170815 - Add in consistency checks bool verbose=true;
    // DX20170815 - Add in consistency checks SYM::CalculatePointGroup(File,a,aflags,WRITE,verbose,cout);
    // DX SYM::CalculateFactorGroup(File,a,aflags,WRITE,verbose,cout);
    //   SYM::CalculateSpaceGroup(File,a,aflags,false,verbose,cout);
    //  SYM::CalculateInequivalentAtoms(File,a,aflags,WRITE,verbose,cout);
    //  SYM::CalculateSitePointGroup(File,a,aflags,WRITE,true,cout);
    // DX20170815 - Add in consistency checks SYM::CalculatePointGroupCrystal(File,a,aflags,WRITE,verbose,cout);
    _kflags kflags; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_FGROUP = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL = true; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL = false; // DX20170815 - Add in consistency checks //DX20171205 - Added pgroupk_xtal
    kflags.KBIN_SYMMETRY_CALCULATE_SGROUP = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_IATOMS = false; // DX20170815 - Add in consistency checks
    kflags.KBIN_SYMMETRY_CALCULATE_AGROUP = false; // DX20170815 - Add in consistency checks
    pflow::PerformFullSymmetry(a, File, aflags, kflags, WRITE, cout); // DX20170815 - Add in consistency checks
  }
} // namespace pflow

// ***************************************************************************
// pflow::PRIM
// ***************************************************************************
namespace pflow {
  xstructure PRIM(istream& input, uint mode) {
    const xstructure a(input, IOAFLOW_AUTO);
    if (mode == 0) {
      return GetPrimitive(a); // DX20210406
    }
    if (mode == 1) {
      return GetPrimitive1(a);
    }
    if (mode == 2) {
      return GetPrimitive2(a);
    }
    if (mode == 3) {
      return GetPrimitive3(a);
    }
    return a;
  }
} // namespace pflow

// ***************************************************************************
bool RequestedAlphabeticLabeling(string& label) {
  if (aurostd::substring2bool(label, ".alphabetic")) {
    aurostd::StringSubstInPlace(label, ".alphabetic", "");
    return true;
  }
  if (aurostd::substring2bool(label, ".alpha")) {
    aurostd::StringSubstInPlace(label, ".alpha", "");
    return true;
  }
  return false;
}

bool AlphabetizePrototypeLabelSpecies(deque<string>& species, deque<string>& species_pp, deque<double>& vvolume, deque<double>& vmass, string& label) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  // DEBUG=true;
  const uint nspeciesHTQC = species.size();
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species.size()=" << species.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " species_pp.size()=" << species_pp.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vvolume.size()=" << vvolume.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " vmass.size()=" << vmass.size() << endl;
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " nspeciesHTQC=" << nspeciesHTQC << endl;
  }
  // if(LDEBUG) cerr << __AFLOW_FUNC__ << " alphabetic=" << alphabetic << endl;
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " label=" << label << endl;
  }
  aurostd::StringSubstInPlace(label, ".alphabetic", "");
  aurostd::StringSubstInPlace(label, ".alpha", "");
  deque<string> rnd_species;
  deque<string> rnd_species_pp;
  deque<string> rnd_label1;
  deque<string> rnd_label2;
  deque<double> rnd_vvolume;
  deque<double> rnd_vmass;
  if (nspeciesHTQC != species.size()) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "nspeciesHTQC!=species.size", _INPUT_ILLEGAL_);
  } // CO20200624
  if (nspeciesHTQC != species_pp.size()) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "nspeciesHTQC!=species_pp.size", _INPUT_ILLEGAL_);
  } // CO20200624
  if (nspeciesHTQC != vvolume.size()) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "nspeciesHTQC!=vvolume.size", _INPUT_ILLEGAL_);
  } // CO20200624
  if (nspeciesHTQC != vmass.size()) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "nspeciesHTQC!=vmass.size", _INPUT_ILLEGAL_);
  } // CO20200624
  for (uint i = 0; i < nspeciesHTQC; i++) {
    rnd_species.push_back(species[i]); // LOAD FROM SPECIES
    rnd_species_pp.push_back(species_pp[i]); // LOAD FROM SPECIES
    rnd_vvolume.push_back(vvolume[i]); // LOAD FROM SPECIES
    rnd_vmass.push_back(vmass[i]); // LOAD FROM SPECIES
    string A = "A";
    A[0] += i;
    rnd_label1.push_back(A);
    rnd_label2.push_back(A); // megatrick thanks C++ to have C inside (SC2011)
  }
  aurostd::sort(rnd_species, rnd_species_pp, rnd_vvolume, rnd_vmass, rnd_label1); // how to fix rnd_atomxX by shufflying rnd_label1
  aurostd::sort(rnd_label1, rnd_label2); // how to go back from rnd_labe1 to rnd_label2
  label = label + ".";
  for (size_t i = 0; i < rnd_label2.size(); i++) {
    label = label + rnd_label2[i]; // FIXING LABEL
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " label=" << label << endl;
  }
  for (uint i = 0; i < nspeciesHTQC; i++) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEFORE species.at(" << i << ")=" << species.at(i) << " species_pp.at(" << i << ")=" << species_pp.at(i) << " vvolume.at(" << i << ")=" << vvolume.at(i) << " vmass.at(" << i
           << ")=" << vmass.at(i) << endl;
    }
  }
  for (uint i = 0; i < nspeciesHTQC; i++) {
    species[i] = rnd_species.at(i);
  }
  for (uint i = 0; i < nspeciesHTQC; i++) {
    species_pp[i] = rnd_species_pp.at(i);
  }
  for (uint i = 0; i < nspeciesHTQC; i++) {
    vvolume[i] = rnd_vvolume.at(i);
  }
  for (uint i = 0; i < nspeciesHTQC; i++) {
    vmass[i] = rnd_vmass.at(i);
  }
  for (uint i = 0; i < nspeciesHTQC; i++) {
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " AFTER species.at(" << i << ")=" << species.at(i) << " species_pp.at(" << i << ")=" << species_pp.at(i) << " vvolume.at(" << i << ")=" << vvolume.at(i) << " vmass.at(" << i
           << ")=" << vmass.at(i) << endl;
    }
  }
  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " END" << endl;
  }
  return true;
}

bool AlphabetizePrototypeLabelSpecies(deque<string>& species, deque<string>& species_pp, string& label) {
  deque<double> vvolume;
  for (size_t i = 0; i < species.size(); i++) {
    vvolume.push_back(double(0));
  }
  deque<double> vmass;
  for (size_t i = 0; i < species.size(); i++) {
    vmass.push_back(double(0));
  }
  return AlphabetizePrototypeLabelSpecies(species, species_pp, vvolume, vmass, label);
}

bool AlphabetizePrototypeLabelSpecies(deque<string>& species, deque<double>& vvolume, string& label) {
  deque<string> species_pp;
  for (size_t i = 0; i < species.size(); i++) {
    species_pp.push_back(species[i]);
  }
  deque<double> vmass;
  for (size_t i = 0; i < species.size(); i++) {
    vmass.push_back(double(0));
  }
  return AlphabetizePrototypeLabelSpecies(species, species_pp, vvolume, vmass, label);
}

bool AlphabetizePrototypeLabelSpecies(deque<string>& species, string& label) {
  deque<string> species_pp;
  for (size_t i = 0; i < species.size(); i++) {
    species_pp.push_back(species[i]);
  }
  deque<double> vvolume;
  for (size_t i = 0; i < species.size(); i++) {
    vvolume.push_back(double(0));
  }
  deque<double> vmass;
  for (size_t i = 0; i < species.size(); i++) {
    vmass.push_back(double(0));
  }
  return AlphabetizePrototypeLabelSpecies(species, species_pp, vvolume, vmass, label);
}

string AlphabetizePrototypeLabelSpeciesArgv(vector<string>& argv) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  //  LDEBUG=true;
  string label = argv.at(2);
  const uint nspeciesHTQC = aflowlib::PrototypeLibrariesSpeciesNumber(label);
  if (LDEBUG) {
    cerr << "AlphabetizePrototypeLabelSpeciesArgv: nspeciesHTQC=" << nspeciesHTQC << endl;
  }
  if (LDEBUG) {
    cerr << "AlphabetizePrototypeLabelSpeciesArgv: label=" << label << endl;
  }
  deque<string> species;
  for (uint i = 0; i < nspeciesHTQC; i++) {
    species.push_back(argv.at(3 + i)); // LOAD FROM ARGV
  }
  AlphabetizePrototypeLabelSpecies(species, label);
  for (uint i = 0; i < nspeciesHTQC; i++) {
    argv.at(3 + i) = species.at(i);
  }
  argv.at(2) = label;
  return label;
}

string AlphabetizePrototypeLabelSpeciesTokens(vector<string>& tokens) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  //  LDEBUG=true;
  string label = tokens.at(0);
  const uint nspeciesHTQC = aflowlib::PrototypeLibrariesSpeciesNumber(label);
  if (LDEBUG) {
    cerr << "AlphabetizePrototypeLabelSpeciesTokens: nspeciesHTQC=" << nspeciesHTQC << endl;
  }
  if (LDEBUG) {
    cerr << "AlphabetizePrototypeLabelSpeciesTokens: label=" << label << endl;
  }
  deque<string> species;
  for (uint i = 0; i < nspeciesHTQC; i++) {
    species.push_back(tokens.at(1 + i)); // LOAD FROM TOKENS
  }
  AlphabetizePrototypeLabelSpecies(species, label);
  for (uint i = 0; i < nspeciesHTQC; i++) {
    tokens.at(1 + i) = species.at(i);
  }
  tokens.at(0) = label;
  return label;
}

// ***************************************************************************
// pflow::PROTO_PARSE_INPUT
// ***************************************************************************
namespace pflow {
  bool PROTO_PARSE_INPUT(const vector<string>& params, vector<vector<string>>& vvstr, vector<vector<double>>& vvnum, bool ignore_label, bool reverse) { // CO20181226
    stringstream message;
    for (size_t i = 0; i < vvstr.size(); i++) {
      vvstr[i].clear();
    }
    vvstr.clear();
    for (size_t i = 0; i < vvnum.size(); i++) {
      vvnum[i].clear();
    }
    vvnum.clear();
    if (params.empty()) {
      return false;
    } // nothing provided...
    if (ignore_label && params.size() == 1) {
      return true;
    } // could just be label, that's ok
    vector<string> tokens;
    bool isfloat1;
    bool isfloat2;
    const bool ignore_first = (ignore_label && !reverse);
    const bool ignore_last = (ignore_label && reverse);
    for (size_t i = (ignore_first ? 1 : 0); i < params.size() - (ignore_last ? 1 : 0); i++) {
      aurostd::string2tokens(params[i], tokens, ",");
      if (tokens.empty()) {
        message << "No inputs found in params[" << i << "]";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      isfloat1 = aurostd::isfloat(tokens[0]);
      if (isfloat1) {
        vvnum.emplace_back(0);
      } else {
        vvstr.emplace_back(0);
      }
      for (size_t j = 0; j < tokens.size(); j++) {
        isfloat2 = aurostd::isfloat(tokens[j]);
        if (isfloat1 != isfloat2) {
          message << "Mixed string/number input at params[" << i << "]=" << params[i];
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
        }
        if (isfloat2) {
          vvnum.back().push_back(aurostd::string2utype<double>(tokens[j]));
        } else {
          vvstr.back().push_back(tokens[j]);
        }
      }
    }
    return true;
  }
} // namespace pflow

// ***************************************************************************
// pflow::PROTO_TEST_INPUT
// ***************************************************************************
namespace pflow {
  bool PROTO_TEST_INPUT(const vector<vector<string>>& vvstr, const vector<vector<double>>& vvnum, uint& nspeciesHTQC, bool patch_nspecies) { // CO20181226
    stringstream message;
    // patch for ICSD/pocc
    uint nspecies = nspeciesHTQC;
    if (patch_nspecies) {
      nspecies = vvstr.size();
    }
    // do some checks
    // check nspecies vs. vvstr
    const bool check_volumes = (!vvnum.empty() && !(vvnum.size() == 1 && vvnum[0].size() == 1)); // either volumes not specified, or global volume specified, i.e., volume=str.atoms.size()*volume_in;
    if (!vvstr.empty() && vvstr.size() != nspecies) {
      message << "Invalid input specification, mismatch between nspecies and species inputs (nspecies==" << nspecies << ", vvstr.size()==" << vvstr.size() << ")" << endl;
      message << "aflow --aflow_proto[=]label*:specieA*[:specieB*][:volumeA*[:volumeB*] | :volume]" << endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
    }
    if (check_volumes) {
      // check nspecies vs. vvnum
      if (vvnum.size() != nspecies) {
        message << "Invalid input specification, mismatch between nspecies and volume inputs (nspecies==" << nspecies << ", vvnum.size()==" << vvnum.size() << ")" << endl;
        message << "aflow --aflow_proto[=]label*:specieA*[:specieB*][:volumeA*[:volumeB*] | :volume]" << endl;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
      }
      // check vvstr vs. vvnum
      if (vvstr.size() != vvnum.size()) {
        message << "Invalid input specification, mismatch between species and volume inputs (vvstr.size()==" << vvstr.size() << ", vvnum.size()==" << vvnum.size() << ")" << endl;
        message << "aflow --aflow_proto[=]label*:specieA*[:specieB*][:volumeA*[:volumeB*] | :volume]" << endl;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
      }
      for (size_t i = 0; i < vvstr.size(); i++) {
        if (vvstr[i].size() != vvnum[i].size()) {
          message << "Invalid input specification, mismatch between combo species and volume inputs (vvstr[" << i << "].size()==" << vvstr[i].size() << ", vvnum[" << i << "].size()==" << vvnum[i].size() << ")" << endl;
          message << "aflow --aflow_proto[=]label*:specieA*[:specieB*][:volumeA*[:volumeB*] | :volume]" << endl;
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
        }
      }
    }
    // check no negative volumes
    for (size_t i = 0; i < vvnum.size(); i++) {
      for (size_t j = 0; j < vvnum[i].size(); j++) {
        if (vvnum[i][j] <= 0) { // also includes 0
          message << "Invalid input specification, negative/zero number input found (vvnum[" << i << "][" << j << "]=" << vvnum[i][j] << ")" << endl;
          message << "aflow --aflow_proto[=]label*:specieA*[:specieB*][:volumeA*[:volumeB*] | :volume]" << endl;
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
        }
      }
    }
    return true;
  }
} // namespace pflow

namespace pflow {
  bool sortPOCCSites(const string& p1, const string& p2) { // CO20181226
    stringstream message;
    vector<string> tokens;
    string designation1;
    string designation2;
    char mode1;
    char mode2;
    string _site1;
    string _site2;
    uint site1;
    uint site2;
    // p1
    aurostd::string2tokens(p1, tokens, "-");
    if (tokens.empty()) {
      message << "No occupants specified for [designation=" << p1 << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    designation1 = tokens[0];
    if (designation1.empty()) {
      message << "No mode specified [designation=" << p1 << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    mode1 = designation1[0];
    _site1 = designation1;
    _site1.erase(_site1.begin());
    if (!aurostd::isfloat(_site1)) {
      message << "Designation ill-written, no site provided [designation=" << p1 << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    site1 = aurostd::string2utype<int>(_site1);
    // p2
    aurostd::string2tokens(p2, tokens, "-");
    if (tokens.empty()) {
      message << "No occupants specified for [designation=" << p2 << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    designation2 = tokens[0];
    if (designation2.empty()) {
      message << "No mode specified [designation=" << p2 << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    mode2 = designation2[0];
    _site2 = designation2;
    _site2.erase(_site2.begin());
    if (!aurostd::isfloat(_site2)) {
      message << "Designation ill-written, no site provided [designation=" << p2 << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    site2 = aurostd::string2utype<int>(_site2);
    // check sites first
    if (site1 != site2) {
      return site1 < site2;
    } // sort 1 then 2 then 3
    // check mode
    if (mode1 != mode2) { // sort P first, as S is a grouped designation
      if (mode1 == 'P') {
        return true;
      } // P < S
      else {
        return false;
      } // S > P
    }
    message << "Double specification for the same site [designation1=" << p1 << ",designation2=" << p2 << "]";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    return true;
  }
  bool sortPOCCOccs(const string& occ1, const string& occ2) { // CO20181226
    const bool LDEBUG = (false || XHOST.DEBUG);
    stringstream message;
    vector<string> tokens;
    string _occupancy1;
    string _occupant1;
    string _occupancy2;
    string _occupant2;
    int occupant1;
    int occupant2;
    double occupancy1;
    double occupancy2;
    // occ1
    aurostd::string2tokens(occ1, tokens, "x");
    if (tokens.size() != 2) {
      message << "Occupancy x Occupant pair ill-defined for [occupation=" << occ1 << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    _occupancy1 = tokens[0];
    _occupant1 = tokens[1];
    occupant1 = (int) _occupant1[0] - 65;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " occupant to int conversion: " << _occupant1 << " == " << occupant1 << endl;
    } //[0] to convert to char
    if (!aurostd::isfloat(_occupancy1)) {
      message << "Occupancy is not a float for [occupation=" << occ1 << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    occupancy1 = aurostd::string2utype<double>(_occupancy1);
    // occ2
    aurostd::string2tokens(occ2, tokens, "x");
    if (tokens.size() != 2) {
      message << "Occupancy x Occupant pair ill-defined for [occupation=" << occ2 << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    _occupancy2 = tokens[0];
    _occupant2 = tokens[1];
    occupant2 = (int) _occupant2[0] - 65;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " occupant to int conversion: " << _occupant2 << " == " << occupant2 << endl;
    } //[0] to convert to char
    if (!aurostd::isfloat(_occupancy2)) {
      message << "Occupancy is not a float for [occupation=" << occ2 << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    occupancy2 = aurostd::string2utype<double>(_occupancy2);
    if (occupant1 != occupant2) {
      return occupant1 < occupant2;
    } // sort A then B then C
    if (occupancy1 != occupancy2) {
      return occupancy1 > occupancy2;
    } // sort most preferred first
    message << "Double specification for the same occupant [occupation1=" << occ1 << ",occupation2=" << occ2 << "]";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    return true;
  }
  bool FIX_PRECISION_POCC(const string& occ, string& new_occ) { // CO20181226
    const bool LDEBUG = (false || XHOST.DEBUG);
    new_occ = "";
    stringstream message;
    vector<string> tokens;
    string _occupancy;
    string _occupant;
    string occupant;
    double occupancy;
    aurostd::string2tokens(occ, tokens, "x");
    if (tokens.size() != 2) {
      message << "Occupancy x Occupant pair ill-defined for [occupation=" << occ << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    _occupancy = tokens[0];
    _occupant = tokens[1];
    occupant = (int) _occupant[0] - 65;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " occupant to int conversion: " << _occupant << " == " << occupant << endl;
    } //[0] to convert to char
    if (!aurostd::isfloat(_occupancy)) {
      message << "Occupancy is not a float for [occupation=" << occ << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    occupancy = aurostd::string2utype<double>(_occupancy);
    const double default_pocc_tol = DEFAULT_POCC_STOICH_TOL; // AFLOWRC_DEFAULT_POCC_STOICH_TOL;  //1e-3
    if (occupancy < default_pocc_tol) {
      message << "Cannot set desired precision (" << default_pocc_tol << ") for [occupation=" << occ << "]" << endl;
      message << "Doping is VERY small, requiring MANY representative ordered structures (not recommended)" << endl;
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, cout, _LOGGER_WARNING_);
      return false;
    }
    new_occ = occ;
    const int prec = (int) ceil(log10(1.0 / default_pocc_tol)); // AFLOWRC_DEFAULT_POCC_STOICH_TOL //use internal default here for consistency
    aurostd::StringSubstInPlace(new_occ, _occupancy, aurostd::utype2string(occupancy, prec));
    return true;
  }
  struct POCCSiteSpecification { // temporary container for sorting, do not put in .h (unless necessary)
    string input_string;
    char mode;
    vector<uint> positions;
  };
  vector<POCCSiteSpecification> poccString2POCCSiteSpecification(const xstructure& xstr, const vector<string> pocc_sites) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) { // CO20230319
      cerr << __AFLOW_FUNC__ << " xstr=" << endl << xstr << endl;
      cerr << __AFLOW_FUNC__ << " pocc_sites=";
      for (size_t i = 0; i < pocc_sites.size(); i++) {
        cerr << pocc_sites[i] << (i < pocc_sites.size() - 1 ? "," : "");
      }
      cerr << endl;
    }
    stringstream message;
    vector<POCCSiteSpecification> vpss;
    if (pocc_sites.empty()) {
      return vpss;
    }
    vector<string> tokens;
    vector<string> tokenstmp;
    string designation;
    char mode;
    string _site;
    uint site;
    POCCSiteSpecification pss;
    uint iatom;
    vector<string> pocc_designations; // quick way to check if there are duplicate specifications
    for (size_t i = 0; i < pocc_sites.size(); i++) {
      const string& ps = pocc_sites[i];
      if (ps.empty()) {
        message << "POCC_site[i=" << i << "] empty";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      aurostd::string2tokens(ps, tokens, "-");
      if (tokens.empty()) {
        message << "No occupants specified for [designation=" << ps << "]";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      designation = tokens[0];
      if (designation.empty()) {
        message << "No mode specified [designation=" << ps << "]";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      if (aurostd::WithinList(pocc_designations, designation)) {
        message << "Duplicate POCC designation: " << designation;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      pocc_designations.push_back(designation);
      mode = designation[0];
      _site = designation;
      _site.erase(_site.begin());
      if (!aurostd::isfloat(_site)) {
        message << "Designation ill-written, no site provided [designation=" << ps << "]";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      site = aurostd::string2utype<int>(_site);
      pss.input_string = ps;
      pss.mode = mode;
      pss.positions.clear();
      pss.positions.push_back(site);
      if (pss.mode == 'S') {
        pss.positions.clear();
        iatom = 0;
        for (size_t ii = 0; ii < xstr.num_each_type.size(); ii++) {
          for (uint jj = 0; jj < (uint) xstr.num_each_type[ii]; jj++) {
            if (ii == site) {
              pss.positions.push_back(iatom);
            }
            iatom++;
          }
        }
      }
      if (pss.positions.empty()) {
        message << "No positions found: " << ps;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      // replace P IFF positions.size()==1
      if (pss.mode != 'P' && pss.positions.size() == 1) {
        //[CO20190629 - need to replace designation + position]pss.input_string[0]='P';
        pss.mode = 'P';
        // CO20190629 START - small bug, need to change pss.input_string to map SPECIES 3 to POSITION X
        tokenstmp.clear();
        for (size_t it = 1; it < tokens.size(); it++) {
          tokenstmp.push_back(tokens[it]);
        } // skip 0, this is being changed
        pss.input_string = "P" + aurostd::utype2string(pss.positions[0]) + "-" + aurostd::joinWDelimiter(tokenstmp, "-");
        // CO20190629 STOP - small bug, need to change pss.input_string to map SPECIES 3 to POSITION X
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " creating new POCCSiteSpecification():" << endl;
        cerr << __AFLOW_FUNC__ << " POCCSiteSpecification.input_string=" << pss.input_string << endl;
        cerr << __AFLOW_FUNC__ << " POCCSiteSpecification.positions=" << aurostd::joinWDelimiter(pss.positions, ",") << endl;
      }
      std::sort(pss.positions.begin(), pss.positions.end());
      vpss.push_back(pss);
    }
    return vpss;
  }
  bool sortPOCCSiteSpecifications(const POCCSiteSpecification& p1, const POCCSiteSpecification& p2) {
    stringstream message;
    if (p1.positions.empty()) {
      message << "No positions found: " << p1.input_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    if (p2.positions.empty()) {
      message << "No positions found: " << p2.input_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    const uint position1 = p1.positions[0];
    const uint position2 = p2.positions[0];
    if (p1.mode != p2.mode) {
      if (p1.mode == 'S' && p2.mode == 'P') {
        if (aurostd::WithinList(p1.positions, position2)) {
          return false;
        } // put p2 first
      } else if (p1.mode == 'P' && p2.mode == 'S') {
        if (aurostd::WithinList(p2.positions, position1)) {
          return true;
        } // put p1 first
      } else {
        message << "Unknown modes: p1.mode=" << p1.mode << ", p2.mode=" << p2.mode;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
    }
    return position1 < position2;
  }
  void FIX_POCC_PARAMS(const xstructure& xstr, string& pocc_params) { // CO20181226
    // this sorts and fixes precision of pocc params
    if (pocc_params.empty()) {
      return;
    }
    aurostd::StringSubstInPlace(pocc_params, "*", "x"); // sometimes it's more intuitive to write 1*A instead of 1xA, let's fix automatically
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      ;
    } // dummy load
    stringstream message;
    vector<string> tokens0;
    vector<string> tokens1;
    vector<string> _tokens2;
    vector<string> tokens2;
    vector<string> tokens2_save;
    // START - do not modify pocc_params inside here
    aurostd::string2tokens(pocc_params, tokens0, ","); // does not modify pocc_params
    string occ;
    string new_occ;
    bool fix_precision_occupancy = true;
    vector<POCCSiteSpecification> vpss;
    for (size_t i = 0; i < tokens0.size(); i++) {
      aurostd::string2tokens(tokens0[i], tokens1, "_");
      vpss = poccString2POCCSiteSpecification(xstr, tokens1); // pack into tmp struct purely for sorting
      // std::sort(tokens1.begin(),tokens1.end(),sortPOCCSites);
      std::sort(vpss.begin(), vpss.end(), sortPOCCSiteSpecifications);
      tokens1.clear();
      for (size_t ip = 0; ip < vpss.size(); ip++) {
        tokens1.push_back(vpss[ip].input_string);
      } // unpack from tmp struct to tokens1
      if (LDEBUG) {
        for (size_t ip = 0; ip < tokens1.size(); ip++) {
          cerr << __AFLOW_FUNC__ << " tokens1[" << ip << "]=" << tokens1[ip] << " (NEW ORDER)" << endl;
        }
      }
      for (size_t j = 0; j < tokens1.size(); j++) {
        // second parse -
        aurostd::string2tokens(tokens1[j], _tokens2, "-");
        if (_tokens2.size() < 2) {
          message << "No occupants specified for [designation=" << tokens1[j] << "]";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
        }
        tokens2.clear();
        tokens2_save.clear();
        fix_precision_occupancy = true;
        for (size_t k = 1; k < _tokens2.size(); k++) { // remove mode part
          tokens2_save.push_back(_tokens2[k]);
          // check that all occupancies are not smaller than desired precision (avoid setting to 0 completely)
          // CO - check that we do not have issues with vacancies, we may also need to check 1-sum(x)
          occ = new_occ = tokens2_save.back();
          fix_precision_occupancy = (fix_precision_occupancy && FIX_PRECISION_POCC(occ, new_occ));
          tokens2.push_back(new_occ);
        }
        if (!fix_precision_occupancy) {
          tokens2.clear();
          for (size_t k = 0; k < tokens2_save.size(); k++) {
            tokens2.push_back(tokens2_save[k]);
          }
        } // avoid setting precision if it's too small
        std::sort(tokens2.begin(), tokens2.end(), sortPOCCOccs);
        tokens2.insert(tokens2.begin(), _tokens2[0]); // add mode part back
        tokens1[j] = aurostd::joinWDelimiter(tokens2, "-");
      }
      tokens0[i] = aurostd::joinWDelimiter(tokens1, "_");
    }
    // STOP - do not modify pocc_params inside here
    pocc_params = aurostd::joinWDelimiter(tokens0, ","); // replace input in place
  }
} // namespace pflow

// ***************************************************************************
// pflow::checkAnionSublattice
// ***************************************************************************
namespace pflow {
  bool checkAnionSublattice(const xstructure& xstr) { // CO20210201
    const bool LDEBUG = (false || XHOST.DEBUG);
    vector<string> vanions;
    aurostd::string2tokens(POCC_ANIONS_LIST, vanions, ",");
    vector<vector<string>> voccupants;
    vector<xvector<double>> vsites; // cpos
    vector<bool> vanions_found;
    uint i = 0;
    uint j = 0;
    bool found = false;
    for (i = 0; i < xstr.atoms.size(); i++) {
      found = false;
      for (j = 0; j < voccupants.size() && !found; j++) {
        if (aurostd::modulus(vsites[j] - xstr.atoms[i].cpos) < _AFLOW_POCC_ZERO_TOL_) {
          if (!aurostd::WithinList(voccupants[j], xstr.atoms[i].cleanname)) {
            voccupants[j].push_back(xstr.atoms[i].cleanname);
            if (aurostd::WithinList(vanions, xstr.atoms[i].cleanname)) {
              vanions_found[j] = true;
            }
            found = true;
          }
        }
      }
      if (!found) {
        vsites.push_back(xstr.atoms[i].cpos);
        voccupants.emplace_back(0);
        voccupants.back().push_back(xstr.atoms[i].cleanname);
        vanions_found.push_back(aurostd::WithinList(vanions, xstr.atoms[i].cleanname));
      }
    }
    if (vsites.size() != voccupants.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vsites.size()!=voccupants.size()", _RUNTIME_ERROR_);
    }
    if (vsites.size() != vanions_found.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vsites.size()!=vanions_found.size()", _RUNTIME_ERROR_);
    }
    if (LDEBUG) {
      for (j = 0; j < voccupants.size(); j++) {
        cerr << __AFLOW_FUNC__ << " site=" << j << " occupants=" << aurostd::joinWDelimiter(voccupants[j], ",") << " [ANIONS_SUBLATTICE=" << vanions_found[j] << "]" << endl;
      }
    }
    for (j = 0; j < voccupants.size(); j++) {
      if (vanions_found[j]) {
        for (i = 0; i < voccupants[j].size(); i++) {
          if (!aurostd::WithinList(vanions, voccupants[j][i])) {
            if (LDEBUG) {
              cerr << __AFLOW_FUNC__ << " found non-anion in anion_sublattice: " << voccupants[j][i] << endl;
            }
            return false;
          }
        }
      }
    }
    return true;
  }
} // namespace pflow

namespace pflow {
  // ***************************************************************************
  // pflow::convertXStr2POCC
  // ***************************************************************************
  // ./aflow --proto=T0009.ABC:Br:Cl:Cs_sv:I:Pb_d:Sm --pocc_params=S0-1xC_S1-0.5xE-0.5xF_S2-0.3333xA-0.3333xB-0.3333xD
  bool convertXStr2POCC(xstructure& xstr, const string& pocc_params, const vector<string>& _vspecies, const vector<double>& vvolumes) { // CO20181226
    if (pocc_params.empty()) {
      return false;
    }
    vector<string> vspecies;
    if (_vspecies.empty()) {
      for (uint i = 0; i < 26; i++) {
        vspecies.push_back(string() + char('A' + i));
      } // put back ABC...
    } else {
      for (size_t i = 0; i < _vspecies.size(); i++) {
        vspecies.push_back(_vspecies[i]);
      }
    }
    const bool LDEBUG = (false || XHOST.DEBUG);
    stringstream message;
    xstructure xstr_orig(xstr);
    // first parse _
    vector<string> tokens1;
    vector<string> tokens2;
    vector<string> tokens3;
    string designation;
    string _occupancy;
    string _occupant;
    string _site;
    double occupancy;
    uint site;
    uint occupant;
    uint iatom;
    char mode;
    aurostd::string2tokens(pocc_params, tokens1, "_");
    _atom atom;
    vector<uint> atoms2remove;
    vector<_atom> atoms2add;
    vector<uint> positions_already_added;
    vector<uint> positions_already_added_tmp;
    for (size_t i = 0; i < tokens1.size(); i++) {
      // second parse -
      aurostd::string2tokens(tokens1[i], tokens2, "-");
      if (tokens2.empty()) {
        message << "No occupants specified for [designation=" << tokens1[i] << "]";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      designation = tokens2[0];
      tokens2.erase(tokens2.begin());
      if (designation.empty()) {
        message << "Designation empty [i=" << i << "]";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      mode = designation[0];
      _site = designation;
      _site.erase(_site.begin());
      if (!aurostd::isfloat(_site)) {
        message << "Designation ill-written, no site provided [designation=" << designation << "]";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      site = aurostd::string2utype<uint>(_site); // this is a hash, so start from 0 //-1
      if (mode == 'P') {
        // site here refers to position (atom)
        if (site >= xstr_orig.atoms.size()) {
          message << "Invalid position index " << site << " >= xstr.atoms.size()=" << xstr_orig.atoms.size();
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
        }
        atoms2remove.push_back(site);
      } else if (mode == 'S') {
        if (site >= xstr_orig.num_each_type.size()) { // site here is type
          message << "Invalid type index " << site << " >= xstr.num_each_type.size()=" << xstr_orig.num_each_type.size();
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
        }
        iatom = 0;
        for (size_t ii = 0; ii < xstr_orig.num_each_type.size(); ii++) {
          for (uint jj = 0; jj < (uint) xstr_orig.num_each_type[ii]; jj++) {
            if (ii == site) {
              atoms2remove.push_back(iatom);
            }
            iatom++;
          }
        }
      } else {
        message << "Unknown POCC designation mode [designation=" << designation << "]";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      positions_already_added_tmp.clear(); // do not ignore at the 'x' level
      for (size_t j = 0; j < tokens2.size(); j++) {
        // third parse x
        aurostd::string2tokens(tokens2[j], tokens3, "x");
        if (tokens3.size() != 2) {
          message << "Occupancy x Occupant pair ill-defined for [designation=" << tokens1[i] << ",occupant=" << j << "]";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
        }
        _occupancy = tokens3[0];
        _occupant = tokens3[1];
        occupant = (int) _occupant[0] - 65;
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " occupant to int conversion: " << _occupant << " == " << occupant << endl;
        } //[0] to convert to char
        if (occupant >= vspecies.size()) {
          message << "Invalid occupant specification: occupant[char=" << _occupant << ",int=" << occupant << "] vs. vspecies.size()=" << vspecies.size();
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
        }
        if (!aurostd::isfloat(_occupancy)) {
          message << "Occupancy is not a float for [designation=" << tokens1[i] << ",occupant=" << j << "]";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
        }
        occupancy = aurostd::string2utype<double>(_occupancy);
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " mode=" << mode << endl;
          cerr << __AFLOW_FUNC__ << " site=" << site << endl;
          cerr << __AFLOW_FUNC__ << " positions_already_added=" << aurostd::joinWDelimiter(positions_already_added, ",") << endl;
        }
        if (mode == 'P') {
          if (!aurostd::WithinList(positions_already_added, site)) { // do NOT over add, each position should be specified exactly once
            atom = xstr_orig.atoms[site];
            positions_already_added_tmp.push_back(site); // do NOT over add, each position should be specified exactly once
            atom.type = occupant;
            atom.name = vspecies[occupant];
            atom.CleanName();
            atom.CleanSpin();
            atom.name_is_given = true;
            atom.partial_occupation_value = occupancy;
            atom.partial_occupation_flag = !aurostd::isequal(atom.partial_occupation_value, 1.0, _AFLOW_POCC_ZERO_TOL_);
            if (atom.partial_occupation_flag) {
              xstr.partial_occupation_flag = true;
            }
            atoms2add.push_back(atom);
          }
        }
        if (mode == 'S') {
          iatom = 0;
          for (size_t ii = 0; ii < xstr_orig.num_each_type.size(); ii++) {
            for (uint jj = 0; jj < (uint) xstr_orig.num_each_type[ii]; jj++) {
              if (ii == site) {
                if (!aurostd::WithinList(positions_already_added, iatom)) { // do NOT over add, each position should be specified exactly once
                  atom = xstr_orig.atoms[iatom];
                  positions_already_added_tmp.push_back(iatom); // do NOT over add, each position should be specified exactly once
                  atom.type = occupant;
                  atom.name = vspecies[occupant];
                  atom.CleanName();
                  atom.CleanSpin();
                  atom.name_is_given = true;
                  atom.partial_occupation_value = occupancy;
                  atom.partial_occupation_flag = !aurostd::isequal(atom.partial_occupation_value, 1.0, _AFLOW_POCC_ZERO_TOL_);
                  if (atom.partial_occupation_flag) {
                    xstr.partial_occupation_flag = true;
                  }
                  atoms2add.push_back(atom);
                }
              }
              iatom++;
            }
          }
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " positions_already_added_tmp=" << aurostd::joinWDelimiter(positions_already_added_tmp, ",") << endl;
        }
      }
      positions_already_added.insert(positions_already_added.end(), positions_already_added_tmp.begin(), positions_already_added_tmp.end()); // insert positions_already_added_tmp into positions_already_added
    }
    xstr.RemoveAtom(atoms2remove);
    std::stable_sort(atoms2add.begin(), atoms2add.end(), sortAtomsTypes); // safe because we do AddAtom() below
    for (size_t i = 0; i < atoms2add.size(); i++) {
      xstr.AddAtom(atoms2add[i], false);
    } // CO20230319 - add by type
    // volumes
    // AddAtom() takes care of setting default volumes/masses
    if (vvolumes.size() == xstr.species_volume.size()) {
      for (size_t i = 0; i < vvolumes.size(); i++) {
        xstr.species_volume[i] = vvolumes[i];
      }
    }
    double volume = 0;
    if (vvolumes.size() == 1) {
      for (size_t i = 0; i < xstr.atoms.size(); i++) {
        volume += vvolumes[0] * xstr.atoms[i].partial_occupation_value;
      }
    } else {
      for (size_t i = 0; i < xstr.atoms.size(); i++) {
        for (size_t j = 0; j < xstr.num_each_type.size(); j++) {
          if (xstr.atoms[i].name == xstr.species[j]) {
            volume += xstr.species_volume[j] * xstr.atoms[i].partial_occupation_value;
          }
        }
      }
    }
    xstr.SetVolume(volume); // CO20190205 - more robust
    xstr.neg_scale = true;
    // patch title
    vector<string> tokens;
    aurostd::string2tokens(xstr.title, tokens, " ");
    if (!tokens.empty()) {
      const string tobereplaced = tokens[0];
      aurostd::string2tokens(tobereplaced, tokens, "/");
      if (tokens.size() == 2) {
        string pp = tokens[0];
        string proto = tokens[1];
        pp = aurostd::joinWDelimiter(xstr.species, ""); // rename
        proto = proto + ":POCC_" + pocc_params;
        aurostd::StringSubstInPlace(xstr.title, tobereplaced, pp + "/" + proto);
      }
    }
    return true;
  }
} // namespace pflow

// ***************************************************************************
// pflow::POccInputs2Xstr
// ***************************************************************************
namespace pflow { // CO20211130
  bool POccInputs2Xstr(const string& pocc_input, aurostd::xoption& pocc_settings, xstructure& xstr, ostream& oss) {
    ofstream FileMESSAGE;
    return POccInputs2Xstr(pocc_input, pocc_settings, xstr, FileMESSAGE, oss);
  } // CO20211130
  bool POccInputs2Xstr(const string& pocc_input, aurostd::xoption& pocc_settings, xstructure& xstr, ofstream& FileMESSAGE, ostream& oss) { // CO20211130
    const bool LDEBUG = (false || XHOST.DEBUG);
    stringstream message;
    // example: Cs_svEuIPb_d:PAW_PBE.AB3C_cP5_221_a_c_b:POCC_S0-1xA_S1-1xC_S2-0.5xB-0.5xD
    // ARUN example: Cs_afEuIPb_d:PAW_PBE.AB3C_cP5_221_a_c_b:POCC_S0-1xA_S1-1xC_S2-0.5xB-0.5xD:ARUN.POCC_1_H0C0
    // convert to: --proto=AB3C_cP5_221_a_c_b:Cs_sv:Eu:I:Pb_d --pocc_params=S0-1xA_S1-1xC_S2-0.5xB-0.5xD
    // arun stuff separate

    if (!aurostd::substring2bool(pocc_input, TAG_TITLE_POCC)) { // use generic
      message << "No TAG_TITLE_POCC found [" << TAG_TITLE_POCC << "], using generic SYSTEM name as title";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_WARNING_); // CO20200404
      return false;
    }
    // Get all the pieces of the default title
    string::size_type loc1 = pocc_input.find(TAG_TITLE_POCC);
    const string elements_prototype_str = pocc_input.substr(0, loc1); // contains elements and prototype
    const string pocc_params_tol_arun_str = pocc_input.substr(loc1 + TAG_TITLE_POCC.length(), string::npos); // pocc_params and ARUN(?)
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " elements_prototype_str=" << elements_prototype_str << endl;
      cerr << __AFLOW_FUNC__ << " pocc_params_tol_arun_str=" << pocc_params_tol_arun_str << endl;
    }
    // parse elements_prototype_str by "."
    // PROBLEM: "." can exist in pp_string (not really important for standard PP, but it exists), as well
    // as proto: .ABC...
    // we will go in loop over "." parses until we get a structure!
    loc1 = elements_prototype_str.find('.');
    string pps;
    string proto;
    vector<string> velements;
    string tmp_str = "";
    string tmp_str2 = "";
    string pocc_params;
    string pocc_tol;
    string pocc_arun;
    string module_arun;
    string::size_type loc2 = 0;
    while (loc1 != string::npos && (loc1 + 1) < elements_prototype_str.length()) {
      pps = elements_prototype_str.substr(0, loc1);
      proto = elements_prototype_str.substr(loc1 + 1, string::npos);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " pps=" << pps << endl;
        cerr << __AFLOW_FUNC__ << " proto=" << proto << endl;
      }

      velements = aurostd::getElements(pps, pp_string, true, false, true); // clean, no sort_elements, pseudopotential string, keep_pp
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " velements=" << aurostd::joinWDelimiter(velements, ",") << endl;
      }

      tmp_str = pocc_params_tol_arun_str;
      pocc_params = tmp_str;
      loc2 = tmp_str.find(TAG_TITLE_POCC_TOL);
      if (loc2 != string::npos && (loc2 + 1) < tmp_str.length()) {
        pocc_params = tmp_str.substr(0, loc2);
        tmp_str = tmp_str.substr(loc2 + 1, string::npos);
        pocc_tol = tmp_str;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " [1]" << endl;
        cerr << __AFLOW_FUNC__ << " proto=" << proto << endl;
        cerr << __AFLOW_FUNC__ << " pocc_params=" << pocc_params << endl;
        cerr << __AFLOW_FUNC__ << " pocc_tol=" << pocc_tol << endl;
        cerr << __AFLOW_FUNC__ << " pocc_arun=" << pocc_arun << endl;
        cerr << __AFLOW_FUNC__ << " module_arun=" << module_arun << endl;
        cerr << __AFLOW_FUNC__ << " tmp_str=" << tmp_str << endl;
      }
      loc2 = tmp_str.find(TAG_TITLE_POCC_ARUN);
      if (loc2 != string::npos && (loc2 + 1) < tmp_str.length()) {
        tmp_str2 = tmp_str.substr(0, loc2);
        if (tmp_str2.find(TAG_TITLE_POCC_TOL) != string::npos) {
          pocc_tol = tmp_str2;
        } else {
          pocc_params = tmp_str2;
        }
        tmp_str = tmp_str.substr(loc2 + 1, string::npos);
        pocc_arun = tmp_str;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " [2]" << endl;
        cerr << __AFLOW_FUNC__ << " proto=" << proto << endl;
        cerr << __AFLOW_FUNC__ << " pocc_params=" << pocc_params << endl;
        cerr << __AFLOW_FUNC__ << " pocc_tol=" << pocc_tol << endl;
        cerr << __AFLOW_FUNC__ << " pocc_arun=" << pocc_arun << endl;
        cerr << __AFLOW_FUNC__ << " module_arun=" << module_arun << endl;
        cerr << __AFLOW_FUNC__ << " tmp_str=" << tmp_str << endl;
      }
      // CO20210315 - might find pocc_params='P0-1xA_P1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF:ARUN.AGL_9_SF_0.95' which is a bad system name, missing ARUN.POCC
      loc2 = tmp_str.find(TAG_TITLE_ARUN);
      if (loc2 != string::npos && (loc2 + 1) < tmp_str.length()) {
        tmp_str2 = tmp_str.substr(0, loc2);
        if (tmp_str2.find(TAG_TITLE_POCC_ARUN) != string::npos) {
          pocc_arun = tmp_str2;
        } else if (tmp_str2.find(TAG_TITLE_POCC_TOL) != string::npos) {
          pocc_tol = tmp_str2;
        } else {
          pocc_params = tmp_str2;
        }
        tmp_str = tmp_str.substr(loc2 + 1, string::npos);
        module_arun = tmp_str;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " [3]" << endl;
        cerr << __AFLOW_FUNC__ << " proto=" << proto << endl;
        cerr << __AFLOW_FUNC__ << " pocc_params=" << pocc_params << endl;
        cerr << __AFLOW_FUNC__ << " pocc_tol=" << pocc_tol << endl;
        cerr << __AFLOW_FUNC__ << " pocc_arun=" << pocc_arun << endl;
        cerr << __AFLOW_FUNC__ << " module_arun=" << module_arun << endl;
        cerr << __AFLOW_FUNC__ << " tmp_str=" << tmp_str << endl;
      }
      loc2 = tmp_str.find(":");
      if (loc2 != string::npos && (loc2 + 1) < tmp_str.length()) {
        tmp_str2 = tmp_str.substr(0, loc2);
        if (tmp_str2.find(TAG_TITLE_POCC_ARUN) != string::npos) {
          pocc_arun = tmp_str2;
        } // more specific, look for first
        else if (tmp_str2.find(TAG_TITLE_ARUN) != string::npos) {
          module_arun = tmp_str2;
        } else if (tmp_str2.find(TAG_TITLE_POCC_TOL) != string::npos) {
          pocc_tol = tmp_str2;
        } else {
          pocc_params = tmp_str.substr(0, loc2);
        }
        tmp_str = pocc_arun.substr(loc2 + 1, string::npos);
        // what's next??
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " [4]" << endl;
        cerr << __AFLOW_FUNC__ << " proto=" << proto << endl;
        cerr << __AFLOW_FUNC__ << " pocc_params=" << pocc_params << endl;
        cerr << __AFLOW_FUNC__ << " pocc_tol=" << pocc_tol << endl;
        cerr << __AFLOW_FUNC__ << " pocc_arun=" << pocc_arun << endl;
        cerr << __AFLOW_FUNC__ << " module_arun=" << module_arun << endl;
        cerr << __AFLOW_FUNC__ << " tmp_str=" << tmp_str << endl;
      }

      aurostd::xoption proto_flags;
      proto_flags.push_attached("PROTO", proto + ":" + aurostd::joinWDelimiter(velements, ":"));
      proto_flags.push_attached("POCC_PARAMS", pocc_params);
      if (!pocc_tol.empty()) {
        tmp_str = pocc_tol;
        if (tmp_str.find(TAG_TOL) != string::npos) {
          aurostd::StringSubstInPlace(tmp_str, TAG_TITLE_POCC_TOL, ""); // remove 'TOL_'
          aurostd::StringSubstInPlace(tmp_str, TAG_TOL + SEP_TAG2, ""); // remove 'TOL_'
          aurostd::StringSubstInPlace(tmp_str, TAG_TOL, ""); // remove 'TOL_'
          aurostd::StringSubstInPlace(tmp_str, SEP_TAG2, ""); // remove 'TOL_'
        }
        proto_flags.push_attached("POCC_TOL", tmp_str);
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " proto_flags.getattachedscheme(\"PROTO\")=" << proto_flags.getattachedscheme("PROTO") << endl;
        cerr << __AFLOW_FUNC__ << " proto_flags.getattachedscheme(\"POCC_PARAMS\")=" << proto_flags.getattachedscheme("POCC_PARAMS") << endl;
        cerr << __AFLOW_FUNC__ << " proto_flags.getattachedscheme(\"POCC_TOL\")=" << proto_flags.getattachedscheme("POCC_TOL") << endl;
      }

      try {
        xstr = pflow::PROTO_LIBRARIES(proto_flags);
        break;
      } catch (aurostd::xerror& excpt) {
        xstr.clear(); // DX20191220 - uppercase to lowercase clear
        loc1 = elements_prototype_str.find('.', loc1 + 1);
        continue;
      }
    }

    if (xstr.atoms.empty()) { // use generic
      message << "Cannot extract identifiable prototype from SYSTEM [" << pocc_input << "], using generic SYSTEM name as title";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_WARNING_); // CO20200404
      return false;
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xstr_found: " << endl;
      cerr << xstr << endl;
    }

    if (xstr.species.size() != xstr.comp_each_type.size()) { // use generic
      message << "Cannot extract composition from prototype [" << proto << "], using generic SYSTEM name as title";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_WARNING_); // CO20200404
      return false;
    }

    pocc_settings.clear();
    pocc_settings.push_attached("PPS", pps);
    pocc_settings.push_attached("PROTO", proto);
    pocc_settings.push_attached("POCC_PARAMS", pocc_params);
    pocc_settings.push_attached("POCC_TOL", pocc_tol);
    pocc_settings.push_attached("POCC_ARUN", pocc_arun);
    pocc_settings.push_attached("MODULE_ARUN", module_arun);

    return true;
  }
} // namespace pflow

// ***************************************************************************
// pflow::PROTO_LIBRARIES
// ***************************************************************************
namespace pflow {
  xstructure PROTO_LIBRARIES(aurostd::xoption vpflow) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    stringstream message;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }

    vector<string> tokens; // CO20181226
    vector<string> params;
    vector<string> vstr;
    vector<string> vstr_orig; // CO20181226
    vector<double> vnum;
    vector<double> vnum_orig; // CO20181226
    vector<vector<string>> vvstr; // CO20181226
    vector<vector<double>> vvnum; // CO20181226
    if (LDEBUG) {
      cerr << XPID << "pflow::PROTO: vpflow.getattachedscheme(\"PROTO\")=" << vpflow.getattachedscheme("PROTO") << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "pflow::PROTO: vpflow.getattachedscheme(\"PROTO_ICSD_AFLOW\")=" << vpflow.getattachedscheme("PROTO_ICSD_AFLOW") << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "pflow::PROTO: vpflow.flag(\"PROTO\")=" << vpflow.flag("PROTO") << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "pflow::PROTO: vpflow.flag(\"PARAMS\")=" << vpflow.flag("PARAMS") << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "pflow::PROTO: vpflow.flag(\"POCC_PARAMS\")=" << vpflow.flag("POCC_PARAMS") << endl; // CO20181226
    }
    if (LDEBUG) {
      cerr << XPID << "pflow::PROTO: vpflow.flag(\"POCC_TOL\")=" << vpflow.flag("POCC_TOL") << endl; // CO20181226
    }
    if (LDEBUG) {
      cerr << XPID << "pflow::PROTO: vpflow.flag(\"PROTO_ICSD_AFLOW\")=" << vpflow.flag("PROTO_ICSD_AFLOW") << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "pflow::PROTO: vpflow.flag(\"PROTO::USE_ANRL_LATTICE_PARAM\")=" << vpflow.flag("PROTO::USE_ANRL_LATTICE_PARAM") << endl; // DX20190227 - add anrl params flag
    }
    aurostd::string2tokens(vpflow.getattachedscheme("PROTO"), params, ":");

    if (params.empty()) {
      init::ErrorOption(vpflow.getattachedscheme("PROTO"), "pflow::PROTO",
                        aurostd::liststring2string(
                            "aflow [options] --proto=label*[:speciesA*[:speciesB*]..[:volumeA*[:volumeB*].. | :volume]] [--params=..... [--hex]]",
                            "                --proto[_icsd]=ICSD_number.{ABC}[:speciesA*[:speciesB*]..[:volumeA*[:volumeB*].. | :volume]]", "                --proto_icsd=label_ICSD_number",
                            "options = [--server=xxxxxxx]  [--vasp | --itc | --qe | --abinit | --aims | --cif | --abccar] [--params=... | --hex]]", "To get the list of prototypes --protos or --protos_icsd ."));
    }

    // CO20181226 parsing label, species, volume, and pocc input
    const string labels_raw = params.at(0);
    aurostd::string2tokens(labels_raw, tokens, ",");
    if (tokens.size() != 1) {
      message << "Too many labels provided, --proto can only handle one label at at time: labels.size()==" << tokens.size();
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    string label; // params.at(0);
    if (!tokens.empty()) {
      label = tokens[0];
    }

    PROTO_PARSE_INPUT(params, vvstr, vvnum, true, false); // CO20181226 - parse into vectors of strings (species) and numbers (volumes), skip first (label), no reverse (label)
    // ignore any multi-specifications
    vstr.clear();
    vstr_orig.clear();
    for (size_t i = 0; i < vvstr.size(); i++) {
      if (!vvstr[i].empty()) {
        vstr.push_back(vvstr[i][0]);
        vstr_orig.push_back(vvstr[i][0]);
      }
    } // SAVE ORIG
    vnum.clear();
    vnum_orig.clear();
    for (size_t i = 0; i < vvnum.size(); i++) {
      if (!vvnum[i].empty()) {
        vnum.push_back(vvnum[i][0]);
        vnum_orig.push_back(vvnum[i][0]);
      }
    } // SAVE ORIG

    string parameters = vpflow.getattachedscheme("PARAMS");

    // HE20250515
    //  Check if the value in --proto is in the database
    const std::string uid = anrl::getPrototypeUID(label);
    const anrl::ProtoData pd = anrl::ProtoData::get();

    // DX20190227 START
    //  default without parameters is to use the atomic volume scaling method (i.e., a=-1)
    //  --use_anrl_lattice_param forces the use of the original lattice parameter defined in ANRL (e.g., a=5.4)
    if (!parameters.empty() && vpflow.flag("PROTO::USE_ANRL_LATTICE_PARAM")) {
      message << "Cannot have both --params=<...> and --use_anrl_lattice_param at the same time; use one or the other.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    } else if (parameters.empty() && vpflow.flag("PROTO::USE_ANRL_LATTICE_PARAM")) {
      parameters = "use_anrl_lattice_param";
    }
    // DX20190227 END

    const string pocc_parameters_raw = vpflow.getattachedscheme("POCC_PARAMS"); // CO20181226
    const string pocc_tol = vpflow.getattachedscheme("POCC_TOL"); // CO20181226
    aurostd::string2tokens(pocc_parameters_raw, tokens, ",");
    if (!(tokens.empty() || tokens.size() == 1)) {
      message << "Too many sets of pocc_parameters provided, --proto can only handle one set of pocc_parameters at at time: pocc_parameters.size()==" << tokens.size();
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    string pocc_parameters;
    if (!tokens.empty()) {
      pocc_parameters = tokens[0];
    }
    const bool pocc = (!pocc_parameters.empty());

    //  cerr << params.at(0) << " " << params.size() << endl;  for(size_t i=0;i<params.size();i++) cerr << params.at(i) << " ";  cerr << endl;

    xstructure str;
    deque<string> atomX;
    deque<double> volumeX;
    uint nspecies = 0; // CO20191110 - NOTE: nspeciesHTQC is what is expected from label (proto), nspecies is REAL input number of species separated by colons, not commas (different for pocc)
    bool alphabetic = true; // DEFAULT

    // ***************************************************************************
    // MODE LIBRARY FROM HTQC OR ICSD
    if (uid.empty()) {
      nspecies = aflowlib::PrototypeLibrariesSpeciesNumber(label); // CO20181226 - recently revamped, now nspecies should be exact!
    } else {
      label = static_cast<string>(pd.content[uid]["label"]);
      nspecies = static_cast<uint>(pd.content[uid]["number_of_species"]);
      if (parameters.empty()) {
        parameters = static_cast<string>(pd.content[uid]["parameter_values"]);
      }
    }
    alphabetic = RequestedAlphabeticLabeling(label);

    bool found = false;
    double vol = 0.0;

    // CO20181226 START
    PROTO_TEST_INPUT(vvstr, vvnum, nspecies, pocc); // test if inputs are correct in number and type (not negative, etc.)  //CO20191110 - note nspecies changes here if pocc, so below we use nspeciesHTQC to fetch correct proto

    // CO20191110 - moved down from above so it prints updated values from pocc
    //  DEBUG=true;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " params.size()=" << params.size() << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vstr.size()=" << vstr.size() << endl; // CO20181226
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vnum.size()=" << vnum.size() << endl; // CO20181226
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " pocc=" << pocc << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " nspecies=" << nspecies << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " alphabetic=" << alphabetic << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " label=" << label << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " parameters=" << parameters << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " pocc_parameters=" << pocc_parameters << endl; // CO20181226
    }

    // get a vector of species and volume that fits exactly what is expected for xproto
    if (pocc) {
      vstr.clear();
      for (size_t i = 0; i < nspecies && i < vstr_orig.size(); i++) {
        vstr.push_back(vstr_orig[i]);
      }
      vnum.clear();
      for (size_t i = 0; i < nspecies && i < vnum_orig.size(); i++) {
        vnum.push_back(vnum_orig[i]);
      }
    }

    const int mode = 0;
    if (vstr.empty()) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " label" << endl;
      }
      str = aflowlib::PrototypeLibraries(cerr, label, parameters, mode);
      found = true; // good for binary, ternary etc etc...
    }

    atomX.clear();
    for (size_t i = 0; i < nspecies && i < vstr.size(); i++) {
      atomX.push_back(vstr[i]);
    }

    if (!found && vstr.size() == nspecies && vnum.empty()) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " label:species:species..." << endl;
      }
      if (LDEBUG) {
        for (size_t i = 0; i < atomX.size(); i++) {
          cerr << __AFLOW_FUNC__ << " atomX[" << i << "]=" << atomX[i] << endl;
        }
      }
      str = aflowlib::PrototypeLibraries(cerr, label, parameters, atomX, mode);
      found = true;
    }

    if (!found && vstr.size() == nspecies && (vnum.size() == 1 || vstr.size() == vnum.size())) {
      volumeX.clear();
      if (vnum.size() == 1) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " label:species:species...:volume" << endl;
        }
        vol = vnum[0];
        for (size_t i = 0; i < nspecies && i < atomX.size(); i++) {
          volumeX.push_back(0);
        }
      } else {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " label:species:species...:volume:volume..." << endl;
        }
        vol = -1.0;
        for (size_t i = 0; i < nspecies && i < vnum.size(); i++) {
          volumeX.push_back(vnum[i]);
        }
      }
      if (LDEBUG) {
        for (size_t i = 0; i < volumeX.size(); i++) {
          cerr << __AFLOW_FUNC__ << " volumeX[" << i << "]=" << volumeX[i] << endl;
        }
        cerr << __AFLOW_FUNC__ << " vol=" << vol << endl;
      }
      str = aflowlib::PrototypeLibraries(cerr, label, parameters, atomX, volumeX, vol, mode);
      found = true;
    }

    // CO20181226 STOP

    if (!found) {
      message << "Unknown label+input specification (label=" << label << ",nspecies=" << nspecies << ",ninput=" << params.size() - 1 << ")"; // params.size()-1 means skip label
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_); // CO20180801
    }

    // CO20181226 - try to POCC structure
    if (pocc) {
      pflow::FIX_POCC_PARAMS(str, pocc_parameters);
      convertXStr2POCC(str, pocc_parameters, vstr_orig, vnum_orig);
      setPOCCTOL(str, pocc_tol);
      if (!pflow::checkAnionSublattice(str)) { // CO20210201
        if (!XHOST.vflag_control.flag("FORCE_POCC")) {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Found non-anion in anion sublattice. Please check (and run with --force_pocc).", _VALUE_ILLEGAL_);
        }
      }
    }

    if (LDEBUG) {
      for (size_t i = 0; i < str.species.size(); i++) {
        if (!str.species[i].empty()) {
          cerr << "DEBUG specie(" << i << ")=" << str.species[i] << endl;
        }
      }
    }
    //  if(str.AlphabeticSpecie(0,1)==false) str.SwapSpecie(0,1);

    if (!pocc) { // CO20181226
      // now checking QUANTUM ESPRESSO
      if (vpflow.flag("PROTO::QE")) {
        str.xstructure2qe();
      }

      // now checking ABINIT
      if (vpflow.flag("PROTO::ABINIT")) {
        str.xstructure2abinit();
      }

      // now checking AIMS
      if (vpflow.flag("PROTO::AIMS")) {
        str.xstructure2aims();
      }

      // now checking ELK //DX20200313
      if (vpflow.flag("PROTO::ELK")) {
        str.xstructure2elk();
      }

      // now checking LMP //SD20240111
      if (vpflow.flag("PROTO::LMP")) {
        str.xstructure2lmp();
      }

      // now checking ITC
      if (vpflow.flag("PROTO::ITC")) {
        str.xstructure2itc(); // CO20220613
      }
    }

    // DX20190123 - add CIF/ABCCAR - START
    //  now checking CIF
    if (vpflow.flag("PROTO::CIF")) {
      str.xstructure2cif();
    }

    // now checking ABCCAR
    if (vpflow.flag("PROTO::ABCCAR")) {
      str.xstructure2abccar();
    }
    // DX20190123 - add CIF/ABCCAR - END

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::PROTO_AFLOW
// ***************************************************************************
// struct _AVASP_PROTO {
//   vector<string> ucell;
//   vector<int> vkppra;
//   vector<double> vpressure;
//   aurostd::xoptions vparams;
// };

namespace pflow {
  bool PROTO_AFLOW(aurostd::xoption vpflow, bool flag_REVERSE) { // too many options
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }

    // check prototypes
    _AVASP_PROTO PARAMS;
    PARAMS.ucell.clear();

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"PROTO_AFLOW\")=" << vpflow.getattachedscheme("PROTO_AFLOW") << endl;
    }
    aurostd::string2tokens(vpflow.getattachedscheme("PROTO_AFLOW"), PARAMS.ucell, ":");
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.ucell.size()=" << PARAMS.ucell.size() << endl;
    }

    // check usage
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK USAGE" << endl;
    }
    if (vpflow.flag("PROTO_AFLOW::USAGE") || PARAMS.ucell.empty()) {
      init::MessageOption(vpflow.getattachedscheme("PROTO_AFLOW"), "pflow::PROTO_AFLOW",
                          aurostd::liststring2string(
                              "aflow [options] --aflow_proto[|_icsd]=label*:speciesA*[:speciesB*][:volumeA*[:volumeB*]|:volume] [--params=... [--hex]]", "       options:", "                --usage",
                              "                --potential=pot_LDA | pot_GGA | potpaw_LDA | potpaw_GGA | potpaw_PBE | potpaw_LDA_KIN | potpaw_PBE_KIN ", "                --potential_complete",
                              "                --module=[APL|QHA|AAPL]", "                --apl_supercell=NxNxN", "                --usage", "                --missing", "                --noautopp",
                              "                --bader (default: DEFAULT_VASP_FORCE_OPTION_BADER)", "                --spin_remove_relax_1 (default: DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1)",
                              "                --spin_remove_relax_2 (default: DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2)",
                              "                --kscheme=[M | G] (default: DEFAULT_KSCHEME in .aflow.rc) --kscheme_static=[M | G] (default: DEFAULT_KSCHEME_STATIC in .aflow.rc)",
                              "                --kppra=NNNN (default: DEFAULT_KPPRA in .aflow.rc) --kppra_static=NNNN (default: DEFAULT_KPPRA_STATIC in .aflow.rc) --bands_grid=NNNN (default: "
                              "DEFAULT_BANDS_GRID in .aflow.rc)",
                              "                --enmax_multiply=NNNN (default:VASP_PREC_ENMAX_XXXX in .aflow.rc)", "                --pressure=0,1,2 (kB) (default:0.0)",
                              "                --potim=XXX (default 0.05)  (VASP) ", "                --relax_type=[ALL | IONS | CELL_SHAPE | CELL_VOLUME | IONS_CELL_VOLUME] ",
                              "                --relax_mode=[ENERGY | FORCES | ENERGY_FORCES | FORCES_ENERGY] (default: DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME in .aflow.rc) (VASP) ",
                              "                --relax_count=XX (default: DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT in .aflow.rc) (VASP) ", "                --run_relax_static",
                              "                --run_relax_static_bands", "                --run_static", "                --run_static_bands",
                              "                --precision=[(LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED] (default: DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME in .aflow.rc) (VASP) ",
                              "                --algorithm=[(NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED] (default: DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME in .aflow.rc) (VASP) ",
                              "                --metagga=[TPSS | RTPSS | M06L | MBJL | SCAN | MS0 | MS1 | MS2 | NONE] (default: DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME in .aflow.rc) (VASP) ",
                              "                --ivdw=[number_for_VASP_see_manual_for_IVDW | 0] (default: DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME in .aflow.rc) (VASP) ",
                              "                --type=[METAL | INSULATOR | SEMICONDUCTOR | DEFAULT] (default: DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME in .aflow.rc) (VASP) ",
                              "                --convert_unit_cell= [SPRIM, SCONV, NIGGLI, MINK, INCELL, COMPACT, WS, CART, FRAC, PRES] ", "                --volume_plus_equal=XXX ",
                              "                --volume_multiply_equal=XXX ", "                --volume_preserved ", "                --ediffg=XXX  (default: DEFAULT_VASP_PREC_EDIFFG in .aflow.rc) (VASP) ",
                              "                --ldau2", "                --noldau2", "                --neglect_nomix", "                --stdout", "                --qe", "                --abinit",
                              "                --aims", "                --list", "                --params=....  { check aflow --readme=anrl }", "                --hex          { check aflow --readme=anrl }"));
      return true; // CO20200624 - the option was expressed successfully
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"PROTO_AFLOW::USAGE\")=" << vpflow.flag("PROTO_AFLOW::USAGE") << endl;
    }

    // label
    string string_LABEL = PARAMS.ucell.at(0);

    //    if(LDEBUG) {    cerr << __AFLOW_FUNC__ << " " << string_LABEL << " " << PARAMS.ucell.size() << endl;     for(size_t i=0;i<PARAMS.ucell.size();i++) cerr << PARAMS.ucell.at(i) << " "; cerr << endl;}

    // check potential
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK POTENTIAL" << endl;
    }

    string string_POTENTIAL = _AVASP_PSEUDOPOTENTIAL_AUTO_;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"PROTO_AFLOW::POTENTIAL\")=" << vpflow.getattachedscheme("PROTO_AFLOW::POTENTIAL") << endl;
    }
    string params_potential = vpflow.getattachedscheme("PROTO_AFLOW::POTENTIAL");
    if (params_potential == _AVASP_PSEUDOPOTENTIAL_AUTO_) {
      string_POTENTIAL = _AVASP_PSEUDOPOTENTIAL_AUTO_;
    }
    aurostd::StringSubstInPlace(params_potential, "pot", "");
    aurostd::StringSubstInPlace(params_potential, "POT", "");
    aurostd::StringSubstInPlace(params_potential, "_", "");
    aurostd::StringSubstInPlace(params_potential, "paw", "PAW");
    if (params_potential == "LDA") {
      string_POTENTIAL = DEFAULT_VASP_POTCAR_DIR_POT_LDA;
    }
    if (params_potential == "GGA") {
      string_POTENTIAL = DEFAULT_VASP_POTCAR_DIR_POT_GGA;
    }
    if (params_potential == "PBE") {
      string_POTENTIAL = DEFAULT_VASP_POTCAR_DIR_POT_PBE;
    }
    if (params_potential == "PAWLDA") {
      string_POTENTIAL = DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA;
    }
    if (params_potential == "PAWGGA") {
      string_POTENTIAL = DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA;
    }
    if (params_potential == "PAWPBE") {
      string_POTENTIAL = DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE;
    }
    if (params_potential == "PAWLDAKIN") {
      string_POTENTIAL = DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN;
    }
    if (params_potential == "PAWPBEKIN") {
      string_POTENTIAL = DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " string_POTENTIAL=" << string_POTENTIAL << endl;
    }

    //   PARAMS.vparams.flag("AFLOWING_STRING_POTENTIAL",vpflow.flag("PROTO_AFLOW::POTENTIAL"));
    // if(PARAMS.vparams.flag("AFLOWING_STRING_POTENTIAL")) PARAMS.vparams.push_attached("AFLOWING_STRING_POTENTIAL",vpflow.getattachedscheme("PROTO_AFLOW::POTENTIAL"));

    // check potential_type  //CO20191110
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK POTENTIAL_TYPE" << endl; // CO20191110
    }
    if (vpflow.flag("PROTO_AFLOW::POTENTIAL_TYPE")) {
      string_POTENTIAL += _AVASP_PSEUDOPOTENTIAL_DELIMITER_ + _AVASP_PSEUDOPOTENTIAL_POTENTIAL_TYPE_; // CO20191110
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"PROTO_AFLOW::POTENTIAL_TYPE\")=" << vpflow.flag("PROTO_AFLOW::POTENTIAL_TYPE") << endl; // CO20191110
    }

    // check potential_complete
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK POTENTIAL_COMPLETE" << endl;
    }
    if (vpflow.flag("PROTO_AFLOW::POTENTIAL_COMPLETE")) {
      string_POTENTIAL += _AVASP_PSEUDOPOTENTIAL_DELIMITER_ + _AVASP_PSEUDOPOTENTIAL_POTENTIAL_COMPLETE_;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"PROTO_AFLOW::POTENTIAL_COMPLETE\")=" << vpflow.flag("PROTO_AFLOW::POTENTIAL_COMPLETE") << endl;
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " string_POTENTIAL=" << string_POTENTIAL << endl;
    }
    PARAMS.vparams.push_attached("AFLOWIN_STRING::POTENTIAL", string_POTENTIAL);

    // reverse
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK REVERSE" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::REVERSE", flag_REVERSE);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::REVERSE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::REVERSE") << endl;
    }

    // check missing
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK MISSING" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::MISSING", vpflow.flag("PROTO_AFLOW::MISSING"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::MISSING\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::MISSING") << endl;
    }

    // check noautopp
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK NOAUTOPP" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::NOAUTOPP", vpflow.flag("PROTO_AFLOW::NOAUTOPP"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::NOAUTOPP\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::NOAUTOPP") << endl;
    }

    // check bader //CO20181226
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK BADER" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::BADER", vpflow.flag("PROTO_AFLOW::BADER"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::BADER\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::BADER") << endl;
    }

    // check spin_remove_relax_1 //CO20181226
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK SPIN_REMOVE_RELAX_1" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::SPIN_REMOVE_RELAX_1", vpflow.flag("PROTO_AFLOW::SPIN_REMOVE_RELAX_1"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::SPIN_REMOVE_RELAX_1\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::SPIN_REMOVE_RELAX_1") << endl;
    }

    // check spin_remove_relax_2 //CO20181226
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK SPIN_REMOVE_RELAX_2" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::SPIN_REMOVE_RELAX_2", vpflow.flag("PROTO_AFLOW::SPIN_REMOVE_RELAX_2"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::SPIN_REMOVE_RELAX_2\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::SPIN_REMOVE_RELAX_2") << endl;
    }

    // check ldau
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK LDAU" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::AUTOLDAU", vpflow.flag("PROTO_AFLOW::LDAU"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::AUTOLDAU\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::AUTOLDAU") << endl;
    }

    // check noldau
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK NOLDAU" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::AUTONOLDAU", vpflow.flag("PROTO_AFLOW::NOLDAU"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::AUTONOLDAU\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::AUTONOLDAU") << endl;
    }

    // check neglect_nomix
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK NEGLECT_NOMIX" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::NEGLECT_NOMIX", vpflow.flag("PROTO_AFLOW::NEGLECT_NOMIX"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::NEGLECT_NOMIX\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::NEGLECT_NOMIX") << endl;
    }

    // check stdout
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK STDOUT" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::STDOUT", vpflow.flag("PROTO_AFLOW::STDOUT"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::STDOUT\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::STDOUT") << endl;
    }

    // check module
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK MODULE" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::MODULE", vpflow.flag("AFLOW::MODULE"));
    if (vpflow.flag("AFLOW::MODULE")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::MODULE", vpflow.getattachedscheme("AFLOW::MODULE"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::MODULE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::MODULE") << endl;
    }

    // check apl_supercell
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK APL_SUPERCELL" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::APL_SUPERCELL", vpflow.flag("PROTO_AFLOW::APL_SUPERCELL"));
    if (vpflow.flag("PROTO_AFLOW::APL_SUPERCELL")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::APL_SUPERCELL", vpflow.getattachedscheme("PROTO_AFLOW::APL_SUPERCELL"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::APL_SUPERCELL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::APL_SUPERCELL") << endl;
    }

    // check potim
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK POTIM" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::POTIM", vpflow.flag("PROTO_AFLOW::POTIM"));
    if (vpflow.flag("PROTO_AFLOW::POTIM")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::POTIM", vpflow.getattachedscheme("PROTO_AFLOW::POTIM"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::POTIM\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::POTIM") << endl;
    }

    // check precision
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK PRECISION" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::PRECISION", vpflow.flag("PROTO_AFLOW::PRECISION"));
    if (vpflow.flag("PROTO_AFLOW::PRECISION")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::PRECISION", vpflow.getattachedscheme("PROTO_AFLOW::PRECISION"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::PRECISION\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::PRECISION") << endl;
    }

    // check algorithm
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK ALGORITHM" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::ALGORITHM", vpflow.flag("PROTO_AFLOW::ALGORITHM"));
    if (vpflow.flag("PROTO_AFLOW::ALGORITHM")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::ALGORITHM", vpflow.getattachedscheme("PROTO_AFLOW::ALGORITHM"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::ALGORITHM\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ALGORITHM") << endl;
    }

    // check metagga
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK METAGGA" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::METAGGA", vpflow.flag("PROTO_AFLOW::METAGGA"));
    if (vpflow.flag("PROTO_AFLOW::METAGGA")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::METAGGA", vpflow.getattachedscheme("PROTO_AFLOW::METAGGA"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::METAGGA\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::METAGGA") << endl;
    }

    // check ivdw
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK IVDW" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::IVDW", vpflow.flag("PROTO_AFLOW::IVDW"));
    if (vpflow.flag("PROTO_AFLOW::IVDW")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::IVDW", vpflow.getattachedscheme("PROTO_AFLOW::IVDW"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::IVDW\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::IVDW") << endl;
    }

    // check relax_type
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK RELAX_TYPE" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::RELAX_TYPE", vpflow.flag("PROTO_AFLOW::RELAX_TYPE"));
    if (vpflow.flag("PROTO_AFLOW::RELAX_TYPE")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::RELAX_TYPE", vpflow.getattachedscheme("PROTO_AFLOW::RELAX_TYPE"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RELAX_TYPE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RELAX_TYPE") << endl;
    }

    // check relax_mode
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK RELAX_MODE" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::RELAX_MODE", vpflow.flag("PROTO_AFLOW::RELAX_MODE"));
    if (vpflow.flag("PROTO_AFLOW::RELAX_MODE")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::RELAX_MODE", vpflow.getattachedscheme("PROTO_AFLOW::RELAX_MODE"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RELAX_MODE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RELAX_MODE") << endl;
    }

    // check relax_count  //CO20181226
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK RELAX_COUNT" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::RELAX_COUNT", vpflow.flag("PROTO_AFLOW::RELAX_COUNT"));
    if (vpflow.flag("PROTO_AFLOW::RELAX_COUNT")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::RELAX_COUNT", vpflow.getattachedscheme("PROTO_AFLOW::RELAX_COUNT"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RELAX_COUNT\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RELAX_COUNT") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::RELAX_COUNT\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::RELAX_COUNT") << endl;
    }

    // check run_relax_static //CO20181226
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK RUN_RELAX_STATIC" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_RELAX_STATIC", vpflow.flag("PROTO_AFLOW::RUN_RELAX_STATIC"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RUN_RELAX_STATIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_RELAX_STATIC") << endl;
    }

    // check run_relax_static_bands //CO20181226
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK RUN_RELAX_STATIC_BANDS" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_RELAX_STATIC_BANDS", vpflow.flag("PROTO_AFLOW::RUN_RELAX_STATIC_BANDS"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RUN_RELAX_STATIC_BANDS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_RELAX_STATIC_BANDS") << endl;
    }

    // check run_relax_static_dielectric
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK RUN_RELAX_STATIC_DIELECTRIC" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_RELAX_STATIC_DIELECTRIC", vpflow.flag("PROTO_AFLOW::RUN_RELAX_STATIC_DIELECTRIC"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RUN_RELAX_STATIC_DIELECTRIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_RELAX_STATIC_DIELECTRIC") << endl;
    }

    // check run_static //CO20181226
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK RUN_STATIC" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_STATIC", vpflow.flag("PROTO_AFLOW::RUN_STATIC"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RUN_STATIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_STATIC") << endl;
    }

    // check run_static_bands //CO20181226
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK RUN_STATIC_BANDS" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_STATIC_BANDS", vpflow.flag("PROTO_AFLOW::RUN_STATIC_BANDS"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RUN_STATIC_BANDS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_STATIC_BANDS") << endl;
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK RUN_STATIC_DIELECTRIC" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_STATIC_DIELECTRIC", vpflow.flag("PROTO_AFLOW::RUN_STATIC_DIELECTRIC"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RUN_STATIC_DIELECTRIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_STATIC_DIELECTRIC") << endl;
    }

    // check type
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK TYPE" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::TYPE", vpflow.flag("PROTO_AFLOW::TYPE"));
    if (vpflow.flag("PROTO_AFLOW::TYPE")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::TYPE", vpflow.getattachedscheme("PROTO_AFLOW::TYPE"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::TYPE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::TYPE") << endl;
    }

    // check CONVERT_UNIT_CELL
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK CONVERT_UNIT_CELL" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::CONVERT_UNIT_CELL", vpflow.flag("PROTO_AFLOW::CONVERT_UNIT_CELL"));
    if (vpflow.flag("PROTO_AFLOW::CONVERT_UNIT_CELL")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::CONVERT_UNIT_CELL", vpflow.getattachedscheme("PROTO_AFLOW::CONVERT_UNIT_CELL"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::CONVERT_UNIT_CELL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::CONVERT_UNIT_CELL") << endl;
    }

    // check volume_plus_equal
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK VOLUME_PLUS_EQUAL" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL", vpflow.flag("PROTO_AFLOW::VOLUME_PLUS_EQUAL"));
    if (vpflow.flag("PROTO_AFLOW::VOLUME_PLUS_EQUAL")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL", vpflow.getattachedscheme("PROTO_AFLOW::VOLUME_PLUS_EQUAL"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::VOLUME_PLUS_EQUAL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL") << endl;
    }

    // check volume_multiply_equal
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK VOLUME_MULTIPLY_EQUAL" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL", vpflow.flag("PROTO_AFLOW::VOLUME_MULTIPLY_EQUAL"));
    if (vpflow.flag("PROTO_AFLOW::VOLUME_MULTIPLY_EQUAL")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL", vpflow.getattachedscheme("PROTO_AFLOW::VOLUME_MULTIPLY_EQUAL"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL") << endl;
    }

    // check volume_preserved //CO20181226
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK VOLUME_PRESERVED" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::VOLUME_PRESERVED", vpflow.flag("PROTO_AFLOW::VOLUME_PRESERVED"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::VOLUME_PRESERVED\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::VOLUME_PRESERVED") << endl;
    }

    // check ediffg
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK EDIFFG" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::EDIFFG", vpflow.flag("PROTO_AFLOW::EDIFFG"));
    if (vpflow.flag("PROTO_AFLOW::EDIFFG")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::EDIFFG", vpflow.getattachedscheme("PROTO_AFLOW::EDIFFG"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::EDIFFG\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::EDIFFG") << endl;
    }

    // check kscheme  //CO20181226
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK KSCHEME" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::KSCHEME", vpflow.flag("PROTO_AFLOW::KSCHEME"));
    if (vpflow.flag("PROTO_AFLOW::KSCHEME")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::KSCHEME", vpflow.getattachedscheme("PROTO_AFLOW::KSCHEME"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::KSCHEME\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::KSCHEME") << endl;
    }

    // check kscheme_static //CO20181226
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK KSCHEME_STATIC" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::KSCHEME_STATIC", vpflow.flag("PROTO_AFLOW::KSCHEME_STATIC"));
    if (vpflow.flag("PROTO_AFLOW::KSCHEME_STATIC")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::KSCHEME_STATIC", vpflow.getattachedscheme("PROTO_AFLOW::KSCHEME_STATIC"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::KSCHEME_STATIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::KSCHEME_STATIC") << endl;
    }

    // check kscheme_dielectric
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK KSCHEME_DIELECTRIC" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::KSCHEME_DIELECTRIC", vpflow.flag("PROTO_AFLOW::KSCHEME_DIELECTRIC"));
    if (vpflow.flag("PROTO_AFLOW::KSCHEME_DIELECTRIC")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::KSCHEME_DIELECTRIC", vpflow.getattachedscheme("PROTO_AFLOW::KSCHEME_DIELECTRIC"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::KSCHEME_DIELECTRIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::KSCHEME_DIELECTRIC") << endl;
    }

    // check kppra
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK KPPRA" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::KPPRA", vpflow.flag("PROTO_AFLOW::KPPRA"));
    if (vpflow.flag("PROTO_AFLOW::KPPRA")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::KPPRA", vpflow.getattachedscheme("PROTO_AFLOW::KPPRA"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::KPPRA\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::KPPRA") << endl;
    }

    // check kppra_static
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK KPPRA_STATIC" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::KPPRA_STATIC", vpflow.flag("PROTO_AFLOW::KPPRA_STATIC"));
    if (vpflow.flag("PROTO_AFLOW::KPPRA_STATIC")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::KPPRA_STATIC", vpflow.getattachedscheme("PROTO_AFLOW::KPPRA_STATIC"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::KPPRA_STATIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::KPPRA_STATIC") << endl;
    }

    // check kppra_dielectric
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK KPPRA_DIELECTRIC" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::KPPRA_DIELECTRIC", vpflow.flag("PROTO_AFLOW::KPPRA_DIELECTRIC"));
    if (vpflow.flag("PROTO_AFLOW::KPPRA_DIELECTRIC")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::KPPRA_DIELECTRIC", vpflow.getattachedscheme("PROTO_AFLOW::KPPRA_DIELECTRIC"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::KPPRA_DIELECTRIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::KPPRA_DIELECTRIC") << endl;
    }

    // check bands_grid
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK BANDS_GRID" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::BANDS_GRID", vpflow.flag("PROTO_AFLOW::BANDS_GRID"));
    if (vpflow.flag("PROTO_AFLOW::BANDS_GRID")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::BANDS_GRID", vpflow.getattachedscheme("PROTO_AFLOW::BANDS_GRID"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::BANDS_GRID\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::BANDS_GRID") << endl;
    }

    // check enmax_multiply
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK ENMAX_MULTIPLY" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::ENMAX_MULTIPLY", vpflow.flag("PROTO_AFLOW::ENMAX_MULTIPLY"));
    if (vpflow.flag("PROTO_AFLOW::ENMAX_MULTIPLY")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::ENMAX_MULTIPLY", vpflow.getattachedscheme("PROTO_AFLOW::ENMAX_MULTIPLY"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::ENMAX_MULTIPLY\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ENMAX_MULTIPLY") << endl;
    }

    // check ABINIT QE VASP AIMS
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK VASP/ABCCAR/ITC/QE/ABINIT/AIMS/CIF/ELK/LMP" << endl; // DX20190123 - add CIF //DX20200313 - add ELK  //CO20220613
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::VASP", true); // default
    PARAMS.vparams.flag("AFLOWIN_FLAG::ITC", vpflow.flag("PROTO_AFLOW::ITC")); // CO20220613
    PARAMS.vparams.flag("AFLOWIN_FLAG::QE", vpflow.flag("PROTO_AFLOW::QE"));
    PARAMS.vparams.flag("AFLOWIN_FLAG::ABCCAR", vpflow.flag("PROTO_AFLOW::ABCCAR")); // DX20190123 - add ABCCAR
    PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT", vpflow.flag("PROTO_AFLOW::ABINIT"));
    PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS", vpflow.flag("PROTO_AFLOW::AIMS"));
    PARAMS.vparams.flag("AFLOWIN_FLAG::CIF", vpflow.flag("PROTO_AFLOW::CIF")); // DX20190123 - add CIF
    PARAMS.vparams.flag("AFLOWIN_FLAG::ELK", vpflow.flag("PROTO_AFLOW::ELK")); // DX20200313 - add ELK
    PARAMS.vparams.flag("AFLOWIN_FLAG::LMP", vpflow.flag("PROTO_AFLOW::LMP")); // SD20240111 - add LMP
    if (PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS")) {
      PARAMS.vparams.flag("AFLOWIN_FLAG::QE", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::VASP", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::CIF", false); // DX20190123 - add CIF
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABCCAR", false); // DX20190123 - add ABCCAR
      PARAMS.vparams.flag("AFLOWIN_FLAG::ELK", false); // DX20200313 - add ELK
      PARAMS.vparams.flag("AFLOWIN_FLAG::LMP", false); // SD20240111 - add LMP
    }
    if (PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT")) {
      PARAMS.vparams.flag("AFLOWIN_FLAG::QE", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::VASP", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::CIF", false); // DX20190123 - add CIF
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABCCAR", false); // DX20190123 - add ABCCAR
      PARAMS.vparams.flag("AFLOWIN_FLAG::ELK", false); // DX20200313 - add ELK
      PARAMS.vparams.flag("AFLOWIN_FLAG::LMP", false); // SD20240111 - add LMP
    }
    // CO20220613 - add ITC - START
    if (PARAMS.vparams.flag("AFLOWIN_FLAG::ITC")) { // CO20220613
      // CO20220613 - --itc can be added with another flag
    }
    // DX20190123 - add ABCCAR - START
    if (PARAMS.vparams.flag("AFLOWIN_FLAG::ABCCAR")) {
      PARAMS.vparams.flag("AFLOWIN_FLAG::QE", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::VASP", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::CIF", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::ELK", false); // DX20200313 - add ELK
      PARAMS.vparams.flag("AFLOWIN_FLAG::LMP", false); // SD20240111 - add LMP
    }
    // DX20190123 - add ABCCAR - END
    // DX20190123 - add CIF - START
    if (PARAMS.vparams.flag("AFLOWIN_FLAG::CIF")) {
      PARAMS.vparams.flag("AFLOWIN_FLAG::QE", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::VASP", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABCCAR", false); // DX20190123 - add ABCCAR
      PARAMS.vparams.flag("AFLOWIN_FLAG::ELK", false); // DX201200313 - add ELK
      PARAMS.vparams.flag("AFLOWIN_FLAG::LMP", false); // SD20240111 - add LMP
    }
    // DX20190123 - add CIF - END
    if (PARAMS.vparams.flag("AFLOWIN_FLAG::QE")) {
      PARAMS.vparams.flag("AFLOWIN_FLAG::VASP", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::CIF", false); // DX20190123 - add CIF
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABCCAR", false); // DX20190123 - add ABCCAR
      PARAMS.vparams.flag("AFLOWIN_FLAG::ELK", false); // DX201200313 - add ELK
      PARAMS.vparams.flag("AFLOWIN_FLAG::LMP", false); // SD20240111 - add LMP
    }
    if (PARAMS.vparams.flag("AFLOWIN_FLAG::VASP")) {
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::QE", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::CIF", false); // DX20190123 - add CIF
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABCCAR", false); // DX20190123 - add ABCCAR
      PARAMS.vparams.flag("AFLOWIN_FLAG::ELK", false); // DX201200313 - add ELK
      PARAMS.vparams.flag("AFLOWIN_FLAG::LMP", false); // SD20240111 - add LMP
    }
    // DX20200313 - add ELK - START
    if (PARAMS.vparams.flag("AFLOWIN_FLAG::ELK")) {
      PARAMS.vparams.flag("AFLOWIN_FLAG::QE", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::VASP", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABCCAR", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::CIF", false);
    }
    // DX20200313 - add ELK - END
    // SD20240111 - add LMP - START
    if (PARAMS.vparams.flag("AFLOWIN_FLAG::LMP")) {
      PARAMS.vparams.flag("AFLOWIN_FLAG::QE", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::VASP", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::ABCCAR", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::CIF", false);
      PARAMS.vparams.flag("AFLOWIN_FLAG::ELK", false);
    }
    // SD20240111 - add LMP - END
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::AIMS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::ABINIT\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::QE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::QE") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::VASP\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::VASP") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::ITC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ITC") << endl; // CO20220613
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::CIF\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::CIF") << endl; // DX20190123 - add CIF
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::ABCCAR\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ABCCAR") << endl; // DX20190123 - add ABCCAR
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::ELK\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ELK") << endl; // DX20200313 - add ELK
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::LMP\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::LMP") << endl; // SD20240111 - add LMP
    }

    // check stdout
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK STDOUT" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::STDOUT", vpflow.flag("PROTO_AFLOW::STDOUT"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::STDOUT\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::STDOUT") << endl;
    }

    // check list
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK LIST" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::LIST", vpflow.flag("PROTO_AFLOW::LIST"));
    if (PARAMS.vparams.flag("AFLOWIN_FLAG::LIST")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::LIST_VCMD", vpflow.getattachedscheme("PROTO_AFLOW::LIST_VCMD"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::LIST\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::LIST") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme((\"AFLOWIN_FLAG::LIST_VCMD\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::LIST_VCMD") << endl;
    }

    // DX20180118 - Add ANRL functionality - START
    //  check params
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK PARAMS" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::PARAMS", vpflow.flag("PARAMS"));
    if (vpflow.flag("PARAMS")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::PARAMS", vpflow.getattachedscheme("PARAMS"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::PARAMS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::PARAMS") << endl;
    }
    // DX20180118 - Add ANRL functionality - END

    // DX20190227 - add anrl lattice parameter flag - START
    PARAMS.vparams.flag("AFLOWIN_FLAG::USE_ANRL_LATTICE_PARAM", vpflow.flag("PROTO::USE_ANRL_LATTICE_PARAM"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::USE_ANRL_LATTICE_PARAM\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::USE_ANRL_LATTICE_PARAM") << endl;
    }
    // DX20190227 - add anrl lattice parameter flag - END

    // check htqc
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK HTQC_ICSD" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::HTQC_ICSD", vpflow.flag("PROTO_AFLOW::HTQC") || (!aurostd::substring2bool(string_LABEL, "_ICSD_") && aurostd::substring2bool(string_LABEL, "ICSD_")));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::HTQC_ICSD\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::HTQC_ICSD") << endl;
    }

    // check pressure
    vector<string> params_pressure;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"PROTO_AFLOW::PRESSURE\")=" << vpflow.getattachedscheme("PROTO_AFLOW::PRESSURE") << endl;
    }
    aurostd::string2tokens(vpflow.getattachedscheme("PROTO_AFLOW::PRESSURE"), params_pressure, ",");
    PARAMS.vpressure.clear();
    if (!params_pressure.empty()) {
      for (size_t i = 0; i < params_pressure.size(); i++) {
        PARAMS.vpressure.push_back(aurostd::string2utype<double>(params_pressure[i]));
      }
    }

    // CO20181226 - check pocc_params
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK POCC_PARAMS" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::POCC_PARAMS", vpflow.flag("POCC_PARAMS"));
    if (vpflow.flag("POCC_PARAMS")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::POCC_PARAMS", vpflow.getattachedscheme("POCC_PARAMS"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::POCC_PARAMS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::POCC_PARAMS") << endl;
    }

    // CO20181226 - check pocc_tol
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK POCC_TOL" << endl;
    }
    PARAMS.vparams.flag("AFLOWIN_FLAG::POCC_TOL", vpflow.flag("POCC_TOL"));
    if (vpflow.flag("POCC_TOL")) {
      PARAMS.vparams.push_attached("AFLOWIN_FLAG::POCC_TOL", vpflow.getattachedscheme("POCC_TOL"));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::POCC_TOL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::POCC_TOL") << endl;
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " --- " << endl;
    }

    const uint nspeciesHTQC = aflowlib::PrototypeLibrariesSpeciesNumber(string_LABEL);
    const bool alphabetic = RequestedAlphabeticLabeling(string_LABEL);
    PARAMS.vparams.push_attached("AFLOWIN_STRING::LABEL", string_LABEL);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " nspeciesHTQC=" << nspeciesHTQC << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " alphabetic=" << alphabetic << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vpressure.size()=" << PARAMS.vpressure.size() << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vkppra.size()=" << PARAMS.vkppra.size() << " " << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme((\"AFLOWIN_STRING::LABEL\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_STRING::LABEL") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme((\"AFLOWIN_STRING::POTENTIAL\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_STRING::POTENTIAL") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::MODULE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::MODULE") << endl; // CO20180214
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::MODULE\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::MODULE") << endl; // CO20180214
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::APL_SUPERCELL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::APL_SUPERCELL") << endl; // CO20180214
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::APL_SUPERCELL\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::APL_SUPERCELL") << endl; // CO20180214
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::POTIM\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::POTIM") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::POTIM\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::POTIM") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::PRECISION\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::PRECISION") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::PRECISION\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::PRECISION") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::ALGORITHM\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ALGORITHM") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::ALGORITHM\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::ALGORITHM") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::METAGGA\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::METAGGA") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::METAGGA\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::METAGGA") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::IVDW\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::IVDW") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::IVDW\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::IVDW") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RELAX_TYPE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RELAX_TYPE") << endl; // CO20180214
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::RELAX_TYPE\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::RELAX_TYPE") << endl; // CO20180214
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RELAX_MODE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RELAX_MODE") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::RELAX_MODE\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::RELAX_MODE") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RELAX_COUNT\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RELAX_COUNT") << endl; // CO20181226
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::RELAX_COUNT\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::RELAX_COUNT") << endl; // CO20181226
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RUN_RELAX_STATIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_RELAX_STATIC") << endl; // CO20181226
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RUN_RELAX_STATIC_BANDS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_RELAX_STATIC_BANDS") << endl; // CO20181226
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RUN_RELAX_STATIC_DIELECTRIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_RELAX_STATIC_DIELECTRIC") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RUN_STATIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_STATIC") << endl; // CO20181226
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RUN_STATIC_BANDS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_STATIC_BANDS") << endl; // CO20181226
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::RUN_STATIC_DIELECTRIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::RUN_STATIC_DIELECTRIC") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::TYPE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::TYPE") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::TYPE\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::TYPE") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::CONVERT_UNIT_CELL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::CONVERT_UNIT_CELL") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::CONVERT_UNIT_CELL\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::CONVERT_UNIT_CELL") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::VOLUME_PLUS_EQUAL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::VOLUME_PLUS_EQUAL\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::VOLUME_PLUS_EQUAL") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::VOLUME_MULTIPLY_EQUAL") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::VOLUME_PRESERVED\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::VOLUME_PRESERVED") << endl; // CO20180214
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::EDIFFG\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::EDIFFG") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::EDIFFG\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::EDIFFG") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::KPPRA\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::KPPRA") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::KPPRA\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::KPPRA") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::KPPRA_STATIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::KPPRA_STATIC") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::KPPRA_STATIC\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::KPPRA_STATIC") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::KPPRA_DIELECTRIC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::KPPRA_DIELECTRIC") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::KPPRA_DIELECTRIC\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::KPPRA_DIELECTRIC") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::BANDS_GRID\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::BANDS_GRID") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::BANDS_GRID\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::BANDS_GRID") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::ENMAX_MULTIPLY\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ENMAX_MULTIPLY") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme(\"AFLOWIN_FLAG::ENMAX_MULTIPLY\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::ENMAX_MULTIPLY") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::USAGE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::USAGE") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::REVERSE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::REVERSE") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::MISSING\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::MISSING") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::NOAUTOPP\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::NOAUTOPP") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::AUTOLDAU\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::AUTOLDAU") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::AUTONOLDAU\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::AUTONOLDAU") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::STDOUT\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::STDOUT") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::LIST\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::LIST") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::PARAMS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::PARAMS") << endl; // DX20180118 - Added ANRL functionality
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::POCC_PARAMS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::POCC_PARAMS") << endl; // CO20181226
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::POCC_TOL\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::POCC_TOL") << endl; // CO20181226
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.getattachedscheme((\"AFLOWIN_FLAG::LIST_VCMD\")=" << PARAMS.vparams.getattachedscheme("AFLOWIN_FLAG::LIST_VCMD") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::HTQC_ICSD\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::HTQC_ICSD") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::NEGLECT_NOMIX\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::NEGLECT_NOMIX") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::VASP\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::VASP") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::ITC\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ITC") << endl; // CO20220613
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::ABINIT\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ABINIT") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::AIMS\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::AIMS") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::QE\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::QE") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::CIF\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::CIF") << endl; // DX20190123 - add CIF
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::ABCCAR\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ABCCAR") << endl; // DX20190123 - add ABCCAR
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::ELK\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::ELK") << endl; // DX202003133 - add ELK
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARAMS.vparams.flag(\"AFLOWIN_FLAG::LMP\")=" << PARAMS.vparams.flag("AFLOWIN_FLAG::LMP") << endl; // SD202401111 - add LMP
    }

    if (PARAMS.ucell.size() == 1 && aurostd::substring2bool(PARAMS.vparams.getattachedscheme("AFLOWIN_STRING::LABEL"), "_ICSD_") && !PARAMS.vparams.flag("AFLOWIN_FLAG::HTQC_ICSD")) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " running AVASP_MakePrototypeICSD_AFLOWIN" << endl;
      }
      //  return AVASP_MakePrototypeICSD_AFLOWIN(PARAMS.ucell,false);
      return AVASP_MakePrototypeICSD_AFLOWIN((&PARAMS), false);
    } else {
      // if(alphabetic==true && (nspeciesHTQC==3 || nspeciesHTQC==4)) AlphabetizePrototypeLabelSpeciesArgv(argv); // it is going to fail if you want a big list of protos
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " running AVASP_MakePrototype_AFLOWIN" << endl;
      }
      return AVASP_MakePrototype_AFLOWIN((&PARAMS));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }

    return false;
  }
} // namespace pflow

namespace pflow {
  vector<string> GENERATE_CERAMICS(const vector<string>& _vnonmetals, const vector<string>& _vmetals, uint metal_arity) { // CO20200731
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " nonmetals=" << aurostd::joinWDelimiter(_vnonmetals, ",") << endl;
      cerr << __AFLOW_FUNC__ << " metals=" << aurostd::joinWDelimiter(_vmetals, ",") << endl;
      cerr << __AFLOW_FUNC__ << " metal_arity=" << metal_arity << endl;
    }
    if (_vnonmetals.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vnonmetals.size()==0", _INPUT_MISSING_);
    } // CO20200624
    if (_vmetals.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vmetals.size()==0", _INPUT_MISSING_);
    } // CO20200624
    vector<string> vnonmetals(_vnonmetals);
    vector<string> vmetals(_vmetals);
    std::sort(vnonmetals.begin(), vnonmetals.end()); // SORT FIRST
    std::sort(vmetals.begin(), vmetals.end()); // SORT FIRST

    uint i = 0;
    uint j = 0;

    // get vmix
    vector<string> vmix = vnonmetals;
    vmix.insert(vmix.end(), vmetals.begin(), vmetals.end());
    std::sort(vmix.begin(), vmix.end());

    // get vind_nm
    vector<uint> vind_nm;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vmix=" << aurostd::joinWDelimiter(vmix, ",") << endl;
    }
    std::vector<string>::iterator it;
    for (i = 0; i < vnonmetals.size(); i++) {
      it = std::find(vmix.begin(), vmix.end(), vnonmetals[i]);
      vind_nm.push_back(it - vmix.begin());
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vind_nm=" << aurostd::joinWDelimiter(vind_nm, ",") << endl;
    }

    // get vvinbetween
    vector<vector<string>> vvinbetween;
    vvinbetween.resize(vnonmetals.size() + 1); // pre==vvinbetween[0], post==vvinbetween[-1]
    for (i = 0, j = 0; i < vmix.size(); i++) {
      if (aurostd::WithinList(vind_nm, i)) {
        j++;
        continue;
      }
      vvinbetween[j].push_back(vmix[i]);
    }
    if (LDEBUG) {
      for (i = 0; i < vvinbetween.size(); i++) {
        cerr << __AFLOW_FUNC__ << " vvinbetween[i=" << i << "]=" << aurostd::joinWDelimiter(vvinbetween[i], ",") << endl;
      }
    }
    // get vNinbetween
    vector<uint> vNinbetween;
    vNinbetween.assign(vvinbetween.size(), 0); // assign all 0
    uint sum_species = 0; // colons
    // think of it as drops of water into separate cups, start with first and fill until you have no more water
    sum_species = 0;
    for (i = 0; i < vvinbetween.size() && sum_species < metal_arity; i++) {
      for (j = 0; j < vvinbetween[i].size() && sum_species < metal_arity; j++) {
        vNinbetween[i]++;
        sum_species++;
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vNinbetween=" << aurostd::joinWDelimiter(vNinbetween, ",") << endl;
    }

    // get vNinbetween_std
    vector<uint> vNinbetween_std;
    vNinbetween_std.assign(vNinbetween.size(), 0);
    vNinbetween_std.back() = metal_arity;

    uint k = 0; // index of nonmetals
    vector<string> vcommands;
    vector<string> vtmp;
    string command = "";
    string tmp;
    while (vcommands.size() < 1e4) { // not simple as number of combinations since we could have 2 or more nonmetals
      command = "";

      for (i = 0, k = 0; i < vNinbetween.size(); i++) {
        tmp = aurostd::joinWDelimiter(vvinbetween[i], ",");
        vtmp.clear();
        if (vNinbetween[i] > 0) {
          if (!command.empty()) {
            command += ":";
          }
          for (j = 0; j < vNinbetween[i]; j++) {
            vtmp.push_back(tmp);
          }
          command += aurostd::joinWDelimiter(vtmp, ":");
        }
        if (k < vnonmetals.size()) {
          if (!command.empty()) {
            command += ":";
          }
          command += vnonmetals[k++];
        }
      }

      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " command=" << command << endl;
      }
      vcommands.push_back(command);
      if (vNinbetween == vNinbetween_std) {
        break;
      }

      // adjust vNinbetween
      for (i = 0; i < vNinbetween.size() - 1; i++) {
        if (vNinbetween[i] == 0) {
          continue;
        }
        vNinbetween[i]--;
        vNinbetween[i + 1]++;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " vNinbetween=" << aurostd::joinWDelimiter(vNinbetween, ",") << endl;
      }
    }

    return vcommands;
  }
  vector<string> GENERATE_CERAMICS(const aurostd::xoption& vpflow) { // CO20200731
    vector<string> vnonmetals;
    vector<string> vmetals;
    const string nonmetals = vpflow.getattachedscheme("GENERATE_CERAMICS::NON_METALS");
    aurostd::string2tokens(nonmetals, vnonmetals, ",");
    const string metals = vpflow.getattachedscheme("GENERATE_CERAMICS::METALS");
    aurostd::string2tokens(metals, vmetals, ",");
    const int metal_arity = vpflow.getattachedutype<uint>("GENERATE_CERAMICS::METAL_ARITY");
    return GENERATE_CERAMICS(vnonmetals, vmetals, metal_arity);
  }
  string GENERATE_CERAMICS_PRINT(const aurostd::xoption& vpflow) { // CO20200731
    return aurostd::joinWDelimiter(GENERATE_CERAMICS(vpflow), "\n");
  }
} // namespace pflow

// ***************************************************************************
// pflow::PSEUDOPOTENTIALS_CHECK
// ***************************************************************************

namespace pflow {
  bool PSEUDOPOTENTIALS_CHECK(const aurostd::xoption& vpflow, const string& file, ostream& oss) { // too many options
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }

    // check usage
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " CHECK USAGE" << endl;
    }
    if (vpflow.flag("PSEUDOPOTENTIALS_CHECK::USAGE") || file.empty()) {
      init::MessageOption(vpflow.getattachedscheme("PSEUDOPOTENTIALS_CHECK"), "pflow::PSEUDOPOTENTIALS_CHECK",
                          aurostd::liststring2string("aflow [options] --pseudopotentials_check=[POTCAR|OUTCAR]["
                                                     "|.bz2|.gz|.xz] | --pp_check= | --ppk=",
                                                     "       options:", "                --usage"));
      return true; // CO20200624 - the option was expressed successfully
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"PSEUDOPOTENTIALS_CHECK::USAGE\")=" << vpflow.flag("PSEUDOPOTENTIALS_CHECK::USAGE") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " file=" << file << endl;
    }
    // now start
    if (aurostd::substring2bool(file, "POTCAR")) {
      XHOST.DEBUG = false;
      XHOST.PSEUDOPOTENTIAL_GENERATOR = false;
      const xPOTCAR xPOT(file);
      oss << xPOT << endl;
    }
    if (aurostd::substring2bool(file, "OUTCAR")) {
      XHOST.DEBUG = false;
      XHOST.PSEUDOPOTENTIAL_GENERATOR = false;
      xOUTCAR xOUT(file);
      const xPOTCAR xPOT;
      oss << " FROM OUTCAR" << endl;
      oss << " vTITEL.size()=" << xOUT.vTITEL.size() << ": ";
      for (size_t i = 0; i < xOUT.vTITEL.size(); i++) {
        oss << xOUT.vTITEL[i] << " ";
      }
      oss << endl;
      oss << " SYSTEM=" << xOUT.SYSTEM << endl;
      oss << " pp_type=" << xOUT.pp_type << endl;
      oss << " species.size()=" << xOUT.species.size() << ": ";
      for (size_t i = 0; i < xOUT.species.size(); i++) {
        oss << xOUT.species[i] << " ";
      }
      oss << endl;
      oss << " species_pp.size()=" << xOUT.species_pp.size() << ": ";
      for (size_t i = 0; i < xOUT.species_pp.size(); i++) {
        oss << xOUT.species_pp[i] << " ";
      }
      oss << endl;
      oss << " species_Z.size()=" << xOUT.species_Z.size() << ": ";
      for (size_t i = 0; i < xOUT.species_Z.size(); i++) {
        oss << xOUT.species_Z[i] << " ";
      }
      oss << endl;
      oss << " species_pp_type.size()=" << xOUT.species_pp_type.size() << ": ";
      for (size_t i = 0; i < xOUT.species_pp_type.size(); i++) {
        oss << xOUT.species_pp_type[i] << " ";
      }
      oss << endl;
      oss << " species_pp_version.size()=" << xOUT.species_pp_version.size() << ": ";
      for (size_t i = 0; i < xOUT.species_pp_version.size(); i++) {
        oss << xOUT.species_pp_version[i] << " ";
      }
      oss << endl;
      oss << " species_pp_AUID.size()=" << xOUT.species_pp_AUID.size() << ": ";
      for (size_t i = 0; i < xOUT.species_pp_AUID.size(); i++) {
        oss << xOUT.species_pp_AUID[i] << " ";
      }
      oss << endl;
      oss << " species_pp_AUID_collisions.size()=" << xOUT.species_pp_AUID_collisions.size() << ": ";
      for (size_t i = 0; i < xOUT.species_pp_AUID_collisions.size(); i++) {
        oss << xOUT.species_pp_AUID_collisions[i] << " ";
      }
      oss << endl;
    }
    //
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return true;
  }
} // namespace pflow

// ***************************************************************************
// pflow::PYTHON_MODULES //ME20211103
// ***************************************************************************

namespace pflow {
  void PYTHON_MODULES(const string& modules, ostream& oss) {
    ofstream FileMESSAGE;
    PYTHON_MODULES(modules, FileMESSAGE, oss);
  }

  void PYTHON_MODULES(const string& modules, ofstream& FileMESSAGE, ostream& oss) {
    vector<string> vmodules;
    aurostd::string2tokens(modules, vmodules, ",");
    PYTHON_MODULES(vmodules, FileMESSAGE, oss);
  }

  void PYTHON_MODULES(const vector<string>& vmodules_in, ostream& oss) {
    ofstream FileMESSAGE;
    PYTHON_MODULES(vmodules_in, FileMESSAGE, oss);
  }

  void PYTHON_MODULES(const vector<string>& vmodules_in, ofstream& FileMESSAGE, ostream& oss) {
    string directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");
    if (directory.empty()) {
      directory = ".";
    }

    vector<string> vmodules;
    vector<string> vskip;
    const vector<string> vavailable{"aflow_cce", "aflow_chull", "aflow_chull_plotter", "aflow_sym", "aflow_xtal_finder"};

    // Get available modules first
    if (!vmodules_in.empty()) {
      for (size_t i = 0; i < vmodules_in.size(); i++) {
        if (aurostd::WithinList(vavailable, vmodules_in[i]) || (vmodules_in[i] == "aflow_xtalfinder")) { // Need to account for inconsistency between aflow command and publication
          vmodules.push_back(vmodules_in[i]);
        } else {
          vskip.push_back(vmodules_in[i]); // For the warning message
        }
      }
    } else {
      vmodules = vavailable;
    }

    stringstream message;
    if (!vskip.empty()) {
      message << "Could not find modules " << aurostd::joinWDelimiter(vskip, ", ") << "."
              << " Available modules: " << aurostd::joinWDelimiter(vavailable, ", ") << ".";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, directory, FileMESSAGE, oss, _LOGGER_WARNING_);
    }
    if (!vmodules.empty()) {
      // Dependencies
      if (aurostd::WithinList(vmodules, "aflow_chull_plotter") && !aurostd::WithinList(vmodules, "aflow_chull")) {
        vmodules.emplace_back("aflow_chull");
      }
      for (size_t i = 0; i < vmodules.size(); i++) {
        const string& mod = vmodules[i];
        const string moddir = aurostd::CleanFileName(directory + "/" + mod);
        if (!aurostd::FileExist(moddir)) {
          aurostd::DirectoryMake(moddir);
        }
        // Only include cpps here to avoid multiple definition errors
        if (mod == "aflow_cce") {
          aurostd::EmbData::save_to_file("aflow_cce_python.py", "SCRIPTS", moddir + "/__init__.py");
        } else if (mod == "aflow_chull") {
          aurostd::EmbData::save_to_file("aflow_chull_python.py", "SCRIPTS", moddir + "/__init__.py");
        } else if (mod == "aflow_chull_plotter") {
          aurostd::EmbData::save_to_file("aflow_chull_jupyter_plotter.py", "SCRIPTS", moddir + "/__init__.py");
        } else if (mod == "aflow_sym") {
          aurostd::EmbData::save_to_file("aflow_sym_python.py", "SCRIPTS", moddir + "/__init__.py");
        } else if ((mod == "aflow_xtal_finder") || (mod == "aflow_xtalfinder")) {
          aurostd::EmbData::save_to_file("aflow_xtalfinder_python.py", "SCRIPTS", moddir + "/__init__.py");
        }
      }
      message << "Successfully installed modules " << aurostd::joinWDelimiter(vmodules, ", ") << ".";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, directory, FileMESSAGE, oss, _LOGGER_NOTICE_);
    } else {
      message << "No modules left to write.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::STATDIEL // CAMILO
// ***************************************************************************
namespace pflow {
  void STATDIEL(vector<string>& argv) { // loop pflow::STATDIEL
    string outcar;
    xvector<double> real(3);
    xvector<double> imag(3);
    if (argv.size() != 3) { // user control lines - aflow specific.
      init::ErrorOption("", "pflow::STATDIEL", "aflow --statdiel OUTCAR*");
    }
    outcar = argv.at(2);
    KBIN::GetStatDiel(outcar, real, imag);
    printf("%11.6f %11.6f %11.6f ", real(1), real(2), real(3));
    printf("%11.6f %11.6f %11.6f \n", imag(1), imag(2), imag(3));
  } // loop pflow::STATDIEL
} // namespace pflow

// ***************************************************************************
// pflow::QMVASP
// ***************************************************************************
// CO20180703 - revamped for vpflow, and fixed INCAR search issue
namespace pflow {
  bool QMVASP(aurostd::xoption& vpflow) {
    // CO20180703
    const bool LDEBUG = (false || XHOST.DEBUG);

    if (vpflow.flag()) {
      ;
    } // keep busy

    // fix directory
    string directory;
    if (XHOST.vflag_control.flag("DIRECTORY")) {
      directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");
    }
    if (directory.empty()) {
      directory = ".";
    }
    if (!aurostd::IsDirectory(directory)) {
      cerr << __AFLOW_FUNC__ << " invalid directory input" << endl;
      return false;
    }

    // organize runs, e.g., .static, .relax2
    const vector<string> vruns{".relax1", ".relax2", ".relax3", ".relax4", ".static", ".bands"};

    // organize CARs, e.g., CONCAR, OUTCAR
    const vector<string> vCARs{"CONTCAR", "OUTCAR", "INCAR", "vasprun.xml"}; // look in aflow_kvasp, KBIN::VASP_Analyze()

    // remove aflow.qmvasp.out
    for (size_t k = 0; k < XHOST.vext.size(); k++) {
      aurostd::RemoveFile(directory + "/" + DEFAULT_AFLOW_QMVASP_OUT + XHOST.vext[k]);
    }

    // go through runs sequentially to fill in aflow.qmvasp.out
    bool found_run = false; // i.e., found COMPLETE run, all files must be .static or all .relax2, no mix match
    string ext = "";
    string vCAR = "";
    string vCAR_compressed = "";
    _xvasp xvasp;
    xvasp.Directory = directory;
    for (size_t j = 0; j < vruns.size(); j++) {
      ext = vruns[j];
      aurostd::StringSubstInPlace(ext, ".", "");
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " checking run=" << aurostd::toupper(ext) << endl;
      }
      found_run = true;
      for (size_t i = 0; i < vCARs.size() && found_run == true; i++) {
        vCAR = directory + "/" + vCARs[i] + vruns[j];
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " looking for " << vCAR << endl;
        }
        if (!aurostd::CompressFileExist(vCAR, vCAR_compressed)) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " " << directory + "/" + vCAR << " NOT found" << endl;
          }
          found_run = false;
          break;
        }
      }
      if (!found_run) {
        continue;
      }
      for (size_t i = 0; i < vCARs.size(); i++) {
        if (aurostd::FileExist(directory + "/" + vCARs[i])) {
          cerr << __AFLOW_FUNC__ << directory + "/" + vCARs[i] << " already exists, cannot overwrite with relax1/relax2/static variant. ";
          cerr << "Do not want to overwrite. Please delete: " << directory + "/" + vCARs[i] << endl;
          return false;
        }
      }
      // now write!
      for (size_t i = 0; i < vCARs.size(); i++) {
        vCAR = directory + "/" + vCARs[i] + vruns[j];
        aurostd::CompressFileExist(vCAR, vCAR_compressed); // get compressed variant
        aurostd::DecompressFile(vCAR_compressed, directory + "/" + vCARs[i], true);
      }
      cout << __AFLOW_FUNC__ << " Performing dir=" << directory << " run=" << aurostd::toupper(ext) << endl;
      xvasp.AnalyzeLabel = ext;
      xvasp.str = xstructure(directory + "/CONTCAR", IOAFLOW_AUTO);
      KBIN::VASP_Analyze(xvasp, true);
      if (!aurostd::FileExist(directory + "/" + DEFAULT_AFLOW_QMVASP_OUT)) {
        cerr << __AFLOW_FUNC__ << directory + "/" + DEFAULT_AFLOW_QMVASP_OUT << " was not successfully generated" << endl;
        return false;
      }
      for (size_t i = 0; i < vCARs.size(); i++) {
        aurostd::RemoveFile(directory + "/" + vCARs[i]);
      }
    }

    // zip in style of vCAR
    for (size_t i = 0; i < vCARs.size(); i++) {
      for (size_t j = 0; j < vruns.size(); j++) {
        for (size_t k = 0; k < XHOST.vext.size(); k++) {
          if (aurostd::FileExist(directory + "/" + vCARs[i] + vruns[j] + XHOST.vext[k])) {
            if (!XHOST.vzip.at(k).empty()) {
              aurostd::execute(XHOST.command(XHOST.vzip.at(k)) + " -9qf " + xvasp.Directory + "/" + DEFAULT_AFLOW_QMVASP_OUT);
              return true;
            }
          }
        }
      }
    }
    return true;
  }
} // namespace pflow

// ***************************************************************************
// pflow::RAYTRACE
// ***************************************************************************
namespace pflow {
  void RAYTRACE(vector<string> argv) {
    // cout << aflow::Banner("BANNER_TINY") << endl;
    pflow::RayTraceManager(argv); // Manage ray tracing related activities.
  }
} // namespace pflow

// ***************************************************************************
// pflow::RASMOL
// ***************************************************************************
namespace pflow {
  void RASMOL(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    const xvector<int> ijk(3); // default
    ijk[1] = 1;
    ijk[2] = 1;
    ijk[3] = 1; // default
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");

    if (tokens.size() > 3) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --rasmol[=n1[,n2[,n3]]] < POSCAR");
    }

    if (!tokens.empty()) {
      ijk[1] = aurostd::string2utype<int>(tokens[0]);
    }
    if (tokens.size() >= 2) {
      ijk[2] = aurostd::string2utype<int>(tokens[1]);
    }
    if (tokens.size() >= 3) {
      ijk[3] = aurostd::string2utype<int>(tokens[2]);
    }

    ofstream FileOUTPUT;
    const string FileOUTPUTName = "aflow.rasmol.xyz." + XHOST.ostrPID.str() + "." + XHOST.ostrTID.str(); // CO20200502 - threadID
    FileOUTPUT.open(FileOUTPUTName.c_str(), std::ios::out);
    if (FileOUTPUT) {
      xvector<int> _ijk(vabs(ijk));
      if (max(_ijk) == 0) {
        _ijk.set(1);
      }
      ostringstream aus_exec;
      const xstructure a(input, IOAFLOW_AUTO);
      PrintXYZ(a, _ijk, FileOUTPUT);
      FileOUTPUT.close();
      //   cerr << _ijk << endl;
      aus_exec << XHOST.command("rasmol") << " -xyz " << FileOUTPUTName << endl; // " && rm -f " << FileOUTPUTName << endl;
      aurostd::execute(aus_exec);
      // aurostd::StringstreamClean(aus_exec);
      aurostd::RemoveFile(FileOUTPUTName);
    } else {
      FileOUTPUT.close();
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::RSM
// ***************************************************************************
namespace pflow {
  void RSM(vector<string> argv, istream& input) {
    int i;
    int j;
    int Ntypes = 0;
    string strtmp = "";
    string AtomSymbol = "";
    bool Z_flag = false;
    int Z[20]; // to store the z values following --z option

    // FINDING -Z
    int iargv = 0;
    Ntypes = 0;
    for (i = 1; i < (int) argv.size(); i++) {
      strtmp = argv[i];
      if (strtmp == "--Z" or strtmp == "--z") {
        Z_flag = true;
        iargv = i + 1;
        break;
      }
    }
    Ntypes = (int) argv.size() - iargv;
    cerr << "Ntypes = " << Ntypes << endl;
    for (j = 0; j < Ntypes; j++) {
      AtomSymbol = argv.at(j + iargv);
      cerr << AtomSymbol << "(";
      if (AtomSymbol[0] < '0' || AtomSymbol[0] > '9') {
        Z[j] = (int) GetAtomNumber(AtomSymbol);
      } else {
        Z[j] = aurostd::string2utype<int>(AtomSymbol); // index in Z starts from ZERO
      }
      cerr << Z[j] << ")  ";
    }
    cerr << endl;

    xstructure a(input, IOAFLOW_AUTO);

    if (Ntypes < 1) {
      Z_flag = false;
      Ntypes = a.num_each_type.size();
    }
    if ((int) a.num_each_type.size() != Ntypes) {
      cerr << "STOP: inconsistency in number of types in POSCAR and --z" << endl << "a.num_each_type.size() = " << a.num_each_type.size() << endl << "Ntypes = " << Ntypes << endl;
      abort();
    }
    if (!Z_flag) {
      for (i = 0; i < Ntypes; i++) {
        Z[i] = i + 1;
      }
    }

    // fixing atoms.type
    int k;
    int ptr = 0;
    for (i = 0; i < Ntypes; i++) {
      for (j = 0; j < a.num_each_type.at(i); j++) {
        ptr = 0;
        for (k = 0; k < i + 1; k++) {
          ptr = ptr + a.num_each_type.at(k);
        }
        ptr = ptr - a.num_each_type.at(i) + j; //        ptr=sum(species(1:i))-species(i)+j;
        a.atoms.at(ptr).type = Z[i];
      }
    }

    //[CO20220623 - not used]int num_atoms=0;
    //[CO20220623 - not used]for(i=0;i<(int) a.num_each_type.size();i++)
    //[CO20220623 - not used]  num_atoms+=a.num_each_type.at(i);

    PrintRSM(a, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::RBANAL
// ***************************************************************************
namespace pflow {
  void RBANAL(vector<string> argv) {
    // cout << aflow::Banner("BANNER_TINY") << endl;
    // Read in input files.
    const int nim = atoi(argv.at(2).c_str());
    const string path_flag(argv.at(3));
    PrintRBAnal(nim, path_flag, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::RBDIST
// ***************************************************************************
namespace pflow {
  void RBDIST(vector<string> argv) {
    xstructure strA;
    xstructure strB;
    ifstream infileA(argv.at(2).c_str());
    ifstream infileB(argv.at(3).c_str());
    infileA >> strA;
    infileB >> strB;
    const string path_flag(argv.at(4).c_str());
    double totdist;
    aurostd::matrix<double> cm(2, 3, 999999); // CO20200404 pflow::matrix()->aurostd::matrix()
    xstructure diffstr;
    pflow::RBPoscarDisp(strA, strB, diffstr, totdist, cm, path_flag);
    PrintRBPoscarDisp(diffstr, totdist, cm, path_flag, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::RDF
// ***************************************************************************
namespace pflow {
  void RDF(const string& options, istream& input, bool raw_counts) { // CO220627 - rewritten
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    } // CO2020062
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() > 4) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --rdf[=rmax[,nbins[,sigma[,window_gaussian]]]] [--raw_counts] < POSCAR");
    }

    xstructure a(input, IOAFLOW_AUTO);
    if (a.species.empty()) {
      a.species = aurostd::vector2deque(pflow::getFakeElements(a.num_each_type.size()));
    }
    double rmax = 5.0;
    int nbins = 25;
    double sigma = 0;
    int window_gaussian = 0;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
    }
    if (!tokens.empty()) {
      rmax = aurostd::string2utype<double>(tokens[0]);
    }
    if (tokens.size() >= 2) {
      nbins = aurostd::string2utype<int>(tokens[1]);
    }
    if (tokens.size() >= 3) {
      sigma = aurostd::string2utype<double>(tokens[2]);
    }
    if (tokens.size() >= 4) {
      window_gaussian = aurostd::string2utype<int>(tokens[3]);
    }

    aurostd::xmatrix<double> rdf_all;
    pflow::GetRDF(a, rdf_all, rmax, nbins, raw_counts, sigma, window_gaussian);
    // to be implemented/patched later: GetRDFShells()
    // this function is supposed to count the number of wiggles in the rdf above/below horizontal asymptote (constant)
    cout << aflow::Banner("BANNER_TINY") << endl;
    PrintRDF(a, rmax, nbins, rdf_all, cout);

    // stringstream gout;
    // xoption plotoptions;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::RDFCMP
// ***************************************************************************
namespace pflow {
  void RDFCMP(const string& options) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 6) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --rdfcmp=rmax,nbins,sigma,nshmax,POSCAR1,POSCAR2");
    }
    const double rmax = (double) aurostd::string2utype<double>(tokens.at(0));
    const int nbins = (int) aurostd::string2utype<int>(tokens.at(1));
    const int smooth_width = (int) aurostd::string2utype<int>(tokens.at(2));
    const int nsh = (int) aurostd::string2utype<int>(tokens.at(3));

    if (!aurostd::FileExist(tokens.at(4))) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "file not found: " + tokens.at(4), _FILE_CORRUPT_); // CO20200624
    }
    if (!aurostd::FileExist(tokens.at(5))) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "file not found: " + tokens.at(5), _FILE_CORRUPT_); // CO20200624
    }
    const xstructure strA(tokens.at(4), IOAFLOW_AUTO);
    const xstructure strB(tokens.at(5), IOAFLOW_AUTO);
    // Get rdfs
    const aurostd::matrix<double> rdf_all_A; // CO20200404 pflow::matrix()->aurostd::matrix()
    const aurostd::matrix<double> rdf_all_B; // CO20200404 pflow::matrix()->aurostd::matrix()
    const aurostd::matrix<double> rdf_all_A_sm = pflow::GetSmoothRDF(rdf_all_A, smooth_width); // CO20200404 pflow::matrix()->aurostd::matrix()
    const aurostd::matrix<double> rdf_all_B_sm = pflow::GetSmoothRDF(rdf_all_B, smooth_width); // CO20200404 pflow::matrix()->aurostd::matrix()
    // Get shells
    aurostd::matrix<double> rdfsh_all_A; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> rdfsh_loc_A; // Radial location of rdf shells. //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> rdfsh_all_B; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> rdfsh_loc_B; // Radial location of rdf shells. //CO20200404 pflow::matrix()->aurostd::matrix()
    pflow::GetRDFShells(strA, rmax, nbins, smooth_width, rdf_all_A_sm, rdfsh_all_A, rdfsh_loc_A);
    pflow::GetRDFShells(strB, rmax, nbins, smooth_width, rdf_all_B_sm, rdfsh_all_B, rdfsh_loc_B);
    vector<int> best_match;
    aurostd::matrix<double> rms_mat; // CO20200404 pflow::matrix()->aurostd::matrix()
    pflow::CmpRDFShells(strA, strB, rdfsh_all_A, rdfsh_all_B, nsh, best_match, rms_mat);
    PrintRDFCmp(strA, strB, rmax, nbins, smooth_width, nsh, rdfsh_all_A, rdfsh_all_B, best_match, rms_mat, cout);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::RMATOM
// ***************************************************************************
namespace pflow {
  xstructure RMATOM(istream& input, const int& iatom) {
    xstructure a(input, IOAFLOW_AUTO);
    a.RemoveAtom(iatom);
    return a;
  }
} // namespace pflow

// ***************************************************************************
// pflow::RMCOPIES
// ***************************************************************************
namespace pflow {
  xstructure RMCOPIES(istream& input) {
    xstructure a(input, IOAFLOW_AUTO);
    a.RemoveCopies();
    return a;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SCALE
// ***************************************************************************
namespace pflow {
  xstructure SCALE(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 1) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --scale=s < POSCAR");
    }
    // move on
    const xstructure a(input, IOAFLOW_AUTO);
    xstructure b;
    b = a;
    double scale = 0.0; // some defaults
    if (!tokens.empty()) {
      scale = aurostd::string2utype<double>(tokens[0]);
    }

    if (LDEBUG) {
      cerr << "DEBUG  b.scale=" << b.scale << endl;
    }
    b = ReScale(b, scale);
    if (LDEBUG) {
      cerr << "DEBUG  b.scale=" << b.scale << endl;
    }
    b.neg_scale = false;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::INFLATE_LATTICE
// ***************************************************************************
namespace pflow {
  xstructure INFLATE_LATTICE(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 1) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --inflate_lattice=coefficient | --ilattice=coefficient < POSCAR");
    }
    // move on
    const xstructure a(input, IOAFLOW_AUTO);
    xstructure b;
    b = a;
    double coefficient = 0.0; // some defaults
    if (!tokens.empty()) {
      coefficient = aurostd::string2utype<double>(tokens[0]);
    }
    if (LDEBUG) {
      cerr << "DEBUG  b.scale=" << b.scale << endl;
    }
    b = InflateLattice(b, coefficient);
    if (LDEBUG) {
      cerr << "DEBUG  b.scale=" << b.scale << endl;
    }
    b.neg_scale = false;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::INFLATE_VOLUME
// ***************************************************************************
namespace pflow {
  xstructure INFLATE_VOLUME(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 1) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --inflate_volume=coefficient | --ivolume=coefficient < POSCAR");
    }
    // move on
    const xstructure a(input, IOAFLOW_AUTO);
    xstructure b;
    b = a;
    double coefficient = 0.0; // some defaults
    if (!tokens.empty()) {
      coefficient = aurostd::string2utype<double>(tokens[0]);
    }
    if (LDEBUG) {
      cerr << "DEBUG  b.scale=" << b.scale << endl;
    }
    b = InflateVolume(b, coefficient);
    if (LDEBUG) {
      cerr << "DEBUG  b.scale=" << b.scale << endl;
    }
    b.neg_scale = false;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SD
// ***************************************************************************
namespace pflow {
  xstructure SD(vector<string> argv, istream& input) {
    const xstructure a(input, IOAFLOW_AUTO);
    xstructure b;
    b = a;
    // Read in input file.
    // Set a sd vector
    const int num_types = a.num_each_type.size();
    vector<string> sd(num_types);
    for (int i = 0; i < num_types; i++) {
      sd[i] = "TTT";
    }
    //  vector<string> sd(argv.size()-2);
    for (int i = 0; i < min(((int) argv.size() - 2), num_types); i++) {
      sd[i] = string(argv.at(i + 2));
    }
    b.isd = true;
    b = SetSDTypes(b, sd);
    return b;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SETCM
// ***************************************************************************
namespace pflow {
  xstructure SETCM(istream& input, const xvector<double>& newcm) {
    xstructure str(input, IOAFLOW_AUTO);
    const xvector<double> oldcm(GetMom1(str));
    str = ShiftCPos(str, (newcm - oldcm) / str.scale);
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SETORIGIN
// ***************************************************************************
namespace pflow {
  xstructure SETORIGIN(istream& input, const xvector<double>& neworigin) {
    xstructure str(input, IOAFLOW_AUTO);
    // xvector<double> oldorigin(GetMom1(str));
    xvector<double> oldorigin(3);
    oldorigin.clear();
    //  str=ShiftCPos(str,(neworigin-oldorigin)/str.scale);
    if (str.coord_flag == _COORDS_CARTESIAN_) {
      str = ShiftCPos(str, -(neworigin - oldorigin));
    }
    if (str.coord_flag == _COORDS_FRACTIONAL_) {
      str = ShiftFPos(str, -(neworigin - oldorigin));
    }
    //    str=pflow::SetOrigin(str,neworigin);
    return str;
  }
} // namespace pflow

namespace pflow {
  xstructure SETORIGIN(istream& input, const int& natom) {
    xstructure str(input, IOAFLOW_AUTO);
    if (natom < 0) {
      return str;
    }
    if (natom < (int) str.atoms.size()) {
      if (str.coord_flag == _COORDS_CARTESIAN_) {
        str = ShiftCPos(str, -str.atoms.at((uint) natom).cpos);
      }
      if (str.coord_flag == _COORDS_FRACTIONAL_) {
        str = ShiftFPos(str, -str.atoms.at((uint) natom).fpos);
      }
    }
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SEWALD
// ***************************************************************************
namespace pflow {
  void SEWALD(vector<string> argv, istream& input) {
    // cout << aflow::Banner("BANNER_TINY") << endl;
    double Ks = atof(argv.at(2).c_str());
    const double SUMTOL = 1.0e-16;
    xstructure str(input, IOAFLOW_AUTO);
    str = GetNiggliStr(str);
    double epoint = 0.0;
    double ereal = 0.0;
    double erecip = 0.0;
    double eewald = 0.0;
    ereal = pflow::ScreenedESEner(str, Ks, SUMTOL);
    eewald = ereal;
    pflow::PrintEwald(str, epoint, ereal, erecip, eewald, Ks, SUMTOL, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::SG
// ***************************************************************************
namespace pflow {
  string SG(aurostd::xoption& vpflow, istream& input, string mode, string print) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    string flag_name = "SG::" + mode; // DX20170926
    stringstream message; // DX20200103
    if (print == "LABEL" || print == "NUMBER") {
      flag_name += "_" + print;
    }
    const string options = vpflow.getattachedscheme(flag_name); // DX20170926
    if (LDEBUG) {
      cerr << XPID << "pflow::SG: mode=" << mode << endl;
    }

    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() == 1) {
      if (tokens[0] == "usage" || tokens[0] == "USAGE") {
        init::MessageOption(options, __AFLOW_FUNC__,
                            aurostd::liststring2string("aflow --aflowSG[=<tolerance_value>|=tight|=loose] [--mag|--magnetic|--magmom=[m1,m2,...|INCAR|OUTCAR]] < POSCAR  default: "
                                                       "(minimum_interatomic_distance)/100.0",
                                                       "aflow --platonSG[_label,_number][=EQUAL| EXACT][,ang,d1,d2,d3] < POSCAR  default:" + string("EQUAL=") + aurostd::utype2string<int>(DEFAULT_PLATON_P_EQUAL) +
                                                           "," + string("EXACT=") + aurostd::utype2string<int>(DEFAULT_PLATON_P_EXACT) + "," + aurostd::utype2string(DEFAULT_PLATON_P_ANG, 5) + "," +
                                                           aurostd::utype2string(DEFAULT_PLATON_P_D1, 5) + "," + aurostd::utype2string(DEFAULT_PLATON_P_D2, 5) + "," + aurostd::utype2string(DEFAULT_PLATON_P_D3, 5) + "",
                                                       "aflow --findsymSG[_label,_number][=tolerance] < POSCAR   default:" + aurostd::utype2string(DEFAULT_FINDSYM_TOL, 5)));
        return ""; // CO20200624 - the option was expressed successfully
      }
    }
    // move on

    if (LDEBUG) {
      cerr << XPID << "pflow::SG: tokens.size()=" << tokens.size() << endl;
    }

    // check usage
    //   if(LDEBUG) cerr << XPID << "pflow::SG: vpflow.getattachedscheme(\"SG::USAGE\")=" << vpflow.flag("SG::USAGE") << endl;

    if (input.peek() == EOF) {
      message << "File is empty. Check POSCAR.";
      throw aurostd::xerror(__AFLOW_FILE__, flag_name, message, _INPUT_MISSING_);
    }
    xstructure a(input, IOAFLOW_AUTO);
    // DX20180527 - use pwd - START
    if (a.directory.empty()) {
      a.directory = aurostd::getPWD();
    }
    // DX20180527 - use pwd - END
    //    cerr << a << endl;
    //  AFLOW ENGINE RHT
    if (mode == "AFLOW" || mode == "aflow") { // RHT
      if (LDEBUG) {
        cerr << XPID << "pflow::SG: aflow" << endl;
      }
      // DX START
      a.ReScale(1.0);
      // DX20170921 - MAGNETIC SYMMETRY - START
      if (vpflow.flag("SG::MAGNETIC")) {
        const string magmom_info = vpflow.getattachedscheme("SG::MAGNETIC");
        ProcessAndAddSpinToXstructure(a, magmom_info); // DX20191108 - condensed into a single function
      }
      // DX20170921 - MAGNETIC SYMMETRY - END
      // DX END

      // get tolerance
      double tolerance = pflow::getSymmetryTolerance(a, vpflow.getattachedscheme("SG::TOLERANCE")); // DX20200820 - consolidated setting tolerance into a function

      bool tolerance_spectrum_analysis = false;
      vector<double> tolerance_spectrum;
      // DX20200817 - SPACEGROUP SPECTRUM - START
      if (vpflow.flag("SG::TOLERANCE_SPECTRUM")) {
        tolerance_spectrum_analysis = true;
        tolerance_spectrum = pflow::getSymmetryToleranceSpectrum(vpflow.getattachedscheme("SG::TOLERANCE_SPECTRUM"));
      } else if (vpflow.flag("SG::TOLERANCE") && vpflow.flag("SG::TOLERANCE_SPECTRUM")) {
        message << "pflow::SG::ERROR: Cannot specify a single tolerance value and perform the tolerance spectrum at the same time. Please choose one or the other.";
        throw aurostd::xerror(__AFLOW_FILE__, flag_name, message, _INPUT_ILLEGAL_);
      }
      // DX20200817 - SPACEGROUP SPECTRUM - END
      // DX20170926 - NO SCAN - START
      if (vpflow.flag("SG::NO_SCAN")) {
        a.sym_eps_no_scan = true; // DX20210406
      }
      // DX20170926 - NO SCAN - END
      if (!tolerance_spectrum_analysis) {
        const uint sgroup = a.SpaceGroup_ITC(tolerance, a.sym_eps_no_scan);
        a.spacegroup = GetSpaceGroupName(sgroup, a.directory) + " #" + aurostd::utype2string(sgroup); // DX20190319 - put directory name
      } else {
        // perform space group analysis through a range of tolerances //DX20200820
        for (size_t i = 0; i < tolerance_spectrum.size(); i++) {
          const uint sgroup = a.SpaceGroup_ITC(tolerance_spectrum[i], a.sym_eps_no_scan);
          try {
            GetSpaceGroupName(sgroup, a.directory);
            cout << "tol=" << tolerance_spectrum[i] << ": " << GetSpaceGroupName(sgroup, a.directory) << " #" << aurostd::utype2string(sgroup) << endl;
          } catch (aurostd::xerror& excpt) {
            cout << "tol=" << tolerance_spectrum[i] << ": -" << endl;
          }
        }
        return "";
      }
      //  return a.spacegroup; //RHT
    }

    if (mode == "PLATON" || mode == "platon") {
      if (LDEBUG) {
        cerr << XPID << "pflow::SG: platon" << endl;
      }
      // SIMPLE CALCULATION
      // output << a.platon2sg() << endl;
      // PERFECT CALCULATION
      bool Platon_EQUAL = DEFAULT_PLATON_P_EQUAL;
      bool Platon_EXACT = DEFAULT_PLATON_P_EXACT;
      double Platon_ang = DEFAULT_PLATON_P_ANG;
      double Platon_d1 = DEFAULT_PLATON_P_D1;
      double Platon_d2 = DEFAULT_PLATON_P_D2;
      double Platon_d3 = DEFAULT_PLATON_P_D3;
      // DX20170926 - FLAGS - START
      if (vpflow.flag("SG::TOLERANCE")) {
        const string tolerance_string = vpflow.getattachedscheme("SG::TOLERANCE");
        vector<string> tol_tokens;
        aurostd::string2tokens(tolerance_string, tol_tokens, ",");
        //     if(tokens.size()==0) Platon_ang=1.0e-2;
        // if(tokens.size()==1) Platon_ang=aurostd::string2utype<double>(tokens.at(0));
        // Read in input file.
        if ((!tol_tokens.empty() && tol_tokens[0] == "EQUAL") || (tol_tokens.size() >= 2 && tol_tokens[1] == "EQUAL")) {
          Platon_EQUAL = true;
        }
        if ((!tol_tokens.empty() && tol_tokens[0] == "EXACT") || (tol_tokens.size() >= 2 && tol_tokens[1] == "EXACT")) {
          Platon_EXACT = true;
        }
        if (tol_tokens.size() >= 3) {
          Platon_ang = aurostd::string2utype<double>(tol_tokens.at(tol_tokens.size() - 4));
          Platon_d1 = aurostd::string2utype<double>(tol_tokens.at(tol_tokens.size() - 3));
          Platon_d2 = aurostd::string2utype<double>(tol_tokens.at(tol_tokens.size() - 2));
          Platon_d3 = aurostd::string2utype<double>(tol_tokens.at(tol_tokens.size() - 1));
        }
      }
      // DX20170926 - FLAGS - END
      a.platon2sg(Platon_EQUAL, Platon_EXACT, Platon_ang, Platon_d1, Platon_d2, Platon_d3);
      //  return a.spacegroup;
    }

    if (mode == "FINDSYM" || mode == "findsym") {
      if (LDEBUG) {
        cerr << XPID << "pflow::SG: findsym" << endl;
      }
      // SIMPLE CALCULATION
      double tolerance = DEFAULT_FINDSYM_TOL;
      if (vpflow.flag("SG::TOLERANCE")) {
        const string tolerance_string = vpflow.getattachedscheme("SG::TOLERANCE");
        vector<string> tol_tokens;
        aurostd::string2tokens(tolerance_string, tol_tokens, ",");
        if (tol_tokens.empty()) {
          tolerance = DEFAULT_FINDSYM_TOL;
        }
        if (tol_tokens.size() == 1) {
          tolerance = aurostd::string2utype<double>(tol_tokens[0]);
        }
      }
      a.findsym2sg(tolerance);
      //   return a.spacegroup;
    }
    // what to print

    //   vector<string> tokens;
    if (print == "LABEL" || print == "label") {
      aurostd::string2tokens(a.spacegroup, tokens, " ");
      return tokens[0];
    }
    if (print == "NUMBER" || print == "number") {
      aurostd::string2tokens(a.spacegroup, tokens, " ");
      return tokens[1];
    }
    return a.spacegroup;

    //  return NOSG;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SHELL
// ***************************************************************************
namespace pflow {
  void SHELL(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 5) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --shell=ns,r1,r2,name,dens < POSCAR");
    }

    int ns = 4; // some defaults
    double r1 = 1.8;
    double r2 = 2.2; // some defaults
    string name; // some defaults
    const int dens = 20; // some defaults
    if (!tokens.empty()) {
      ns = aurostd::string2utype<int>(tokens[0]);
    }
    if (tokens.size() >= 2) {
      r1 = aurostd::string2utype<double>(tokens[1]);
    }
    if (tokens.size() >= 3) {
      r2 = aurostd::string2utype<double>(tokens[2]);
    }
    if (tokens.size() >= 4) {
      name = tokens[3];
    }
    if (tokens.size() >= 5) {
      ns = aurostd::string2utype<int>(tokens[4]);
    }

    cout << aflow::Banner("BANNER_TINY") << endl;
    const xstructure a(input, IOAFLOW_AUTO);

    PrintShell(a, ns, r1, r2, name, dens, cout);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::SHIFT
// ***************************************************************************
namespace pflow {
  xstructure SHIFT(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 3 && tokens.size() != 4) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --shift=Sx,Sy,Sz[,cCdD] < POSCAR");
    }
    // cout << aflow::Banner("BANNER_TINY") << endl;
    // Get shift vector
    const xvector<double> shift(3);
    for (int ic = 1; ic <= 3; ic++) {
      shift[ic] = aurostd::string2utype<double>(tokens.at(ic - 1));
    }
    // argv(4)=c cart,d for direct shift
    bool flag = false;
    if (tokens.size() == 3) {
      cerr << "WARNING - pflow::SHIFT: 4th argument can be [cCfFdD], defaulting to [cC]=Cartesian shift" << endl;
    } else {
      if (tokens.at(3) == "c" || tokens.at(3) == "C") {
        flag = false;
      } else {
        if (tokens.at(3) == "d" || tokens.at(3) == "D" || tokens.at(3) == "f" || tokens.at(3) == "F") {
          flag = true;
        } else {
          cerr << "WARNING - pflow::SHIFT: 4th argument can be [cCfFdD], Defaulting to [cC]=Cartesian shift" << endl;
          flag = false;
        } // else
      } // else
    } // else
    // Read in input file.
    xstructure str(input, IOAFLOW_AUTO);
    str = ShiftPos(str, shift, flag);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SG
// ***************************************************************************
namespace pflow {
  void SG(istream& input) {
    xstructure a(input, IOAFLOW_AUTO);
    a.CalculateSymmetryFactorGroup(true);
    spacegroup::SpaceGroupInitialize();
    cerr << a.pgroup.size() << "," << a.fgroup.size() << endl;
    cerr.flush();
    xmatrix<double> I(3, 3);
    xmatrix<double> U(3, 3);
    xmatrix<double> A(3, 3);
    xvector<double> tau(3);
    I = identity(I);
    for (size_t i = 0; i < a.fgroup.size(); i++) {
      U = a.fgroup[i].Uf;
      tau = a.fgroup[i].ftau;
      A = I - U;
      if (std::abs(det(A)) > 0.01) {
        cerr << inverse(A) * tau << endl;
      }
    }
    //  spacegroup::SpaceGroupNumberStructure(a);
  }
} // namespace pflow

// ***************************************************************************
// pflow::SGROUP
// ***************************************************************************
namespace pflow {
  void SGROUP(_aflags& aflags, istream& input, double radius) {
    cout << aflow::Banner("BANNER_TINY") << endl;
    aflags.QUIET = true;
    xstructure a(input, IOAFLOW_AUTO);
    const bool WRITE = true;
    ofstream File("/dev/null");
    SYM::CalculatePointGroup(File, a, aflags, WRITE, true, cout);
    SYM::CalculatePointGroupKLattice(File, a, aflags, WRITE, true, cout);
    SYM::CalculateFactorGroup(File, a, aflags, WRITE, true, cout);
    SYM::CalculatePointGroupCrystal(File, a, aflags, WRITE, true, cout);
    a.sgroup_radius = radius;
    SYM::CalculateSpaceGroup(File, a, aflags, WRITE, true, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::SPECIES
// ***************************************************************************
namespace pflow {
  string SPECIES(istream& input) {
    xstructure a(input, IOAFLOW_AUTO);
    string strout = a.SpeciesString() + "\n";
    return strout;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SPLINE
// ***************************************************************************
namespace pflow {
  void SPLINE(vector<string> argv) {
    // Read in input data.
    vector<double> x;
    vector<double> y;
    double p;
    double q;
    char dum[500];
    while (cin >> p >> q) {
      x.push_back(p);
      y.push_back(q);
      cin.getline(dum, 500); // Get anything else on line.
    }
    const int npts = atoi(argv.at(2).c_str());
    pflow::PrintSpline(x, y, npts, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::SUMPDOS
// ***************************************************************************
namespace pflow {
  void SUMPDOS(vector<string> argv) {
    cerr << "# WARNING: THIS REQUIRES THAT YOU HAVE PDOS IN DOSCAR FILE - THIS IS OBTAINED BY RUNNING WITH LORBIT=2 - SEE aflow --h" << endl;
    ifstream SumPDOSParams_infile(argv.at(2).c_str());
    aurostd::InFileExistCheck("convasp.cc", argv.at(2), SumPDOSParams_infile);
    ifstream PDOS_infile(argv.at(3).c_str());
    aurostd::InFileExistCheck("convasp.cc", argv.at(3), PDOS_infile);
    pflow::pdosdata pdd;
    aurostd::matrix<aurostd::matrix<double>> allpdos; // CO20200404 pflow::matrix()->aurostd::matrix()
    pflow::ReadSumDOSParams(SumPDOSParams_infile, pdd);
    //      projdata prd;
    //      pdd.PrintParams(cout,prd.LMnames);
    pflow::ReadInPDOSData(allpdos, pdd, PDOS_infile);
    SumPDOS(allpdos, pdd);
  }
} // namespace pflow

// ***************************************************************************
// pflow::SUPERCELL
// ***************************************************************************
namespace pflow {
  xstructure SUPERCELL(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    const xstructure str(input, IOAFLOW_AUTO);
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    const xmatrix<double> msc(3, 3);

    if (tokens.size() != 9 && tokens.size() != 3 && tokens.size() != 1) {
      init::ErrorOption(options, __AFLOW_FUNC__,
                        aurostd::liststring2string("aflow --supercell=a11,a12,a13,a21,a22,a23,a31,a32,a33 < POSCAR", "aflow --supercell=a11,a22,a33 < POSCAR", "aflow --supercell=file < POSCAR"));
    }

    if (tokens.size() == 9) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " 9 entries" << endl;
      }
      msc(1, 1) = aurostd::string2utype<double>(tokens[0]);
      msc(1, 2) = aurostd::string2utype<double>(tokens[1]);
      msc(1, 3) = aurostd::string2utype<double>(tokens[2]);
      msc(2, 1) = aurostd::string2utype<double>(tokens[3]);
      msc(2, 2) = aurostd::string2utype<double>(tokens[4]);
      msc(2, 3) = aurostd::string2utype<double>(tokens[5]);
      msc(3, 1) = aurostd::string2utype<double>(tokens[6]);
      msc(3, 2) = aurostd::string2utype<double>(tokens[7]);
      msc(3, 3) = aurostd::string2utype<double>(tokens[8]);
      if (std::abs(det(msc)) < 0.01) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "singular supercell matrix", _INPUT_ILLEGAL_); // CO20200624
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " END" << endl;
      }
      return GetSuperCell(str, msc);
    }

    if (tokens.size() == 3) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " 3 entries" << endl;
      }
      msc(1, 1) = aurostd::string2utype<double>(tokens[0]);
      msc(2, 2) = aurostd::string2utype<double>(tokens[1]);
      msc(3, 3) = aurostd::string2utype<double>(tokens[2]);
      if (std::abs(det(msc)) < 0.01) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "singular supercell matrix", _INPUT_ILLEGAL_); // CO20200624
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " END" << endl;
      }
      return GetSuperCell(str, msc);
    }

    if (tokens.size() == 1) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " 1 entries" << endl;
      }
      ifstream infile(tokens[0].c_str());
      aurostd::InFileExistCheck(__AFLOW_FUNC__, tokens[0].c_str(), infile);
      for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
          infile >> msc(i, j);
        }
      }
      if (std::abs(det(msc)) < 0.01) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "singular supercell matrix", _INPUT_ILLEGAL_); // CO20200624
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " END" << endl;
      }
      return GetSuperCell(str, msc);
    }
    return str;
  }
} // namespace pflow

// ***************************************************************************
// pflow::SUPERCELLSTRLIST
// ***************************************************************************
namespace pflow {
  void SUPERCELLSTRLIST(const string& options) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    const xmatrix<double> msc(3, 3);

    string infile_name;

    if (tokens.size() != 10 && tokens.size() != 4 && tokens.size() != 2) {
      init::ErrorOption(
          options, __AFLOW_FUNC__,
          aurostd::liststring2string("aflow --supercell_strlist=a11,a12,a13,a21,a22,a23,a31,a32,a33,strlist", "aflow --supercell_strlist=a11,a22,a33,strlist", "aflow --supercell_strlist=file,strlist"));
    }

    if (tokens.size() == 10) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " 10 entries" << endl;
      }
      msc(1, 1) = aurostd::string2utype<double>(tokens[0]);
      msc(1, 2) = aurostd::string2utype<double>(tokens[1]);
      msc(1, 3) = aurostd::string2utype<double>(tokens[2]);
      msc(2, 1) = aurostd::string2utype<double>(tokens[3]);
      msc(2, 2) = aurostd::string2utype<double>(tokens[4]);
      msc(2, 3) = aurostd::string2utype<double>(tokens[5]);
      msc(3, 1) = aurostd::string2utype<double>(tokens[6]);
      msc(3, 2) = aurostd::string2utype<double>(tokens[7]);
      msc(3, 3) = aurostd::string2utype<double>(tokens[8]);
      infile_name = tokens[9];
      if (std::abs(det(msc)) < 0.01) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "singular supercell matrix", _INPUT_ILLEGAL_); // CO20200624
      }
    }
    if (tokens.size() == 4) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " 4 entries" << endl;
      }
      msc(1, 1) = aurostd::string2utype<double>(tokens[0]);
      msc(2, 2) = aurostd::string2utype<double>(tokens[1]);
      msc(3, 3) = aurostd::string2utype<double>(tokens[2]);
      infile_name = tokens[3];
      if (std::abs(det(msc)) < 0.01) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "singular supercell matrix", _INPUT_ILLEGAL_); // CO20200624
      }
    }
    if (tokens.size() == 2) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " 2 entries" << endl;
      }
      ifstream infile(tokens[0].c_str());
      aurostd::InFileExistCheck(__AFLOW_FUNC__, tokens[0].c_str(), infile);
      for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
          infile >> msc(i, j);
        }
      }
      infile_name = tokens[1];
      if (std::abs(det(msc)) < 0.01) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "singular supercell matrix", _INPUT_ILLEGAL_); // CO20200624
      }
    }
    //    ifstream infile(infile_name.c_str());
    // aurostd::InFileExistCheck("aflow",infile_name,infile,cerr);
    aurostd::matrix<double> mmsc(3, 3); // CO20200404 pflow::matrix()->aurostd::matrix()
    // if(LDEBUG) cerr << __AFLOW_FUNC__ << " msc=" << msc << endl;
    mmsc = aurostd::xmatrix2matrix(msc); // CO20200404 pflow::matrix()->aurostd::matrix()
    //  pflow::Mout(mmsc,cout);
    ifstream list_inf(infile_name.c_str());
    aurostd::InFileExistCheck(__AFLOW_FUNC__, infile_name, list_inf);
    vector<xstructure> vstr;
    vector<xstructure> vstr_sc;
    pflow::ReadInStrVec(vstr, list_inf);
    for (size_t i = 0; i < vstr.size(); i++) {
      vstr_sc.push_back(vstr[i]);
    }
    pflow::SuperCellStrVec(vstr_sc, mmsc);
    pflow::PrintStrVec(vstr_sc, cout);
    cerr << __AFLOW_FUNC__ << " vstr_sc.size()=" << vstr_sc.size() << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::xstrSWAP
// ***************************************************************************
namespace pflow {
  xstructure xstrSWAP(vector<string> argv, istream& input) {
    xstructure a(input, IOAFLOW_AUTO);
    const int speciesA = aurostd::string2utype<int>(argv.at(2));
    const int speciesB = aurostd::string2utype<int>(argv.at(3));
    if (speciesA < 0) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "speciesA<0 (speciesA=" + aurostd::utype2string(speciesA) + ")", _INPUT_ILLEGAL_); // CO20200624
    }
    if (speciesA >= (int) a.num_each_type.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "speciesA>=num_each_type.size() (speciesA=" + aurostd::utype2string(speciesA) + ")", _INPUT_ILLEGAL_); // CO20200624
    }
    if (speciesB < 0) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "speciesB<0 (speciesB=" + aurostd::utype2string(speciesB) + ")", _INPUT_ILLEGAL_); // CO20200624
    }
    if (speciesB >= (int) a.num_each_type.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "speciesB>=num_each_type.size() (speciesB=" + aurostd::utype2string(speciesB) + ")", _INPUT_ILLEGAL_); // CO20200624
    }
    a.SpeciesSwap(speciesA, speciesB);
    return a;
  }
} // namespace pflow

// ***************************************************************************
// pflow::VOLUME
// ***************************************************************************
namespace pflow {
  xstructure VOLUME(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 2) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --volume[|*|+]=x < POSCAR");
    }
    xstructure a(input, IOAFLOW_AUTO);
    if (tokens.at(0) == "VOLUME::EQUAL") {
      a = SetVolume(a, aurostd::string2utype<double>(tokens.at(1)));
    }
    if (tokens.at(0) == "VOLUME::MULTIPLY_EQUAL") {
      a = SetVolume(a, a.Volume() * aurostd::string2utype<double>(tokens.at(1)));
    }
    if (tokens.at(0) == "VOLUME::PLUS_EQUAL") {
      a = SetVolume(a, a.Volume() + aurostd::string2utype<double>(tokens.at(1)));
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return a;
  }
} // namespace pflow

// DX20180807 - put wyccar options into pflow - START
//  ***************************************************************************
//  pflow::WyckoffPositions()
//  ***************************************************************************
namespace pflow {
  string WyckoffPositions(aurostd::xoption& vpflow, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }

    const string options = vpflow.getattachedscheme("WYCKOFF_POSITIONS");
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (vpflow.flag("WYCKOFF_POSITIONS::USAGE")) {
      init::MessageOption(options, __AFLOW_FUNC__,
                          aurostd::liststring2string("aflow --Wyckoff|--Wyckoff_positions|--wyckoff|--wyckoff_positions|--wyccar[=<tolerance_value>|=tight|=loose] [--no_scan] [--setting=1|2|aflow] "
                                                     "[--magmom=[m1,m2,...|INCAR|OUTCAR]] < file  default: tolerance=(minimum_interatomic_distance)/100.0, setting=1"));
      return ""; // CO20200624 - the option was expressed successfully
    }
    // DX20190201 START
    const bool wyccar = vpflow.flag("WYCKOFF_POSITIONS::PRINT_WYCCAR"); // DX20210525
    const bool letters_only = vpflow.flag("WYCKOFF_POSITIONS::PRINT_LETTERS_ONLY");
    const bool site_symmetries_only = vpflow.flag("WYCKOFF_POSITIONS::PRINT_SITE_SYMMETRIES_ONLY");
    const bool multiplicities_only = vpflow.flag("WYCKOFF_POSITIONS::PRINT_MULTIPLICITIES_ONLY");

    // DX20190201 END
    xstructure str(input, IOAFLOW_AUTO);
    str.ReScale(1.0);

    // get tolerance
    double tolerance = pflow::getSymmetryTolerance(str, vpflow.getattachedscheme("WYCKOFF_POSITIONS::TOLERANCE")); // DX20200820 - consolidated setting tolerance into a function

    // ---------------------------------------------------------------------------
    // get space group setting
    const uint setting = pflow::getSpaceGroupSetting(vpflow.getattachedscheme("WYCKOFF_POSITIONS::SETTING")); // DX20210421 - consolidated space group setting into function

    // get magnetic moment
    if (vpflow.flag("WYCKOFF_POSITIONS::MAGNETIC")) {
      const string magmom_info = vpflow.getattachedscheme("WYCKOFF_POSITIONS::MAGNETIC");
      ProcessAndAddSpinToXstructure(str, magmom_info); // DX20191108 - condensed into a single function
    }

    // tolerance scan
    if (vpflow.flag("WYCKOFF_POSITIONS::NO_SCAN")) {
      str.sym_eps_no_scan = true; // DX20210406
    }
    // DX20190201 START
    const uint space_group_number = str.SpaceGroup_ITC(tolerance, -1, setting, str.sym_eps_no_scan);

    if (wyccar) {
      stringstream ss_output;
      str.iomode = IOVASP_WYCKCAR;
      ss_output << str;
      return ss_output.str();
    } else if (letters_only) {
      const string letters = SYM::ExtractWyckoffLettersString(str.wyccar_ITC);
      string sym_info = "sg=" + aurostd::utype2string<uint>(space_group_number) + ": " + letters;
      return sym_info;
    } else if (site_symmetries_only) {
      const string site_symmetries = SYM::ExtractWyckoffSiteSymmetriesString(str.wyccar_ITC);
      string sym_info = "sg=" + aurostd::utype2string<uint>(space_group_number) + ": " + site_symmetries;
      return sym_info;
    } else if (multiplicities_only) {
      const string multiplicities = SYM::ExtractWyckoffMultiplicitiesString(str.wyccar_ITC);
      string sym_info = "sg=" + aurostd::utype2string<uint>(space_group_number) + ": " + multiplicities;
      return sym_info;
    }

    const bool already_calculated = true;
    // ---------------------------------------------------------------------------
    // file type //DX20210525 - string to filetype
    filetype ftype = txt_ft;
    if (XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      ftype = txt_ft;
    } else if (XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
      ftype = json_ft;
    }

    return PrintWyckoffData(str, ftype, already_calculated); // no need to pass in setting info, etc. already calculated
  }
} // namespace pflow

// ***************************************************************************
// pflow::WYCKOFF
// ***************************************************************************
namespace pflow {
  xstructure WYCKOFF(vector<string> argv, istream& input) {
    xstructure a(input, IOAFLOW_AUTO);
    const int sg = aurostd::string2utype<int>(argv.at(2));
    a = WyckoffPOSITIONS(sg, a);
    cerr << a.spacegroup << endl;
    cerr << a.spacegroupnumber << endl;
    cerr << a.spacegroupoption << endl;
    return a;
  }
} // namespace pflow

// ***************************************************************************
// pflow::XRAY
// ***************************************************************************
namespace pflow {
  void XRAY(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 1) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --xray=lambda < POSCAR");
    }
    double l = 0.0;
    if (!tokens.empty()) {
      l = aurostd::string2utype<double>(tokens[0]);
    }

    const xstructure a(input, IOAFLOW_AUTO);
    cout << aflow::Banner("BANNER_TINY") << endl;
    PrintXray(a, l, cout);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
  void XRAY_PEAKS(const aurostd::xoption& vpflow, istream& input) { // CO20190520
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }

    const double lambda = aurostd::string2utype<double>(vpflow.getattachedscheme("XRAY_PEAKS"));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " lambda=" << lambda << endl;
    }
    if (lambda <= 0.0 || aurostd::isequal(lambda, 0.0)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "lambda <= 0", _VALUE_ILLEGAL_);
    }

    cout << aflow::Banner("BANNER_TINY") << endl;

    const xstructure a(input, IOAFLOW_AUTO);

    vector<double> v_twotheta;
    vector<double> v_intensity;
    vector<double> v_intensity_smooth;
    vector<uint> peak_indices = GetXrayPeaks(a, v_twotheta, v_intensity, v_intensity_smooth, lambda);
    if (v_twotheta.size() != v_intensity.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "v_twotheta.size()!=v_intensity.size()", _VALUE_ILLEGAL_);
    }
    if (v_twotheta.size() != v_intensity_smooth.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "v_twotheta.size()!=v_intensity_smooth.size()", _VALUE_ILLEGAL_);
    }

    // get amplitude
    vector<double> v_amplitude;
    vector<double> v_peaks_amplitude;
    double intmax = 1e-8;
    for (size_t i = 0; i < v_intensity.size(); i++) {
      if (v_intensity[i] > intmax) {
        intmax = v_intensity[i];
      }
    }
    for (size_t i = 0; i < v_intensity.size(); i++) {
      v_amplitude.push_back(100 * v_intensity[i] / intmax);
    }

    // match peak_indices to peaks
    vector<double> v_peaks_twotheta;
    vector<double> v_peaks_intensity;
    for (size_t i = 0; i < peak_indices.size(); i++) {
      v_peaks_twotheta.push_back(v_twotheta[peak_indices[i]]);
      v_peaks_intensity.push_back(v_intensity[peak_indices[i]]);
      v_peaks_amplitude.push_back(v_amplitude[peak_indices[i]]);
    }
    if (v_peaks_twotheta.size() != v_peaks_intensity.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "v_peaks_twotheta.size()!=v_peaks_intensity.size()", _VALUE_ILLEGAL_);
    }

    cout << "X-Ray Peaks:" << endl;
    cout << "Two-Theta=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(v_peaks_twotheta, 5, false), ",") << endl; // no roff
    cout << "Intensity=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(v_peaks_intensity, 10, false), ",") << endl; // no roff
    cout << "Amplitude=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(v_peaks_amplitude, 5, false), ",") << endl; // no roff

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
  void READ_XRAY_DATA(const string& filename, vector<double>& v_twotheta, vector<double>& v_intensity) { // CO20190620
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }

    if (filename.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No filename provided", _FILE_ERROR_);
    }
    if (!aurostd::FileExist(filename)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "File does not exist: ", _FILE_NOT_FOUND_);
    }
    if (aurostd::FileEmpty(filename)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "File empty: ", _FILE_ERROR_);
    }

    vector<string> data_file_lines;
    aurostd::file2vectorstring(filename, data_file_lines);
    string line;
    string::size_type loc;
    vector<double> tokens;
    for (size_t i = 0; i < data_file_lines.size(); i++) {
      line = data_file_lines[i];
      loc = line.find("#");
      line = line.substr(0, loc);
      aurostd::RemoveControlCodeCharactersFromString(line, line);
      line = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(line);
      if (!line.empty()) {
        aurostd::string2tokens<double>(line, tokens, ","); // csv style
        // if(tokens.size()!=2){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Line[i="+aurostd::utype2string(i)+"] has "+aurostd::utype2string(tokens.size())+" tokens (expected 2): line=\""+line+"\"",_FILE_WRONG_FORMAT_);}
        if (tokens.size() != 2) {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Line[i=" + aurostd::utype2string(i) + "] has " + aurostd::utype2string(tokens.size()) + " tokens (expected 2)", _FILE_WRONG_FORMAT_);
        }
        v_twotheta.push_back(tokens[0]);
        v_intensity.push_back(tokens[1]);
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " v_twotheta.size()=" << v_twotheta.size() << endl;
      cerr << __AFLOW_FUNC__ << " v_intensity.size()=" << v_intensity.size() << endl;
    }
  }
#define XRAY_DATA_PLOT_FILE "aflow_xray_data_plot_file.txt"
#define XRAY_DATA_PEAKS_FILE "aflow_xray_data_peaks_file.txt"
  void PRINT_XRAY_DATA_PLOT(const aurostd::xoption& vpflow, istream& input) {
    const xstructure str(input, IOAFLOW_AUTO);
    return PRINT_XRAY_DATA_PLOT(vpflow, str);
  } // CO20190520
  void PRINT_XRAY_DATA_PLOT(const aurostd::xoption& vpflow, const xstructure& str) {
    const double lambda = aurostd::string2utype<double>(vpflow.getattachedscheme("PRINT_XRAY_DATA_PLOT"));
    const string directory = vpflow.getattachedscheme("PRINT_XRAY_DATA_PLOT::DIRECTORY");
    return PRINT_XRAY_DATA_PLOT(str, lambda, directory);
  } // don't use XHOST.vflag_control.getattachedscheme("DIRECTORY") for directory, might interfere with LIB2RAW //CO20190520
  void PRINT_XRAY_DATA_PLOT(istream& input, double lambda, const string& directory) {
    const xstructure str(input, IOAFLOW_AUTO);
    return PRINT_XRAY_DATA_PLOT(str, lambda, directory);
  } // CO20190520
  void PRINT_XRAY_DATA_PLOT(const xstructure& str, double lambda, const string& directory) { // CO20190520
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " lambda=" << lambda << endl;
    }
    if (lambda <= 0.0 || aurostd::isequal(lambda, 0.0)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "lambda <= 0", _VALUE_ILLEGAL_);
    }
    vector<double> v_twotheta;
    vector<double> v_intensity;
    GetXray2ThetaIntensity(str, v_twotheta, v_intensity, lambda); // v_amplitude can be grabbed later
    if (v_twotheta.size() != v_intensity.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "v_twotheta.size()!=v_intensity.size()", _INDEX_MISMATCH_);
    }
    return PRINT_XRAY_DATA_PLOT(v_twotheta, v_intensity, directory);
  }
  void PRINT_XRAY_DATA_PLOT(const aurostd::xoption& vpflow, const string& directory) { // CO20190520
    // assume a file input from vpflow
    const bool LDEBUG = (false || XHOST.DEBUG);

    const string filename = vpflow.getattachedscheme("PLOT_XRAY_FILE");
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " filename=" << filename << endl;
    }

    return PRINT_XRAY_DATA_PLOT(filename, directory);
  }
  void PRINT_XRAY_DATA_PLOT(const string& filename, const string& directory) { // CO20190520
    vector<double> v_twotheta;
    vector<double> v_intensity;
    READ_XRAY_DATA(filename, v_twotheta, v_intensity);
    return PRINT_XRAY_DATA_PLOT(v_twotheta, v_intensity, directory);
  }
  void PRINT_XRAY_DATA_PLOT(const vector<double>& v_twotheta, const vector<double>& v_intensity, const string& _directory) { // CO20190620
    const bool LDEBUG = (false || XHOST.DEBUG);

    string directory = _directory;
    if (directory.empty()) {
      directory = ".";
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " directory=" << directory << endl;
    }

    cout << aflow::Banner("BANNER_TINY") << endl;

    vector<double> v_intensity_smooth;
    vector<uint> peak_indices = GetXrayPeaks(v_twotheta, v_intensity, v_intensity_smooth);
    if (v_twotheta.size() != v_intensity_smooth.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "v_twotheta.size()!=v_intensity_smooth.size()", _VALUE_ILLEGAL_);
    }

    // get amplitude
    vector<double> v_amplitude;
    vector<double> v_peaks_amplitude;
    double intmax = 1e-8;
    for (size_t i = 0; i < v_intensity.size(); i++) {
      if (v_intensity[i] > intmax) {
        intmax = v_intensity[i];
      }
    }
    for (size_t i = 0; i < v_intensity.size(); i++) {
      v_amplitude.push_back(100 * v_intensity[i] / intmax);
    }

    // match peak_indices to peaks
    vector<double> v_peaks_twotheta;
    vector<double> v_peaks_intensity;
    for (size_t i = 0; i < peak_indices.size(); i++) {
      v_peaks_twotheta.push_back(v_twotheta[peak_indices[i]]);
      v_peaks_intensity.push_back(v_intensity[peak_indices[i]]);
      v_peaks_amplitude.push_back(v_amplitude[peak_indices[i]]);
    }
    if (v_peaks_twotheta.size() != v_peaks_intensity.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "v_peaks_twotheta.size()!=v_peaks_intensity.size()", _INDEX_MISMATCH_);
    }

    stringstream data_file_ss;
    string data_file;

    // PLOT FOR TWO-THETA, INTENSITY, AMPLITUDE, AND SMOOTH
    data_file_ss.str("");
    data_file_ss << aurostd::PaddedPOST("#Two-theta(degrees)", 20) << " ";
    data_file_ss << aurostd::PaddedPOST("Intensity", 20) << " ";
    data_file_ss << aurostd::PaddedPOST("Amplitude", 20) << " ";
    data_file_ss << aurostd::PaddedPOST("Intensity_Smoothed", 20) << " ";
    data_file_ss << endl;
    for (size_t i = 0; i < v_twotheta.size(); i++) {
      data_file_ss << aurostd::PaddedPOST<double>(v_twotheta[i], 20) << " ";
      data_file_ss << aurostd::PaddedPOST<double>(v_intensity[i], 20) << " ";
      data_file_ss << aurostd::PaddedPOST<double>(v_amplitude[i], 20) << " ";
      data_file_ss << aurostd::PaddedPOST<double>(v_intensity_smooth[i], 20) << " ";
      data_file_ss << endl;
    }
    data_file = directory + "/" + XRAY_DATA_PLOT_FILE;
    aurostd::StringSubstInPlace(data_file, "//", "/");
    aurostd::stringstream2file(data_file_ss, data_file);

    // PLOT FOR TWO-THETA, INTENSITY, AMPLITUDE, AND SMOOTH
    data_file_ss.str("");
    data_file_ss << aurostd::PaddedPOST("#Two-theta(degrees)", 20) << " ";
    data_file_ss << aurostd::PaddedPOST("Intensity", 20) << " ";
    data_file_ss << aurostd::PaddedPOST("Amplitude", 20) << " ";
    data_file_ss << endl;
    for (size_t i = 0; i < v_peaks_twotheta.size(); i++) {
      data_file_ss << aurostd::PaddedPOST<double>(v_peaks_twotheta[i], 20) << " ";
      data_file_ss << aurostd::PaddedPOST<double>(v_peaks_intensity[i], 20) << " ";
      data_file_ss << aurostd::PaddedPOST<double>(v_peaks_amplitude[i], 20) << " ";
      data_file_ss << endl;
    }
    data_file = directory + "/" + XRAY_DATA_PEAKS_FILE;
    aurostd::StringSubstInPlace(data_file, "//", "/");
    aurostd::stringstream2file(data_file_ss, data_file);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
  void PLOT_XRAY(const aurostd::xoption& vpflow, istream& input) {
    const xstructure str(input, IOAFLOW_AUTO);
    return PLOT_XRAY(vpflow, str);
  } // CO20190520
  void PLOT_XRAY(const aurostd::xoption& vpflow, const xstructure& str) {
    const bool force_generic_title = vpflow.flag("PLOT_XRAY::FORCE_GENERIC_TITLE");
    const double lambda = aurostd::string2utype<double>(vpflow.getattachedscheme("PLOT_XRAY"));
    const string directory = vpflow.getattachedscheme("PLOT_XRAY::DIRECTORY");
    const bool keep_gp = vpflow.flag("PLOT_XRAY::KEEP_GP");
    return PLOT_XRAY(str, lambda, directory, keep_gp, force_generic_title);
  } // don't use XHOST.vflag_control.getattachedscheme("DIRECTORY") for directory, might interfere with LIB2RAW //CO20190520
  void PLOT_XRAY(istream& input, double lambda, const string& directory, bool keep_gp, bool force_generic_title) {
    const xstructure str(input, IOAFLOW_AUTO);
    return PLOT_XRAY(str, lambda, directory, keep_gp, force_generic_title);
  } // CO20190520
  void PLOT_XRAY(const xstructure& str, double lambda, const string& directory, bool keep_gp, bool force_generic_title) { // CO20190520
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " lambda=" << lambda << endl;
    }
    if (lambda <= 0.0 || aurostd::isequal(lambda, 0.0)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "lambda <= 0", _VALUE_ILLEGAL_);
    }
    vector<double> v_twotheta;
    vector<double> v_intensity;
    GetXray2ThetaIntensity(str, v_twotheta, v_intensity, lambda); // v_amplitude can be grabbed later
    if (v_twotheta.size() != v_intensity.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "v_twotheta.size()!=v_intensity.size()", _INDEX_MISMATCH_);
    }

    string title = aurostd::fixStringLatex(str.title, true, false); // double_back_slash==true (gnuplot), not symmetry sting
    if (false || force_generic_title) {
      title.clear();
    } // force generic
    if (title.empty()) {
      title = getGenericTitleXStructure(str, true);
    } // latex
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " title=\"" << title << "\"" << endl;
    }

    return PLOT_XRAY(v_twotheta, v_intensity, title, directory, keep_gp);
  }
  void PLOT_XRAY(const aurostd::xoption& vpflow, const string& title, const string& directory, bool keep_gp) { // CO20190520
    // assume a file input from vpflow
    const bool LDEBUG = (false || XHOST.DEBUG);

    const string filename = vpflow.getattachedscheme("PLOT_XRAY_FILE");
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " filename=" << filename << endl;
    }

    return PLOT_XRAY(filename, title, directory, keep_gp);
  }
  void PLOT_XRAY(const string& filename, const string& title, const string& directory, bool keep_gp) { // CO20190520
    vector<double> v_twotheta;
    vector<double> v_intensity;
    READ_XRAY_DATA(filename, v_twotheta, v_intensity);
    return PLOT_XRAY(v_twotheta, v_intensity, title, directory, keep_gp);
  }
  void PLOT_XRAY(const vector<double>& v_twotheta, const vector<double>& v_intensity, const string& _title, const string& _directory, bool keep_gp) { // CO20190620
    const bool LDEBUG = (false || XHOST.DEBUG);
    stringstream message;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }

    if (!aurostd::IsCommandAvailable("gnuplot")) {
      message << "gnuplot is not available, please put in path" << endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    if (!aurostd::IsCommandAvailable("pdflatex")) {
      message << "pdflatex is not available, please put in path" << endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }

    string directory = _directory;
    if (directory.empty()) {
      directory = ".";
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " directory=" << directory << endl;
    }

    const string curdir = aurostd::getPWD();

    const string PLOT_tmp_dir = aurostd::TmpDirectoryCreate("XRAY_PLOT");
    fs::current_path(PLOT_tmp_dir.c_str());
    PRINT_XRAY_DATA_PLOT(v_twotheta, v_intensity, ".");

    const bool plot_intensity = true; // else plot amplitude

    vector<string> data_file_lines;
    aurostd::file2vectorstring(XRAY_DATA_PLOT_FILE, data_file_lines);
    vector<double> x;
    vector<double> y1;
    vector<double> y2;
    vector<double> tokens;
    string::size_type loc;
    string line;
    double y_max = 0.0; // this is really all we are after
    for (size_t i = 0, fl_size_i = data_file_lines.size(); i < fl_size_i; i++) {
      line = data_file_lines[i];
      loc = line.find("#");
      line = line.substr(0, loc);
      line = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(line);
      if (!line.empty()) {
        aurostd::string2tokens<double>(line, tokens, " ");
        x.push_back(tokens[0]);
        y1.push_back(tokens[(plot_intensity ? 1 : 2)]);
        y2.push_back(tokens[3]);
        if (y1.back() > y_max) {
          y_max = y1.back();
        }
        if (y2.back() > y_max) {
          y_max = y2.back();
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " x=" << x.back() << " y1=" << y1.back() << " y2=" << y2.back() << endl;
        }
      }
    }

    const int exp = (int) floor(log10(y_max));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " y_max=" << y_max << endl;
      cerr << __AFLOW_FUNC__ << " exp=" << exp << endl;
    }

    string title = aurostd::fixStringLatex(_title, true, false); // double_back_slash==true (gnuplot), not symmetry sting
    if (title.empty()) {
      title = "X-Ray Plot";
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " title=\"" << title << "\"" << endl;
    }

    stringstream plot_file_ss;
    const string xray_tex = "aflow_xray_plot.tex";
    string xray_pdf = xray_tex;
    aurostd::StringSubstInPlace(xray_pdf, ".tex", ".pdf");

    plot_file_ss << "set terminal epslatex color standalone" << endl;
    plot_file_ss << "set output \"" + xray_tex + "\"" << endl;
    plot_file_ss << endl;
    // plot_file_ss << "unset key" << endl;
    plot_file_ss << "set yrange [0:*]" << endl;
    plot_file_ss << "set xrange [0:180]" << endl;
    plot_file_ss << "set xlabel \"2$\\\\theta$ (degrees)\"" << endl;
    plot_file_ss << "set ylabel \"" << (plot_intensity ? string("intensity") : string("amplitude")) << R"( $\\left[\\times 10^{)" << exp << "}\\\\right]$" << "\"" << endl;
    // plot_file_ss << "set ylabel \"amplitude\"" << endl;
    plot_file_ss << "set title \"" << title << "\"" << endl;
    plot_file_ss << endl;

    bool plot_smooth = true;
    const bool plot_peaks = true;
    plot_smooth = (plot_smooth && plot_intensity); // smooth is based on intensity, not amplitude

    plot_file_ss << "plot \\" << endl;
    plot_file_ss << "\"" + string(XRAY_DATA_PLOT_FILE) + "\" using 1:(10**(-" << exp << ")*$" << (plot_intensity ? 2 : 3) << R"() with lines lc rgb "blue" lw 4 title ")"
                 << (plot_intensity ? string("intensity") : string("amplitude")) << "\"" << (plot_smooth || plot_peaks ? string(", \\") : string("")) << endl;
    if (plot_smooth) {
      plot_file_ss << "\"" + string(XRAY_DATA_PLOT_FILE) + "\" using 1:(10**(-" << exp << ")*$" << 4 << R"() with lines lc rgb "red" lw 4 title "smoothed")" << (plot_peaks ? string(", \\") : string("")) << endl;
    } // can't plot this if we're plotting amplitude
    if (plot_peaks) {
      plot_file_ss << "\"" + string(XRAY_DATA_PEAKS_FILE) + "\" using 1:(10**(-" << exp << ")*$" << (plot_intensity ? 2 : 3) << R"() lt 2 lc rgb "#006400" ps 2 lw 4 title "peaks")" << endl;
    } // green - #006400

    const string plot_file = "aflow_xray_plot.gp";
    aurostd::stringstream2file(plot_file_ss, plot_file);
    aurostd::execute2string("gnuplot " + plot_file); // 2string so it does not go into output, put to cout if you need
    aurostd::execute2string("pdflatex " + xray_tex); // 2string so it does not go into output, put to cout if you need

    vector<string> files2move;

    ostream& oss = cout;
    ofstream FileMESSAGE;
    _aflags aflags;
    aflags.Directory = directory;

    if (false || keep_gp) {
      if (aurostd::FileExist(XRAY_DATA_PLOT_FILE)) {
        files2move.push_back(PLOT_tmp_dir + "/" + XRAY_DATA_PLOT_FILE);
      } else {
        message << XRAY_DATA_PLOT_FILE << " was not created";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
      }
      if (aurostd::FileExist(XRAY_DATA_PEAKS_FILE)) {
        files2move.push_back(PLOT_tmp_dir + "/" + XRAY_DATA_PEAKS_FILE);
      } else {
        message << XRAY_DATA_PEAKS_FILE << " was not created";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
      }
      if (aurostd::FileExist(plot_file)) {
        files2move.push_back(PLOT_tmp_dir + "/" + plot_file);
      } else {
        message << plot_file << " was not created";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
      }
    }
    if (aurostd::FileExist(xray_pdf)) {
      files2move.push_back(PLOT_tmp_dir + "/" + xray_pdf);
    } else {
      message << xray_pdf << " was not created";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
    }

    fs::current_path(curdir.c_str());
    aurostd::file2directory(files2move, directory);
    aurostd::RemoveDirectory(PLOT_tmp_dir);

    string destination = directory + "/" + xray_pdf;
    aurostd::StringSubstInPlace(destination, "//", "/");
    message << "Done. See " << destination << "." << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::XYZ
// ***************************************************************************
namespace pflow {
  void XYZ(const string& options, istream& input) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    const xvector<int> ijk(3); // default
    ijk[1] = 1;
    ijk[2] = 1;
    ijk[3] = 1; // default
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");

    if (tokens.size() > 3) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --xyz[=n1[,n2[,n3]]] < POSCAR");
    }

    if (!tokens.empty()) {
      ijk[1] = aurostd::string2utype<int>(tokens[0]);
    }
    if (tokens.size() >= 2) {
      ijk[2] = aurostd::string2utype<int>(tokens[1]);
    }
    if (tokens.size() >= 3) {
      ijk[3] = aurostd::string2utype<int>(tokens[2]);
    }

    const xstructure a(input, IOAFLOW_AUTO);
    PrintXYZ(a, ijk, cout);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::XYZINSPHERE
// ***************************************************************************
namespace pflow {
  void XYZINSPHERE(istream& input, double radius) {
    cout << aflow::Banner("BANNER_TINY") << endl;
    const xstructure a(input, IOAFLOW_AUTO);
    PrintXYZInSphere(a, radius, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::XYZWS
// ***************************************************************************
namespace pflow {
  void XYZWS(istream& input) {
    cout << aflow::Banner("BANNER_TINY") << endl;
    const xstructure a(input, IOAFLOW_AUTO);
    PrintXYZws(a, cout);
  }
} // namespace pflow

// ***************************************************************************
// pflow::ZVAL
// ***************************************************************************
namespace pflow {
  void ZVAL(const string& options) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.empty() || tokens.size() > 2) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow " + tokens.at(0) + "[=directory]");
    }
    string Directory;
    if (tokens.size() == 1) {
      Directory = "./";
    }
    if (tokens.size() == 2) {
      Directory = tokens[1];
    }

    vector<double> vZVAL;
    vector<double> sZVAL;
    vector<double> vPOMASS;
    vector<double> sPOMASS;
    if (tokens.at(0) == "ZVAL") {
      cout << "Total ZVAL (from PP) = " << GetZVAL(Directory, vZVAL) << endl;
    }
    if (tokens.at(0) == "ZVAL::CELL") {
      cout << "Total ZVAL_CELL (from PP and xstructure) = " << GetCellAtomZVAL(Directory, vZVAL, sZVAL, "CELL") << endl;
    }
    if (tokens.at(0) == "ZVAL::ATOM") {
      cout << "Total ZVAL_ATOM (from PP and xstructure) = " << GetCellAtomZVAL(Directory, vZVAL, sZVAL, "ATOM") << endl;
    }
    if (tokens.at(0) == "POMASS") {
      cout << "Total POMASS (from PP) = " << GetPOMASS(Directory, vPOMASS) << endl;
    }
    if (tokens.at(0) == "POMASS::CELL") {
      cout << "Total POMASS_CELL (from PP and xstructure) = " << GetCellAtomPOMASS(Directory, vPOMASS, sPOMASS, "CELL") << endl;
    }
    if (tokens.at(0) == "POMASS::ATOM") {
      cout << "Total POMASS_ATOM (from PP and xstructure) = " << GetCellAtomPOMASS(Directory, vPOMASS, sPOMASS, "ATOM") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::OPARAMETER
// ***************************************************************************

bool AreAtomsEquivalent(const deque<_atom>& atoms1, const deque<_atom>& atoms2, bool oparameter_check, bool poccupation_check) {
  double epsilon = 0.0001;
  double epsilon_pocc = 0.05;
  if (atoms1.size() != atoms2.size()) {
    return false; // cant be equivalent if different numbers
  }
  for (size_t iat1 = 0; iat1 < atoms1.size(); iat1++) {
    bool found = false;
    for (size_t iat2 = 0; iat2 < atoms2.size() && found == false; iat2++) {
      if (atoms1[iat1].type == atoms2[iat2].type && aurostd::isequal(atoms1[iat1].fpos, atoms2[iat2].fpos, epsilon)) {
        if (oparameter_check == false && poccupation_check == false) {
          found = true;
        }
        if (oparameter_check == true && poccupation_check == false) {
          if (atoms1[iat1].order_parameter_value == atoms2[iat2].order_parameter_value) {
            found = true;
          }
        }
        //  atoms1[iat1].order_parameter_atom==true && atoms2[iat2].order_parameter_atom==true &&
        if (oparameter_check == false && poccupation_check == true) {
          if (aurostd::isequal(atoms1[iat1].partial_occupation_value, atoms2[iat2].partial_occupation_value, epsilon_pocc)) {
            found = true;
          }
        }
        //  atoms1[iat1].partial_occupation_atom==true && atoms2[iat2].partial_occupation_atom==true &&
        if (oparameter_check == true && poccupation_check == true) {
          if (atoms1[iat1].order_parameter_value == atoms2[iat2].order_parameter_value && aurostd::isequal(atoms1[iat1].partial_occupation_value, atoms2[iat2].partial_occupation_value, epsilon_pocc)) {
            found = true;
          }
        }
        //  atoms1[iat1].order_parameter_atom==true && atoms2[iat2].order_parameter_atom==true &&
        //  atoms1[iat1].partial_occupation_atom==true && atoms2[iat2].partial_occupation_atom==true &&
      }
    }
    if (found == false) {
      return false;
    }
  }
  return true; // survived everything
}

bool AreAtomsEquivalent(const deque<_atom>& atoms1, const deque<_atom>& atoms2) {
  return AreAtomsEquivalent(atoms1, atoms2, false, false); // no order check
}

bool AreAtomsOrderEquivalent(const deque<_atom>& atoms1, const deque<_atom>& atoms2) {
  return AreAtomsEquivalent(atoms1, atoms2, true, false); // no order check
}

bool AreAtomsOccupationEquivalent(const deque<_atom>& atoms1, const deque<_atom>& atoms2) {
  return AreAtomsEquivalent(atoms1, atoms2, false, true); // no order check
}

bool AreAtomsOccupationOrderEquivalent(const deque<_atom>& atoms1, const deque<_atom>& atoms2) {
  return AreAtomsEquivalent(atoms1, atoms2, true, true); // no order check
}

class _sort_xstructure_order_parameter_orbit { // sorting through reference
public:
  bool operator()(const xstructure& str1, const xstructure& str2) const { return (bool) (str1.order_parameter_orbit < str2.order_parameter_orbit); }
};

class _rsort_xstructure_order_parameter_orbit { // sorting through reference
public:
  bool operator()(const xstructure& str1, const xstructure& str2) const { return (bool) (str1.order_parameter_orbit > str2.order_parameter_orbit); }
};

int order_parameter_sum(const xstructure& str) {
  int _sum = 0;
  for (size_t iat = 0; iat < str.order_parameter_atoms.size(); iat++) {
    _sum += (int) str.atoms.at(str.order_parameter_atoms[iat]).order_parameter_value;
  }
  return _sum;
}

// ***************************************************************************
// helpIndividualOption
// ***************************************************************************
void helpIndividualOption(vector<string>& argv) {
  // "aflow --help option"
  // to print out the help information of given option
  // use the information stored in README_AFLOW_PFLOW.txt

  stringstream help_file;
  help_file << aurostd::EmbData::get_content("README_AFLOW_PFLOW.TXT", "README");

  bool help_found = false;
  string option = argv.at(2);

  const string seperation_line_pattern = "******";

  string line;

  // to get standard form of option --option+" "
  option.append(" ");
  option.insert(0, "--");

  bool output_flag = false;
  bool begin_section_flag = false;
  while (getline(help_file, line)) {
    // get the begining of the help section of an option
    // it begins with "aflow --option_name"
    // it may have one or several space between the two words

    size_t pos = line.find(option);
    const string program_name = "aflow";
    const string option_prefix = "--";
    const size_t pos2 = line.find(program_name);
    const size_t pos3 = line.find(option_prefix);
    string::iterator str_itr;

    if ((pos2 < pos3) && (pos3 != string::npos)) {
      for (str_itr = line.begin() + pos2 + program_name.size(); str_itr < line.begin() + pos3; str_itr++) {
        if (*str_itr != ' ') {
          begin_section_flag = false;
          break;
        }
      }
      begin_section_flag = true;
    } else {
      begin_section_flag = false;
    }

    if (begin_section_flag && (pos != string::npos) && (!output_flag)) {
      // find the line with --option
      // output this line and lines below
      // until next line begining with "aflow --" pattern

      const int last_space_pos = line.find_first_not_of(" ");
      line.erase(line.begin(), line.begin() + last_space_pos);
      cout << line << endl;
      output_flag = true;
      help_found = true;
      continue;
    } else if (begin_section_flag && output_flag) {
      // reset the flag and output
      // the other part of help information if found

      pos = line.find(option);
      if ((pos != string::npos)) {
        const int last_space_pos = line.find_first_not_of(" ");
        line.erase(line.begin(), line.begin() + last_space_pos);
        cout << line << endl;
        output_flag = true;
      } else {
        output_flag = false;
      }
    } else if (!begin_section_flag && output_flag) {
      // output the part of help section not
      // beginning with "aflow --option" pattern

      if (line.find(seperation_line_pattern) == string::npos) {
        // do not output seperation line
        const int last_space_pos = line.find_first_not_of(" ");
        line.erase(line.begin(), line.begin() + last_space_pos);
        cout << "    " << line << endl;
      }
    }
  }

  if (!help_found) {
    // not supported option
    const size_t pos = option.find_last_of(" ");
    const string option_raw = option.substr(2, pos - 2);
    cerr << "No option \"" << option_raw << "\" is found." << endl;
    cerr << "Try \"aflow --help\" to see all supported options." << endl;
    return;
  }
  return;
}

// ***************************************************************************

// ----------------------------------------------------------------------------
// get_itemized_vector_string_from_input
// Stefano Curtarolo
// namespace aurostd {
//   bool get_itemized_vector_string_from_input(vector<string> &argv,const string& s0,vector<string>& tokens,const string& delimiter) {// =":") {
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2flag(argv,s0neq)) {aurostd::args2string(argv,s0neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s0neq,EMPTY_WORDING);}
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return false;
//     return true;
//   }
//   bool get_itemized_vector_string_from_input(vector<string> &argv,const string& s0,const string& s1,vector<string>& tokens,const string& delimiter) {// =":") {
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2attachedflag(argv,s1equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2flag(argv,s0neq)) {aurostd::args2string(argv,s0neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s0neq,EMPTY_WORDING);}
//     if(aurostd::args2flag(argv,s1neq)) {aurostd::args2string(argv,s1neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s1neq,EMPTY_WORDING);}
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return false;
//     return true;
//   }
//   bool get_itemized_vector_string_from_input(vector<string> &argv,const string& s0,const string& s1,const string& s2,vector<string>& tokens,const string& delimiter) {// =":") {
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
//     string s2neq=s2,s2equ;aurostd::StringSubst(s2neq,"=","");s2equ=s2neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2attachedflag(argv,s1equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2attachedflag(argv,s2equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s2equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2flag(argv,s0neq)) {aurostd::args2string(argv,s0neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s0neq,EMPTY_WORDING);}
//     if(aurostd::args2flag(argv,s1neq)) {aurostd::args2string(argv,s1neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s1neq,EMPTY_WORDING);}
//     if(aurostd::args2flag(argv,s2neq)) {aurostd::args2string(argv,s2neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s2neq,EMPTY_WORDING);}
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return false;
//     return true;
//   }
//   bool get_itemized_vector_string_from_input(vector<string> &argv,const string& s0,const string& s1,const string& s2,const string& s3,vector<string>& tokens,const string& delimiter) {// =":") {
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
//     string s2neq=s2,s2equ;aurostd::StringSubst(s2neq,"=","");s2equ=s2neq+"=";
//     string s3neq=s3,s3equ;aurostd::StringSubst(s3neq,"=","");s3equ=s3neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2attachedflag(argv,s1equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2attachedflag(argv,s2equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s2equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2attachedflag(argv,s3equ)) {aurostd::string2tokens(aurostd::args2attachedstring(argv,s3equ,EMPTY_WORDING),tokens,delimiter);}
//     if(aurostd::args2flag(argv,s0neq)) {aurostd::args2string(argv,s0neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s0neq,EMPTY_WORDING);}
//     if(aurostd::args2flag(argv,s1neq)) {aurostd::args2string(argv,s1neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s1neq,EMPTY_WORDING);}
//     if(aurostd::args2flag(argv,s2neq)) {aurostd::args2string(argv,s2neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s2neq,EMPTY_WORDING);}
//     if(aurostd::args2flag(argv,s3neq)) {aurostd::args2string(argv,s3neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s3neq,EMPTY_WORDING);}
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return false;
//     return true;
//   }
// }

// // ----------------------------------------------------------------------------
// // getproto_itemized_vector_string_from_input
// // Stefano Curtarolo
// namespace aurostd {
//   bool getproto_itemized_vector_string_from_input(vector<string> &argv,const string& s0,vector<string>& tokens,const string& delimiter) {// =":") {
//     string ss;tokens.clear();vector<string> stokens;
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) ss=aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s0neq)) ss=aurostd::args2string(argv,s0neq,EMPTY_WORDING);
//     if(aurostd::substring2bool(ss,"./")) aurostd::StringSubst(ss,"./","");
//     if(ss!="") {
//       if(!aurostd::substring2bool(ss,"/")) { return get_itemized_vector_string_from_input(argv,s0,tokens,delimiter);
//       } else {
// 	aurostd::string2tokens(ss,stokens,"/");
// 	tokens.push_back(stokens.at(1));
// 	KBIN::VASP_SplitAlloyPseudoPotentials(stokens.at(0),stokens);
// 	for(size_t i=0;i<stokens.size();i++) tokens.push_back(stokens.at(i));
//       }
//     }
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return false;
//     return true;
//   }
//   bool getproto_itemized_vector_string_from_input(vector<string> &argv,const string& s0,const string& s1,vector<string>& tokens,const string& delimiter) {// =":") {
//     string ss;tokens.clear();vector<string> stokens;
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) ss=aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s0neq)) ss=aurostd::args2string(argv,s0neq,EMPTY_WORDING);
//     if(aurostd::args2attachedflag(argv,s1equ)) ss=aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s1neq)) ss=aurostd::args2string(argv,s1neq,EMPTY_WORDING);
//     if(aurostd::substring2bool(ss,"./")) aurostd::StringSubst(ss,"./","");
//     if(ss!="") {
//       if(!aurostd::substring2bool(ss,"/")) { return get_itemized_vector_string_from_input(argv,s0,s1,tokens,delimiter);
//       } else {
// 	aurostd::string2tokens(ss,stokens,"/");
// 	tokens.push_back(stokens.at(1));
// 	KBIN::VASP_SplitAlloyPseudoPotentials(stokens.at(0),stokens);
// 	for(size_t i=0;i<stokens.size();i++) tokens.push_back(stokens.at(i));
//       }
//     }
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return false;
//     return true;
//   }
//   bool getproto_itemized_vector_string_from_input(vector<string> &argv,const string& s0,const string& s1,const string& s2,vector<string>& tokens,const string& delimiter) {// =":") {
//     string ss;tokens.clear();vector<string> stokens;
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
//     string s2neq=s2,s2equ;aurostd::StringSubst(s2neq,"=","");s2equ=s2neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) ss=aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s0neq)) ss=aurostd::args2string(argv,s0neq,EMPTY_WORDING);
//     if(aurostd::args2attachedflag(argv,s1equ)) ss=aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s1neq)) ss=aurostd::args2string(argv,s1neq,EMPTY_WORDING);
//     if(aurostd::args2attachedflag(argv,s2equ)) ss=aurostd::args2attachedstring(argv,s2equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s2neq)) ss=aurostd::args2string(argv,s2neq,EMPTY_WORDING);
//     if(aurostd::substring2bool(ss,"./")) aurostd::StringSubst(ss,"./","");
//     if(ss!="") {
//       if(!aurostd::substring2bool(ss,"/")) { return get_itemized_vector_string_from_input(argv,s0,s1,s2,tokens,delimiter);
//       } else {
// 	aurostd::string2tokens(ss,stokens,"/");
// 	tokens.push_back(stokens.at(1));
// 	KBIN::VASP_SplitAlloyPseudoPotentials(stokens.at(0),stokens);
// 	for(size_t i=0;i<stokens.size();i++) tokens.push_back(stokens.at(i));
//       }
//     }
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return false;
//     return true;
//   }
//   bool getproto_itemized_vector_string_from_input(vector<string> &argv,const string& s0,const string& s1,const string& s2,const string& s3,vector<string>& tokens,const string& delimiter) {// =":") {
//     string ss;tokens.clear();vector<string> stokens;
//     string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//     string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
//     string s2neq=s2,s2equ;aurostd::StringSubst(s2neq,"=","");s2equ=s2neq+"=";
//     string s3neq=s3,s3equ;aurostd::StringSubst(s3neq,"=","");s3equ=s3neq+"=";
//     if(aurostd::args2attachedflag(argv,s0equ)) ss=aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s0neq)) ss=aurostd::args2string(argv,s0neq,EMPTY_WORDING);
//     if(aurostd::args2attachedflag(argv,s1equ)) ss=aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s1neq)) ss=aurostd::args2string(argv,s1neq,EMPTY_WORDING);
//     if(aurostd::args2attachedflag(argv,s2equ)) ss=aurostd::args2attachedstring(argv,s2equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s2neq)) ss=aurostd::args2string(argv,s2neq,EMPTY_WORDING);
//     if(aurostd::args2attachedflag(argv,s3equ)) ss=aurostd::args2attachedstring(argv,s3equ,EMPTY_WORDING);
//     if(aurostd::args2flag(argv,s3neq)) ss=aurostd::args2string(argv,s3neq,EMPTY_WORDING);
//     if(aurostd::substring2bool(ss,"./")) aurostd::StringSubst(ss,"./","");
//     if(ss!="") {
//       if(!aurostd::substring2bool(ss,"/")) { return get_itemized_vector_string_from_input(argv,s0,s1,s2,s3,tokens,delimiter);
//       } else {
// 	aurostd::string2tokens(ss,stokens,"/");
// 	tokens.push_back(stokens.at(1));
// 	KBIN::VASP_SplitAlloyPseudoPotentials(stokens.at(0),stokens);
// 	for(size_t i=0;i<stokens.size();i++) tokens.push_back(stokens.at(i));
//       }
//     }
//     if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);aurostd::string2tokens(s0equ,tokens,delimiter);}
//     if(tokens.size()==0) return false;
//     return true;
//   }
// }

// ***************************************************************************
// RICHARD stuff on XRD
namespace pflow {
  double GetAtomicPlaneDist(const string& options, istream& cin) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() != 3) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --xrd_dist=h,k,l < POSCAR");
    }
    // move on
    double h = 0;
    double k = 0;
    double l = 0;
    if (!tokens.empty()) {
      h = aurostd::string2utype<double>(tokens[0]);
    }
    if (tokens.size() >= 2) {
      k = aurostd::string2utype<double>(tokens[1]);
    }
    if (tokens.size() >= 3) {
      l = aurostd::string2utype<double>(tokens[2]);
    }

    _aflags aflags;
    aflags.QUIET = true;
    const xstructure a(cin, IOVASP_AUTO);
    const double tol = 1e-6;
    double dist;
    //  a = GetStandardConventional(a);

    const xvector<double> origin;
    origin(1) = 0;
    origin(2) = 0;
    origin(3) = 0;

    xvector<double> A(3);
    xvector<double> B(3);
    xvector<double> C(3);

    // NO ZEROS
    if (h > tol && k > tol && l > tol) {
      A = (1 / h) * a.lattice(1);
      B = (1 / k) * a.lattice(2);
      C = (1 / l) * a.lattice(3);
    }
    // ONE ZERO
    if (h < tol && k > tol && l > tol) { // 0 k l
      B = (1 / k) * a.lattice(2);
      A = a.lattice(1) + B;
      C = (1 / l) * a.lattice(3);
    }
    if (h > tol && k < tol && l > tol) { // h 0 l
      A = (1 / h) * a.lattice(1);
      B = a.lattice(2) + A;
      C = (1 / l) * a.lattice(3);
    }
    if (h > tol && k > tol && l < tol) { // h k 0
      A = (1 / h) * a.lattice(1);
      B = (1 / k) * a.lattice(2);
      C = a.lattice(3) + B;
    }
    // TWO ZEROS
    if (h < tol && k < tol && l > tol) { // 0 0 l
      C = (1 / l) * a.lattice(3);
      A = a.lattice(1) + C;
      B = a.lattice(2) + C;
    }
    if (h < tol && k > tol && l < tol) { // 0 k 0
      B = (1 / k) * a.lattice(2);
      A = a.lattice(1) + B;
      C = a.lattice(3) + B;
    }
    if (h > tol && k < tol && l < tol) { // h 0 0
      A = (1 / h) * a.lattice(1);
      B = a.lattice(2) + A;
      C = a.lattice(3) + A;
    }
    // GET PLANE NORMAL VECTOR & NORMAL ALONG A
    const xvector<double> n = 1 / aurostd::modulus(aurostd::vector_product(A - B, A - C)) * aurostd::vector_product(A - B, A - C);
    const xvector<double> nA = 1 / aurostd::modulus(A) * A;
    dist = aurostd::modulus(A) * aurostd::scalar_product(nA, n);
    cout << setprecision(6) << dist << endl;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return dist;
  }
} // namespace pflow

namespace pflow {
  int whereischar(string str, char c) {
    int out = -1;
    for (size_t i = 0; i < str.size(); i++) {
      if (str[i] == c) {
        out = i;
      }
    }
    return out;
  }
} // namespace pflow

namespace pflow {
  bool havechar(string str_in, char c) {
    bool contains = false;
    for (uint i = 0; i < str_in.length(); i++) {
      if (str_in[i] == c) {
        contains = true;
      }
    }
    return contains;
  }
} // namespace pflow

namespace pflow {
  void cleanupstring(string& str) {
    uint i = str.size() - 1;
    while (str[i] == ' ') {
      i--;
    }
    str.erase(i + 1, str.size() - 1);
    i = 0;
    while (str[i] == ' ') {
      i++;
    }
    str.erase(0, i);
  }
} // namespace pflow

namespace pflow {
  double frac2dbl(string str) {
    double out;
    bool neg = false;
    cleanupstring(str);
    if (str[0] == '-') {
      neg = true;
      str.erase(str.begin(), str.begin() + 1);
    }
    if (str[0] == '+') {
      str.erase(str.begin(), str.begin() + 1);
    }
    if (havechar(str, '/')) {
      double numerator;
      double denominator;
      const uint slash = whereischar(str, '/');
      char num[256];
      char den[256];
      for (uint i = 0; i < slash; i++) {
        num[i] = str[i];
      }
      for (size_t i = 0; i < (str.size() - slash); i++) {
        den[i] = str[i + slash + 1];
      }
      numerator = atof(num);
      denominator = atof(den);
      if (neg == true) {
        out = -numerator / denominator;
      } else {
        out = numerator / denominator;
      }
    } else {
      out = atof(str.c_str());
    }
    return out;
  }
} // namespace pflow

#endif //

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2024              *
// *                                                                        *
// **************************************************************************
