// ***************************************************************************
// *                                                                         *
// *         aflow - Automatic FLOW for materials discovery project          *
// *             Stefano Curtarolo - Duke University - 2003-2023             *
// *                                                                         *
// ***************************************************************************
//
//  Copyright 2003-2023 - Stefano Curtarolo - AFLOW.ORG consortium
//
//  This file is part of AFLOW software.
//
//  AFLOW is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ***************************************************************************

#include "aflow.h"

#include "config.h"

#include <cstddef>
#include <cstdio>
#include <deque>
#include <fstream>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <curl/curl.h>

#include <modules/HULL/aflow_nhull.h>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_argv.h"
#include "AUROSTD/aurostd_time.h"
#include "AUROSTD/aurostd_xcomplex.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xrandom.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "aflow_init.h"
#include "aflow_xhost.h"
#include "aflowlib/aflowlib_libraries_scrubber.h"
#include "aflowlib/aflowlib_web_outreach.h"
#include "extern/SQLITE/sqlite3.h"
#include "flow/aflow_ivasp.h"
#include "flow/aflow_pflow.h"
#include "interfaces/aflow_pthreads.h"
#include "modules/POCC/aflow_pocc.h"  //CO20200624
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"
#include "test/aflow_test.h"
#include "test/aflow_unit_test.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::iostream;
using std::istringstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::string;
using std::stringstream;
using std::vector;

using aurostd::xmatrix;
using aurostd::xvector;

int main(int _argc, char** _argv) {
  ostream& oss = cout;  // CO20180419
  try {
    const bool LDEBUG = false; // true;
    int return_code = 0;  // ME20200901
    if (LDEBUG) {
      cerr << "AFLOW-MAIN [1]" << endl;
    }
    std::vector<string> argv(aurostd::get_arguments_from_input(_argc, _argv));
    if (LDEBUG) {
      cerr << "AFLOW-MAIN [2]" << endl;
    }
    std::vector<string> cmds;

    // init curl library (only needs to be called once) HE20250602
    curl_global_init(CURL_GLOBAL_SSL);

    // MACHINE
    // ME20200724
    const int code = init::InitMachine(false, argv, cmds, cerr);
    if (code >= 0) {
      return code;
    }
    if (LDEBUG || XHOST.DEBUG) {
      cerr << "AFLOW-MAIN [3]" << endl;
    }

    // aurostd::TmpDirectoryCreate("test");
    // cerr << args2flag(argv,"--aaa|--bbb |--ccc") << endl;
    // CHECK USERS MACHINES - DEBUG

    // initialize_templates_never_call_this_procedure(1);

    // INITIALIZE ***************************************************
    // INIT LOOK UP TABLES
    atoms_initialize();
    xelement::Initialize();
    // spacegroup::SpaceGroupInitialize(); only if necessary
    // INFORMATION **************************************************
    AFLOW_PTHREADS::FLAG = AFLOW_PTHREADS::Check_Threads(argv, !XHOST.QUIET);

    bool Arun = false;
    if (!Arun && aurostd::args2flag(argv, cmds, "--pocc_old2new|--pocc_o2n")) {
      Arun = true;
      pocc::poccOld2New();
    } // CO20200624
    if (!Arun && aurostd::args2flag(argv, cmds, "--prx|--prx=")) {
      Arun = true;
      PERFORM_PRX(cout);
    }
    if (!Arun && aurostd::args2flag(argv, cmds, "--test_getpp")) {
      if (KBIN::VASP_PseudoPotential_CleanName_TEST()) {
        return 0;
      }
      return 1;
    }

    if (!Arun && aurostd::args2flag(argv, cmds, "--test_xmatrix")) { // CO20190911
      const bool LDEBUG = true;// true;
      xmatrix<double> mat;
      mat(1, 1) = 5;
      mat(1, 2) = 9;
      mat(1, 3) = 12;
      mat(2, 1) = 7;
      mat(2, 2) = 10;
      mat(2, 3) = 13;
      mat(3, 1) = 8;
      mat(3, 2) = 11;
      mat(3, 3) = 14;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " mat=" << endl;
        cerr << mat << endl;
      }
      // getxmat()
      xmatrix<double> submat;
      mat.getxmatInPlace(submat, 2, 3, 2, 3);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " submat=" << endl;
        cerr << submat << endl;
      }
      // setmat()
      mat.setmat(submat, 1, 1); // do nothing
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " replacing with submat at 1,1" << endl;
        cerr << __AFLOW_FUNC__ << " mat=" << endl;
        cerr << mat << endl;
      }
      xvector<double> xv;
      xv(1) = 2;
      xv(2) = 3;
      xv(3) = 4;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " xv=" << xv << endl;
      }
      mat.setmat(xv, 1, false); // row
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " replacing with xv at row=1" << endl;
        cerr << __AFLOW_FUNC__ << " mat=" << endl;
        cerr << mat << endl;
      }
      mat.setmat(xv, 2, true); // col
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " replacing with xv at col=2" << endl;
        cerr << __AFLOW_FUNC__ << " mat=" << endl;
        cerr << mat << endl;
      }
      // setrow()
      mat.setrow(xv, 2);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " replacing with xv at row=2" << endl;
        cerr << __AFLOW_FUNC__ << " mat=" << endl;
        cerr << mat << endl;
      }
      // setcol()
      mat.setcol(xv, 3);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " replacing with xv at col=3" << endl;
        cerr << __AFLOW_FUNC__ << " mat=" << endl;
        cerr << mat << endl;
      }
      return 1;
    }
    if (!Arun && aurostd::args2flag(argv, cmds, "--test_stefano")) {
      uint y = 2017;
      uint m = 11;
      m += 1;
      for (uint i = 0; i < 200; i++) {
        if (m == 0) {
          cout << "mv \"unknown.pdf\" stefano_" << y << m << ".pdf" << endl;
        } else {
          if (m < 10) {
            cout << "mv \"unknown(" << i << ").pdf\" stefano_" << y << "0" << m << ".pdf" << endl;
          } else {
            cout << "mv \"unknown(" << i << ").pdf\" stefano_" << y << m << ".pdf" << endl;
          }
        }
        m--;
        if (m == 0) {
          y--;
          m += 12;
        }
      }
      return 0; // CO20180419
    }

    // ME20220128
    if (!Arun && aurostd::args2attachedflag(argv, cmds, "--unit_test")) {
      aurostd::xoption unittest_options;
      unittest_options.args2addattachedscheme(argv, cmds, "UNIT_TESTS", "--unit_test=", "all");
      vector<string> tests;
      aurostd::string2tokens(unittest_options.getattachedscheme("UNIT_TESTS"), tests, ",");
      unittest::UnitTest ut(std::cout);
      return (ut.runTestSuites(tests) ? 0 : 1);
    }
    if (!Arun && aurostd::args2flag(argv, cmds, "--test")) {
      if (XHOST.vext.size() != XHOST.vcat.size()) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "XHOST.vext.size()!=XHOST.vcat.size(), aborting.", _RUNTIME_ERROR_);
      }

      for (size_t iext = 0; iext < XHOST.vext.size(); iext++) {
        cout << "\"" << XHOST.vext[iext] << "\"" << " " << "\"" << XHOST.vcat[iext] << "\"" << endl;
      }

      const int dim = 7;
      cout << "dimm=" << dim << endl;

      xmatrix<double> m(dim, dim);
      xmatrix<double> mi(dim, dim);
      for (int i = 1; i <= dim; i++) {
        for (int j = 1; j <= dim; j++) {
          m(i, j) = aurostd::ran0();
        }
      }
      cout << "m=" << endl << m << endl;
      mi = inverse(m);
      cout << "mi=" << endl << mi << endl;
      cout << "mi*m=" << endl << det(mi * m) << endl;

      // CO how to create 64bit string from binary file
      // string b64String;
      // aurostd::bin2base64("aflow_logo.pdf",b64String);
      // cout << b64String << endl;

      const string test = "2.730747137  -2.730747137-12.397646334";
      vector<string> _tokens;
      aurostd::string2tokens(test, _tokens, "-");
      for (size_t i = 0; i < _tokens.size(); i++) {
        cerr << _tokens[i] << endl;
      }
      return 0; // CO20180419
      // CO START 20170614 - some SQLITE tests
      // http://zetcode.com/db/sqlitec/ - more tests here
      // this will create test.db file
      sqlite3* db;
      char* err_msg = nullptr;
      int rc = sqlite3_open("test.db", &db);
      if (rc != SQLITE_OK) {
        fprintf(stderr, "Cannot open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return 1;
      }
      const char* sql =
          "DROP TABLE IF EXISTS Cars;"
          "CREATE TABLE Cars(Id INT, Name TEXT, Price INT);"
          "INSERT INTO Cars VALUES(1, 'Audi', 52642);"
          "INSERT INTO Cars VALUES(2, 'Mercedes', 57127);"
          "INSERT INTO Cars VALUES(3, 'Skoda', 9000);"
          "INSERT INTO Cars VALUES(4, 'Volvo', 29000);"
          "INSERT INTO Cars VALUES(5, 'Bentley', 350000);"
          "INSERT INTO Cars VALUES(6, 'Citroen', 21000);"
          "INSERT INTO Cars VALUES(7, 'Hummer', 41400);"
          "INSERT INTO Cars VALUES(8, 'Volkswagen', 21600);";
      rc = sqlite3_exec(db, sql, nullptr, nullptr, &err_msg);
      if (rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);
        sqlite3_close(db);
        return 1;
      }
      sqlite3_close(db);
      // return 0;

      // MORE TESTS
      // printf("%s\n,sqlite3_libversion()");
      //     sqlite3 *db;
      // sqlite3_stmt *res;
      // int rc = sqlite3_open(":memory:", &db);
      // if(rc != SQLITE_OK) {
      //     fprintf(stderr, "Cannot open database: %s\n", sqlite3_errmsg(db));
      //     sqlite3_close(db);
      //     return 1;
      // }
      // rc = sqlite3_prepare_v2(db, "SELECT SQLITE_VERSION()", -1, &res, 0);
      // if(rc != SQLITE_OK) {
      //     fprintf(stderr, "Failed to fetch data: %s\n", sqlite3_errmsg(db));
      //     sqlite3_close(db);
      //     return 1;
      // }
      // rc = sqlite3_step(res);
      // if(rc == SQLITE_ROW) {
      //     printf("%s\n", sqlite3_column_text(res, 0));
      // }
      // sqlite3_finalize(res);
      // sqlite3_close(db);
      // return 0;

      // quick easy test
      cerr << sqlite3_libversion() << endl;
      // CO END 20170614 - some SQLITE tests
      aurostd::xcomplex<double> x(123.0, 456.0);
      cout << x.re << "," << x.im << " - " << x.real() << "," << x.imag() << " - " << x << endl;
      x.re = 111;
      x.im = 222;
      cout << x.re << "," << x.im << " - " << x.real() << "," << x.imag() << " - " << x << endl;
      cout << aurostd::PaddedPOST("EMIN= -30.0", 10) << endl;
      ;
      stringstream for_corey;
      for_corey << "scatter/use mapped color={draw=black,fill=mapped color,solid}";
      const string corey = for_corey.str();
      cout << corey << endl;
      stringstream aus;
      aus << "************************   00000  MESSAGE KPOINTS KSHIFT=[" << 1 << " " << 2 << " " << 3 << "]" << " ************************ " << endl;
      cout << aus.str() << endl;
      return 0; // CO20180419
    }

    if (!Arun && aurostd::args2attachedflag(argv, "--bin2base64=")) {
      string b64String;
      aurostd::bin2base64(aurostd::args2attachedstring(argv, "--bin2base64=", ""), b64String);
      cout << b64String << endl;
      return 0; // CO20180419
    }

    if (!Arun && aurostd::args2flag(argv, cmds,
                                    "--test=POTCAR|--test=POTCAR.relax1" + DEFAULT_KZIP_EXT + "|--test=POTCAR.relax2" + DEFAULT_KZIP_EXT + "|--test=POTCAR.static" + DEFAULT_KZIP_EXT + "|--test=POTCAR.bands" +
                                        DEFAULT_KZIP_EXT + "")) {
      XHOST.DEBUG = true;
      xPOTCAR(aurostd::args2attachedstring(argv, "--test=", ""));
      return 0;
    } // CO20180419
    if (!Arun && aurostd::args2flag(argv, cmds,
                                    "--test=DOSCAR|--test=DOSCAR.relax1" + DEFAULT_KZIP_EXT + "|--test=DOSCAR.relax2" + DEFAULT_KZIP_EXT + "|--test=DOSCAR.static" + DEFAULT_KZIP_EXT + "|--test=DOSCAR.bands" +
                                        DEFAULT_KZIP_EXT + "")) {
      XHOST.DEBUG = true;
      xDOSCAR(aurostd::args2attachedstring(argv, "--test=", ""));
      return 0;
    } // CO20180419
    if (!Arun && aurostd::args2flag(argv, cmds,
                                    "--test=EIGENVAL|--test=EIGENVAL.relax1" + DEFAULT_KZIP_EXT + "|--test=EIGENVAL.relax2" + DEFAULT_KZIP_EXT + "|--test=EIGENVAL.static" + DEFAULT_KZIP_EXT +
                                        "|--test=EIGENVAL.bands" + DEFAULT_KZIP_EXT + "")) {
      XHOST.DEBUG = true;
      xEIGENVAL(aurostd::args2attachedstring(argv, "--test=", ""));
      return 0;
    } // CO20180419
    if (!Arun && aurostd::args2flag(argv, cmds,
                                    "--test=OUTCAR|--test=OUTCAR.relax1" + DEFAULT_KZIP_EXT + "|--test=OUTCAR.relax2" + DEFAULT_KZIP_EXT + "|--test=OUTCAR.static" + DEFAULT_KZIP_EXT + "|--test=OUTCAR.bands" +
                                        DEFAULT_KZIP_EXT + "")) {
      XHOST.DEBUG = true;
      xOUTCAR(aurostd::args2attachedstring(argv, "--test=", ""));
      return 0;
    } // CO20180419
    if (!Arun && aurostd::args2flag(argv, cmds,
                                    "--test=KPOINTS|--test=KPOINTS.relax1" + DEFAULT_KZIP_EXT + "|--test=KPOINTS.relax2" + DEFAULT_KZIP_EXT + "|--test=KPOINTS.static" + DEFAULT_KZIP_EXT +
                                        "|--test=KPOINTS.bands" + DEFAULT_KZIP_EXT + "")) {
      XHOST.DEBUG = true;
      xKPOINTS(aurostd::args2attachedstring(argv, "--test=", ""));
      return 0;
    }  // CO20180419
    if (!Arun && aurostd::args2flag(argv, cmds,
                                    "--test=vasprun|--test=vasprun.xml.relax1" + DEFAULT_KZIP_EXT + "|--test=vasprun.xml.relax2" + DEFAULT_KZIP_EXT + "|--test=vasprun.xml.static" + DEFAULT_KZIP_EXT +
                                        "|--test=vasprun.xml.bands" + DEFAULT_KZIP_EXT + "")) {
      XHOST.DEBUG = true;
      xVASPRUNXML(aurostd::args2attachedstring(argv, "--test=", ""));
      return 0;
    } // CO20180419
    if (!Arun && aurostd::args2flag(argv, cmds,
                                    "--test=IBZKPT|--test=IBZKPT.relax1" + DEFAULT_KZIP_EXT + "|--test=IBZKPT.relax2" + DEFAULT_KZIP_EXT + "|--test=IBZKPT.static" + DEFAULT_KZIP_EXT + "|--test=IBZKPT.bands" +
                                        DEFAULT_KZIP_EXT + "")) {
      XHOST.DEBUG = true;
      xIBZKPT(aurostd::args2attachedstring(argv, "--test=", ""));
      return 0;
    } // CO20180419
    if (!Arun && aurostd::args2flag(argv, cmds,
                                    "--test=CHGCAR|--test=CHGCAR.relax1" + DEFAULT_KZIP_EXT + "|--test=CHGCAR.relax2" + DEFAULT_KZIP_EXT + "|--test=CHGCAR.static" + DEFAULT_KZIP_EXT + "|--test=CHGCAR.bands" +
                                        DEFAULT_KZIP_EXT + "")) {
      XHOST.DEBUG = true;
      xCHGCAR(aurostd::args2attachedstring(argv, "--test=", ""));
      return 0;
    } // CO20180419

    // SCRUB things
    if (!Arun && aurostd::args2flag(argv, cmds, "--scrub=POTCAR")) {
      XHOST.DEBUG = false;
      XHOST.PSEUDOPOTENTIAL_GENERATOR = true;
      vector<string> vfile(aurostd::args2vectorstring(argv, "--FILE|--file|--F|--f", "./"));
      for (size_t ifile = 0; ifile < vfile.size(); ifile++) {
        cerr << "PROCESSING = " << vfile[ifile] << endl;
        const xPOTCAR xPOT(vfile[ifile]);
      }
      return 0; // CO20180419
    }
    if (!Arun && aurostd::args2flag(argv, cmds, "--scrub=OUTCAR")) {
      XHOST.DEBUG = false;
      XHOST.PSEUDOPOTENTIAL_GENERATOR = false;
      vector<string> vfile(aurostd::args2vectorstring(argv, "--FILE|--file|--F|--f", "./"));
      for (size_t ifile = 0; ifile < vfile.size(); ifile++) {
        cerr << "PROCESSING = " << vfile[ifile] << endl;
        const xOUTCAR xOUT(vfile[ifile]);
        // cout << xOUT << endl;
      }
      return 0; // CO20180419
    }

    if (!Arun && (aurostd::args2flag(argv, cmds, "--scrub") || aurostd::args2attachedflag(argv, "--scrub="))) {
      //  XHOST.DEBUG=true;
      aflowlib::LIB2SCRUB(aurostd::args2attachedstring(argv, "--scrub=", "ALL"), true);
      return 0; // CO20180419
    }
    if (!Arun && (aurostd::args2attachedflag(argv, "--lib2auid="))) {
      //  XHOST.DEBUG=true;
      aflowlib::LIB2AUID(aurostd::args2attachedstring(argv, "--lib2auid=", "ALL"), false, true); // no test, act
      return 0; // CO20180419
    }
    if (!Arun && (aurostd::args2flag(argv, cmds, "--mosfet") || aurostd::args2attachedflag(argv, "--mosfet="))) {
      //  XHOST.DEBUG=true;
      aflowlib::MOSFET(aurostd::args2attachedutype<int>(argv, "--mosfet=", 0), true);
      return 0; // CO20180419
    }
    if (!Arun && (aurostd::args2flag(argv, cmds, "--mail2scan") || aurostd::args2attachedflag(argv, "--mail2scan="))) {
      //  XHOST.DEBUG=true;
      aflowlib::MAIL2SCAN(aurostd::args2attachedstring(argv, "--mail2scan=", "/var/mail/auro"), true);
      return 0; // CO20180419
    }

    if (!Arun && aurostd::args2flag(argv, cmds, "--test_proto1")) {
      vector<xstructure> vstr;
      aflowlib::_aflowlib_entry data;
      vector<aflowlib::_aflowlib_entry> vdata;
      for (size_t i = 4; i < BRAVAIS_LATTICES.size(); i++) { // HE20220420 switch to global bravais lattices list
        data.clear();
        data.url2aflowlib("materials.duke.edu:AFLOWDATA/ICSD_WEB/" + BRAVAIS_LATTICES[i] + "/?format=text", cout, false);
        vdata.push_back(data);
        cout << "AFLOWLIB " << BRAVAIS_LATTICES[i] << "=" << data.vaflowlib_entries.size() << endl;
        for (size_t j = 0; j < data.vaflowlib_entries.size(); j++) {
          aflowlib::_aflowlib_entry dataj;
          dataj.url2aflowlib("materials.duke.edu:AFLOWDATA/ICSD_WEB/" + BRAVAIS_LATTICES[i] + "/" + data.vaflowlib_entries[j], cout, true);
          aurostd::StringSubstInPlace(dataj.aurl, "aflowlib", "materials");
          if (!dataj.aurl.empty()) {
            const xstructure str(dataj.aurl, "CONTCAR.relax.vasp", IOAFLOW_AUTO);
            xEIGENVAL xEIGENVAL;
            xEIGENVAL.GetPropertiesUrlFile(dataj.aurl, "EIGENVAL.bands" + DEFAULT_KZIP_EXT + "", false);
            xOUTCAR xOUTCAR;
            xOUTCAR.GetPropertiesUrlFile(dataj.aurl, "OUTCAR.static" + DEFAULT_KZIP_EXT + "", false);
            xDOSCAR xDOSCAR;
            xDOSCAR.GetPropertiesUrlFile(dataj.aurl, "DOSCAR.static" + DEFAULT_KZIP_EXT + "", false);
            // if(aurostd::args2flag(argv,cmds,"--vasp")) aurostd::url2stringstream(dataj.aurl+"/CONTCAR.relax.vasp",stream);
            // if(aurostd::args2flag(argv,cmds,"--qe")) aurostd::url2stringstream(dataj.aurl+"/CONTCAR.relax.qe",stream);
            // if(aurostd::args2flag(argv,cmds,"--abinit")) aurostd::url2stringstream(dataj.aurl+"/CONTCAR.relax.abinit",stream);
            // if(aurostd::args2flag(argv,cmds,"--aims")) aurostd::url2stringstream(dataj.aurl+"/CONTCAR.relax.aims",stream);
            vstr.push_back(str);
            cerr << "vstr.size()=" << vstr.size() << "  "
                 << "str.atoms.size()=" << str.atoms.size() << "  "
                 << "OUTCAR.static" + DEFAULT_KZIP_EXT + ".size()=" << xOUTCAR.vcontent.size() << "  "
                 << "DOSCAR.static" + DEFAULT_KZIP_EXT + ".size()=" << xDOSCAR.vcontent.size() << "  "
                 << "EIGENVAL.static" + DEFAULT_KZIP_EXT + ".size()=" << xEIGENVAL.vcontent.size() << "  " << endl;
            //	  cerr << str << endl;
          }
        }
      }
      return 0; // CO20180419
    }

    if (!Arun && aurostd::args2flag(argv, cmds, "--testJ")) {
      Arun = true;
      PERFORM_TESTJ(cout);
    }
    if (!Arun && aurostd::args2flag(argv, cmds, "--test3")) {
      Arun = true;
      PERFORM_TEST3(cout);
    }
    if (!Arun && XHOST.vflag_control.flag("MACHINE")) {
      // ME20200724
      const int code = init::InitMachine(false, argv, cmds, cerr);
      if (code >= 0) {
        return code;
      }
    }

    // **************************************************************
    // INTERCEPT AFLOW
    if (!Arun && XHOST.vflag_control.flag("SWITCH_AFLOW")) {
      Arun = true;
      AFLOW_main(argv);
    }
    // DX
    if (!Arun && XHOST.vflag_control.flag("AFLOWIN_SYM")) {
      Arun = true;
      AFLOW_main(argv);
    }
    // DX
    if (!Arun && (XHOST.vflag_aflow.flag("CLEAN") || XHOST.vflag_aflow.flag("XCLEAN") || XHOST.AFLOW_RUNDIRflag || XHOST.AFLOW_MULTIflag || XHOST.AFLOW_RUNXflag)) {
      Arun = true;
      AFLOW_main(argv);
    }

    // **************************************************************
    // INTERCEPT aconvasp/aqe by title
    if (!Arun && aurostd::substring2bool(XHOST.progname, "aconvasp", "convasp")) {
      Arun = true;
      pflow::main(argv, cmds);
    }
    if (!Arun && aurostd::substring2bool(XHOST.progname, "aqe")) {
      Arun = true;
      pflow::main(argv, cmds);
    }

    // **************************************************************
    // intercept commands
    if (!Arun && XHOST.vflag_control.flag("MULTI=SH")) {
      Arun = true;
      AFLOW_PTHREADS::MULTI_sh(argv);
      return 0;
    }  // CO20180419
    if (!Arun && XHOST.vflag_control.flag("MONITOR")) {
      Arun = true;
      AFLOW_monitor(argv);
      return 0;
    } // CO20180419
    if (!Arun && XHOST.vflag_control.flag("MONITOR_VASP")) {
      Arun = true;
      AFLOW_monitor_VASP();
      return 0;
    } // CO20180419
    if (!Arun && XHOST.vflag_control.flag("GETTEMP")) {
      Arun = true;
      AFLOW_getTEMP(argv);
      return 0;
    } // CO20180419

    // **************************************************************
    // INTERCEPT HELP
    // ME20200921 - Restructured to make web processing easier
    stringstream banner_message;
    if (XHOST.vflag_control.flag("AFLOW_HELP")) {
      banner_message << aflow::Banner("BANNER_BIG") << endl << aflow::Intro_HELP("aflow") << aflow::Banner("BANNER_BIG") << endl;
    } else if (XHOST.vflag_control.flag("AFLOW_EXCEPTIONS")) {
      banner_message << aflow::Banner("BANNER_BIG") << endl << aflow::Banner("EXCEPTIONS") << endl;
    } else if (XHOST.vflag_control.flag("README_AFLOW_LICENSE_GPL3")) {
      banner_message << aflow::License_Preamble_aflow() << endl;
      banner_message << " " << endl;
      banner_message << aurostd::EmbData::get_content("README_AFLOW_LICENSE_GPL3.TXT", "README") << endl;
      banner_message << " " << endl;
      banner_message << "*************************************************************************** " << endl;
    } else if (XHOST.vflag_control.flag("README_AFLOW")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_AFLOW_VERSIONS_HISTORY")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_VERSIONS_HISTORY.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_AFLOW_PFLOW")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_PFLOW.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_FROZSL")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_FROZSL.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_APL")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_APL.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_QHA")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_APL.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_AAPL")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_APL.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_AGL")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_AGL.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_AEL")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_AEL.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_ANRL")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_ANRL.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_COMPARE")) { // CO20190401
      banner_message << aurostd::EmbData::get_content("README_AFLOW_COMPARE.TXT", "README") << endl; // CO20190401
    } else if (XHOST.vflag_control.flag("README_GFA")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_GFA.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_SYMMETRY")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_SYM.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_CCE")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_CCE.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_NHULL")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_NHULL.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_PARTIAL_OCCUPATION")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_POCC.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_SCRIPTING")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_SCRIPTING.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_EXCEPTIONS")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_EXCEPTIONS.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_XAFLOW")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_XAFLOW.TXT", "README") << endl;
    } else if (XHOST.vflag_control.flag("README_AFLOWRC")) {
      banner_message << aurostd::EmbData::get_content("README_AFLOW_AFLOWRC.TXT", "README") << endl;
    }

    if (!banner_message.str().empty()) {
      std::cout << (XHOST.vflag_control.flag("WWW") ? aurostd::text2html(banner_message.str()) : banner_message.str()) << std::endl;
      return 0;
    }

    // **************************************************************
    // PHP-WEB AND CURRICULUM AND HIGH-THROUGHPUT STUFF
    if (ProcessPhpLatexCv()) {
      return 0;
    }
    //  ProcessSecurityOptions(argv,cmds); OLD STUFF AFLOW SECURITY

    // **************************************************************
    const bool VVERSION = aurostd::args2flag(argv, cmds, "-v|--version");
    // ME20200921 - Added web mode
    if (!Arun && VVERSION) {
      // look for version IMMEDIATELY //CO20180419
      Arun = true;
      cout << (XHOST.vflag_control.flag("WWW") ? aurostd::text2html(aflow::Banner("AFLOW_VERSION")) : aflow::Banner("AFLOW_VERSION"));
      return 0;
    }
    if (!Arun && XHOST.TEST) {
      Arun = true;
      cerr << "test" << endl;
      return 0;
    } // CO20180419

    if (!Arun && XHOST.argv.size() == 1 && (aurostd::substring2bool(XHOST.progname, "aflow") || aurostd::substring2bool(XHOST.progname, "aflowd"))) {
      //   Arun=true;AFLOW_main(argv);
      Arun = true;
      //    cout << "******************************************************************************************************" << endl;
      //  cout << aflow::Banner("BANNER_TINY") << endl;
      cout << aflow::Banner("BANNER_BIG") << endl;
      cout << aflow::Intro_aflow("aflow") << endl;
      cout << pflow::Intro_pflow("aflow") << endl;
      cout << aflow::Intro_sflow("aflow") << endl;
      cout << aflow::Intro_HELP("aflow") << endl;
      cout << aflow::Banner("BANNER_BIG") << endl;
      //    cout << "XHOST.argv.size()=" << XHOST.argv.size()<< endl;
      //    cout << "******************************************************************************************************" << endl;
    }

    // **************************************************************
    // LAST RESOURCE PFLOW
    if (!Arun) {
      Arun = true;
      return_code = pflow::main(argv, cmds);   // ME20200901 - use pflow::main return code for database handling
    }
    // **************************************************************
    // END
    curl_global_cleanup();
    return (Arun ? return_code : 1); // Arun==true is 1, so flip because return 0 is normal  //CO20190629 - more explicit return 0//ME20200901 - use return_code
  } catch (aurostd::xerror& excpt) {
    pflow::logger(excpt.whereFileName(), excpt.whereFunction(), excpt.buildMessageString(), oss, _LOGGER_ERROR_);
    return excpt.whatCode();
  }
}

// ***************************************************************************
// AFLOW_main
// ***************************************************************************
int AFLOW_main(vector<string>& argv) {
  if (!XHOST.QUIET) {
    cout << aflow::Banner("INTRODUCTION");// << endl;
  }
  KBIN::KBIN_Main(argv);
  // if(!XHOST.QUIET) cout << "MMMMM  AFLOW VERSION " << string(AFLOW_VERSION) << "  " << endl;
  return 0; // 1; //CO20180419 - return 0 is normal
}

// ***************************************************************************
// aflow::License_aflow
// ***************************************************************************
namespace aflow {
  string License_Preamble_aflow() {
    //(C) 2003-2023 Stefano Curtarolo, MIT-Duke University stefano@duke.edu
    stringstream strstream;
    strstream << endl;
    strstream << "***************************************************************************" << endl;
    strstream << "*                                                                         *" << endl;
    strstream << "*                    AFLOW - Duke University 2003-2024                    *" << endl; // CO20200502 - SC -> AFLOW consortium
    strstream << "*                                                                         *" << endl;
    strstream << "***************************************************************************" << endl;
    strstream << "Copyright 2003-2023 - AFLOW.ORG consortium" << endl;  // CO20200502 - SC -> AFLOW consortium
    strstream << endl;
    strstream << "This file is part of AFLOW software." << endl;
    strstream << endl;
    strstream << "AFLOW is free software: you can redistribute it and/or modify" << endl;
    strstream << "it under the terms of the GNU General Public License as published by" << endl;
    strstream << "the Free Software Foundation, either version 3 of the License, or" << endl;
    strstream << "(at your option) any later version." << endl;
    strstream << endl;
    strstream << "This program is distributed in the hope that it will be useful," << endl;
    strstream << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
    strstream << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
    strstream << "GNU General Public License for more details." << endl;
    strstream << endl;
    strstream << "You should have received a copy of the GNU General Public License" << endl;
    strstream << "along with this program.  If not, see <http://www.gnu.org/licenses/>." << endl;
    strstream << endl;
    strstream << "***************************************************************************" << endl;
    strstream << endl;
    return strstream.str();
  }
} // namespace aflow

// ***************************************************************************
// aflow::Intro_aflow
// ***************************************************************************
// patched by CO20200106 to avoid long string construction (patching for indents)
namespace aflow {
  string Intro_aflow(string x) {
    //(C) 2003-2023 Stefano Curtarolo, MIT-Duke University stefano@duke.edu
    stringstream strstream;
    const string tab = "  ";
    string xspaces;
    for (size_t i = 0; i < x.size(); i++) {
      xspaces += " ";
    } // spaces size of x
    strstream << endl;
    strstream << "******* BEGIN INFORMATION MODE *********************************************************************" << endl;
    strstream << tab << x << " -h|--help|--readme_aflow   CHECK" << endl;
    strstream << tab << x << " --version|-v                 CHECK" << endl;
    strstream << tab << x << " --machine                      CHECK" << endl;
    strstream << "******* END INFORMATION MODE ***********************************************************************" << endl;
    strstream << endl;
    strstream << "******* BEGIN RUNNING MODE *************************************************************************" << endl;
    strstream << tab << x << " --run|--run=multi|--run=N   COMMANDS for running" << endl;
    strstream << tab << x << " --clean|-c                    COMMAND for cleaning  " << endl;
    strstream << endl;
    strstream << " MODIFIERS" << endl;
    strstream << tab << " --DIRECTORY[=| ]dir|--D[=| ]dir|--d[=| ]dir" << endl;
    strstream << tab << " --FILE[=| ]file|--F[=| ]file|--f[=| ]file " << endl;
    strstream << tab << " --quiet|-quiet|-q " << endl;
    strstream << tab << " --loop " << endl;
    strstream << tab << " --sort|-sort" << endl;
    strstream << tab << " --reverse|-rsort" << endl;
    strstream << tab << " --random|-rnd" << endl;
    strstream << tab << " --force|-force" << endl;
    strstream << tab << " --mem=XX|--maxmem=XX" << endl;
    strstream << tab << " --readme= xaflow|aflow|aconvasp|aflowrc|scripting|apl|agl|ael|anrl|compare|gfa|symmetry|nhull|errors|exceptions|frozsl CHECK !!!!" << endl;
    strstream << tab << " --np=NUMBER|--npmax" << endl;
    strstream << tab << " --generate_aflowin_from_vasp" << endl;
    strstream << tab << " --generate_vasp_from_aflowin|--generate" << endl;
    strstream << tab << " --use_aflow.in=XXX" << endl;
    strstream << tab << " --use_LOCK=XXX" << endl;
    strstream << tab << " --use_tmpfs=XXXX " << endl;
    strstream << endl;
    strstream << " MODIFIERS MPI/SERIAL PARAMETERS" << endl;
    strstream << tab << " --mpi|-nompi|--serial " << endl;
    strstream << endl;
    strstream << " HOST ORIENTED OPTION" << endl;
    strstream << tab << " --machine=beta|beta_openmpi|qrats|qflow|x|conrad|eos|materials|habana|aflowlib|ranger|kraken" << endl;
    strstream << tab << "           marylou|parsons|jellium|ohad|host1" << endl;
    strstream << tab << "           raptor --np=N|diamond --np=N" << endl;
    strstream << tab << " --machine_name=XXXX " << endl;
    strstream << "******* END RUNNING MODE ***************************************************************************" << endl;
    strstream << endl;
    // --readme=htresources
    return strstream.str();
  }
} // namespace aflow

// ***************************************************************************
// aflow::Intro_sflow
// ***************************************************************************
// patched by CO20200106 to avoid long string construction (patching for indents)
namespace aflow {
  string Intro_sflow(string x) {
    stringstream strstream;
    const string tab = "  ";
    // string xspaces="";for(size_t i=0;i<x.size();i++){xspaces+=" ";} //spaces size of x
    strstream << "******* BEGIN SCRIPTING MODE ***********************************************************************" << endl;
    strstream << " AFLOW SCRIPTING COMMANDS" << endl;
    strstream << tab << x << " --justafter=string" << endl;
    strstream << tab << x << " --justbefore=string" << endl;
    strstream << tab << x << " --justbetween=string_from[,string_to]" << endl;
    strstream << tab << x << " --qsub=N,file" << endl;
    strstream << tab << x << " --qdel=aaa,nnn:mmm,aaa,bbb,ccc" << endl;
    strstream << tab << x << " --bsub=N,file" << endl;
    strstream << tab << x << " --bkill=aaa,nnn:mmm,aaa,bbb,ccc" << endl;
    strstream << tab << x << " --sbatch=N,file" << endl;
    strstream << tab << x << " --scancel=aaa,nnn:mmm,aaa,bbb,ccc" << endl;
    strstream << tab << x << " --kill=aaa,nnn:mmm,aaa,bbb,ccc" << endl;
    strstream << tab << x << " --multi=sh [--np=NUMBER|npmax|nothing] [--F[ILE]] file" << endl;
    strstream << tab << x << " --getTEMP [--runstat|--runbar|--refresh=X|--warning_beep=T|--warning_halt=T|--mem=XX]" << endl;
    strstream << tab << x << " --monitor [--mem=XX]" << endl;
    strstream << "******* END SCRIPTING MODE *************************************************************************" << endl;
    strstream << endl;
    return strstream.str();
  }
} // namespace aflow

// ***************************************************************************
// aflow::Intro_HELP
// ***************************************************************************
// patched by CO20200106 to avoid long string construction (patching for indents)
namespace aflow {
  string Intro_HELP(string x) {
    stringstream strstream;
    const string tab = "  ";
    string xspaces;
    for (size_t i = 0; i < x.size(); i++) {
      xspaces += " ";
    } // xspacess size of x
    strstream << "******* BEGIN HELP MODE ****************************************************************************" << endl;
    strstream << " AFLOW HELP AVAILABLE HELPS" << endl;
    strstream << tab << x << " --license" << endl;
    strstream << tab << xspaces << " " << tab << "License information." << endl;
    strstream << tab << x << " --help" << endl;
    strstream << tab << xspaces << " " << tab << "This help." << endl;
    strstream << tab << x << " --readme" << endl;
    strstream << tab << xspaces << " " << tab << "The list of all the commands available." << endl;
    strstream << tab << x << " --readme=xaflow" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the installation of aflow." << endl;
    strstream << tab << x << " --readme=aflow|--readme=run|--readme_aflow" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"running machinery\"." << endl;
    strstream << tab << x << " --readme=aflowrc" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the installation of aflow." << endl;
    strstream << tab << x << " --readme=pflow|--readme=processor|--readme=aconvasp|--readme_aconvasp" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"processing machinery\"." << endl;
    strstream << tab << x << " --readme=scripting|--readme_scripting" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"scripting\" operations." << endl;
    strstream << tab << x << " --readme=apl|--readme_apl" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"aflow-harmonic-phonon-library\"." << endl;
    strstream << tab << x << " --readme=qha|--readme_qha|--readme=qha3p|--readme_qha3p|--readme=scqha|--readme_scqha" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"aflow-quasi-harmonic-library\"." << endl;
    strstream << tab << x << " --readme=aapl|--readme_aapl" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"aflow-anharmonic-phonon-library (AFLOW-AAPL)\"." << endl;
    strstream << tab << x << " --readme=agl|--readme_agl" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"aflow-gibbs-library (AFLOW-AGL)\"." << endl;
    strstream << tab << x << " --readme=ael|--readme_ael" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"aflow-elastic-library (AFLOW-AEL)\"." << endl;
    strstream << tab << x << " --readme=prototypes|--readme_prototypes|--readme=anrl|--readme_anrl" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"aflow library of prototypes\"." << endl;
    strstream << tab << x << " --readme=xtalfinder|--readme_xtalfinder|--readme=compare|--readme_compare" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"aflow-crystal-finder (AFLOW-XtalFinder) code\"." << endl;
    strstream << tab << x << " --readme=gfa|--readme_gfa" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"glass-forming-ability code\"." << endl;
    strstream << tab << x << " --readme=symmetry|--readme_symmetry" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"symmetry library (AFLOW-SYM)\"." << endl;
    strstream << tab << x << " --readme=nhull|--readme_nhull" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"convex hull library (AFLOW-hull)\"." << endl;
    strstream << tab << x << " --readme=errors|--readme=exceptions|--readme_errors|--readme_exceptions" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for exception handling in AFLOW." << endl;
    strstream << tab << x << " --readme=partial_occupation|--readme=pocc|--readme_pocc" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"partial occupation library\"." << endl;
    strstream << tab << x << " --readme=frozsl|--readme_frozsl" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"frozsl\" add ons." << endl;
    strstream << "******* END HELP MODE ******************************************************************************" << endl;
    strstream << endl;
    return strstream.str();
  }
} // namespace aflow

// ***************************************************************************
// aflow::Banner
// ***************************************************************************
namespace aflow {
  string Banner(string type) {
    stringstream oss;
    if (type == "VERSION" || type == "AFLOW_VERSION") {
      oss << "" << endl;
      oss << "AFLOW VERSION " << string(AFLOW_VERSION) << " Automatic-FLOW [(C) " << XHOST.Copyright_Years << " aflow.org consortium]" << endl;
      oss << "New versions are available here: <http://" << XHOST.AFLOW_MATERIALS_SERVER << "/AFLOW/>" << endl;
      oss << "" << endl;
      oss << "AFLOW is free software: you can redistribute it and/or modify it under the terms of the" << endl;
      oss << "GNU General Public License as published by the Free Software Foundation, either version 3" << endl;
      oss << "of the License, or (at your option) any later version." << endl;
      oss << "" << endl;
      oss << "This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;" << endl;
      oss << "without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl;
      oss << "See the GNU General Public License for more details." << endl;
      oss << "" << endl;
      oss << "You should have received a copy of the GNU General Public License along with this program." << endl;
      oss << "If not, see <http://www.gnu.org/licenses/>." << endl;
      oss << "" << endl;
      oss << "AFLOW V" << string(AFLOW_VERSION) << " [" << XHOST.hostname << "] [" << XHOST.machine_type << "] [" << XHOST.CPU_Cores << "] [" << XHOST.Find_Parameters << "]" << endl;
      //     oss << endl;
      return oss.str();
    }
    if (type == "INTRODUCTION") {
      stringstream oss;
      oss << "MMMMM  AFLOW VERSION " << string(AFLOW_VERSION) << " Automatic-FLOW for materials discovery - [" << TODAY << "] -" << aflow_get_time_string() << endl;
      //    oss << "MMMMM  AFLOW VERSION " << aurostd::PaddedPOST(string(AFLOW_VERSION),5) <<" - BUILT ["<<TODAY<<"] - (C) " << XHOST.Copyright_Years << "    " << endl;
      oss << "MMMMM  AFLOW.org consortium - High-Throughput ab-initio Computing Project - (C) " << XHOST.Copyright_Years << "    " << endl;
      oss << "MMMMM  ";// << endl;
      oss << "MMMMM  AFLOW is free software: you can redistribute it and/or modify it under the terms of the" << endl;
      oss << "MMMMM  GNU General Public License as published by the Free Software Foundation, either version 3" << endl;
      oss << "MMMMM  of the License, or (at your option) any later version." << endl;
      oss << "MMMMM  " << endl;
      oss << "MMMMM  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;" << endl;
      oss << "MMMMM  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl;
      oss << "MMMMM  See the GNU General Public License for more details." << endl;
      oss << "MMMMM  " << endl;
      oss << "MMMMM  You should have received a copy of the GNU General Public License along with this program." << endl;
      oss << "MMMMM  If not, see <http://www.gnu.org/licenses/>." << endl;
      oss << "MMMMM " << endl;
      return oss.str();
    }
    if (type == "BANNER_NORMAL") {
      oss << "****************************************************************************************************" << endl;
      oss << "MMMMM  AFLOW V" << string(AFLOW_VERSION) << " Automatic-FLOW [" << TODAY << "] - " << endl; // << aflow_get_time_string() << endl; //CO
      oss << "****************************************************************************************************" << endl;
      return oss.str();
    }
    if (type == "BANNER_BIG") {
      oss << "****************************************************************************************************" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*                          aflow - Automatic-FLOW for materials discovery                          *" << endl;
      oss << "*                aflow.org consortium - High-Throughput ab-initio Computing Project                *" << endl;
      oss << "*"
          << aurostd::PaddedCENTER(string("version " + string(AFLOW_VERSION) + " - g++/gcc " + aurostd::utype2string(__GNUC__) + "." + aurostd::utype2string(__GNUC_MINOR__) + "." +
                                          aurostd::utype2string(__GNUC_PATCHLEVEL__) + " - built [" + string(TODAY) + "] - (C) " + XHOST.Copyright_Years),
                                   100)
          << "*" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "****************************************************************************************************" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*    AFLOW is free software: you can redistribute it and/or modify it under the terms of the       *" << endl;
      oss << "*    GNU General Public License as published by the Free Software Foundation, either version 3     *" << endl;
      oss << "*    of the License, or (at your option) any later version.                                        *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;     *" << endl;
      oss << "*    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.     *" << endl;
      oss << "*    See the GNU General Public License for more details.                                          *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*    You should have received a copy of the GNU General Public License along with this program.    *" << endl;
      oss << "*    If not, see <http://www.gnu.org/licenses/>.                                                   *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "****************************************************************************************************" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*     Use of AFLOW software and repositories welcomes references to the following publications:    *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*  Divilov et al.   High Entropy Alloys & Mater. 3, 178 (2025) 10.1007/s44210-025-00058-2 (AFLOW4) *" << endl;
      oss << "*  Eckert et al.    npj Comput. Mater. 11, 40 (2025)   10.1038/s41524-025-01529-1      (Soliquidy) *" << endl;
      oss << "*  Eckert et al.    Comp. Mat. Sci. 240, 112988 (2024) 10.1016/j.commatsci.2024.112988 (proto4)    *" << endl;
      oss << "*  Divilov et al.   Acta Mater. 266, 119667 (2024)     10.1016/j.actamat.2024.119667   (spinodal)  *" << endl;
      oss << "*  Divilov et al.   Nature 625, 66-73 (2024)           10.1038/s41586-023-06786-y      (DEED)      *" << endl;
      oss << "*  Friedrich et al. J. Chem. Phys. 160, 042501 (2024)  10.1063/5.0184917               (CCE2)      *" << endl;
      oss << "*  Esters et al.    Comp. Mat. Sci. 216, 111808 (2023) 10.1016/j.commatsci.2022.111808 (AFLOW.org) *" << endl;
      oss << "*  Oses et al.      Comp. Mat. Sci. 217, 111889 (2023) 10.1016/j.commatsci.2022.111889 (aflow++)   *" << endl;
      oss << "*  Friedrich et al. npj Comput. Mater. 5, 59 (2019)  10.1038/s41524-019-0192-1       (CCE)         *" << endl;
      oss << "*  Hicks et al.     Comp. Mat. Sci. 161, S1 (2019)   10.1016/j.commatsci.2018.10.043 (ANRL proto2) *" << endl;
      oss << "*  Oses et al.      J. Chem. Inf. Model. (2018)      10.1021/acs.jcim.8b00393        (AFLOW-CHULL) *" << endl;
      oss << "*  Gossett et al.   Comp. Mat. Sci. 152, 134 (2018)  10.1016/j.commatsci.2018.03.075 (AFLOW-ML)    *" << endl;
      oss << "*  Hicks et al.     Acta Cryst. A74, 184-203 (2018)  10.1107/S2053273318003066       (AFLOW-SYM)   *" << endl;
      oss << "*  MBNardelli et al Comp. Mat. Sci. 143, 462 (2018)  10.1016/j.commatsci.2017.11.034 (PAOFLOW)     *" << endl;
      oss << "*  Rose et al.      Comp. Mat. Sci. 137, 362 (2017)  10.1016/j.commatsci.2017.04.036 (AFLUX lang)  *" << endl;
      oss << "*  Supka et al.     Comp. Mat. Sci. 136, 76 (2017)   10.1016/j.commatsci.2017.03.055 (AFLOWpi)     *" << endl;
      oss << "*  Plata et al.     npj Comput. Mater. 3, 45 (2017)  10.1038/s41524-017-0046-7       (AAPL kappa)  *" << endl;
      oss << "*  Toher et al.     Phys. Rev.Mater.1, 015401 (2017) 10.1103/PhysRevMaterials.1.015401 (AEL elast) *" << endl;
      oss << "*  Mehl et al.      Comp. Mat. Sci. 136, S1 (2017)   10.1016/j.commatsci.2017.01.017 (ANRL proto1) *" << endl;
      oss << "*  Calderon et al.  Comp. Mat. Sci. 108A, 233 (2015) 10.1016/j.commatsci.2015.07.019 (standard)    *" << endl;
      oss << "*  Toher et al.     Phys. Rev. B 90, 174107 (2014)   10.1103/PhysRevB.90.174107      (AGL Gibbs)   *" << endl;
      oss << "*  Taylor et al.    Comp. Mat. Sci. 93, 178 (2014)   10.1016/j.commatsci.2014.05.014 (REST-API)    *" << endl;
      oss << "*  Curtarolo et al. Comp. Mat. Sci. 58, 227 (2012)   10.1016/j.commatsci.2012.02.002 (AFLOW.org)   *" << endl;
      oss << "*  Curtarolo et al. Comp. Mat. Sci. 58, 218 (2012)   10.1016/j.commatsci.2012.02.005 (AFLOW C++)   *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "****************************************************************************************************" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*" << aurostd::PaddedCENTER(string("aflow/aflow.org - CONTRIBUTORS"), 100) << "*" << endl;
      oss << "*  2000-2019 Stefano Curtarolo (aflow); 2002-2004 Dane Morgan (convasp); 2007-2011 Wahyu Setyawan  *" << endl;
      oss << "*  (--rsm --edos --kband --icsd*); 2008-2011 Roman Chepulskyy (--edos --kband  surfaces);          *" << endl;
      oss << "*  2008 Gus Hart (lattice reductions - prototypes); 2009-2011, Ohad Levy (prototypes);             *" << endl;
      oss << "*  2009-2010, Michal Jahnatek (APL); 2010-2013 Shidong Wang (cluster expansion); 2010-2013         *" << endl;
      oss << "*  Richard Taylor (surfaces, apennsy); 2010-2013 Junkai Xue (prototyper); 2010-2013 Kesong Yang    *" << endl;
      oss << "*  (findsym, frozsl, plotband/dos); 2013-2019 Cormac Toher (AGL Debye-Gruneisen, AEL elastic);     *" << endl;
      oss << "*  2013-2019 Frisco Rose (API, Aflux); 2013-2018 Pinku Nath (Quasi-harmonic approximation);        *" << endl;
      oss << "*  2013-2017 Jose J. Plata (AAPL, thermal cond.); 2014-2019 David Hicks (symmetry, structure       *" << endl;
      oss << "*  comparison, prototypes); 2014-2019 Corey Oses (Egap, bader, chull, APL, pocc); 2018-2019 Marco  *" << endl;
      oss << "*  Esters (AAPL, thermal cond.); 2016-2019 Denise Ford (GFA); 2018-2019 Rico Friedrich (CCE);      *" << endl;
      oss << "*  2021-2023 Simon Divilov (QCA, fitting); 2021-2023 Hagen Eckert (GFA, entry loader, JSON);       *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "****************************************************************************************************" << endl;
      oss << "*"
          << aurostd::PaddedCENTER(string("version " + string(AFLOW_VERSION) + " - g++/gcc " + aurostd::utype2string(__GNUC__) + "." + aurostd::utype2string(__GNUC_MINOR__) + "." +
                                          aurostd::utype2string(__GNUC_PATCHLEVEL__) + " - built [" + string(TODAY) + "] - (C) " + XHOST.Copyright_Years),
                                   100)
          << "*" << endl;
      oss << "****************************************************************************************************" << endl;
      return oss.str();
    }
    if (type == "BANNER_TINY") {
      // called 63 times
      oss << "AFLOW VERSION " << AFLOW_VERSION << ":  [aflow.org consortium - 2003-2023] ";
      return oss.str();
    }
    if (type == "EXCEPTIONS") {
      oss << "List of AFLOW exceptions with error codes. See README_AFLOW_EXCEPTIONS.TXT for more information" << endl;
      oss << endl;
      oss << "----------------------------------------------------------------------------------------------------" << endl;
      oss << "Error Code      Error Type             Error                              Name of Constant          " << endl;
      oss << "----------------------------------------------------------------------------------------------------" << endl;
      oss << "        1       N/A                    Generic error                      _GENERIC_ERROR_           " << endl;
      oss << "        2                              Illegal error code                 _ILLEGAL_CODE_            " << endl;
      oss << "       10       Input Error            generic                            _INPUT_ERROR_             " << endl;
      oss << "       11                              unknown flag                       _INPUT_UNKNOWN_           " << endl;
      oss << "       12                              missing flag                       _INPUT_MISSING_           " << endl;
      oss << "       13                              input ambiguous                    _INPUT_AMBIGUOUS_         " << endl;
      oss << "       14                              illegal parameter                  _INPUT_ILLEGAL_           " << endl;
      oss << "       15                              number of parameters               _INPUT_NUMBER_            " << endl;
      oss << "       20       File Error             generic                            _FILE_ERROR_              " << endl;
      oss << "       21                              file not found                     _FILE_NOT_FOUND_          " << endl;
      oss << "       22                              wrong format                       _FILE_WRONG_FORMAT_       " << endl;
      oss << "       23                              file corrupt                       _FILE_CORRUPT_            " << endl;
      oss << "       30       Value Error            generic                            _VALUE_ERROR_             " << endl;
      oss << "       31                              illegal value                      _VALUE_ILLEGAL_           " << endl;
      oss << "       32                              out of range                       _VALUE_RANGE_             " << endl;
      oss << "       40       Index Error            generic                            _INDEX_ERROR_             " << endl;
      oss << "       41                              illegal value                      _INDEX_ILLEGAL_           " << endl;
      oss << "       42                              out of bounds                      _INDEX_BOUNDS_            " << endl;
      oss << "       43                              mismatch                           _INDEX_MISMATCH_          " << endl;
      oss << "       50       Runtime Error          generic                            _RUNTIME_ERROR_           " << endl;
      oss << "       51                              not initialized                    _RUNTIME_INIT_            " << endl;
      oss << "       52                              SQL error                          _RUNTIME_SQL_             " << endl;
      oss << "       53                              busy                               _RUNTIME_BUSY_            " << endl;
      oss << "       54                              external command not found         _RUNTIME_EXTERNAL_MISS_   " << endl;  // CO20200531
      oss << "       55                              external command failed            _RUNTIME_EXTERNAL_FAIL_   " << endl;  // CO20200531
      oss << "       60       Allocation Error       generic                            _ALLOC_ERROR_             " << endl;
      oss << "       61                              could not allocate memory          _ALLOC_ALLOCATE_          " << endl;
      oss << "       62                              insufficient memory                _ALLOC_INSUFFICIENT_      " << endl;
      oss << "----------------------------------------------------------------------------------------------------" << endl;
      return oss.str();
    }
    cerr << XPID << "aflow::Banner type=" << type << " not found..." << endl;
    oss << "aflow::Banner type=" << type << " not found..." << endl;
    return oss.str();
  }
} // namespace aflow

// ***************************************************************************
// *                                                                         *
// *         aflow - Automatic FLOW for materials discovery project          *
// *                                                                         *
// ***************************************************************************

// Update Mon Mar  7 14:05:40 EST 2011 (by WSETYAWAN):
// ICSD prototypes of compounds used in aflow are generated using the following
// command:
// xzcat /common/NIST/$1ary.icsd.xz | pflow --icsd_nobrokenbasis |
// pflow --icsd_nopartialocc | pflow --icsd2proto > README_LIBRARY_ICSD$1.TXT

// #include "AFLOW3_AURO/aflow_auro.cpp"

// ***************************************************************************

// ***************************************************************************
// *                                                                         *
// *         aflow - Automatic FLOW for materials discovery project          *
// *                                                                         *
// ***************************************************************************
