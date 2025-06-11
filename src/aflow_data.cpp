// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

#include <iostream>                 // OUTSIDE XLIBS.H
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <cmath>
//#include <ctime>
#include <sys/time.h>
#include <cstdlib>  //CO20180515 - system()
//#include <time.h>

using std::string;
using std::ifstream;
using std::ostream;
using std::ofstream;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::ios_base;
using std::istream;
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::vector;

typedef unsigned uint;
#define TRUE 1
#define FALSE 0

// LINKED #include "aflow_data_htqc.cpp"  // created automatically
extern std::string Library_HTQC;

// LINKED #include "aflow_data_calculated.cpp"  // created automatically
extern std::string vLIBS;
extern std::string vAURL,vAUID,vLOOP;

// [OBSOLETE] #include "aflow_data_binary.cpp"  // created automatically
// LINKED #include "aflow_data_libraries.cpp"  // created automatically
extern std::string aflowlib_icsd;
extern std::string aflowlib_lib0;
extern std::string aflowlib_lib1;
extern std::string aflowlib_lib2;
extern std::string aflowlib_lib3;
extern std::string aflowlib_lib4;
extern std::string aflowlib_lib5;
extern std::string aflowlib_lib6;
extern std::string aflowlib_lib7;
extern std::string aflowlib_lib8;
extern std::string aflowlib_lib9;

// LINKED #include "aflow_data_stokes.cpp"  // created automatically
extern std::string FINDSYM_data_space_txt;
extern std::string FINDSYM_data_wyckoff_txt;
extern std::string FROZSL_data_space_txt;
extern std::string FROZSL_data_wyckoff_txt;
extern std::string FROZSL_data_images_txt;
extern std::string FROZSL_data_irreps_txt;
extern std::string FROZSL_data_isotropy_txt;
extern std::string FROZSL_data_little_txt;
extern std::string FROZSL_symmetry2_dat;
extern std::string FROZSL_const_dat;
extern std::string FROZSL_phvaspsetup_AFLOW;
extern std::string FROZSL_phvaspsetup_POSCAR;

// LINKED #include "aflow_data_nist.cpp"  // created automatically
extern std::string ElectronStoppingPower_txt;
extern std::string PhotonCrossSection_txt;
extern std::string PhotonStoppingPower_txt;
extern std::string ICSD_List_txt;

// LINKED #include "aflow_data_readme.cpp"  // created automatically
extern std::string README_AFLOW_LICENSE_GPL3_TXT;
extern std::string README_AFLOW_TXT;
extern std::string README_AFLOW_VERSIONS_HISTORY_TXT;
extern std::string README_AFLOW_PFLOW_TXT;
extern std::string README_AFLOW_APENNSY_TXT;
extern std::string README_AFLOW_SCRIPTING_TXT;
extern std::string README_AFLOW_FROZSL_TXT;
extern std::string README_AFLOW_POCC_TXT;
extern std::string README_AFLOW_APL_TXT;
extern std::string README_AFLOW_AGL_TXT;
extern std::string README_AFLOW_AEL_TXT;
extern std::string README_AFLOW_ANRL_TXT;
extern std::string README_AFLOW_COMPARE_TXT;
extern std::string README_AFLOW_GFA_TXT;  //CO20190401
extern std::string README_AFLOW_SYM_TXT;
extern std::string README_AFLOW_CCE_TXT;  //CO20190620
extern std::string README_AFLOW_CHULL_TXT;  //CO20190620
extern std::string README_AFLOW_EXCEPTIONS_TXT; //ME20180705
extern std::string README_PROTO_TXT;
extern std::string README_AFLOW_XAFLOW_TXT;
extern std::string README_AFLOW_AFLOWRC_TXT;

// LINKED #include "aflow_data_latex.cpp"  // created automatically
extern std::string f144468a7ccc2d3a72ba44000715efdb;
// [OBSOLETE] extern std::string d0f1b0e47f178ae627a388d3bf65d2d2;
// [OBSOLETE] extern std::string decf00ca3ad2fe494eea8e543e929068;

// LINKED #include "aflow_data_aflow_potcars.cpp"  // created automatically
extern bool AFLOW_PSEUDOPOTENTIALS;
extern std::string AFLOW_PSEUDOPOTENTIALS_TXT;
extern std::string AFLOW_PSEUDOPOTENTIALS_LIST_TXT;

// #include "aflow_data_aurostd.cpp"  // created automatically
#include "aflow_data_extra.cpp"

#ifndef defined_Library_ICSD
std::string Library_ICSD="";
#endif

std::string conststring="";

namespace aurostd {
  string StringSubst(string &strstring, const string &strfind, const string &strreplace);
  uint string2tokens(const string& str,vector<string>& tokens,const string& delimiters = " ");
  uint string2vectorstring(const string& stringIN,vector<string> &vstringout);
  bool substring2bool(const string& strstream, const string& strsub1, bool CLEAN);
  bool substring2bool(const string& strstream, const string& strsub1);
  string PaddedPOST(string input,int depth);
  string PaddedCENTER(string input,int depth);
  template<typename utype> string utype2string(const utype& from);
  template<typename utype> string utype2string(const utype& from,int precision);
  string CleanFileName(string fileIN);
  bool FileExist(const string& _FileName);
  bool FileEmpty(const string& _FileName);
  bool string2file(const string& StringOUTPUT,const string& FileNameOUTPUT);
  uint file2string(string _FileNameIN,string& StringIN);
  int FileSize(const string& _FileName);
  //  string RemoveWhiteSpaces(const string& s);
  // uint string2vectorstring(const string& stringIN,vector<string> &vstringout,bool consecutive,bool trim_edges);
  inline bool isBase64(unsigned char c);
  std::string base64Decoder(std::string const& encoded_string);
};

namespace aurostd {
  bool string2file(const string& StringOUTPUT,const string& FileNameOUTPUT) {
    ofstream FileOUTPUT;
    FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
    FileOUTPUT << StringOUTPUT;    
    FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
    return TRUE;  // return FALSE if something got messed up
  }

  bool execute(string command) {
    system(command.c_str());
    return TRUE;
  }

  long double get_useconds(void) {timeval tim;gettimeofday(&tim,NULL);return tim.tv_usec;}

  string TmpFileCreate(string identifier) {
    return "/tmp/_aflow_"+identifier+".a"+AFLOW_VERSION+"."+aurostd::utype2string(uint((double) get_useconds()))+".tmp";
  }
};


using aurostd::utype2string;
using aurostd::string2vectorstring;


namespace aflow {
  string Banner(string type) {
    stringstream oss;
    if(type=="BANNER_NORMAL") {
      oss << "****************************************************************************************************" << endl;
      oss << "MMMMM  AFLOW VERSION " << string(AFLOW_VERSION) << " Automatic-Flow [" << TODAY << "] - " << endl;// << aflow_get_time_string() << endl;
      oss << "****************************************************************************************************" << endl;
      return oss.str();
    }
    if(type=="BANNER_BIG") {
      oss << "****************************************************************************************************" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*                          aflow - Automatic-FLOW for materials discovery                          *" << endl;
      oss << "*                aflow.org consortium - High-Throughput ab-initio Computing Project                *" << endl;
      //      oss << "*                     VERSION "<<aurostd::PaddedPOST(string(AFLOW_VERSION),5) <<" - BUILT ["<<TODAY<<"] - Copyright " << XHOST.Copyright_Years << "                     *" << endl;
      oss << "*" << aurostd::PaddedCENTER(string("version "+string(AFLOW_VERSION)+" - g++/gcc "+aurostd::utype2string(__GNUC__)+"."+aurostd::utype2string(__GNUC_MINOR__)+"."+aurostd::utype2string(__GNUC_PATCHLEVEL__)+" - built ["+string(TODAY)+"] - (C) " +string("2003-2023")),100) << "*" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "****************************************************************************************************" << endl;
      return oss.str();
    }
    if(type=="BANNER_TINY") {
      oss << "AFLOW VERSION "<<string(AFLOW_VERSION)<<":  [Stefano Curtarolo - 2003-2023] ";
      return oss.str();
    }
    cerr << "aflow::Banner type=" << type << " not found..." << endl;
    oss << "aflow::Banner type=" << type << " not found..." << endl;
    return oss.str();
  }
} // namespace aflow

#define AFLOW_LIBRARY_DIRECTORIES          string("/common/AFLOW/LIBS/,/common/VASP,/home/aflow/common/AFLOW/LIBS/,/fslhome/glh43/src/,/usr/local/bin/,/fslhome/fslcollab8/group/bin/,/home/auro/work/AFLOW3/,~/common/AFLOW/LIBS/,./,/nics/a/proj/aflow/common/AFLOW/LIBS/,/home/users/aflow/common/AFLOW/LIBS,/share/apps/AFLOW3/VASP,/home/junkai/PROTO_DATABASE/,/projects/kyang-group/common/LIBS,/somewhere/")  // first is default, tokenized with ","

int main(int _argc,char **_argv) {
  vector<string> vtemp;
  bool LDEBUG=FALSE;

  vector<string> vAFLOW_LIBRARY_DIRECTORIES;
  aurostd::string2tokens(string(AFLOW_LIBRARY_DIRECTORIES),vAFLOW_LIBRARY_DIRECTORIES,",");// vAFLOW_LIBRARY_DIRECTORIES;

  if(_argc==1) {
    cout << aflow::Banner("BANNER_BIG");
  }
  for(int i=1;i<_argc;i++) {
    string argvi(_argv[i]);

    //     if(_argc==1) {
    //       cout << aflow::Banner("BANNER_BIG");
    //       cout << " V=" << string(AFLOW_VERSION) << " -> " << argvi << endl;
    //       cout << "********************************************************************************" << endl;
    // }
    //  cerr << "i=" << i << " _argv[i]=[" << argvi << "]" << endl;
    if(argvi=="--avail") {
      cout << aflow::Banner("BANNER_BIG");
      cout << " V=" << string(AFLOW_VERSION) << " -> " << argvi << endl;
      cout << "********************************************************************************" << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_LICENSE_GPL3_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_LICENSE_GPL3_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_LICENSE_GPL3_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_VERSIONS_HISTORY_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_VERSIONS_HISTORY_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_VERSIONS_HISTORY_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_PFLOW_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_PFLOW_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_PFLOW_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_ACONVASP_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_PFLOW_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_PFLOW_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_APENNSY_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_APENNSY_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_APENNSY_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_SCRIPTING_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_SCRIPTING_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_SCRIPTING_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_FROZSL_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_FROZSL_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_FROZSL_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_POCC_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_POCC_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_POCC_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_APL_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_APL_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_APL_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_AGL_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_AGL_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_AGL_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_AEL_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_AEL_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_AEL_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_ANRL_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_ANRL_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_ANRL_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_COMPARE_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_COMPARE_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_COMPARE_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_GFA_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_GFA_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_GFA_TXT,vtemp) << endl; //CO20190401
      cout << aurostd::PaddedPOST("README_AFLOW_SYM_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_SYM_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_SYM_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_CCE_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_CCE_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_CCE_TXT,vtemp) << endl; //CO20190620
      cout << aurostd::PaddedPOST("README_AFLOW_CHULL_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_CHULL_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_CHULL_TXT,vtemp) << endl; //CO20190620
      cout << aurostd::PaddedPOST("README_AFLOW_EXCEPTIONS_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_EXCEPTIONS_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_EXCEPTIONS_TXT,vtemp) << endl; //ME20180705
      cout << aurostd::PaddedPOST("README_PROTO_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_PROTO_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_PROTO_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_XAFLOW_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_XAFLOW_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_XAFLOW_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("README_AFLOW_AFLOWRC_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(README_AFLOW_AFLOWRC_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(README_AFLOW_AFLOWRC_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("Library_HTQC",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(Library_HTQC.size()),10) << "lines=" << aurostd::string2vectorstring(Library_HTQC,vtemp) << endl;
      cout << aurostd::PaddedPOST("Library_ICSD",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(Library_ICSD.size()),10) << "lines=" << aurostd::string2vectorstring(Library_ICSD,vtemp) << endl;
      cout << aurostd::PaddedPOST("aflowlib_lib0",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(aflowlib_lib0.size()),10) << "lines=" << aurostd::string2vectorstring(aflowlib_lib0,vtemp) << endl;
      cout << aurostd::PaddedPOST("aflowlib_lib1",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(aflowlib_lib1.size()),10) << "lines=" << aurostd::string2vectorstring(aflowlib_lib1,vtemp) << endl;
      cout << aurostd::PaddedPOST("aflowlib_lib2",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(aflowlib_lib2.size()),10) << "lines=" << aurostd::string2vectorstring(aflowlib_lib2,vtemp) << endl;
      cout << aurostd::PaddedPOST("aflowlib_lib3",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(aflowlib_lib3.size()),10) << "lines=" << aurostd::string2vectorstring(aflowlib_lib3,vtemp) << endl;
      cout << aurostd::PaddedPOST("aflowlib_lib4",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(aflowlib_lib4.size()),10) << "lines=" << aurostd::string2vectorstring(aflowlib_lib4,vtemp) << endl;
      cout << aurostd::PaddedPOST("aflowlib_lib5",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(aflowlib_lib5.size()),10) << "lines=" << aurostd::string2vectorstring(aflowlib_lib5,vtemp) << endl;
      cout << aurostd::PaddedPOST("aflowlib_lib6",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(aflowlib_lib6.size()),10) << "lines=" << aurostd::string2vectorstring(aflowlib_lib6,vtemp) << endl;
      cout << aurostd::PaddedPOST("aflowlib_lib7",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(aflowlib_lib7.size()),10) << "lines=" << aurostd::string2vectorstring(aflowlib_lib7,vtemp) << endl;
      cout << aurostd::PaddedPOST("aflowlib_lib8",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(aflowlib_lib8.size()),10) << "lines=" << aurostd::string2vectorstring(aflowlib_lib8,vtemp) << endl;
      cout << aurostd::PaddedPOST("aflowlib_lib9",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(aflowlib_lib9.size()),10) << "lines=" << aurostd::string2vectorstring(aflowlib_lib9,vtemp) << endl;
      cout << aurostd::PaddedPOST("aflowlib_icsd",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(aflowlib_icsd.size()),10) << "lines=" << aurostd::string2vectorstring(aflowlib_icsd,vtemp) << endl;
      uint vLIBS_size=vLIBS.size(),vLIBS_lines=aurostd::string2vectorstring(vLIBS,vtemp);
      cout << aurostd::PaddedPOST("vLIBS",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(vLIBS_size),10) << "lines=" << vLIBS_lines << endl;
      cout << aurostd::PaddedPOST("vAUID",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(vLIBS_size/3),10) << "lines=" << vLIBS_lines/3 << endl;
      cout << aurostd::PaddedPOST("vAURL",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(vLIBS_size/3),10) << "lines=" << vLIBS_lines/3 << endl;
      cout << aurostd::PaddedPOST("vLOOP",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(vLIBS_size/3),10) << "lines=" << vLIBS_lines/3 << endl;
      cout << aurostd::PaddedPOST("FINDSYM_data_space_txt",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(FINDSYM_data_space_txt.size()),10) << "lines=" << aurostd::string2vectorstring(FINDSYM_data_space_txt,vtemp) << endl;
      cout << aurostd::PaddedPOST("FINDSYM_data_wyckoff_txt",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(FINDSYM_data_wyckoff_txt.size()),10) << "lines=" << aurostd::string2vectorstring(FINDSYM_data_wyckoff_txt,vtemp) << endl;
      cout << aurostd::PaddedPOST("FROZSL_data_space_txt",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(FROZSL_data_space_txt.size()),10) << "lines=" << aurostd::string2vectorstring(FROZSL_data_space_txt,vtemp) << endl;
      cout << aurostd::PaddedPOST("FROZSL_data_wyckoff_txt",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(FROZSL_data_wyckoff_txt.size()),10) << "lines=" << aurostd::string2vectorstring(FROZSL_data_wyckoff_txt,vtemp) << endl;
      cout << aurostd::PaddedPOST("FROZSL_data_images_txt",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(FROZSL_data_images_txt.size()),10) << "lines=" << aurostd::string2vectorstring(FROZSL_data_images_txt,vtemp) << endl;
      cout << aurostd::PaddedPOST("FROZSL_data_irreps_txt",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(FROZSL_data_irreps_txt.size()),10) << "lines=" << aurostd::string2vectorstring(FROZSL_data_irreps_txt,vtemp) << endl;
      cout << aurostd::PaddedPOST("FROZSL_data_isotropy_txt",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(FROZSL_data_isotropy_txt.size()),10) << "lines=" << aurostd::string2vectorstring(FROZSL_data_isotropy_txt,vtemp) << endl;
      cout << aurostd::PaddedPOST("FROZSL_data_little_txt",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(FROZSL_data_little_txt.size()),10) << "lines=" << aurostd::string2vectorstring(FROZSL_data_little_txt,vtemp) << endl;
      cout << aurostd::PaddedPOST("FROZSL_symmetry2_dat",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(FROZSL_symmetry2_dat.size()),10) << "lines=" << aurostd::string2vectorstring(FROZSL_symmetry2_dat,vtemp) << endl;
      cout << aurostd::PaddedPOST("FROZSL_const_dat",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(FROZSL_const_dat.size()),10) << "lines=" << aurostd::string2vectorstring(FROZSL_const_dat,vtemp) << endl;
      cout << aurostd::PaddedPOST("FROZSL_phvaspsetup_AFLOW",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(FROZSL_phvaspsetup_AFLOW.size()),10) << "lines=" << aurostd::string2vectorstring(FROZSL_phvaspsetup_AFLOW,vtemp) << endl;
      cout << aurostd::PaddedPOST("FROZSL_phvaspsetup_POSCAR",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(FROZSL_phvaspsetup_POSCAR.size()),10) << "lines=" << aurostd::string2vectorstring(FROZSL_phvaspsetup_POSCAR,vtemp) << endl;
      cout << aurostd::PaddedPOST("ElectronStoppingPower_txt",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(ElectronStoppingPower_txt.size()),10) << "lines=" << aurostd::string2vectorstring(ElectronStoppingPower_txt,vtemp) << endl;
      cout << aurostd::PaddedPOST("PhotonCrossSection_txt",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(PhotonCrossSection_txt.size()),10) << "lines=" << aurostd::string2vectorstring(PhotonCrossSection_txt,vtemp) << endl;
      cout << aurostd::PaddedPOST("PhotonStoppingPower_txt",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(PhotonStoppingPower_txt.size()),10) << "lines=" << aurostd::string2vectorstring(PhotonStoppingPower_txt,vtemp) << endl;
      cout << aurostd::PaddedPOST("ICSD_List_txt",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(ICSD_List_txt.size()),10) << "lines=" << aurostd::string2vectorstring(ICSD_List_txt,vtemp) << endl;
      cout << aurostd::PaddedPOST("AFLOW_PSEUDOPOTENTIALS",40) << " size=1" << endl;
      cout << aurostd::PaddedPOST("AFLOW_PSEUDOPOTENTIALS_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(AFLOW_PSEUDOPOTENTIALS_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(AFLOW_PSEUDOPOTENTIALS_TXT,vtemp) << "       BASE64" << endl;
      cout << aurostd::PaddedPOST("AFLOW_PSEUDOPOTENTIALS_LIST_TXT",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(AFLOW_PSEUDOPOTENTIALS_LIST_TXT.size()),10) << "lines=" << aurostd::string2vectorstring(AFLOW_PSEUDOPOTENTIALS_LIST_TXT,vtemp) << endl;
      cout << aurostd::PaddedPOST("f144468a7ccc2d3a72ba44000715efdb",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(f144468a7ccc2d3a72ba44000715efdb.size()),10) << "lines=" << aurostd::string2vectorstring(f144468a7ccc2d3a72ba44000715efdb,vtemp) << endl;
// [OBSOLETE]       cout << aurostd::PaddedPOST("d0f1b0e47f178ae627a388d3bf65d2d2",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(d0f1b0e47f178ae627a388d3bf65d2d2.size()),10) << "lines=" << aurostd::string2vectorstring(d0f1b0e47f178ae627a388d3bf65d2d2,vtemp) << endl;
// [OBSOLETE]       cout << aurostd::PaddedPOST("decf00ca3ad2fe494eea8e543e929068",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(decf00ca3ad2fe494eea8e543e929068.size()),10) << "lines=" << aurostd::string2vectorstring(decf00ca3ad2fe494eea8e543e929068,vtemp) << endl;

#ifdef ICSD
      cout << aurostd::PaddedPOST("amir_natan_icsd_1_ary",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(icsd_1_ary.size()),10) << "lines=" << aurostd::string2vectorstring(icsd_1_ary,vtemp) << endl;
      cout << aurostd::PaddedPOST("amir_natan_icsd_2_ary",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(icsd_2_ary.size()),10) << "lines=" << aurostd::string2vectorstring(icsd_2_ary,vtemp) << endl;
      cout << aurostd::PaddedPOST("amir_natan_icsd_3_ary",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(icsd_3_ary.size()),10) << "lines=" << aurostd::string2vectorstring(icsd_3_ary,vtemp) << endl;
      cout << aurostd::PaddedPOST("amir_natan_icsd_4_ary",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(icsd_4_ary.size()),10) << "lines=" << aurostd::string2vectorstring(icsd_4_ary,vtemp) << endl;
      cout << aurostd::PaddedPOST("amir_natan_icsd_5_ary",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(icsd_5_ary.size()),10) << "lines=" << aurostd::string2vectorstring(icsd_5_ary,vtemp) << endl;
      cout << aurostd::PaddedPOST("amir_natan_icsd_6_ary",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(icsd_6_ary.size()),10) << "lines=" << aurostd::string2vectorstring(icsd_6_ary,vtemp) << endl;
      cout << aurostd::PaddedPOST("amir_natan_icsd_7_ary",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(icsd_7_ary.size()),10) << "lines=" << aurostd::string2vectorstring(icsd_7_ary,vtemp) << endl;
      cout << aurostd::PaddedPOST("amir_natan_icsd_8_ary",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(icsd_8_ary.size()),10) << "lines=" << aurostd::string2vectorstring(icsd_8_ary,vtemp) << endl;
      cout << aurostd::PaddedPOST("amir_natan_icsd_9_ary",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(icsd_9_ary.size()),10) << "lines=" << aurostd::string2vectorstring(icsd_9_ary,vtemp) << endl;
      cout << aurostd::PaddedPOST("amir_natan_icsd_10_ary",40) << " size=" << aurostd::PaddedPOST(aurostd::utype2string(icsd_10_ary.size()),10) << "lines=" << aurostd::string2vectorstring(icsd_10_ary,vtemp) << endl;
#endif
      cout << "********************************************************************************" << endl;
    }

    string *pstr=NULL;
    bool found=FALSE;
    bool BASE64=TRUE;
    //   (*pstr)=" ";
    if(!found && argvi=="Library_HTQC") {found=TRUE;pstr=&Library_HTQC;}// << endl;
    if(!found && argvi=="Library_ICSD") {found=TRUE;pstr=&Library_ICSD;} // << endl;
    if(!found && argvi=="aflowlib_lib0") {found=TRUE;pstr=&aflowlib_lib0;}//  << endl;
    if(!found && argvi=="aflowlib_lib1") {found=TRUE;pstr=&aflowlib_lib1;}//  << endl;
    if(!found && argvi=="aflowlib_lib2") {found=TRUE;pstr=&aflowlib_lib2;}//  << endl;
    if(!found && argvi=="aflowlib_lib3") {found=TRUE;pstr=&aflowlib_lib3;}//  << endl;
    if(!found && argvi=="aflowlib_lib4") {found=TRUE;pstr=&aflowlib_lib4;}//  << endl;
    if(!found && argvi=="aflowlib_lib5") {found=TRUE;pstr=&aflowlib_lib5;}//  << endl;
    if(!found && argvi=="aflowlib_lib6") {found=TRUE;pstr=&aflowlib_lib6;}//  << endl;
    if(!found && argvi=="aflowlib_lib7") {found=TRUE;pstr=&aflowlib_lib7;}//  << endl;
    if(!found && argvi=="aflowlib_lib8") {found=TRUE;pstr=&aflowlib_lib8;}//  << endl;
    if(!found && argvi=="aflowlib_lib9") {found=TRUE;pstr=&aflowlib_lib9;}//  << endl;
    if(!found && argvi=="aflowlib_icsd") {found=TRUE;pstr=&aflowlib_icsd;}//  << endl;
    if(!found && argvi=="FINDSYM_data_space_txt") {found=TRUE;pstr=&FINDSYM_data_space_txt;}//  << endl;
    if(!found && argvi=="FINDSYM_data_wyckoff_txt") {found=TRUE;pstr=&FINDSYM_data_wyckoff_txt;}//  << endl;
    if(!found && argvi=="FROZSL_data_space_txt") {found=TRUE;pstr=&FROZSL_data_space_txt;}//  << endl;
    if(!found && argvi=="FROZSL_data_wyckoff_txt") {found=TRUE;pstr=&FROZSL_data_wyckoff_txt;}//  << endl;
    if(!found && argvi=="FROZSL_data_images_txt") {found=TRUE;pstr=&FROZSL_data_images_txt;}//  << endl;
    if(!found && argvi=="FROZSL_data_irreps_txt") {found=TRUE;pstr=&FROZSL_data_irreps_txt;}// << endl;
    if(!found && argvi=="FROZSL_data_isotropy_txt") {found=TRUE;pstr=&FROZSL_data_isotropy_txt;}//  << endl;
    if(!found && argvi=="FROZSL_data_little_txt") {found=TRUE;pstr=&FROZSL_data_little_txt;}//  << endl;
    if(!found && argvi=="FROZSL_symmetry2_dat") {found=TRUE;pstr=&FROZSL_symmetry2_dat;}//  << endl;
    if(!found && argvi=="FROZSL_const_dat") {found=TRUE;pstr=&FROZSL_const_dat;}//  << endl;
    if(!found && argvi=="FROZSL_phvaspsetup_AFLOW") { //   \#/_ASCII_23_/g' | perl -p -e 's/\"/_ASCII_22_/g' | perl -p -e 's/\\/_ASCII_5C_
      aurostd::StringSubst(FROZSL_phvaspsetup_AFLOW,"_ASCII_22_","\"");
      aurostd::StringSubst(FROZSL_phvaspsetup_AFLOW,"_ASCII_5C_","\\");
      aurostd::StringSubst(FROZSL_phvaspsetup_AFLOW,"_ASCII_23_","#");
      found=TRUE;pstr=&FROZSL_phvaspsetup_AFLOW;}
    if(!found && argvi=="FROZSL_phvaspsetup_POSCAR") { //   \#/_ASCII_23_/g' | perl -p -e 's/\"/_ASCII_22_/g' | perl -p -e 's/\\/_ASCII_5C_
      aurostd::StringSubst(FROZSL_phvaspsetup_POSCAR,"_ASCII_22_","\"");
      aurostd::StringSubst(FROZSL_phvaspsetup_POSCAR,"_ASCII_5C_","\\");
      aurostd::StringSubst(FROZSL_phvaspsetup_POSCAR,"_ASCII_23_","#");
      found=TRUE;pstr=&FROZSL_phvaspsetup_POSCAR;}
    if(!found && argvi=="README_AFLOW_LICENSE_GPL3_TXT") {found=TRUE;pstr=&README_AFLOW_LICENSE_GPL3_TXT;} // << endl;
    if(!found && argvi=="README_AFLOW_TXT") {found=TRUE;pstr=&README_AFLOW_TXT;} // << endl;
    if(!found && argvi=="README_AFLOW_VERSIONS_HISTORY_TXT") {found=TRUE;pstr=&README_AFLOW_VERSIONS_HISTORY_TXT;} // << endl;
    if(!found && argvi=="README_AFLOW_PFLOW_TXT") {found=TRUE;pstr=&README_AFLOW_PFLOW_TXT;} // << endl;
    if(!found && argvi=="README_AFLOW_ACONVASP_TXT") {found=TRUE;pstr=&README_AFLOW_PFLOW_TXT;} // << endl;
    if(!found && argvi=="README_AFLOW_APENNSY_TXT") {found=TRUE;pstr=&README_AFLOW_APENNSY_TXT;}//  << endl;
    if(!found && argvi=="README_AFLOW_SCRIPTING_TXT") {found=TRUE;pstr=&README_AFLOW_SCRIPTING_TXT;}//  << endl;
    if(!found && argvi=="README_AFLOW_FROZSL_TXT") {found=TRUE;pstr=&README_AFLOW_FROZSL_TXT;} // << endl;
    if(!found && argvi=="README_AFLOW_POCC_TXT") {found=TRUE;pstr=&README_AFLOW_POCC_TXT;} // << endl;
    if(!found && argvi=="README_AFLOW_APL_TXT") {found=TRUE;pstr=&README_AFLOW_APL_TXT;} // << endl;
    if(!found && argvi=="README_AFLOW_AGL_TXT") {found=TRUE;pstr=&README_AFLOW_AGL_TXT;} // << endl;
    if(!found && argvi=="README_AFLOW_AEL_TXT") {found=TRUE;pstr=&README_AFLOW_AEL_TXT;} // << endl;
    if(!found && argvi=="README_AFLOW_ANRL_TXT") {found=TRUE;pstr=&README_AFLOW_ANRL_TXT;} // << endl;
    if(!found && argvi=="README_AFLOW_COMPARE_TXT") {found=TRUE;pstr=&README_AFLOW_COMPARE_TXT;} // << endl;
    if(!found && argvi=="README_AFLOW_GFA_TXT") {found=TRUE;pstr=&README_AFLOW_GFA_TXT;} // << endl;  //CO20190401
    if(!found && argvi=="README_AFLOW_SYM_TXT") {found=TRUE;pstr=&README_AFLOW_SYM_TXT;} // << endl;
    if(!found && argvi=="README_AFLOW_CCE_TXT") {found=TRUE;pstr=&README_AFLOW_CCE_TXT;} // << endl;  //CO20190620
    if(!found && argvi=="README_AFLOW_CHULL_TXT") {found=TRUE;pstr=&README_AFLOW_CHULL_TXT;} // << endl;  //CO20190620
    if(!found && argvi=="README_AFLOW_EXCEPTIONS_TXT") {found=TRUE;pstr=&README_AFLOW_EXCEPTIONS_TXT;} //ME20180705
    if(!found && argvi=="README_PROTO_TXT") {found=TRUE;pstr=&README_PROTO_TXT;} // << endl;
    if(!found && argvi=="README_AFLOW_XAFLOW_TXT") {found=TRUE;pstr=&README_AFLOW_XAFLOW_TXT;} // << endl;
    if(!found && argvi=="README_AFLOW_AFLOWRC_TXT") {found=TRUE;pstr=&README_AFLOW_AFLOWRC_TXT;} // << endl;
    if(!found && argvi=="ElectronStoppingPower_txt") {found=TRUE;pstr=&ElectronStoppingPower_txt;} // << endl;
    if(!found && argvi=="PhotonCrossSection_txt") {found=TRUE;pstr=&PhotonCrossSection_txt;} // << endl;
    if(!found && argvi=="PhotonStoppingPower_txt") {found=TRUE;pstr=&PhotonStoppingPower_txt;} // << endl;
    if(!found && argvi=="ICSD_List_txt") {found=TRUE;pstr=&ICSD_List_txt;} // << endl;
    //    if(!found && argvi=="AFLOW_PSEUDOPOTENTIALS") {found=TRUE;if(AFLOW_PSEUDOPOTENTIALS) {conststring="1";} else {conststring="0";} pstr=&conststring;}
    if(!found && argvi=="AFLOW_PSEUDOPOTENTIALS") {found=TRUE;cout << AFLOW_PSEUDOPOTENTIALS << endl;return 1;}


    if(!found && argvi=="AFLOW_PSEUDOPOTENTIALS_LIST_TXT") {found=TRUE;pstr=&AFLOW_PSEUDOPOTENTIALS_LIST_TXT;} // << endl;
    if(!found && aurostd::substring2bool(argvi,"AFLOW_PSEUDOPOTENTIALS_LIST_TXT:")) {
      if(LDEBUG) {cerr << "AFLOW_DATA: AFLOW_PSEUDOPOTENTIALS_LIST_TXT: " << endl;}
      found=TRUE;
      vector<string> tokens;
      aurostd::string2tokens(argvi,tokens,":");
      string query="";
      if(tokens.size()>0) {query=tokens.at(1);}
      if(LDEBUG) {cerr << "AFLOW_DATA: AFLOW_PSEUDOPOTENTIALS_LIST_TXT: " << tokens.at(0) << " & " << tokens.at(0) << endl;}
      aurostd::string2vectorstring(AFLOW_PSEUDOPOTENTIALS_LIST_TXT,tokens);
      if(LDEBUG) {cerr << "AFLOW_DATA: AFLOW_PSEUDOPOTENTIALS_LIST_TXT: " << tokens.size() << endl;}
      for(uint i=0;i<tokens.size();i++) {
        if(aurostd::substring2bool(tokens.at(i),query)) {conststring+=tokens.at(i);}
      }
      pstr=&conststring;
    }

    // vAUID vAURL vLOOP vLIBS
    // [OBSOLETE] if(!found && argvi=="vAUID") {found=TRUE;pstr=&vAUID;}//  << endl;
    // [OBSOLETE] if(!found && argvi=="vAURL") {found=TRUE;pstr=&vAURL;}
    // [OBSOLETE] if(!found && argvi=="vLOOP") {found=TRUE;pstr=&vLOOP;}
    if(!found && argvi=="vLIBS") {found=TRUE;pstr=&vLIBS;}
    if(!found && (argvi=="vAURL" || argvi=="vAUID" || argvi=="vLOOP")) {
      found=TRUE;
      vector<string> vtemp;
      aurostd::string2vectorstring(vLIBS,vtemp);
      for(uint i=0;i<vtemp.size();) {
        vAURL+=vtemp.at(i++)+"\n";
        vAUID+=vtemp.at(i++)+"\n";
        vLOOP+=vtemp.at(i++)+"\n";
      }
      if(argvi=="vAURL") pstr=&vAURL;
      if(argvi=="vAUID") pstr=&vAUID;
      if(argvi=="vLOOP") pstr=&vLOOP;
    } 

    // AFLOW_PSEUDOPOTENTIALS_TXT
    if(!BASE64) { //  no BASE64
      if(!found && argvi=="AFLOW_PSEUDOPOTENTIALS_TXT") {found=TRUE;pstr=&AFLOW_PSEUDOPOTENTIALS_TXT;} // << endl;
      if(!found && aurostd::substring2bool(argvi,"AFLOW_PSEUDOPOTENTIALS_TXT:")) {
        if(LDEBUG) {cerr << "AFLOW_DATA: AFLOW_PSEUDOPOTENTIALS_TXT: " << endl;}
        found=TRUE;
        vector<string> tokens;
        aurostd::string2tokens(argvi,tokens,":");
        string query="";
        if(tokens.size()>0) {query=tokens.at(1);}
        if(LDEBUG) {cerr << "AFLOW_DATA: AFLOW_PSEUDOPOTENTIALS_TXT: " << tokens.at(0) << " & " << tokens.at(0) << endl;}
        aurostd::string2vectorstring(AFLOW_PSEUDOPOTENTIALS_TXT,tokens);
        if(LDEBUG) {cerr << "AFLOW_DATA: AFLOW_PSEUDOPOTENTIALS_TXT: " << tokens.size() << endl;}
        uint istart=0,istop=0;
        for(uint i=0;i<tokens.size();i++) {
          if(aurostd::substring2bool(tokens.at(i),"START."+query)) {istart=i;}
          if(aurostd::substring2bool(tokens.at(i),"STOP."+query)) {istop=i;}
        }
        for(uint i=istart+1;i<istop;i++) {
          for(uint j=0;j<tokens.at(i).size()-1;j++) {conststring+=tokens.at(i).at(j);}
          if(i<istop-1) {conststring+=string("\n");}
        }
        pstr=&conststring;
      }
    }
    if(BASE64) { // with BASE64
      string tmpfile=aurostd::TmpFileCreate("aflow_data");
      if(!found && argvi=="AFLOW_PSEUDOPOTENTIALS_TXT") {  // [BASE64]
        found=TRUE;  
        aurostd::execute("rm -f "+tmpfile+"*");
        aurostd::string2file(aurostd::base64Decoder(AFLOW_PSEUDOPOTENTIALS_TXT),tmpfile+".xz");
        aurostd::execute("xz -d "+tmpfile+".xz");
        aurostd::file2string(""+tmpfile+"",conststring);
        aurostd::execute("rm -f "+tmpfile+"*");
        pstr=&conststring;  
      }  // [BASE64]
      if(!found && aurostd::substring2bool(argvi,"AFLOW_PSEUDOPOTENTIALS_TXT:")) {  // [BASE64]
        if(LDEBUG) {cerr << "AFLOW_DATA: AFLOW_PSEUDOPOTENTIALS_TXT: " << endl;}
        found=TRUE;
        vector<string> tokens;
        aurostd::string2tokens(argvi,tokens,":");
        string query="",strtmp="";
        if(tokens.size()>0) {query=tokens.at(1);}
        if(LDEBUG) {cerr << "AFLOW_DATA: AFLOW_PSEUDOPOTENTIALS_TXT: " << tokens.at(0) << " & " << tokens.at(0) << endl;}
        aurostd::execute("rm -f "+tmpfile+"*");
        aurostd::string2file(aurostd::base64Decoder(AFLOW_PSEUDOPOTENTIALS_TXT),tmpfile+".xz");
        aurostd::execute("xz -d "+tmpfile+".xz");
        aurostd::file2string(""+tmpfile+"",strtmp);
        aurostd::execute("rm -f "+tmpfile+"*");
        aurostd::string2vectorstring(strtmp,tokens);
        if(LDEBUG) {cerr << "AFLOW_DATA: AFLOW_PSEUDOPOTENTIALS_TXT: " << tokens.size() << endl;}
        uint istart=0,istop=0;
        for(uint i=0;i<tokens.size();i++) {
          if(aurostd::substring2bool(tokens.at(i),"START."+query)) {istart=i;}
          if(aurostd::substring2bool(tokens.at(i),"STOP."+query)) {istop=i;}
        }
        for(uint i=istart+1;i<istop;i++) {
          for(uint j=0;j<tokens.at(i).size()-1;j++) {conststring+=tokens.at(i).at(j);}
          if(i<istop-1) {conststring+=string("\n");}
        }
        pstr=&conststring;
      }  // [BASE64]
    }

    //      std::string base64Decoder(std::string const& encoded_string);

    if(!found && argvi=="f144468a7ccc2d3a72ba44000715efdb") {found=TRUE;pstr=&f144468a7ccc2d3a72ba44000715efdb;} // << endl;
// [OBSOLETE]     if(!found && argvi=="d0f1b0e47f178ae627a388d3bf65d2d2") {found=TRUE;pstr=&d0f1b0e47f178ae627a388d3bf65d2d2;} // << endl;
// [OBSOLETE]     if(!found && argvi=="decf00ca3ad2fe494eea8e543e929068") {found=TRUE;pstr=&decf00ca3ad2fe494eea8e543e929068;} // << endl;
#ifdef ICSD
    if(!found && argvi=="amir_natan_icsd_1_ary") {found=TRUE;pstr=&icsd_1_ary;} // << endl;
    if(!found && argvi=="amir_natan_icsd_2_ary") {found=TRUE;pstr=&icsd_2_ary;} // << endl;
    if(!found && argvi=="amir_natan_icsd_3_ary") {found=TRUE;pstr=&icsd_3_ary;} // << endl;
    if(!found && argvi=="amir_natan_icsd_4_ary") {found=TRUE;pstr=&icsd_4_ary;} // << endl;
    if(!found && argvi=="amir_natan_icsd_5_ary") {found=TRUE;pstr=&icsd_5_ary;} // << endl;
    if(!found && argvi=="amir_natan_icsd_6_ary") {found=TRUE;pstr=&icsd_6_ary;} // << endl;
    if(!found && argvi=="amir_natan_icsd_7_ary") {found=TRUE;pstr=&icsd_7_ary;} // << endl;
    if(!found && argvi=="amir_natan_icsd_8_ary") {found=TRUE;pstr=&icsd_8_ary;} // << endl;
    if(!found && argvi=="amir_natan_icsd_9_ary") {found=TRUE;pstr=&icsd_9_ary;} // << endl;
    if(!found && argvi=="amir_natan_icsd_10_ary") {found=TRUE;pstr=&icsd_10_ary;} // << endl;
#endif

    if(0) {
      if((found && (*pstr).length()==0) || pstr==NULL) {
        if(argvi=="Library_ICSD" || argvi=="aflowlib_lib0" || argvi=="aflowlib_lib1" || argvi=="aflowlib_lib2" || argvi=="aflowlib_lib3" || argvi=="aflowlib_lib4" || argvi=="aflowlib_lib5" || argvi=="aflowlib_lib6" || argvi=="aflowlib_lib7" || argvi=="aflowlib_lib8" || argvi=="aflowlib_lib9" || argvi=="aflowlib_icsd" ) { string FileLibrary;
          string *pstr;
          if(argvi=="Library_ICSD")  {pstr=&Library_ICSD;}
          if(argvi=="aflowlib_icsd") {pstr=&aflowlib_icsd;}
          if(argvi=="aflowlib_lib0") {pstr=&aflowlib_lib0;}
          if(argvi=="aflowlib_lib1") {pstr=&aflowlib_lib1;}
          if(argvi=="aflowlib_lib2") {pstr=&aflowlib_lib2;}
          if(argvi=="aflowlib_lib3") {pstr=&aflowlib_lib3;}
          if(argvi=="aflowlib_lib4") {pstr=&aflowlib_lib4;}
          if(argvi=="aflowlib_lib5") {pstr=&aflowlib_lib5;}
          if(argvi=="aflowlib_lib6") {pstr=&aflowlib_lib6;}
          if(argvi=="aflowlib_lib7") {pstr=&aflowlib_lib7;}
          if(argvi=="aflowlib_lib8") {pstr=&aflowlib_lib8;}
          if(argvi=="aflowlib_lib9") {pstr=&aflowlib_lib9;}

          // check if available
          if((*pstr).empty()) { // find and LOAD
            string str2search=argvi;
            aurostd::StringSubst(str2search,"Library_ICSD","aflow_library_icsd");
            (*pstr)="";
            for(uint j=0;j<vAFLOW_LIBRARY_DIRECTORIES.size() && (*pstr).empty();j++) { // cycle through possible directories
              FileLibrary=aurostd::CleanFileName(vAFLOW_LIBRARY_DIRECTORIES.at(j)+"/"+str2search+".dat");
              if(LDEBUG) {cerr << "00000  AFLOW LIBRARY  (" << j << ")  FileLibrary=" << FileLibrary << endl;}
              if(aurostd::FileExist(FileLibrary) && !aurostd::FileEmpty(FileLibrary)) {
                if(LDEBUG) {cerr << "00000  AFLOW LIBRARY  (" << j << ")  found=" <<  FileLibrary << endl;}
                if(LDEBUG) {cerr << "loading... ";}
                if(LDEBUG) {cerr.flush();}
                //       if(grep=="") {
                aurostd::file2string(FileLibrary,(*pstr));
                //   } else {(*pstr)=aurostd::execute2string("cat "+FileLibrary+" | grep -E '"+grep+"'");}
                if(LDEBUG) {cerr << "length=" << (*pstr).size();} // << " " << endl;
                if(LDEBUG) {cerr.flush();}
              }
            } // cycle through possible directories
            if((*pstr).empty()) { if(LDEBUG) {cerr << "WARNING - init::InitGlobalObject: " << argvi << " not found! " << endl;}
              found=FALSE;
            }
            cout << (*pstr) << endl;
            return 1;
          }
        }
      }
    }

    if(found) {
      if((*pstr).length()>10 && (string(_argv[1])=="--size" || string(_argv[1])=="size")) {cout << (*pstr).size() << endl;return 1;}
      if((*pstr).length()>10 && (string(_argv[1])=="--lines" || string(_argv[1])=="lines")) {cout << aurostd::string2vectorstring((*pstr),vtemp) << endl;return 1;}
      if((*pstr).length()>10) {cout << (*pstr) << endl;return 1;}
    }
  }

  return 1;
}


#define COMMENT_NEGLECT_2 string("//")
#define COMMENT_NEGLECT_3 string("!")


namespace aurostd {

  // ***************************************************************************
  // Function StringSubst
  // ***************************************************************************
  // Substitute strings here and there
  // Stefano Curtarolo
  string StringSubst(string &strstring, const string &strfind, const string &strreplace) {
    if(strfind.empty()) {return strstring;}
    string::size_type pos=0;
    while((pos=strstring.find(strfind, pos))!=string::npos) {
      strstring.erase(pos, strfind.length());
      strstring.insert(pos, strreplace);
      pos+=strreplace.length();
    }
    return strstring;
  }
  // ***************************************************************************
  // Function string2tokens and string2tokensAdd
  // ***************************************************************************
  // Finds string2tokens to split strings in tokens
  // Stefano Curtarolo
  // void string2tokens(const string& str,vector<string>& tokens,const string& delimiters=" ") {
  uint string2tokens(const string& str,std::vector<string>& tokens,const string& delimiters) {
    //  cerr << "string2tokens str='" <<  str << "'" << endl;
    //    cerr << "string2tokens delimiters='" <<  delimiters << "'" << endl;
    tokens.clear();// clear in the case there was something already in!
    string::size_type lastPos=str.find_first_not_of(delimiters,0);  // Skip delimiters at beginning.
    string::size_type pos=str.find_first_of(delimiters,lastPos);    // Find first "non-delimiter".
    while (string::npos!=pos || string::npos!=lastPos) {
      tokens.push_back(str.substr(lastPos,pos-lastPos));            // Found a token, add it to the vector.
      lastPos=str.find_first_not_of(delimiters,pos);                // Skip delimiters.  Note the "not_of"
      pos=str.find_first_of(delimiters,lastPos);                    // Find next "non-delimiter"
    }
    return tokens.size();
  }

  // ***************************************************************************
  // Function string2vectorstring
  // ***************************************************************************
  // take itring into a vector strings - Stefano Curtarolo
  uint string2vectorstring(const string& stringIN,vector<string> &vstringout) {
    return aurostd::string2tokens(stringIN,vstringout,"\n");
  }

  // ***************************************************************************
  // Function RemoveWhiteSpaces
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Morgan / Curtarolo
  string RemoveWhiteSpaces(const string& s) {
    if(s.size()==0) {return s;} // nothing to do
    string ss;
    for (uint i=0;i<s.size();i++) {
      if(s[i]!=' ' && s[i]!='\t') { ss+=s[i]; }
    }
    return ss;
  }
  string RemoveWhiteSpaces(const string& s, const char toogle) {
    if(s.size()==0) {return s;} // nothing to do
    string ss;
    bool copy=TRUE;
    for (uint i=0;i<s.size();i++) {
      if(s[i]==toogle) { copy=!copy; }
      if(copy) {
        if(s[i]!=' ' && s[i]!='\t') { ss+=s[i]; }
      }
      if(!copy) { ss+=s[i]; }
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveWhiteSpacesFromTheBack
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Morgan / Curtarolo
  string RemoveWhiteSpacesFromTheBack(const string& s) {
    if(s.size()==0) {return s;} // nothing to do
    string ss=s;
    while(ss[ss.size()-1]==' ' || ss[ss.size()-1]=='\t') {
      ss.erase(ss.size()-1,1);
      if(ss.size()==0) {
        return ss;
      } // nothing to do
    }
    return ss;
  }

  // ***************************************************************************
  // Function SubStringsPresent
  // ***************************************************************************
  bool substring2bool(const string& strstream, const string& strsub1, bool CLEAN) {
    //  if(LDEBUG) cerr << "bool aurostd::substring2bool(const string& strstream, const string& strsub1, bool CLEAN)" << endl;
    //  if(LDEBUG) cerr << "DEBUG substring_present3: (BEGIN) " << strsub1 << " " << CLEAN << endl;
    string _strstream(strstream),_strline,_strsub1(strsub1);
    string strout="";
    string::size_type idxS1;
    if(CLEAN==TRUE) {_strstream=aurostd::RemoveWhiteSpaces(_strstream,'"');}
    if(_strstream.find(_strsub1)==string::npos) {return (bool) FALSE;}
    vector<string> tokens;
    aurostd::string2tokens(_strstream,tokens,"\n");
    for(uint i=0;i<tokens.size();i++) {
      _strline=tokens.at(i);
      if(_strline.find("#")!=string::npos  
          && _strline.find("#1")==string::npos && _strline.find("#2")==string::npos && _strline.find("#3")==string::npos // avoid space groups
          && _strline.find("#4")==string::npos && _strline.find("#5")==string::npos && _strline.find("#6")==string::npos // avoid space groups
          && _strline.find("#7")==string::npos && _strline.find("#8")==string::npos && _strline.find("#9")==string::npos // avoid space groups
        ) {_strline=_strline.substr(0,_strline.find("#"));} // avoid space groups
      if(_strline.find(COMMENT_NEGLECT_2)!=string::npos) {_strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_2));}
      if(_strline.find(COMMENT_NEGLECT_3)!=string::npos) {_strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_3));}
      idxS1=_strline.find(_strsub1);
      if(idxS1!=string::npos) {
        strout=_strline.substr(_strline.find(_strsub1)+_strsub1.length());
        strout=aurostd::RemoveWhiteSpacesFromTheBack(strout);
        //   if(LDEBUG) cerr << "DEBUG substring_present3: (END) " << strsub1 << " " << CLEAN << endl;
        return (bool) TRUE;
      }
    }
    return (bool) FALSE;
  }

  bool substring2bool(const string& strstream, const string& strsub1) {
    //  if(LDEBUG)   cerr << "bool aurostd::substring2bool(const string& strstream, const string& strsub1)" << endl;
    return (bool) aurostd::substring2bool(strstream,strsub1,FALSE);
  }

  template<class utype>
    string _PaddedPOST(utype input,int depth,string ch) {
      stringstream aus("");
      aus << input;
      string strout=aus.str();
      for(int i=0;i<depth-(int) aus.str().length();i++) {strout+=ch;}
      return strout;
    }

  string PaddedPOST(string input,int depth) {
    return _PaddedPOST(input,depth," ");
  }

  // Add PRE AND POST characters to pad so that string is in the center
  template<class utype>
    string _PaddedCENTER(utype input,int depth,string ch) {
      stringstream aus("");
      int pre=(depth-(int) input.length())/2;
      int post=depth-pre-(int) input.length();
      for(int i=1;i<pre;i++) {aus << ch;} 
      aus << input;
      for(int i=1;i<post;i++) {aus << ch;}
      return aus.str();
    }
  string PaddedCENTER(string input,int depth) {
    return _PaddedCENTER(input,depth," ");
  }

  // convert whatever into a string !
  template<typename typeTo, typename typeFrom> typeTo stream2stream(typeFrom from,int precision) {
    std::stringstream temp;temp.precision(precision);temp << from;typeTo to=typeTo();temp >> to;return to;
  }

  template<typename typeTo, typename typeFrom> typeTo stream2stream(typeFrom from) {
    return (typeTo) stream2stream<typeTo>(from,20);
  }


  template<typename utype> string utype2string(const utype& from) {
    return (string) stream2stream<string>(from);
  }
  template<typename utype> string utype2string(const utype& from,int precision) {
    return (string) stream2stream<string>(from,precision);
  }

  // ***************************************************************************
  // Function CleanFileName
  // ***************************************************************************
  // Stefano Curtarolo
  // cleans file names from obvious things
  string CleanFileName(string fileIN) {
    bool LDEBUG=(FALSE);
    string fileOUT=fileIN;
    if(LDEBUG) {cerr << "aurostd::CleanFileName: " << fileOUT << endl;}
    if(aurostd::substring2bool(fileOUT,"~/")) { aurostd::StringSubst(fileOUT,"//","/"); }
    aurostd::StringSubst(fileOUT,"/./","/");
    if(LDEBUG) {cerr << "aurostd::CleanFileName: " << fileOUT << endl;}
    return fileOUT;
  }

  // ***************************************************************************
  // Function FileExist
  // ***************************************************************************
  // Stefano Curtarolo
  // return a simple bool, nothing else, from a string (which is supposed to be
  // a file name)
  bool FileExist(const string& _FileName) {
    string FileName(CleanFileName(_FileName));
    //   cerr << FileName << endl;
    bool exist=FALSE;
    ifstream FileStream;
    FileStream.open(FileName.c_str(),std::ios::in);
    FileStream.clear();
    FileStream.close();
    if(FileStream.good()) {
      exist=TRUE;
    } else {
      exist=FALSE;
    }
    return exist;
  }

  // ***************************************************************************
  // Function FileEmpty && FileNotEmpty
  // ***************************************************************************
  // Stefano Curtarolo - jan 08
  // returns in bytes the size of a file
  bool FileEmpty(const string& _FileName) {
    string FileName(CleanFileName(_FileName));
    if(FileExist(FileName)==FALSE) {return TRUE;}  // does not exist hence empty
    // it exists
    if(1) {
      int i=0;
      ifstream FileStream;
      FileStream.open(FileName.c_str(),std::ios::in);
      char c;
      while (FileStream.get(c)&&i<256) {i++;};
      // count no more that 16... it is not worth to count more
      FileStream.close();
      if(i>0) {
        return FALSE;
      } else {
        return TRUE;
      }
    }
    if(0) {
      if(aurostd::FileSize(FileName)<=1) {return TRUE;}
      else return FALSE;
    }
    return FALSE;
  }
  // ***************************************************************************
  // Function file2string bz2file2string gzfile2string efile2string
  // ***************************************************************************
  // write file to string - Stefano Curtarolo
  uint file2string(string _FileNameIN,string& StringIN) {
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    if(!FileExist(FileNameIN)) {
      // cerr << "ERROR - file2string:  file=" << FileNameIN << " not present !" << endl;
      return 0;}   
    ifstream FileIN;
    FileIN.open(FileNameIN.c_str(),std::ios::in);
    char c;while (FileIN.get(c)) {StringIN+=c;}
    FileIN.clear();FileIN.close();
    return StringIN.length(); // return 0 if something got messed up
  }
  // ***************************************************************************
  // Function FileSize
  // ***************************************************************************
  // Stefano Curtarolo - jan 08
  // returns in bytes the size of a file
  int FileSize(const string& _FileName) {
    string FileName(CleanFileName(_FileName));
    ifstream FileStream;
    int sizeout;
    FileStream.open(FileName.c_str(),std::ios::in);
    if(!FileStream.good()) {
      sizeout=0;
    } else {
      string FileString;FileString="";
      char c;
      while (FileStream.get(c)) {FileString+=c;}
      sizeout=FileString.length();
    }
    FileStream.close();
    return sizeout;
  }
  // ***************************************************************************
  // aurostd::isBase64(unsigned char c)
  // ***************************************************************************
  // determines if char is base64
  // http://www.adp-gmbh.ch/cpp/common/base64.html
  //static inline bool isBase64(unsigned char c) {
  inline bool isBase64(unsigned char c) {
    return (isalnum(c) || (c=='+') || (c=='/'));
  }

  // ***************************************************************************
  // aurostd::base64Decoder(std::string const& encoded_string)
  // ***************************************************************************
  // decodes base64 to bytes
  // http://www.adp-gmbh.ch/cpp/common/base64.html
  const std::string base64_chars=
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz"
    "0123456789+/";
  std::string base64Decoder(std::string const& encoded_string) {
    int in_len=encoded_string.size();
    int i=0,j=0,in_=0;
    unsigned char char_array_4[4], char_array_3[3];
    std::string ret;
    while (in_len-- && ( encoded_string[in_] != '=') && isBase64(encoded_string[in_])) {
      char_array_4[i++]=encoded_string[in_]; in_++;
      if (i==4) {
        for (i=0; i <4; i++) { char_array_4[i]=base64_chars.find(char_array_4[i]); } 
        char_array_3[0]=(char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
        char_array_3[1]=((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
        char_array_3[2]=((char_array_4[2] & 0x3) << 6) + char_array_4[3]; 
        for (i=0; (i<3); i++) { ret += char_array_3[i]; }
        i=0;
      }
    }
    if (i) {
      for (j=i; j <4; j++) { char_array_4[j]=0; }
      for (j=0; j <4; j++) { char_array_4[j]=base64_chars.find(char_array_4[j]);  }
      char_array_3[0]=(char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
      char_array_3[1]=((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
      char_array_3[2]=((char_array_4[2] & 0x3) << 6) + char_array_4[3];
      for (j=0; (j<i - 1); j++) { ret += char_array_3[j]; }
    }

    return ret;
  }
}


// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************
