// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo

#ifndef _AUROSTD_MAIN_CPP_
#define _AUROSTD_MAIN_CPP_
//#include "aflow.h"
#include "aurostd.h"

#define _CIN_LINE_BUFFER_LENGTH_     16384
#ifndef  CHMOD_BIN
#define  CHMOD_BIN  string("chmod")
#endif

using std::vector;   // for pennsy
using std::deque;   // for pennsy
using std::ostream;
using std::istream;
using std::ofstream;
using std::ifstream;
using std::string;
using std::cerr;

using aurostd::utype2string;
using aurostd::xvector;
using aurostd::xmatrix;
using aurostd::sign;
using aurostd::ran0;

#define COMMENT_NEGLECT_1 string("#")
//#define COMMENT_NEGLECT_2 string("// ")
#define COMMENT_NEGLECT_2 string("//")
#define COMMENT_NEGLECT_3 string("!")

#ifdef AFLOW_MULTITHREADS_ENABLE  //CO+HE20221116
// Global mutex that prevents two xThread instances from executing a system call.
// system calls are not generally thread-safe: https://stackoverflow.com/questions/12025640/how-can-i-know-whether-a-linux-syscall-is-thread-safe
static std::mutex xthread_execute;
#endif


//CO20171215 - moved to xscalar
// ***************************************************************************
// ROUNDOFF for scalars
//namespace aurostd { //DX add roundoff for scalar values
//  template<class utype>
//  utype roundoff(const utype& x, utype _tol_){
//    return ((abs(x)<(utype) _tol_) ? (utype) 0.0 : x);
//  }
//  double _aurostd_initialize_roundoff(const double& x,double y) {return roundoff(x,y);}
//  float _aurostd_initialize_roundoff(const float& x,float y) {return roundoff(x,y);}
//  int _aurostd_initialize_roundoff(const int& x,int y) {return roundoff(x,y);}
//}

// ***************************************************************************
// TIME evolution stuff
namespace aurostd {
  int get_day(void) {time_t t=time(0);struct tm *ptr_now=localtime(&t);return get_day(*ptr_now);} //CO20200624
  int get_day(const tm& tstruct) {return tstruct.tm_mday;} //CO20200624
  int get_month(void) {time_t t=time(0);struct tm *ptr_now=localtime(&t);return get_month(*ptr_now);} //CO20200624
  int get_month(const tm& tstruct) {return tstruct.tm_mon+1;}  //CO20200624
  int get_year(void) {time_t t=time(0);struct tm *ptr_now=localtime(&t);return get_year(*ptr_now);}
  int get_year(const tm& tstruct) {return tstruct.tm_year+1900;} //CO20200624
  void get_offset_utc(int& offset_hours,int& offset_mins) {time_t t=time(0);struct tm *ptr_now=localtime(&t);return get_offset_utc(*ptr_now,offset_hours,offset_mins);}  //CO20210601: https://codereview.stackexchange.com/questions/175353/getting-current-timezone
  void get_offset_utc(const tm& _tstruct_inp,int& offset_hours,int& offset_mins) {  //CO20210601
    //https://codereview.stackexchange.com/questions/175353/getting-current-timezone
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    char buffer[30];
    tm tstruct_inp=_tstruct_inp; //mktime modifies tstruct, make copy
    time_t t_inp=std::mktime(&tstruct_inp);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " ///////////////////////////////////////////////" << endl;
      cerr << __AFLOW_FUNC__ << " LOOKING AT: tstruct_inp" << endl;
      cerr << __AFLOW_FUNC__ << " tstruct_inp.tm_sec=" << tstruct_inp.tm_sec << endl;
      cerr << __AFLOW_FUNC__ << " tstruct_inp.tm_min=" << tstruct_inp.tm_min << endl;
      cerr << __AFLOW_FUNC__ << " tstruct_inp.tm_hour=" << tstruct_inp.tm_hour << endl;
      cerr << __AFLOW_FUNC__ << " tstruct_inp.tm_mday=" << tstruct_inp.tm_mday << endl;
      cerr << __AFLOW_FUNC__ << " tstruct_inp.tm_mon=" << tstruct_inp.tm_mon << endl;
      cerr << __AFLOW_FUNC__ << " tstruct_inp.tm_year=" << tstruct_inp.tm_year << endl;
      cerr << __AFLOW_FUNC__ << " tstruct_inp.tm_wday=" << tstruct_inp.tm_wday << endl;
      cerr << __AFLOW_FUNC__ << " tstruct_inp.tm_yday=" << tstruct_inp.tm_yday << endl;
      cerr << __AFLOW_FUNC__ << " tstruct_inp.tm_isdst=" << tstruct_inp.tm_isdst << endl;
      cerr << __AFLOW_FUNC__ << " mktime(tstruct_inp)=" << t_inp << endl;
      strftime(buffer,30,"%F %T %Z",&tstruct_inp);cerr << __AFLOW_FUNC__ << " tstruct_inp=" << buffer << endl;  //%Y:%m:%d %H:%M:%S
      cerr << __AFLOW_FUNC__ << " ///////////////////////////////////////////////" << endl;
    }
    //
    time_t t_local=t_inp;
    bool fix_utc_2_now=false; //this is good for debugging different time zones
    if(fix_utc_2_now){t_local=time(0);} //struct tm *tstruct_now=localtime(&t_now);
    struct tm *ptr_tstruct_gmt=std::gmtime(&t_local); //get gmt wrt to local (now vs. input)
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " ///////////////////////////////////////////////" << endl;
      cerr << __AFLOW_FUNC__ << " LOOKING AT: ptr_tstruct_gmt (BEFORE DST CHANGE)" << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_sec=" << ptr_tstruct_gmt->tm_sec << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_min=" << ptr_tstruct_gmt->tm_min << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_hour=" << ptr_tstruct_gmt->tm_hour << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_mday=" << ptr_tstruct_gmt->tm_mday << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_mon=" << ptr_tstruct_gmt->tm_mon << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_year=" << ptr_tstruct_gmt->tm_year << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_wday=" << ptr_tstruct_gmt->tm_wday << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_yday=" << ptr_tstruct_gmt->tm_yday << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_isdst=" << ptr_tstruct_gmt->tm_isdst << endl;
      tm tstruct_tmp=*ptr_tstruct_gmt;  //mktime modifies tstruct
      cerr << __AFLOW_FUNC__ << " mktime(ptr_tstruct_gmt)=" << std::mktime(&tstruct_tmp) << endl;
      strftime(buffer,30,"%F %T %Z",ptr_tstruct_gmt);cerr << __AFLOW_FUNC__ << " tstruct_gmt=" << buffer << endl;  //%Y:%m:%d %H:%M:%S
      cerr << __AFLOW_FUNC__ << " ///////////////////////////////////////////////" << endl;
    }
    //NB: before the following DST change, the %Z of ptr_struct_gmt is GMT, after it is EST (or EDT)
    ptr_tstruct_gmt->tm_isdst=-1; //VERY IMPORTANT, forces mktime to figure out dst
    time_t t_gmt=std::mktime(ptr_tstruct_gmt);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " ///////////////////////////////////////////////" << endl;
      cerr << __AFLOW_FUNC__ << " LOOKING AT: ptr_tstruct_gmt (AFTER DST CHANGE)" << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_sec=" << ptr_tstruct_gmt->tm_sec << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_min=" << ptr_tstruct_gmt->tm_min << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_hour=" << ptr_tstruct_gmt->tm_hour << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_mday=" << ptr_tstruct_gmt->tm_mday << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_mon=" << ptr_tstruct_gmt->tm_mon << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_year=" << ptr_tstruct_gmt->tm_year << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_wday=" << ptr_tstruct_gmt->tm_wday << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_yday=" << ptr_tstruct_gmt->tm_yday << endl;
      cerr << __AFLOW_FUNC__ << " ptr_tstruct_gmt->tm_isdst=" << ptr_tstruct_gmt->tm_isdst << endl;
      cerr << __AFLOW_FUNC__ << " mktime(ptr_tstruct_gmt)=" << t_gmt << endl;
      strftime(buffer,30,"%F %T %Z",ptr_tstruct_gmt);cerr << __AFLOW_FUNC__ << " tstruct_gmt=" << buffer << endl;  //%Y:%m:%d %H:%M:%S
      cerr << __AFLOW_FUNC__ << " ///////////////////////////////////////////////" << endl;
    }
    //
    long int t_diff=static_cast<long int>(t_inp-t_gmt); //flip to get right sign
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " t_diff=" << t_diff << endl;}
    double offset=(double)t_diff/3600.0;
    offset_hours=(int)std::floor(offset);
    offset_mins=(int)(round((offset-(double)offset_hours)*60.0));
  }
  long int get_date(void) {time_t t=time(0);struct tm *ptr_now=localtime(&t);return get_date(*ptr_now);}  //CO20200624
  long int get_date(const tm& tstruct) {return aurostd::get_year(tstruct)*10000+aurostd::get_month(tstruct)*100+aurostd::get_day(tstruct);}  //CO20200624
  int get_hour(void) {time_t t=time(0);struct tm *ptr_now=localtime(&t);return get_hour(*ptr_now);}  //CO20200624
  int get_hour(const tm& tstruct) {return tstruct.tm_hour;}  //CO20200624
  int get_min(void) {time_t t=time(0);struct tm *ptr_now=localtime(&t);return get_min(*ptr_now);}  //CO20200624
  int get_min(const tm& tstruct) {return tstruct.tm_min;}  //CO20200624
  int get_sec(void) {time_t t=time(0);struct tm *ptr_now=localtime(&t);return get_sec(*ptr_now);}  //CO20200624
  int get_sec(const tm& tstruct) {return tstruct.tm_sec;}  //CO20200624
  long double get_seconds(void) {timeval tim;gettimeofday(&tim,NULL);return tim.tv_sec+tim.tv_usec/1e6;}
  long double get_seconds(long double reference_seconds) {return get_seconds()-reference_seconds;}
  long double get_delta_seconds(long double& seconds_begin) {long double out=get_seconds()-seconds_begin;seconds_begin=get_seconds();return out;}
  long double get_mseconds(void) {timeval tim;gettimeofday(&tim,NULL);return tim.tv_usec/1000.0;}
  long double get_mseconds(long double reference_useconds) {return (aurostd::get_useconds()-reference_useconds)/1000.0;}
  long double get_delta_mseconds(long double& useconds_begin) {long double out=(aurostd::get_useconds()-useconds_begin)/1000.0;useconds_begin=aurostd::get_useconds()/1000.0;return out;}
  long double get_useconds(void) {timeval tim;gettimeofday(&tim,NULL);return tim.tv_usec;}
  long double get_useconds(long double reference_useconds) {return aurostd::get_useconds()-reference_useconds;}
  long double get_delta_useconds(long double& useconds_begin) {long double out=aurostd::get_useconds()-useconds_begin;useconds_begin=aurostd::get_useconds();return out;}
  string get_time(void) {time_t t=time(0);struct tm *ptr_now=localtime(&t);return get_time(*ptr_now);} //CO20200624
  string get_time(const tm& tstruct) {int h=get_hour(tstruct),m=get_min(tstruct),s=get_sec(tstruct);return (h<10?"0":"")+aurostd::utype2string(h)+":"+(m<10?"0":"")+aurostd::utype2string(m)+":"+(s<10?"0":"")+aurostd::utype2string(s);} //CO20200624
  string get_datetime(bool include_utc_offset) {time_t t=time(0);struct tm *ptr_now=localtime(&t);return get_datetime(*ptr_now,include_utc_offset);} //CO20200624
  string get_datetime(const tm& tstruct,bool include_utc_offset) {  //CO20200624
    string datetime=utype2string(get_date(tstruct))+"_"+get_time(tstruct);
    if(include_utc_offset){ //CO20210624
      int offset_hours=0,offset_mins=0;
      get_offset_utc(tstruct,offset_hours,offset_mins);
      datetime+="_GMT";
      datetime+=(std::signbit(offset_hours)?string("-"):string("+")); //sign +/-
      bool pad_hours=false; //default AFLOW behavior does not pad hours
      datetime+=aurostd::PaddedNumString(abs(offset_hours),(pad_hours?2:1)); //CO20210624 - PaddedNumString() struggles with negative numbers
      bool print_mins=false;  //default AFLOW behavior does not print mins
      if(print_mins||offset_mins!=0){datetime+=":"+aurostd::PaddedNumString(offset_mins,2);}
    }
    return datetime;
  }
  string get_datetime_formatted(const string& date_delim,bool include_time,const string& date_time_sep,const string& time_delim){time_t t=time(0);struct tm *ptr_now=localtime(&t);return get_datetime_formatted(*ptr_now,date_delim,include_time,date_time_sep,time_delim);}  //CO20171215  //CO20200624
  string get_datetime_formatted(const tm& tstruct,const string& date_delim,bool include_time,const string& date_time_sep,const string& time_delim){  //CO20171215 //CO20200624
    stringstream misc_ss;
    int y=aurostd::get_year(tstruct),b=aurostd::get_month(tstruct),d=aurostd::get_day(tstruct); //CO20200624
    misc_ss << y << date_delim << (b<10?"0":"") << b << date_delim << (d<10?"0":"") << d;
    if(include_time){
      int h=get_hour(tstruct),m=get_min(tstruct),s=get_sec(tstruct);
      misc_ss << date_time_sep << (h<10?"0":"") << h << time_delim << (m<10?"0":"") << m << time_delim << (s<10?"0":"") << s; //CO20200624
    }
    return misc_ss.str();
  }
  bool beep(uint freq,uint duration) {
    return aurostd::execute("beep -f "+aurostd::utype2string<uint>(freq)+" -l "+aurostd::utype2string<uint>(duration));
  }
}

// ***************************************************************************
// get threadID
namespace aurostd {
  unsigned long long int getTID(void){ //CO20200502 - threadID
    //for mac these numbers can be QUITE large, so better to be safe and return unsigned long long int
    //see here: http://elliotth.blogspot.com/2012/04/gettid-on-mac-os.html
    //also for macs: pid!=tid
#ifdef _MACOSX_
#if MAC_OS_X_VERSION_MAX_ALLOWED >= MAC_OS_X_VERSION_10_12
    uint64_t tid64;
    pthread_threadid_np(NULL, &tid64);
    pid_t tid = (pid_t)tid64;
    return (unsigned long long int)tid;
#else
    //////////////////////////////////////////////////////////
#ifdef __GLIBC__
#include <sys/syscall.h>  //CO20200502 - need for gettid()
    pid_t tid = syscall(__NR_gettid);
    return (unsigned long long int)tid;
#else //ONLY if _MACOSX_ AND not __GLIBC__
    return (unsigned long long int)getpid();
#endif
    //////////////////////////////////////////////////////////
#endif  //END _MACOSX_
#else //if NOT _MACOSX_
    //////////////////////////////////////////////////////////
#ifdef __GLIBC__
    return (unsigned long long int)gettid();
#else //for example CYGWIN
    return (unsigned long long int)getpid();
#endif
    //////////////////////////////////////////////////////////
#endif
  }
}

// ***************************************************************************
// FILES creation/destruction
namespace aurostd {
  string TmpStrCreate(const string& _identifier,const string& _tmpdir,bool hidden,bool directory){
    string identifier=_identifier;if(identifier.empty()){identifier="tmp";} //CO20210624
    string tmpdir=_tmpdir;if(tmpdir.empty()){tmpdir=XHOST.tmpfs;} //CO20210315
    string str=tmpdir+"/"+(hidden?".":"")+"_aflow_"+identifier+"."+XHOST.user+".pid"+XHOST.ostrPID.str()+".tid"+XHOST.ostrTID.str()+".a"+AFLOW_VERSION+".rnd"+aurostd::utype2string(uint((double) std::floor((double)100000*aurostd::ran0())))+".u"+aurostd::utype2string(uint((double) aurostd::get_useconds()))+(directory?"_":".")+"tmp"; //CO20200502 - threadID
    str=aurostd::CleanFileName(str);
    return str;
  }
  string TmpFileCreate(const string& identifier,const string& tmpdir,bool hidden) {  //CO20210315
    bool directory=false; //creating a file
    return TmpStrCreate(identifier,tmpdir,hidden,directory);
  }
  string TmpDirectoryCreate(const string& identifier,const string& tmpdir,bool hidden) { //CO20210315
    bool directory=true; //creating a directory
    string dir=TmpStrCreate(identifier,tmpdir,hidden,directory);
    DirectoryMake(dir);
    return dir;
  }
}

// ***************************************************************************
// Function extra operator << for vector
template<class utype>                            // operator <<  vector<>
std::ostream& operator<< (std::ostream& buf,const std::vector<utype>& x) {
  for(uint i=0;i<x.size();i++) {
    buf << x[i] << " ";
  }
  return buf;
}
// ***************************************************************************
// Function extra operator << for deque
template<class utype>                            // operator <<  deque<>
std::ostream& operator<< (std::ostream& buf,const std::deque<utype>& x) {
  for(uint i=0;i<x.size();i++) {
    buf << x[i] << " ";
  }
  return buf;
}
// ***************************************************************************
std::ostream& operator<< (std::ostream& b,const vector<uint>& x) {for(uint i=0;i<x.size();i++) b << x[i] << " "; return b;}
std::ostream& operator<< (std::ostream& b,const deque<uint>& x) {for(uint i=0;i<x.size();i++) b << x[i] << " "; return b;}
std::ostream& operator<< (std::ostream& b,const vector<char>& x) {for(uint i=0;i<x.size();i++) b << x[i] << " "; return b;}
std::ostream& operator<< (std::ostream& b,const deque<char>& x) {for(uint i=0;i<x.size();i++) b << x[i] << " "; return b;}
std::ostream& operator<< (std::ostream& b,const vector<int>& x) {for(uint i=0;i<x.size();i++) b << x[i] << " "; return b;}
std::ostream& operator<< (std::ostream& b,const deque<int>& x) {for(uint i=0;i<x.size();i++) b << x[i] << " "; return b;}
std::ostream& operator<< (std::ostream& b,const vector<long>& x) {for(uint i=0;i<x.size();i++) b << x[i] << " "; return b;}
std::ostream& operator<< (std::ostream& b,const deque<long>& x) {for(uint i=0;i<x.size();i++) b << x[i] << " "; return b;}
std::ostream& operator<< (std::ostream& b,const vector<double>& x) {for(uint i=0;i<x.size();i++) b << x[i] << " "; return b;}
std::ostream& operator<< (std::ostream& b,const deque<double>& x) {for(uint i=0;i<x.size();i++) b << x[i] << " "; return b;}
std::ostream& operator<< (std::ostream& b,const vector<long double>& x) {for(uint i=0;i<x.size();i++) b << x[i] << " "; return b;}
std::ostream& operator<< (std::ostream& b,const deque<long double>& x) {for(uint i=0;i<x.size();i++) b << x[i] << " "; return b;}
std::ostream& operator<< (std::ostream& b,const vector<string>& x) {for(uint i=0;i<x.size();i++) b << x[i] << " "; return b;}
std::ostream& operator<< (std::ostream& b,const deque<string>& x) {for(uint i=0;i<x.size();i++) b << x[i] << " "; return b;}

namespace aurostd {
  // ***************************************************************************
  // Function aswap
  // ***************************************************************************
  // namespace aurostd
  template<class utype> void aswap(utype &a,utype &b) {utype temp=a;a=b;b=temp;}
  void _aurostd_initialize_aswap(bool& x,bool& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(char& x,char& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(int& x,int& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(uint& x,uint& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(float& x,float& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(double& x,double& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(string& x,string& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(long int& x,long int& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(long long int& x,long long int& y) {aswap(x,y);}
  void _aurostd_initialize_aswap(long double& x,long double& y) {aswap(x,y);}
#ifdef _AUROSTD_XCOMPLEX_
  //   void _aurostd_initialize_aswap(xcomplex<float>& x,xcomplex<float>& y) {aswap(x,y);}
  //   void _aurostd_initialize_aswap(xcomplex<double>& x,xcomplex<double>& y) {aswap(x,y);}
  //   void _aurostd_initialize_aswap(xcomplex<long double>& x,xcomplex<long double>& y) {aswap(x,y);}
#endif

  // ***************************************************************************
  // Function max of a vector/deque
  // ***************************************************************************
  //SC
  template<class utype> utype max(const vector<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=vec.at(0);
    for(uint i=0;i<vec.size();i++) if(vec[i]>=out) out=vec[i];
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_max(const vector<bool> vec) { return max(vec);}
  char _aurostd_initialize_max(const vector<char> vec) { return max(vec);}
  string _aurostd_initialize_max(const vector<string> vec) { return max(vec);}
  int _aurostd_initialize_max(const vector<int> vec) { return max(vec);}
  long _aurostd_initialize_max(const vector<long> vec) { return max(vec);}
  uint _aurostd_initialize_max(const vector<uint> vec) { return max(vec);}
  float _aurostd_initialize_max(const vector<float> vec) { return max(vec);}
  double _aurostd_initialize_max(const vector<double> vec) { return max(vec);}
  long double _aurostd_initialize_max(const vector<long double> vec) { return max(vec);}

  template<class utype> utype max(const deque<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=vec.at(0);
    for(uint i=0;i<vec.size();i++) if(vec[i]>=out) out=vec[i];
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_max(const deque<bool> vec) { return max(vec);}
  char _aurostd_initialize_max(const deque<char> vec) { return max(vec);}
  string _aurostd_initialize_max(const deque<string> vec) { return max(vec);}
  int _aurostd_initialize_max(const deque<int> vec) { return max(vec);}
  long _aurostd_initialize_max(const deque<long> vec) { return max(vec);}
  uint _aurostd_initialize_max(const deque<uint> vec) { return max(vec);}
  float _aurostd_initialize_max(const deque<float> vec) { return max(vec);}
  double _aurostd_initialize_max(const deque<double> vec) { return max(vec);}
  long double _aurostd_initialize_max(const deque<long double> vec) { return max(vec);}

  // ***************************************************************************
  // Function max of a vector<vector<>>
  // ***************************************************************************
  template<class utype> utype max(const vector<vector<utype> > mat) {
    if(mat.size()==0) return (utype) 0;
    if(mat.at(0).size()==0) return (utype) 0;
    utype out=mat.at(0).at(0);
    for(uint i=0;i<mat.size();i++)
      for(uint j=0;j<mat[i].size();j++)
        if(mat[i][j]>=out) out=mat[i][j];
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_max(const vector<vector<bool> > mat) { return max(mat);}
  char _aurostd_initialize_max(const vector<vector<char> > mat) { return max(mat);}
  string _aurostd_initialize_max(const vector<vector<string> > mat) { return max(mat);}
  int _aurostd_initialize_max(const vector<vector<int> > mat) { return max(mat);}
  long _aurostd_initialize_max(const vector<vector<long> > mat) { return max(mat);}
  uint _aurostd_initialize_max(const vector<vector<uint> > mat) { return max(mat);}
  float _aurostd_initialize_max(const vector<vector<float> > mat) { return max(mat);}
  double _aurostd_initialize_max(const vector<vector<double> > mat) { return max(mat);}
  long double _aurostd_initialize_max(const vector<vector<long double> > mat) { return max(mat);}

  // ***************************************************************************
  // Function min of a vector/deque
  // ***************************************************************************
  //SC
  template<class utype> utype min(const vector<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=vec.at(0);
    for(uint i=0;i<vec.size();i++) if(vec[i]<=out) out=vec[i];
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_min(const vector<bool> vec) { return min(vec);}
  char _aurostd_initialize_min(const vector<char> vec) { return min(vec);}
  string _aurostd_initialize_min(const vector<string> vec) { return min(vec);}
  int _aurostd_initialize_min(const vector<int> vec) { return min(vec);}
  long _aurostd_initialize_min(const vector<long> vec) { return min(vec);}
  uint _aurostd_initialize_min(const vector<uint> vec) { return min(vec);}
  float _aurostd_initialize_min(const vector<float> vec) { return min(vec);}
  double _aurostd_initialize_min(const vector<double> vec) { return min(vec);}
  long double _aurostd_initialize_min(const vector<long double> vec) { return min(vec);}

  //SC
  template<class utype> utype min(const deque<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=vec.at(0);
    for(uint i=0;i<vec.size();i++) if(vec[i]<=out) out=vec[i];
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_min(const deque<bool> vec) { return min(vec);}
  char _aurostd_initialize_min(const deque<char> vec) { return min(vec);}
  string _aurostd_initialize_min(const deque<string> vec) { return min(vec);}
  int _aurostd_initialize_min(const deque<int> vec) { return min(vec);}
  long _aurostd_initialize_min(const deque<long> vec) { return min(vec);}
  uint _aurostd_initialize_min(const deque<uint> vec) { return min(vec);}
  float _aurostd_initialize_min(const deque<float> vec) { return min(vec);}
  double _aurostd_initialize_min(const deque<double> vec) { return min(vec);}
  long double _aurostd_initialize_min(const deque<long double> vec) { return min(vec);}

  // ***************************************************************************
  // Function min of a vector<vector<>>
  // ***************************************************************************
  template<class utype> utype min(const vector<vector<utype> > mat) {
    if(mat.size()==0) return (utype) 0;
    if(mat.at(0).size()==0) return (utype) 0;
    utype out=mat.at(0).at(0);
    for(uint i=0;i<mat.size();i++)
      for(uint j=0;j<mat[i].size();j++)
        if(mat[i][j]<=out) out=mat[i][j];
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_min(const vector<vector<bool> > mat) { return min(mat);}
  char _aurostd_initialize_min(const vector<vector<char> > mat) { return min(mat);}
  string _aurostd_initialize_min(const vector<vector<string> > mat) { return min(mat);}
  int _aurostd_initialize_min(const vector<vector<int> > mat) { return min(mat);}
  long _aurostd_initialize_min(const vector<vector<long> > mat) { return min(mat);}
  uint _aurostd_initialize_min(const vector<vector<uint> > mat) { return min(mat);}
  float _aurostd_initialize_min(const vector<vector<float> > mat) { return min(mat);}
  double _aurostd_initialize_min(const vector<vector<double> > mat) { return min(mat);}
  long double _aurostd_initialize_min(const vector<vector<long double> > mat) { return min(mat);}

  // ***************************************************************************
  // Function sum of a vector/deque
  // ***************************************************************************
  template<class utype> utype sum(const vector<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=0;
    for(uint i=0;i<vec.size();i++) out+=vec[i];
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_sum(const vector<bool> vec) { return sum(vec);}
  char _aurostd_initialize_sum(const vector<char> vec) { return sum(vec);}
  string _aurostd_initialize_sum(const vector<string> vec) { return sum(vec);}
  int _aurostd_initialize_sum(const vector<int> vec) { return sum(vec);}
  long _aurostd_initialize_sum(const vector<long> vec) { return sum(vec);}
  uint _aurostd_initialize_sum(const vector<uint> vec) { return sum(vec);}
  float _aurostd_initialize_sum(const vector<float> vec) { return sum(vec);}
  double _aurostd_initialize_sum(const vector<double> vec) { return sum(vec);}
  long double _aurostd_initialize_sum(const vector<long double> vec) { return sum(vec);}

  template<class utype> utype sum(const deque<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=0;
    for(uint i=0;i<vec.size();i++) out+=vec[i];
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_sum(const deque<bool> vec) { return sum(vec);}
  char _aurostd_initialize_sum(const deque<char> vec) { return sum(vec);}
  string _aurostd_initialize_sum(const deque<string> vec) { return sum(vec);}
  int _aurostd_initialize_sum(const deque<int> vec) { return sum(vec);}
  long _aurostd_initialize_sum(const deque<long> vec) { return sum(vec);}
  uint _aurostd_initialize_sum(const deque<uint> vec) { return sum(vec);}
  float _aurostd_initialize_sum(const deque<float> vec) { return sum(vec);}
  double _aurostd_initialize_sum(const deque<double> vec) { return sum(vec);}
  long double _aurostd_initialize_sum(const deque<long double> vec) { return sum(vec);}

  // ***************************************************************************
  // Function sum of a vector<vector<>>
  // ***************************************************************************
  template<class utype> utype sum(const vector<vector<utype> > mat) {
    if(mat.size()==0) return (utype) 0;
    if(mat.at(0).size()==0) return (utype) 0;
    utype out=0;
    for(uint i=0;i<mat.size();i++)
      for(uint j=0;j<mat[i].size();j++)
        out+=mat[i][j];
    return out;
  }
  // overload to force compiling
  bool _aurostd_initialize_sum(const vector<vector<bool> > mat) { return sum(mat);}
  char _aurostd_initialize_sum(const vector<vector<char> > mat) { return sum(mat);}
  string _aurostd_initialize_sum(const vector<vector<string> > mat) { return sum(mat);}
  int _aurostd_initialize_sum(const vector<vector<int> > mat) { return sum(mat);}
  long _aurostd_initialize_sum(const vector<vector<long> > mat) { return sum(mat);}
  uint _aurostd_initialize_sum(const vector<vector<uint> > mat) { return sum(mat);}
  float _aurostd_initialize_sum(const vector<vector<float> > mat) { return sum(mat);}
  double _aurostd_initialize_sum(const vector<vector<double> > mat) { return sum(mat);}
  long double _aurostd_initialize_sum(const vector<vector<long double> > mat) { return sum(mat);}

  // ***************************************************************************
  // Function mean of a vector/deque
  // ***************************************************************************
  template<class utype> utype mean(const vector<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=0;
    for(uint i=0;i<vec.size();i++) out+=vec[i];
    return (utype) out/((utype) vec.size());
  }
  // overload to force compiling
  bool _aurostd_initialize_mean(const vector<bool> vec) { return mean(vec);}
  char _aurostd_initialize_mean(const vector<char> vec) { return mean(vec);}
  int _aurostd_initialize_mean(const vector<int> vec) { return mean(vec);}
  long _aurostd_initialize_mean(const vector<long> vec) { return mean(vec);}
  uint _aurostd_initialize_mean(const vector<uint> vec) { return mean(vec);}
  float _aurostd_initialize_mean(const vector<float> vec) { return mean(vec);}
  double _aurostd_initialize_mean(const vector<double> vec) { return mean(vec);}
  long double _aurostd_initialize_mean(const vector<long double> vec) { return mean(vec);}

  template<class utype> utype mean(const deque<utype> vec) {
    if(vec.size()==0) return (utype) 0;
    utype out=0;
    for(uint i=0;i<vec.size();i++) out+=vec[i];
    return (utype) out/((utype) vec.size());
  }
  // overload to force compiling
  bool _aurostd_initialize_mean(const deque<bool> vec) { return mean(vec);}
  char _aurostd_initialize_mean(const deque<char> vec) { return mean(vec);}
  int _aurostd_initialize_mean(const deque<int> vec) { return mean(vec);}
  long _aurostd_initialize_mean(const deque<long> vec) { return mean(vec);}
  uint _aurostd_initialize_mean(const deque<uint> vec) { return mean(vec);}
  float _aurostd_initialize_mean(const deque<float> vec) { return mean(vec);}
  double _aurostd_initialize_mean(const deque<double> vec) { return mean(vec);}
  long double _aurostd_initialize_mean(const deque<long double> vec) { return mean(vec);}

  // ***************************************************************************
  // Function mean of a vector<vector<>>
  // ***************************************************************************
  template<class utype> utype mean(const vector<vector<utype> > mat) {
    if(mat.size()==0) return (utype) 0;
    if(mat.at(0).size()==0) return (utype) 0;
    utype out=0;
    for(uint i=0;i<mat.size();i++)
      for(uint j=0;j<mat[i].size();j++)
        out+=mat[i][j];
    return (utype) out/((utype) mat.size()*mat.at(0).size());
  }
  // overload to force compiling
  bool _aurostd_initialize_mean(const vector<vector<bool> > mat) { return mean(mat);}
  char _aurostd_initialize_mean(const vector<vector<char> > mat) { return mean(mat);}
  int _aurostd_initialize_mean(const vector<vector<int> > mat) { return mean(mat);}
  long _aurostd_initialize_mean(const vector<vector<long> > mat) { return mean(mat);}
  uint _aurostd_initialize_mean(const vector<vector<uint> > mat) { return mean(mat);}
  float _aurostd_initialize_mean(const vector<vector<float> > mat) { return mean(mat);}
  double _aurostd_initialize_mean(const vector<vector<double> > mat) { return mean(mat);}
  long double _aurostd_initialize_mean(const vector<vector<long double> > mat) { return mean(mat);}

  // ***************************************************************************
  // Function reset of a vector/deque
  // ***************************************************************************
  template<class utype> vector<utype> reset(vector<utype>& vec) {
    for(uint i=0;i<vec.size();i++) vec[i]=(utype) 0;
    return vec;
  }
  // overload to force compiling
  vector<bool>  _aurostd_initialize_reset(vector<bool>& vec) { return reset(vec);}
  vector<char>  _aurostd_initialize_reset(vector<char>& vec) { return reset(vec);}
  vector<string>  _aurostd_initialize_reset(vector<string>& vec) { return reset(vec);}
  vector<int>   _aurostd_initialize_reset(vector<int>& vec) { return reset(vec);}
  vector<long>  _aurostd_initialize_reset(vector<long>& vec) { return reset(vec);}
  vector<uint>  _aurostd_initialize_reset(vector<uint>& vec) { return reset(vec);}
  vector<float> _aurostd_initialize_reset(vector<float>& vec) { return reset(vec);}
  vector<double>  _aurostd_initialize_reset(vector<double>& vec) { return reset(vec);}
  vector<long double>  _aurostd_initialize_reset(vector<long double>& vec) { return reset(vec);}

  template<class utype> deque<utype> reset(deque<utype>& vec) {
    for(uint i=0;i<vec.size();i++) vec[i]=(utype) 0;
    return vec;
  }
  // overload to force compiling
  deque<bool>  _aurostd_initialize_reset(deque<bool>& vec) { return reset(vec);}
  deque<char>  _aurostd_initialize_reset(deque<char>& vec) { return reset(vec);}
  deque<string>  _aurostd_initialize_reset(deque<string>& vec) { return reset(vec);}
  deque<int>   _aurostd_initialize_reset(deque<int>& vec) { return reset(vec);}
  deque<long>  _aurostd_initialize_reset(deque<long>& vec) { return reset(vec);}
  deque<uint>  _aurostd_initialize_reset(deque<uint>& vec) { return reset(vec);}
  deque<float> _aurostd_initialize_reset(deque<float>& vec) { return reset(vec);}
  deque<double>  _aurostd_initialize_reset(deque<double>& vec) { return reset(vec);}
  deque<long double>  _aurostd_initialize_reset(deque<long double>& vec) { return reset(vec);}

  // ***************************************************************************
  // Function reset of a vector<vector<>>
  // ***************************************************************************
  template<class utype> vector<vector<utype> > reset(vector<vector<utype> > mat) {
    for(uint i=0;i<mat.size();i++)
      for(uint j=0;j<mat[i].size();j++)
        mat[i][j]=(utype) 0;
    return mat;
  }
  // overload to force compiling
  vector<vector<bool> > _aurostd_initialize_reset(vector<vector<bool> > mat) { return reset(mat);}
  vector<vector<char> > _aurostd_initialize_reset(vector<vector<char> > mat) { return reset(mat);}
  vector<vector<string> > _aurostd_initialize_reset(vector<vector<string> > mat) { return reset(mat);}
  vector<vector<int> > _aurostd_initialize_reset(vector<vector<int> > mat) { return reset(mat);}
  vector<vector<long> > _aurostd_initialize_reset(vector<vector<long> > mat) { return reset(mat);}
  vector<vector<uint> > _aurostd_initialize_reset(vector<vector<uint> > mat) { return reset(mat);}
  vector<vector<float> > _aurostd_initialize_reset(vector<vector<float> > mat) { return reset(mat);}
  vector<vector<double> > _aurostd_initialize_reset(vector<vector<double> > mat) { return reset(mat);}
  vector<vector<long double> > _aurostd_initialize_reset(vector<vector<long double> > mat) { return reset(mat);}

  // ***************************************************************************
  // Function clear of a vector/deque
  // ***************************************************************************
  template<class utype> vector<utype> clear(vector<utype>& vec) {
    for(uint i=0;i<vec.size();i++) vec[i]=(utype) 0;
    return vec;
  }
  // overload to force compiling
  vector<bool>  _aurostd_initialize_clear(vector<bool>& vec) { return clear(vec);}
  vector<char>  _aurostd_initialize_clear(vector<char>& vec) { return clear(vec);}
  vector<string>  _aurostd_initialize_clear(vector<string>& vec) { return clear(vec);}
  vector<int>   _aurostd_initialize_clear(vector<int>& vec) { return clear(vec);}
  vector<long>  _aurostd_initialize_clear(vector<long>& vec) { return clear(vec);}
  vector<uint>  _aurostd_initialize_clear(vector<uint>& vec) { return clear(vec);}
  vector<float> _aurostd_initialize_clear(vector<float>& vec) { return clear(vec);}
  vector<double>  _aurostd_initialize_clear(vector<double>& vec) { return clear(vec);}
  vector<long double>  _aurostd_initialize_clear(vector<long double>& vec) { return clear(vec);}

  template<class utype> deque<utype> clear(deque<utype>& vec) {
    for(uint i=0;i<vec.size();i++) vec[i]=(utype) 0;
    return vec;
  }
  // overload to force compiling
  deque<bool>  _aurostd_initialize_clear(deque<bool>& vec) { return clear(vec);}
  deque<char>  _aurostd_initialize_clear(deque<char>& vec) { return clear(vec);}
  deque<string>  _aurostd_initialize_clear(deque<string>& vec) { return clear(vec);}
  deque<int>   _aurostd_initialize_clear(deque<int>& vec) { return clear(vec);}
  deque<long>  _aurostd_initialize_clear(deque<long>& vec) { return clear(vec);}
  deque<uint>  _aurostd_initialize_clear(deque<uint>& vec) { return clear(vec);}
  deque<float> _aurostd_initialize_clear(deque<float>& vec) { return clear(vec);}
  deque<double>  _aurostd_initialize_clear(deque<double>& vec) { return clear(vec);}
  deque<long double>  _aurostd_initialize_clear(deque<long double>& vec) { return clear(vec);}

  // ***************************************************************************
  // Function clear of a vector<vector<>>
  // ***************************************************************************
  template<class utype> vector<vector<utype> > clear(vector<vector<utype> > mat) {
    for(uint i=0;i<mat.size();i++)
      for(uint j=0;j<mat[i].size();j++)
        mat[i][j]=(utype) 0;
    return mat;
  }
  // overload to force compiling
  vector<vector<bool> > _aurostd_initialize_clear(vector<vector<bool> > mat) { return clear(mat);}
  vector<vector<char> > _aurostd_initialize_clear(vector<vector<char> > mat) { return clear(mat);}
  vector<vector<string> > _aurostd_initialize_clear(vector<vector<string> > mat) { return clear(mat);}
  vector<vector<int> > _aurostd_initialize_clear(vector<vector<int> > mat) { return clear(mat);}
  vector<vector<long> > _aurostd_initialize_clear(vector<vector<long> > mat) { return clear(mat);}
  vector<vector<uint> > _aurostd_initialize_clear(vector<vector<uint> > mat) { return clear(mat);}
  vector<vector<float> > _aurostd_initialize_clear(vector<vector<float> > mat) { return clear(mat);}
  vector<vector<double> > _aurostd_initialize_clear(vector<vector<double> > mat) { return clear(mat);}
  vector<vector<long double> > _aurostd_initialize_clear(vector<vector<long double> > mat) { return clear(mat);}

  // ***************************************************************************
  // Function random_shuffle of a vector/deque
  // ***************************************************************************
  template<class utype> void random_shuffle(vector<utype>& vec) {
    // switch std::random_shuffle to not deprecated std::shuffle //HE20230620
    // https://en.cppreference.com/w/cpp/algorithm/random_shuffle
    std::random_device rd;
    std::shuffle(vec.begin(),vec.end(), rd);
  }
  // overload to force compiling
  void _aurostd_initialize_random_shuffle(vector<bool>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<char>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<string>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<int>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<long>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<uint>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<float>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<double>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(vector<long double>& vec) {random_shuffle(vec);}

  template<class utype> void random_shuffle(deque<utype>& vec) {
    // switch std::random_shuffle to not deprecated std::shuffle //HE20230620
    // https://en.cppreference.com/w/cpp/algorithm/random_shuffle
    std::random_device rd;
    std::shuffle(vec.begin(),vec.end(), rd);
  }
  // overload to force compiling
  void _aurostd_initialize_random_shuffle(deque<bool>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<char>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<string>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<int>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<long>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<uint>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<float>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<double>& vec) {random_shuffle(vec);}
  void _aurostd_initialize_random_shuffle(deque<long double>& vec) {random_shuffle(vec);}

  // ***************************************************************************
  // Function isequal of vector vector
  // ***************************************************************************
  template<class utype> bool identical(vector<utype> vec1,vector<utype> vec2,utype epsilon) {
    if(vec1.size()!=vec2.size()) return FALSE;
    for(uint i=0;i<vec1.size();i++) if(aurostd::abs(vec1[i]-vec2[i])>epsilon) return FALSE;
    return TRUE;
  }
  #define AST_TEMPLATE(utype) template bool identical(vector<utype>, vector<utype>, utype);
    AST_GEN_1(AST_UTYPE_NUM)
  #undef AST_TEMPLATE

  template<class utype> bool identical(deque<utype> vec1,deque<utype> vec2,utype epsilon) {
    if(vec1.size()!=vec2.size()) return FALSE;
    for(uint i=0;i<vec1.size();i++) if(aurostd::abs(vec1[i]-vec2[i])>epsilon) return FALSE;
    return TRUE;
  }
  #define AST_TEMPLATE(utype) template bool identical(deque<utype>, deque<utype>, utype);
    AST_GEN_1(AST_UTYPE_NUM)
  #undef AST_TEMPLATE


  bool identical(vector<int> vec1,vector<int> vec2,int epsilon) {
    if(vec1.size()!=vec2.size()) return FALSE;
    for(uint i=0;i<vec1.size();i++) if(aurostd::abs(vec1[i]-vec2[i])>epsilon) return FALSE;
    return TRUE;
  }

  bool identical(deque<int> vec1,deque<int> vec2,int epsilon) {
    if(vec1.size()!=vec2.size()) return FALSE;
    for(uint i=0;i<vec1.size();i++) if(aurostd::abs(vec1[i]-vec2[i])>epsilon) return FALSE;
    return TRUE;
  }

  //   // overload to force compiling
  //  bool _aurostd_initialize_isequal(vector<bool> v1,vector<bool> v2,bool epsilon) { return isequal(v1,v2,epsilon);}
  //  bool _aurostd_initialize_isequal(vector<char> v1,vector<char> v2,char epsilon) { return isequal(v1,v2,epsilon);}
  //  bool _aurostd_initialize_isequal(vector<int> v1,vector<int> v2,int epsilon) { return isequal(v1,v2,epsilon);}
  //  bool _aurostd_initialize_isequal(vector<long> v1,vector<long> v2,long epsilon) { return isequal(v1,v2,epsilon);}
  //  bool _aurostd_initialize_isequal(vector<uint> v1,vector<uint> v2,uint epsilon) { return isequal(v1,v2,epsilon);}
  //  bool _aurostd_initialize_isequal(vector<float> v1,vector<float> v2,float epsilon) { return isequal(v1,v2,epsilon);}
  //  bool _aurostd_initialize_isequal(vector<double> v1,vector<double> v2,double epsilon) { return isequal(v1,v2,epsilon);}
  //  bool _aurostd_initialize_isequal(vector<long double> v1,vector<long double> v2,long double epsilon) { return isequal(v1,v2,epsilon);}

  // ***************************************************************************
  // Function isequal of vector<vector<>> vector<vector<>>
  // ***************************************************************************
  template<class utype> bool identical(const vector<vector<utype> >& mat1,const vector<vector<utype> >& mat2,utype epsilon) {
    if(mat1.size()!=mat2.size()) return FALSE;
    for(uint i=0;i<mat1.size();i++) {
      if(mat1[i].size()!=mat2[i].size()) return FALSE;
      for(uint j=0;j<mat1[i].size();j++)
        if(aurostd::abs(mat1[i][j]-mat2[i][j])>epsilon) return FALSE;
    }
    return TRUE;
  }

  //   // overload to force compiling
  //   bool _aurostd_initialize_isequal(const vector<vector<bool> >& m1,const vector<vector<bool> >& m2,bool epsilon) { return isequal(m1,m2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<char> >& m1,const vector<vector<char> >& m2,char epsilon) { return isequal(m1,m2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<int> >& m1,const vector<vector<int> >& m2,int epsilon) { return isequal(m1,m2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<long> >& m1,const vector<vector<long> >& m2,long epsilon) { return isequal(m1,m2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<uint> >& m1,const vector<vector<uint> >& m2,uint epsilon) { return isequal(m1,m2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<float> >& m1,const vector<vector<float> >& m2,float epsilon) { return isequal(m1,m2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<double> >& m1,const vector<vector<double> >& m2,double epsilon) { return isequal(m1,m2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<long double> >& m1,const vector<vector<long double> >& m2,long double epsilon) { return isequal(m1,m2,epsilon);}

  // ***************************************************************************
  // Function isequal of vector<vector<vector<>>> vector<vector<vector<>>>
  // ***************************************************************************
  template<class utype> bool identical(const vector<vector<vector<utype> > >& t1,const vector<vector<vector<utype> > >& t2,utype epsilon) {
    if(t1.size()!=t2.size()) return FALSE;
    for(uint i=0;i<t1.size();i++) {
      if(t1[i].size()!=t2[i].size()) return FALSE;
      for(uint j=0;j<t1[i].size();j++) {
        if(t1[i][j].size()!=t2[i][j].size()) return FALSE;
        for(uint k=0;k<t1[i][i].size();k++) 	
          if(aurostd::abs(t1[i][j][k]-t2[i][j][k])>epsilon) return FALSE;
      }
    }
    return TRUE;
  }
  //   // overload to force compiling
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<bool> > >& t1,const vector<vector<vector<bool> > >& t2,bool epsilon) { return isequal(t1,t2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<char> > >& t1,const vector<vector<vector<char> > >& t2,char epsilon) { return isequal(t1,t2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<int> > >& t1,const vector<vector<vector<int> > >& t2,int epsilon) { return isequal(t1,t2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<long> > >& t1,const vector<vector<vector<long> > >& t2,long epsilon) { return isequal(t1,t2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<uint> > >& t1,const vector<vector<vector<uint> > >& t2,uint epsilon) { return isequal(t1,t2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<float> > >& t1,const vector<vector<vector<float> > >& t2,float epsilon) { return isequal(t1,t2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<double> > >& t1,const vector<vector<vector<double> > >& t2,double epsilon) { return isequal(t1,t2,epsilon);}
  //   bool _aurostd_initialize_isequal(const vector<vector<vector<long double> > >& t1,const vector<vector<vector<long double> > >& t2,long double epsilon) { return isequal(t1,t2,epsilon);}

  // ***************************************************************************
  // Function all identical vector (one argument, checks if all entries are equal) //DX20210422
  // ***************************************************************************
  template<class utype> bool identical(const vector<utype>& vec, utype eps) {
    for(uint i=0;i<vec.size();i++){
      if(isdifferent(vec[i],vec[0],eps)){ return false; }
    }
    return true; //includes case when vec is empty
  }

  // ***************************************************************************
  // Function all identical deque (one argument, checks if all entries are equal) //DX20210422
  // ***************************************************************************
  template<class utype> bool identical(const deque<utype>& deq, utype eps) {
    for(uint i=0;i<deq.size();i++){
      if(isdifferent(deq[i],deq[0],eps)){ return false; }
    }
    return true; //includes case when deq is empty
  }

  // ***************************************************************************
  // Function toupper/tolower
  // ***************************************************************************
  string toupper(const string& in) {
    string out(in);
    for(uint i=0;i<out.length();i++)
      out[i]=std::toupper(out[i]);
    return out;
  }

  string tolower(const string& in) {
    string out(in);
    for(uint i=0;i<out.length();i++)
      out[i]=std::tolower(out[i]);
    return out;
  }

  char toupper(const char& in) {
    return std::toupper(in);
  }

  char tolower(const char& in) {
    return std::tolower(in);
  }

  // ***************************************************************************
  // Function getPWD()
  // ***************************************************************************
  // a function to get current directory
  string getPWD(){  //CO20191112
    //old way also needs PATH_LENGTH_MAX, not good for current LONG pocc directories
    //[old way - need to convert char array -> string]const int PATH_LENGTH_MAX=1024;
    //[old way - need to convert char array -> string]char work_dir[PATH_LENGTH_MAX];
    //[old way - need to convert char array -> string]getcwd(work_dir, PATH_LENGTH_MAX); 

    string pwd=aurostd::getenv2string("PWD");
    if(!pwd.empty()){return pwd;}
    return aurostd::execute2string("pwd"); //XHOST.command("pwd") ?
  }

  // ***************************************************************************
  // Function GetNumFields
  // ***************************************************************************
  // Dane Morgan
  int GetNumFields(const string& s) {
    int nf=0;
    int in_a_field=0;
    for(uint i=0;i<s.size();i++) {
      if(!(in_a_field) && s[i]!=' ') {
        in_a_field=1;
        nf++;
      }
      if(in_a_field && s[i]==' ') {
        in_a_field=0;
      }
    }
    return nf;
  }

  // ***************************************************************************
  // Function GetNextVal
  // ***************************************************************************
  // Dane Morgan - Stefano Curtarolo
  string GetNextVal(const string& s, int& id) {
    string ss;
    int i=id;
    while (s[i]==' ') {i++;} // ignore leading spaces.
    while (i<(int) s.size() && s[i]!=' ') { // pull out all text until next space.
      ss+=s[i]; i++;
    }
    id=i;
    return ss;
  }

  // ***************************************************************************
  // Function PaddedNumString
  // ***************************************************************************
  // Dane Morgan - Stefano Curtarolo
  string PaddedNumString(const int num, const int ndigits) {
    ostringstream oss;
    oss << std::setw(ndigits) << std::setfill('0') << num;// << ends;
    return oss.str();
  }

  // ***************************************************************************
  // Function getZeroPadding
  // ***************************************************************************
  // Corey Oses
  int getZeroPadding(double d) {return int(log10(d))+1;}
  int getZeroPadding(int num) {return getZeroPadding((double)num);}
  int getZeroPadding(uint num) {return getZeroPadding((double)num);}
  int getZeroPadding(long int num) {return getZeroPadding((double)num);}
  int getZeroPadding(unsigned long int num) {return getZeroPadding((double)num);}
  int getZeroPadding(long long int num) {return getZeroPadding((double)num);}
  int getZeroPadding(unsigned long long int num) {return getZeroPadding((double)num);}

  // ***************************************************************************
  // Function PaddedPRE
  // ***************************************************************************
  // Add PRE characters to pad
  string PaddedPRE(string input,int depth,string ch) {
    stringstream aus("");
    aus << input;
    string strout="";
    for(int i=0;i<depth-(int) aus.str().length();i++) { strout+=ch; }
    strout+=aus.str();
    return strout;
  } 
  template<class utype> string PaddedPRE(utype input,int depth,string ch) {
    stringstream sss;
    sss << input;
    return PaddedPRE(sss.str(),depth,ch);
  }

  // ***************************************************************************
  // Function PaddedPOST
  // ***************************************************************************
  // Add POST characters to pad
  string PaddedPOST(string input,int depth,string ch) {
    stringstream aus("");
    aus << input;
    string strout=aus.str();
    for(int i=0;i<depth-(int) aus.str().length();i++) { strout+=ch; }
    return strout;
  }
  template<class utype> string PaddedPOST(utype input,int depth,string ch) {
    stringstream sss;
    sss << input;
    return PaddedPOST(sss.str(),depth,ch);
  }

  // ***************************************************************************
  // Function PaddedCENTER
  // ***************************************************************************
  // Add PRE AND POST characters to pad so that string is in the center
  string PaddedCENTER(string input,int depth,string ch) {
    stringstream aus("");
    int pre=(depth-(int) input.length())/2;
    int post=depth-pre-(int) input.length();
    // if(DEBUG) cerr << "aurostd::PaddedCENTER: input.length()=" << input.length() << endl;
    // if(DEBUG) cerr << "aurostd::PaddedCENTER: pre=" << pre << endl;
    // if(DEBUG) cerr << "aurostd::PaddedCENTER: post=" << post << endl;   
    // if(DEBUG) cerr << "aurostd::PaddedCENTER: depth=" << depth << endl;   
    // if(DEBUG) cerr << "aurostd::PaddedCENTER: pre+post+input.length()=" << pre+post+input.length() << endl;   
    for(int i=1;i<pre;i++) { aus << ch; }  
    aus << input;
    for(int i=1;i<post;i++) { aus << ch; }  
    return aus.str();
  }
  template<class utype> string PaddedCENTER(utype input,int depth,string ch) {
    stringstream sss;
    sss << input;
    return PaddedCENTER(sss.str(),depth,ch);
  }

  // ***************************************************************************
  // Function ProgressBar
  // ***************************************************************************
  uint ProgressBar(std::ostream& oss,string prelim,uint j,uint jmax,bool VERBOSE_PERCENTAGE,bool VERBOSE_ROLLER,bool VERBOSE_CURSOR) {
    uint position=0;
    double percentage=double(j)/double(jmax);
    if(j==0) {
      oss << prelim; // position+=prelim.size();
    }
    // VERBOSE PERCENTAGE
    if(!mod<uint>(j,50) || j==jmax-1 || j==jmax) {
      if(VERBOSE_PERCENTAGE) {
        if(j==jmax-1 || j==jmax)  {
          if(0) {oss << "[100.%]";}
          if(1) {oss << "[100.0%]";}
          if(0) {oss << "[100.00%]";}
          position+=8;
        } else {
          if(0) {
            // 99.9%
            oss << "[" << (percentage<0.1?" ":"")
              << mod<uint>(uint(percentage*100),100) << "."
              << mod<uint>(uint(percentage*1000),10)
              << "%]";
            position+=7;
          }
          if(1) {
            // 99.99%
            oss << "[" << (percentage<0.1?" ":"")
              << mod<uint>(uint(percentage*100),100) << "."
              << mod<uint>(uint(percentage*1000),10)
              << mod<uint>(uint(percentage*10000),10)
              << "%]";
            position+=8;
          }
          if(0) {
            // 99.999%
            oss << "[" << (percentage<0.1?" ":"")
              << mod<uint>(uint(percentage*100),100) << "."
              << mod<uint>(uint(percentage*1000),10)
              << mod<uint>(uint(percentage*10000),10)
              << mod<uint>(uint(percentage*100000),10)
              << "%]";
            position+=9;
          }
        }  
        oss << " ";position++;
      }
    }
    // VERBOSE_ROLLER
    if(!mod<uint>(j,50) || j==jmax-1 || j==jmax) {
      if(VERBOSE_ROLLER) {
        if(j==jmax-1 || j==jmax)  {
          oss << "[=]";
        } else {
          if(mod<uint>(j/513,4)==0) oss << "[\\]";
          if(mod<uint>(j/513,4)==1) oss << "[|]";
          if(mod<uint>(j/513,4)==2) oss << "[/]";
          if(mod<uint>(j/513,4)==3) oss << "[-]";
        }
        position+=3;
        oss << " ";position++;
      }
    }
    // VERBOSE CURSOR
    if(j==0 || !mod<uint>(j,50) || j==jmax-1 || j==jmax) {
      if(VERBOSE_CURSOR) {
        if(j==jmax-1 || j==jmax)  {
          oss << "[======================================================================================================]";
          position+=102;
        } else {
          oss << "["; position++;
          for(double k=0;k<percentage*100;k+=1.0) {oss << "=";position++;}
          // if(mod<uint>(j,500)==0)
          {
            if(mod<uint>(j/478,4)==0) {oss << "\\";position++;}
            if(mod<uint>(j/478,4)==1) {oss << "|";position++;}
            if(mod<uint>(j/478,4)==2) {oss << "/";position++;}
            if(mod<uint>(j/478,4)==3) {oss << "-";position++;}
          }
          if(j==0)
          {
            for(double k=0;k<(1.0-percentage)*100.0+0.01;k+=1.0) {oss << " ";position++;}
            oss << "]";position++;
          }
        }
        oss << " ";position++;
      }
      // NOW GO BACK
      for(uint k=0;k<position;k++) {oss << "\b";}
      if(j==jmax-1 || j==jmax) oss << endl;
    }
    return position;
  }

  uint ProgressBar(std::ostream& oss,string prelim,uint j,uint jmax) {
    return ProgressBar(oss,prelim,j,jmax,TRUE,TRUE,TRUE);
  }

  uint ProgressBar(std::ostream& oss,string prelim,double j,bool VERBOSE_PERCENTAGE,bool VERBOSE_ROLLER,bool VERBOSE_CURSOR) {
    return ProgressBar(oss,prelim,uint(double(j*100)),100,VERBOSE_PERCENTAGE,VERBOSE_ROLLER,VERBOSE_CURSOR);
  }

  uint ProgressBar(std::ostream& oss,string prelim,double j) {
    return ProgressBar(oss,prelim,uint(double(j*100)),100,TRUE,TRUE,TRUE);
  }

  // ***************************************************************************
  // Function PercentEncodeASCII //DX20210706
  // ***************************************************************************
  // Converts a single ASCII character into its percent-encoded form
  string PercentEncodeASCII(const char c) {
    stringstream char_percent_encoded;
    char_percent_encoded << "%" << std::hex << (int)c;
    return char_percent_encoded.str();
  }

  // ***************************************************************************
  // Function CleanStringASCII
  // ***************************************************************************
  // Clean a string from ASCII junk
  // Stefano Curtarolo
  string CleanStringASCII(const string& s) {return CleanStringASCII_20190712(s);} //CO20190712
  string CleanStringASCII_20190712(const string& s) { //CO20190712
    string ss=s;
    CleanStringASCII_InPlace(ss);
    return ss;
  }
  string CleanStringASCII_20190101(const string& s) { //CO20190712
    string ss="";
    for(uint i=0;i<s.length();i++) {
      if(s[i]>='A' && s[i]<='Z') ss+=s[i];  // LETTERS
      else if(s[i]>='a' && s[i]<='z') ss+=s[i];  // letters
      else if(s[i]>='0' && s[i]<='9') ss+=s[i];  // numbers
      else if(s[i]=='.' || s[i]=='+' || s[i]=='-' || s[i]=='*' || s[i]=='/') ss+=s[i];  // operations
      else if(s[i]=='_' || s[i]=='#' || s[i]=='&' || s[i]==':' || s[i]==',' || s[i]=='@' || s[i]=='$') ss+=s[i];  // underscore
      else if(s[i]=='=' || s[i]=='|' || s[i]=='\'' || s[i]=='\"' || s[i]==' ') ss+=s[i];  // underscore
    }
    return ss;
  }

  // ***************************************************************************
  // Function CleanStringASCII_InPlace
  // ***************************************************************************
  // Similar to CleanStringASCII, but does NOT create a new string (costly if done MANY times)
  // Corey Oses 20190712
  void CleanStringASCII_InPlace(string& s) {
    //[CO20190712 - slight optimization if we go backwards]for(uint i=0;i<s.length();i++)
    for(uint i=s.length()-1;i<s.length();i--)
    { //CO20200106 - patching for auto-indenting
      //[CO20200624 - not inclusive enough]if(!(
      //[CO20200624 - not inclusive enough]      (s[i]>='A' && s[i]<='Z') || //LETTERS
      //[CO20200624 - not inclusive enough]      (s[i]>='a' && s[i]<='z') || //letters
      //[CO20200624 - not inclusive enough]      (s[i]>='0' && s[i]<='9') || //numbers
      //[CO20200624 - not inclusive enough]      (s[i]=='.' || s[i]=='+' || s[i]=='-' || s[i]=='*' || s[i]=='/') ||  //operations
      //[CO20200624 - not inclusive enough]      (s[i]=='_' || s[i]=='#' || s[i]=='&' || s[i]==':' || s[i]==',' || s[i]=='@' || s[i]=='$') ||  //punctuation1
      //[CO20200624 - not inclusive enough]      (s[i]=='=' || s[i]=='|' || s[i]=='\'' || s[i]=='\"' || s[i]==' ') ||  //punctuation2
      //[CO20200624 - not inclusive enough]      FALSE)
      //[CO20200624 - not inclusive enough]  ){RemoveCharacterInPlace(s,s[i]);}
      //https://stackoverflow.com/questions/48212992/how-to-find-out-if-there-is-any-non-ascii-character-in-a-string-with-a-file-path
      //cerr << s[i] << " " << static_cast<unsigned int>(s[i]) << endl;
      if(static_cast<unsigned int>(s[i])>127){RemoveCharacterInPlace(s,s[i]);}
    }
  }

  // ***************************************************************************
  // Function RemoveTrailingCharacter
  // ***************************************************************************
  // Removes trailing character
  // CO+ME20200825
  string RemoveTrailingCharacter(const string& s,char c){
    string ss=s;
    RemoveTrailingCharacter_InPlace(ss,c);
    return ss;
  }

  // ***************************************************************************
  // Function RemoveTrailingCharacter_InPlace
  // ***************************************************************************
  // Similar to RemoveTrailingCharacter, but does NOT create a new string (costly if done MANY times)
  // CO+ME20200825
  void RemoveTrailingCharacter_InPlace(string& s,char c){
    while(s.size()>0 && s.at(s.size()-1)==c){s=s.substr(0,s.size()-1);}
  }

  //DX20190516 - remove control code characters - START
  // ***************************************************************************
  // Function removeControlCodeCharactersFromString
  // ***************************************************************************
  bool RemoveControlCodeCharactersFromString(const string& in, string& out){  //CO20190620

    // removes control code and backspace characters (e.g., NUL, DEL, etc.)
    // only keep printable characters (i.e., digits, letters, punctuation, and spaces) 
    // and white space characters (i.e., space, newline, tabs, and carrage returns)
    // a boolean indicates if the stringstream contained a control code character
    // string input version

    stringstream ss_in, ss_out;
    ss_in << in;
    bool detected_control_char = RemoveControlCodeCharactersFromStringstream(ss_in, ss_out);
    out = ss_out.str();
    return detected_control_char;
  }

  // ***************************************************************************
  // Function removeControlCodeCharactersFromStringStream
  // ***************************************************************************
  bool RemoveControlCodeCharactersFromStringstream(std::stringstream& ss_in, std::stringstream& ss_out){

    // removes control code and backspace characters (e.g., NUL, DEL, etc.)
    // only keep printable characters (i.e., digits, letters, punctuation, and spaces) 
    // and white space characters (i.e., space, newline, tabs)
    // a boolean indicates if the stringstream contained a control code character
    // stringstream input version
    //
    //ME20190614: We don't want carriage returns either because they mess up string additions.
    // Since they point to the beginning of the string, adding to a string with a carriage
    // return would overwrite instead of append

    bool detected_control_char = false;
    char c;
    //char c1;
    //char c2;

    //stringstream tmp; tmp << ss_in.str();
    while(ss_in.get(c)){
      //[CO20190620 - still doesn't work]if(isprint(c) || isspace(c) || (c != '\r')) {  //ME20190614 //[CO20200106 - close bracket for indenting]}
      if((isprint(c) || isspace(c) || FALSE) && ((c != '\r') || FALSE)) {ss_out << c;}  //CO20190620 - add more cases before FALSE
      else{detected_control_char = true;}
    }

    return detected_control_char;
  }
  //DX20190516 - remove control code characters - END

  //DX20190211 - remove control code characters from file - START
  // ***************************************************************************
  // Function RemoveControlCodeCharactersFromFile
  // ***************************************************************************
  bool RemoveControlCodeCharactersFromFile(const string& directory,const string& filename, bool keep_orig_file){  //CO20210315

    // removes control code and backspace characters (e.g., NUL, DEL, etc.)
    // overwrites file if control characters are detected, otherwise the file is untouched (preserve original timestamp)
    // matches original compression

    if(aurostd::FileExist(directory+"/"+filename)){
      stringstream ss_in, ss_out;
      aurostd::efile2stringstream(directory+"/"+filename, ss_in);
      // if file contains control code characters, then overwrite file
      if(RemoveControlCodeCharactersFromStringstream(ss_in, ss_out)){
        string uncompressed_filename = "";
        //compressed files
        if(IsCompressed(filename,uncompressed_filename)){
          stringstream2file(ss_out,directory+"/"+uncompressed_filename+"_tmp");
          string extension = GetCompressionExtension(filename);
          CompressFile(directory+"/"+uncompressed_filename+"_tmp", extension);
          if(keep_orig_file){ file2file(directory+"/"+filename,directory+"/"+uncompressed_filename+"_old"+extension); } //move original
          file2file(directory+"/"+uncompressed_filename+"_tmp"+extension,directory+"/"+filename); //overwrite
        }
        //uncompressed files
        else{
          stringstream2file(ss_out,directory+"/"+filename+".tmp");
          if(keep_orig_file){ file2file(directory+"/"+filename,directory+"/"+filename+"_old"); } //move original
          file2file(directory+"/"+filename+".tmp",directory+"/"+filename); //overwrite
        }
        return true;
      }
      // file is ok, do not update
      else{ return false; } //signals file is unchanged
    } else {
      string message = "File does not exist: " + directory + "/" + filename;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    return false;
  }
  //DX20190211 - remove control code characters from file - END

  //DX20190125 - remove null bytes - START
  // ***************************************************************************
  // Function isNullbyte
  // ***************************************************************************
  // Deterine if char is a null byte (e.g., ^@)
  bool isNullByte(char c){
    return (c=='\0');
  }

  // ***************************************************************************
  // Function removeNullBytes
  // ***************************************************************************
  // Remove all null bytes in string (e.g., ^@)
  string removeNullBytes(string in){
    string out=in;
    out.erase(remove_if(out.begin(),out.end(),isNullByte), out.end());
    return out;
  }
  //DX20190125 - remove null bytes - END

  //DX20190211 - remove null characters from file - START
  // ***************************************************************************
  // Function RemoveBinaryCharactersFromFile()
  // ***************************************************************************
  // Remove all null bytes in file
  bool RemoveBinaryCharactersFromFile(const string& directory,const string& filename){  //CO20210315
    stringstream aus_exec;
    deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");vext.push_front(""); // cheat for void string
    deque<string> vcmd; aurostd::string2tokens("cat,bzcat,xzcat,gzcat",vcmd,",");
    deque<string> vzip; aurostd::string2tokens("bzip2,xz,gzip",vzip,",");vzip.push_front(""); // cheat for void string
    if(vext.size()!=vcmd.size()) {
      string message = "vext.size()!=vcmd.size()";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }

    for(uint iext=0;iext<vext.size();iext++){ // check filename.EXT
      if(aurostd::FileExist(directory+"/" + filename + vext.at(iext))){
        aus_exec << "cd \"" << directory << "\"" << endl;
        aus_exec << vcmd[iext] << " " << filename << vext[iext] << " | sed \"s/[^[:print:]\\r\\t]//g\" > " << filename << ".tmp && mv " << filename << ".tmp " << filename << endl;
        if(vext[iext]!=""){
          aus_exec << vzip[iext] << " " << filename << endl;
        }
        aurostd::execute(aus_exec);
      }
    }
    return true;
  }
  //DX20190211 - remove binary characters from file - END

  // ***************************************************************************
  // Function CGI_StringClean
  // ***************************************************************************
  // Clean a string from CGI junk
  string CGI_StringClean(const string& stringIN) {
    string stringOUT=stringIN;
    aurostd::StringSubst(stringOUT,"%0D%0A","\n");aurostd::StringSubst(stringOUT,"%0d%0a","\n");   // newlines
    aurostd::StringSubst(stringOUT,"+"," ");    // spaces
    aurostd::StringSubst(stringOUT,"%28","(");aurostd::StringSubst(stringOUT,"%29",")");   // ()
    aurostd::StringSubst(stringOUT,"%5B","[");aurostd::StringSubst(stringOUT,"%5D","]");   // []
    aurostd::StringSubst(stringOUT,"%7B","{");aurostd::StringSubst(stringOUT,"%7D","}");   // brackets (do not write, it screws up indent)
    aurostd::StringSubst(stringOUT,"%2B","+");aurostd::StringSubst(stringOUT,"%2F","/");   //  operations
    aurostd::StringSubst(stringOUT,"%23","#");aurostd::StringSubst(stringOUT,"%21","!");
    aurostd::StringSubst(stringOUT,"%3F","?");aurostd::StringSubst(stringOUT,"%2C",",");
    aurostd::StringSubst(stringOUT,"%3A",":");aurostd::StringSubst(stringOUT,"%3B",";");
    aurostd::StringSubst(stringOUT,"%27","'");aurostd::StringSubst(stringOUT,"%22","\"");
    aurostd::StringSubst(stringOUT,"%60","`");aurostd::StringSubst(stringOUT,"%40","@");
    aurostd::StringSubst(stringOUT,"%24","$");aurostd::StringSubst(stringOUT,"%25","%");
    aurostd::StringSubst(stringOUT,"%5E","^");aurostd::StringSubst(stringOUT,"%26","&");
    aurostd::StringSubst(stringOUT,"%3D","=");aurostd::StringSubst(stringOUT,"%7E","~");
    aurostd::StringSubst(stringOUT,"%5C","\\");aurostd::StringSubst(stringOUT,"%7C","|");
    aurostd::StringSubst(stringOUT,"%3C","<");aurostd::StringSubst(stringOUT,"%3E",">");  // <>
    aurostd::StringSubst(stringOUT,"\n\n","\n");
    return stringOUT;
  }

  // ***************************************************************************
  // Function RemoveWhiteSpaces
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Morgan / Curtarolo
  string RemoveWhiteSpaces(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss;
    for (uint i=0;i<s.size();i++) if(s[i]!=' ' && s[i]!='\t') ss+=s[i];
    return ss;
  }
  string RemoveWhiteSpaces(const string& s, const char toggle) {  //CO20190710
    if(s.size()==0) return s;  // nothing to do
    string ss;
    bool copy=TRUE;
    for (uint i=0;i<s.size();i++) {
      if(s[i]==toggle){copy=!copy;}  //CO20190710
      if(copy){if(s[i]!=' ' && s[i]!='\t'){ss+=s[i];}}
      else{ss+=s[i];}
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveWhiteSpacesFromTheBack
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Morgan / Curtarolo
  string RemoveWhiteSpacesFromTheBack(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    while(ss[ss.size()-1]==' ' || ss[ss.size()-1]=='\t') {
      ss.erase(ss.size()-1,1);
      if(ss.size()==0) return ss;  // nothing to do
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveWhiteSpacesFromTheFront
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Oses
  string RemoveWhiteSpacesFromTheFront(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    while(ss[0]==' ' || ss[ss.size()-1]=='\t') {
      ss.erase(0,1);
      if(ss.size()==0) return ss;  // nothing to do
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveWhiteSpacesFromTheFrontAndBack
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Oses
  string RemoveWhiteSpacesFromTheFrontAndBack(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    ss=RemoveWhiteSpacesFromTheBack(ss);
    ss=RemoveWhiteSpacesFromTheFront(ss);
    return ss;
  }

  // ***************************************************************************
  // Function RemoveSpaces
  // ***************************************************************************
  // Removes all spaces from a string. Morgan / Curtarolo
  string RemoveSpaces(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss;
    for (uint i=0;i<s.size();i++) if(s[i]!=' ') ss+=s[i];
    return ss;
  }
  string RemoveSpaces(const string& s, const char toggle) { //CO20190710
    if(s.size()==0) return s;  // nothing to do
    string ss;
    bool copy=TRUE;
    for (uint i=0;i<s.size();i++) {
      if(s[i]==toggle){copy=!copy;}  //CO20190710
      if(copy){if(s[i]!=' '){ss+=s[i];}}
      else{ss+=s[i];}
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveSpacesFromTheBack
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Morgan / Curtarolo
  string RemoveSpacesFromTheBack(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    while(ss[ss.size()-1]==' ') {
      ss.erase(ss.size()-1,1);
      if(ss.size()==0) return ss;  // nothing to do
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveTabs
  // ***************************************************************************
  // Removes all tabs from a string. Stefano Curtarolo
  string RemoveTabs(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss;
    for (uint i=0;i<s.size();i++) if(s[i]!='\t') ss+=s[i];
    return ss;
  }
  string RemoveTabs(const string& s, const char toggle) { //CO20190710
    if(s.size()==0) return s;  // nothing to do
    string ss;
    bool copy=TRUE;
    for (uint i=0;i<s.size();i++) {
      if(s[i]==toggle){copy=!copy;}  //CO20190710
      if(copy){if(s[i]!='\t'){ss+=s[i];}}
      else{ss+=s[i];}
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveTabsFromTheBack
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string.
  // Dane Morgan / Stefano Curtarolo
  string RemoveTabsFromTheBack(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    while(ss[ss.size()-1]=='\t') {
      ss.erase(ss.size()-1,1);
      if(ss.size()==0) return ss;  // nothing to do
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveComments
  // ***************************************************************************
  // Removes all comments from a string.
  // Stefano Curtarolo
  //   string RemoveComments(const string& s) {
  //     if(s.size()==0) return s;  // nothing to do
  //     string ss;
  //     bool copy=TRUE;
  //     for (uint i=0;i<s.size();i++) {
  //       if(s[i]=='#')  copy=FALSE;
  //       if(s[i]=='\n') copy=TRUE;
  //       if(copy) ss+=s[i];
  //     }
  //     return ss;
  //   }

  //ME20190614 - added vector<string> version of RemoveComments
  vector<string> RemoveComments(const vector<string>& vstrin) { //CO20210315 - cleaned up
    vector<string> vstrout;
    string::size_type loc;
    string line="";
    for (uint i = 0; i < vstrin.size(); i++) {
      line = vstrin[i];
      // COMMENT_NEGLECT_1
      loc = line.find(COMMENT_NEGLECT_1);
      while (loc != string::npos) {
        // Do not remove #[1-9] since it is not a comment (spacegroup)
        if (!((loc > 0) && (loc < line.size()) && (isdigit(line[loc+1])))) {
          line = line.substr(0, loc);
          break;
        }
        loc = line.find(COMMENT_NEGLECT_1, loc + 1);
      }
      // COMMENT_NEGLECT_2
      loc = line.find(COMMENT_NEGLECT_2);
      while (loc != string::npos) {
        // Do not remove :// since it is not a comment (web address)
        if (!((loc > 0) && (loc < line.size()) && (line[loc-1] == ':'))) {
          line = line.substr(0, loc);
          break;
        }
        loc = line.find(COMMENT_NEGLECT_2, loc + 1);
      }
      // COMMENT_NEGLECT_3
      loc = line.find(COMMENT_NEGLECT_3);
      line = line.substr(0, loc);
      if (!line.empty()) vstrout.push_back(line);
    }
    return vstrout;
  }
  deque<string> RemoveComments(const deque<string>& vstrin) {
    return aurostd::vector2deque(RemoveComments(aurostd::deque2vector(vstrin)));
  }

  string RemoveComments(const string& strin) {  //CO20210315 - cleaned up
    vector<string> vlines;
    aurostd::string2vectorstring(strin, vlines);
    vlines = RemoveComments(vlines);
    if(vlines.size()==0){return "";}
    if(vlines.size()==1){return vlines[0];}
    string strout="";
    for (uint i = 0; i < vlines.size(); i++) strout += vlines[i] + '\n';
    return strout;
  }

  //[OBSOLETE]  string RemoveComments(const string &strin) {
  //[OBSOLETE]    string strout=strin;
  //[OBSOLETE]    vector<string> vstrout;aurostd::string2vectorstring(strout,vstrout);
  //[OBSOLETE]    strout.clear();
  //[OBSOLETE]    string::size_type loc;  //CO20180409, don't do find twice (expensive)
  //[OBSOLETE]    for(uint i=0;i<vstrout.size();i++) {
  //[OBSOLETE]      //COMMENT_NEGLECT_1
  //[OBSOLETE]      loc=vstrout[i].find(COMMENT_NEGLECT_1);  //CO20180409
  //[OBSOLETE]      vstrout[i]=vstrout[i].substr(0,loc);  //no NEED TO ask if()..., it will be set to npos anyway
  //[OBSOLETE]
  //[OBSOLETE]      //COMMENT_NEGLECT_2, but not ":"+//COMMENT_NEGLECT_1
  //[OBSOLETE]      loc=vstrout[i].find(COMMENT_NEGLECT_2);  //CO20180409
  //[OBSOLETE]      while(loc!=string::npos){
  //[OBSOLETE]        if(!(loc>0&&loc<vstrout[i].size()&&vstrout[i].at(loc-1)==':')){  //find the NOT case where we are in the range and we find ':' before comment
  //[OBSOLETE]          vstrout[i]=vstrout[i].substr(0,loc);
  //[OBSOLETE]          break;
  //[OBSOLETE]        }
  //[OBSOLETE]        loc=vstrout[i].find(COMMENT_NEGLECT_2,loc+1);
  //[OBSOLETE]      }
  //[OBSOLETE]
  //[OBSOLETE]      //COMMENT_NEGLECT_3
  //[OBSOLETE]      loc=vstrout[i].find(COMMENT_NEGLECT_3);  //CO20180409
  //[OBSOLETE]      vstrout[i]=vstrout[i].substr(0,loc);  //no NEED TO ask if()..., it will be set to npos anyway
  //[OBSOLETE]
  //[OBSOLETE]if(vstrout[i].find(COMMENT_NEGLECT_1)!=string::npos) 
  //[OBSOLETE]  vstrout[i]=vstrout[i].substr(0,vstrout[i].find(COMMENT_NEGLECT_1));  
  //[OBSOLETE]if(vstrout[i].find(COMMENT_NEGLECT_2)!=string::npos && vstrout[i].find(":"+COMMENT_NEGLECT_2)==string::npos)  // look for // but dont touch ://
  //[OBSOLETE]  vstrout[i]=vstrout[i].substr(0,vstrout[i].find(COMMENT_NEGLECT_2));
  //[OBSOLETE]if(vstrout[i].find(COMMENT_NEGLECT_3)!=string::npos)
  //[OBSOLETE]  vstrout[i]=vstrout[i].substr(0,vstrout[i].find(COMMENT_NEGLECT_3));
  //[OBSOLETE]  if(!vstrout[i].empty()) cout << vstrout[i] << endl;
  //[OBSOLETE]      if(!vstrout[i].empty()) strout+=vstrout[i]+"\n";
  //[OBSOLETE]    }  
  //[OBSOLETE]    return strout;
  //[OBSOLETE]  }

  // ***************************************************************************
  // Function RemoveCharacter
  // ***************************************************************************
  // Removes charecters from string
  // Stefano Curtarolo
  string RemoveCharacter(const string& s, const char character) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    RemoveCharacterInPlace(ss,character);
    //[CO20190712 - OBSOLETE with RemoveCharacterInPlace()]for (uint i=0;i<s.size();i++) {
    //[CO20190712 - OBSOLETE with RemoveCharacterInPlace()]  if(s[i]!=character) ss+=s[i];
    //[CO20190712 - OBSOLETE with RemoveCharacterInPlace()]}
    return ss;
  }

  // ***************************************************************************
  // Function RemoveCharacterInPlace
  // ***************************************************************************
  // Similar to RemoveCharacter, but does NOT create a new string (costly if done MANY times)
  // Corey Oses 20190712
  void RemoveCharacterInPlace(string& t, const char character) {t.erase(std::remove(t.begin(), t.end(), character), t.end());}

  // ***************************************************************************
  // Function RemoveCharacterFromTheBack
  // ***************************************************************************
  // Remove character from the back of a string. DX (Hicks)
  string RemoveCharacterFromTheBack(const string& s, const char character) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    if(ss[ss.size()-1]==character) {
      ss.erase(ss.size()-1,1);
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveCharacterFromTheFront
  // ***************************************************************************
  // Removes character from the front of a string. DX (Hicks)
  string RemoveCharacterFromTheFront(const string& s, const char character) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    if(ss[0]==character) {
      ss.erase(0,1);
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveCharacterFromTheFrontAndBack
  // ***************************************************************************
  // Removes character from the front and back of a string. DX (Hicks)
  string RemoveCharacterFromTheFrontAndBack(const string& s, const char character) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    ss=RemoveCharacterFromTheBack(ss,character);
    if(s.size()==0) return s; // cannot remove anything else
    ss=RemoveCharacterFromTheFront(ss,character);
    return ss;
  }

  // ***************************************************************************
  // Function RemoveNumbers
  // ***************************************************************************
  // Removes numbers from string
  // Stefano Curtarolo
  string RemoveNumbers(const string& s) {return RemoveNumbers_20190712(s);} //CO20190712
  string RemoveNumbers_20190712(const string& s) {  //CO20190712 - avoids creating many copies of string
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    RemoveNumbersInPlace(ss);
    return ss;
  }
  string RemoveNumbers_20190101(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    ss=RemoveCharacter(ss,'0');ss=RemoveCharacter(ss,'1');ss=RemoveCharacter(ss,'2');
    ss=RemoveCharacter(ss,'3');ss=RemoveCharacter(ss,'4');ss=RemoveCharacter(ss,'5');
    ss=RemoveCharacter(ss,'6');ss=RemoveCharacter(ss,'7');ss=RemoveCharacter(ss,'8');
    ss=RemoveCharacter(ss,'9');ss=RemoveCharacter(ss,'.');
    return ss;
  }

  // ***************************************************************************
  // Function RemoveNumbersInPlace
  // ***************************************************************************
  // Similar to RemoveNumbers, but does NOT create a new string (costly if done MANY times)
  // Corey Oses 20190712
  void RemoveNumbersInPlace(string& s) { //CO20190712
    RemoveCharacterInPlace(s,'0');RemoveCharacterInPlace(s,'1');RemoveCharacterInPlace(s,'2');
    RemoveCharacterInPlace(s,'3');RemoveCharacterInPlace(s,'4');RemoveCharacterInPlace(s,'5');
    RemoveCharacterInPlace(s,'6');RemoveCharacterInPlace(s,'7');RemoveCharacterInPlace(s,'8');
    RemoveCharacterInPlace(s,'9');RemoveCharacterInPlace(s,'.');
  }

  // ***************************************************************************
  // Function RemoveRounding
  // ***************************************************************************
  // Removes rounding from string
  // Stefano Curtarolo
  string RemoveRounding(const string& s) {
    if(s.size()==0) return s;  // nothing to do
    string ss=s;
    ss=RemoveSubString(ss,"(0)");ss=RemoveSubString(ss,"(1)");ss=RemoveSubString(ss,"(2)");
    ss=RemoveSubString(ss,"(3)");ss=RemoveSubString(ss,"(4)");ss=RemoveSubString(ss,"(5)");
    ss=RemoveSubString(ss,"(6)");ss=RemoveSubString(ss,"(7)");ss=RemoveSubString(ss,"(8)");
    ss=RemoveSubString(ss,"(9)");
    return ss;
  }

  // ***************************************************************************
  // Function RemoveSubStringFirst
  // ***************************************************************************
  // Removes the first substring from string
  // Stefano Curtarolo
  string RemoveSubStringFirst(const string& str_orig, const string& str_rm) {
    string t=str_orig;
    RemoveSubStringFirstInPlace(t,str_rm);  //CO20190712
    //[CO20190712 - moved to RemoveSubStringFirstInPlace()]std::string::size_type i = t.find(str_rm);
    //[CO20190712 - moved to RemoveSubStringFirstInPlace()]if(i != std::string::npos)
    //[CO20190712 - moved to RemoveSubStringFirstInPlace()]  t.erase(i, str_rm.length( ));
    return t;
  }

  // ***************************************************************************
  // Function RemoveSubStringFirstInPlace
  // ***************************************************************************
  // Similar to RemoveSubStringFirst, but does NOT create a new string (costly if done MANY times)
  // Corey Oses 20190712
  void RemoveSubStringFirstInPlace(string& t, const string& str_rm) {
    std::string::size_type i = t.find(str_rm);
    if(i != std::string::npos)
      t.erase(i, str_rm.length( ));
  }

  // ***************************************************************************
  // Function RemoveSubString
  // ***************************************************************************
  // Removes all instances of substring from string
  // Stefano Curtarolo
  string RemoveSubString(const string& str_orig, const string& str_rm) {
    string t=str_orig;
    RemoveSubStringInPlace(t,str_rm); //CO20190712
    //[CO20190712 - moved to RemoveSubStringInPlace()]string::size_type i;
    //[CO20190712 - moved to RemoveSubStringInPlace()]while(t.find(str_rm)!=string::npos) {
    //[CO20190712 - moved to RemoveSubStringInPlace()]  i = t.find(str_rm);
    //[CO20190712 - moved to RemoveSubStringInPlace()]  if(i != std::string::npos) t.erase(i, str_rm.length( ));
    //[CO20190712 - moved to RemoveSubStringInPlace()]};
    return t;
  }

  // ***************************************************************************
  // Function RemoveSubStringInPlace
  // ***************************************************************************
  // Similar to RemoveSubString, but does NOT create a new string (costly if done MANY times)
  // Corey Oses 20190712
  void RemoveSubStringInPlace(string& t, const string& str_rm) {
    string::size_type i = t.find(str_rm); //CO20190712 - fewer operations
    while(i != string::npos) {  //CO20190712 - fewer operations
      t.erase(i, str_rm.length( )); //CO20190712 - fewer operations
      i = t.find(str_rm); //CO20190712 - fewer operations
    }
    //[CO20190712 - OBSOLETE]while(t.find(str_rm)!=string::npos) {
    //[CO20190712 - OBSOLETE]  i = t.find(str_rm);
    //[CO20190712 - OBSOLETE]  if(i != std::string::npos) t.erase(i, str_rm.length( ));
    //[CO20190712 - OBSOLETE]};
  }

  // ***************************************************************************
  // Function VersionString2Double
  // ***************************************************************************
  // 5.1.311 -> 5.0013311
  // 4.2.34 -> 4.002034
  double VersionString2Double(const string& version_str){ //SD20220331 
    vector<string> tokens;
    aurostd::string2tokens(version_str,tokens,".");
    double version=0.0;;
    for (uint i=0;i<tokens.size();i++){version+=aurostd::string2utype<double>(tokens[i])*std::pow(10.0,-3.0*i);}
    return version;
  }

  // ***************************************************************************
  // Function ProcessPIDs
  // ***************************************************************************
  //CO20210315
  vector<string> ProcessPIDs(const string& process,bool user_specific){ //CO20210315
    string output_syscall="";
    return ProcessPIDs(process,output_syscall,user_specific);
  }
  vector<string> ProcessPIDs(const string& process,string& output_syscall,bool user_specific){ //CO20210315
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking for process=" << process << endl;}
    if(0){cerr << __AFLOW_FUNC__ << " ps table:" << endl << aurostd::execute2string("ps aux") << endl;}  //not a good idea to run this all the time

    string command="";
    vector<string> vlines,vtokens,vpids;
    uint i=0,j=0;
    if(aurostd::IsCommandAvailable("pgrep")) {
      string command_pgrep="pgrep -a";  //needed over -l on fossies installs: https://fossies.org/linux/procps-ng/pgrep.1
      if(!aurostd::execute2string(command_pgrep+" test",stderr_fsio).empty()){command_pgrep="pgrep -l";}  //should work on all linux  //the "test" is a dummy to see if the pgrep command works
      if(aurostd::execute2string(command_pgrep+" test",stderr_fsio).empty()){ //the "test" is a dummy to see if the pgrep command works
        command=command_pgrep; //the -a/-l is important, we will need to neglect the subshell call below
        if(user_specific && !XHOST.user.empty()){command+=" -u "+XHOST.user;}
        command+=" -f "+process+" 2> /dev/null";  //the -f is important, will match mpivasp46s in /usr/bin/mpivasp46s
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " running command=\"" << command << "\"" << endl;}
        string output=output_syscall=aurostd::execute2string(command);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " pgrep output:" << endl << "\"" << output << "\"" << endl;}
        if(0){  //before -f and -l
          aurostd::StringSubst(output,"\n"," ");
          output=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(output);
          if(output.empty()){return vpids;}
          aurostd::string2tokens(output,vpids," ");
        }
        aurostd::string2vectorstring(output,vlines);
        for(i=0;i<vlines.size();i++){
          aurostd::string2tokens(vlines[i],vtokens," ");
          if(vtokens.size()<2){continue;}
          const string& pid=vtokens[0];
          string proc=vtokens[1]; //since we split on " ", we need to join columns 11-onward
          for(j=2;j<vtokens.size();j++){proc+=" "+vtokens[j];}
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " proc[i=" << i << "]=\"" << proc << "\"" << endl;}
          if(proc.find(process)==string::npos){continue;}
          if(proc.find(command)!=string::npos){continue;} //ps aux | grep ... always returns itself, neglect  //do a find() instead of == here
          vpids.push_back(pid);
        }
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " vpids=" << aurostd::joinWDelimiter(vpids,",") << endl;
          cerr << __AFLOW_FUNC__ << " vpids.empty()=" << vpids.empty() << endl;
        }
        return vpids;
      }
    }

    if(aurostd::IsCommandAvailable("ps") && aurostd::IsCommandAvailable("grep")) {
      //FR recommends ps aux vs. ps -e
      //tested on linux and mac, PIDs are in second column, process is the last column
      string command_grep="grep "+process;
      command="ps";
      if(user_specific){command+=" ux";}
      else{command+=" aux";}
      command+=" 2>/dev/null | "+command_grep+" 2> /dev/null";
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " running command=\"" << command << "\"" << endl;}
      string output=output_syscall=aurostd::execute2string(command);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " ps/grep output:" << endl << output << endl;}
      aurostd::string2vectorstring(output,vlines);
      for(i=0;i<vlines.size();i++){
        aurostd::string2tokens(vlines[i],vtokens," ");
        if(vtokens.size()<11){continue;}
        const string& pid=vtokens[1];
        string proc=vtokens[10]; //since we split on " ", we need to join columns 11-onward
        for(j=11;j<vtokens.size();j++){proc+=" "+vtokens[j];}
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " proc[i=" << i << "]=\"" << proc << "\"" << endl;}
        if(proc.find(process)==string::npos){continue;}
        if(proc.find(command)!=string::npos){continue;} //ps aux | grep ... always returns itself, neglect  //do a find() instead of == here
        if(proc==command_grep){continue;} //ps aux | grep ... always returns itself, neglect  //do a == instead a find() here
        vpids.push_back(pid);
      }
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " vpids=" << aurostd::joinWDelimiter(vpids,",") << endl;
        cerr << __AFLOW_FUNC__ << " vpids.empty()=" << vpids.empty() << endl;
      }
      return vpids;
    }
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"\"pgrep\"-type command not found",_INPUT_ILLEGAL_);
    return vpids;
  }

  //SD20220329 - overload to allow for only getting the PIDs with a specific PGID
  vector<string> ProcessPIDs(const string& process,const string& pgid,string& output_syscall,bool user_specific){
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(pgid.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"PGID is empty",_INPUT_ILLEGAL_);
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " looking for pgid=" << pgid << endl;
      cerr << __AFLOW_FUNC__ << " looking for process=" << process << endl;
    }
    if(!aurostd::IsCommandAvailable("ps") || !aurostd::IsCommandAvailable("grep")) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"\"pgrep\"-type command not found",_INPUT_ILLEGAL_);
    }
    string ps_opts=" uid,pgid,pid,etime,pcpu,pmem,args"; // user-defined options, since just "u" or "j" might not be good enough
    string command="ps";
    vector<string> vlines,vtokens,vpids;
    uint i=0,j=0;
    aurostd::string2tokens(ps_opts,vtokens,",");
    uint nopts = vtokens.size();
    string command_grep="grep "+process;
    if(user_specific){command+=" xo";}
    else{command+=" axo";}
    command+=ps_opts;
    if(!aurostd::execute2string(command+" > /dev/null",stderr_fsio).empty()) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown options in \"ps\"",_INPUT_ILLEGAL_);
    }
    command+=" 2>/dev/null | "+command_grep+" 2> /dev/null";
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " running command=\"" << command << "\"" << endl;}
    string output=output_syscall=aurostd::execute2string(command);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " ps/grep output:" << endl << output << endl;}
    aurostd::string2vectorstring(output,vlines);
    for(i=0;i<vlines.size();i++){
      aurostd::string2tokens(vlines[i],vtokens," ");
      if(vtokens.size()<nopts){continue;} // set by ps_opts
      const string& pid=vtokens[2]; // set by ps_opts
      if (vtokens[1]==pgid) { // set by ps_opts
        string proc=vtokens[nopts-1]; // set by ps_opts
        for(j=nopts;j<vtokens.size();j++){proc+=" "+vtokens[j];} // set by ps_opts
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " proc[i=" << i << "]=\"" << proc << "\"" << endl;}
        if(proc.find(process)==string::npos){continue;}
        if(proc.find(command)!=string::npos){continue;}
        if(proc==command_grep){continue;}
        vpids.push_back(pid);
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " vpids=" << aurostd::joinWDelimiter(vpids,",") << endl;
      cerr << __AFLOW_FUNC__ << " vpids.empty()=" << vpids.empty() << endl;
    }
    return vpids;
  }

  // ***************************************************************************
  // Function ProcessRunning
  // ***************************************************************************
  //CO20210315
  bool ProcessRunning(const string& process,bool user_specific){return !aurostd::ProcessPIDs(process,user_specific).empty();} //CO20210315

  //SD20220329 - overload to allow for only getting the PIDs with a specific PGID
  bool ProcessRunning(const string& process,const string& pgid,bool user_specific){
    string output_syscall="";
    return !aurostd::ProcessPIDs(process,pgid,output_syscall,user_specific).empty();
  }

  // ***************************************************************************
  // Function ProcessKill
  // ***************************************************************************
  //CO20210315
  //SD20220627 - Typicall signals that we use to kill processes are 9 (SIGKILL), 15 (SIGTERM) and 5 (SIGTRAP)
  //Do not throw an error since only SIGKILL and SIGTERM actually kill the process
  void ProcessKill(const string& process,bool user_specific,uint signal){ //CO20210315
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(signal<1 || signal>64){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"invalid signal specification",_VALUE_ILLEGAL_);
    }
    vector<string> vpids=aurostd::ProcessPIDs(process,user_specific);
    if(vpids.empty()){return;}
    string command="kill";
    //[CO20210315 - does not work, user-specific comes from PID search]if(user_specific && !XHOST.user.empty()){command+=" -u "+XHOST.user;}
    command+=" -"+aurostd::utype2string(signal)+" "+aurostd::joinWDelimiter(vpids," ")+" 2>/dev/null";
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " running command=\"" << command << "\"" << endl;}
    aurostd::execute(command);
    //[SD20220627 - OBSOLETE]string command="";
    //[SD20220627 - OBSOLETE]bool process_killed=(!aurostd::ProcessRunning(process,user_specific));
    //[SD20220627 - OBSOLETE]uint sleep_seconds=5; //2 seconds is too few
    //[SD20220627 - OBSOLETE]if(!process_killed){
    //[SD20220627 - OBSOLETE]  if(aurostd::IsCommandAvailable("killall")) {
    //[SD20220627 - OBSOLETE]    command="killall";
    //[SD20220627 - OBSOLETE]    if(user_specific && !XHOST.user.empty()){command+=" -u "+XHOST.user;}
    //[SD20220627 - OBSOLETE]    command+=" -"+aurostd::utype2string(signal)+" "+process+" 2>/dev/null";
    //[SD20220627 - OBSOLETE]    if(LDEBUG){cerr << __AFLOW_FUNC__ << " running command=\"" << command << "\"" << endl;}
    //[SD20220627 - OBSOLETE]    aurostd::execute(command);
    //[SD20220627 - OBSOLETE]    aurostd::Sleep(sleep_seconds);process_killed=(!aurostd::ProcessRunning(process,user_specific));
    //[SD20220627 - OBSOLETE]  }
    //[SD20220627 - OBSOLETE]}
    //[SD20220627 - OBSOLETE]if(!process_killed){
    //[SD20220627 - OBSOLETE]  if(aurostd::IsCommandAvailable("pkill")) {
    //[SD20220627 - OBSOLETE]    command="pkill";
    //[SD20220627 - OBSOLETE]    if(user_specific && !XHOST.user.empty()){command+=" -u "+XHOST.user;}
    //[SD20220627 - OBSOLETE]    command+=" -"+aurostd::utype2string(signal)+" "+process+" 2>/dev/null";
    //[SD20220627 - OBSOLETE]    if(LDEBUG){cerr << __AFLOW_FUNC__ << " running command=\"" << command << "\"" << endl;}
    //[SD20220627 - OBSOLETE]    aurostd::execute(command);
    //[SD20220627 - OBSOLETE]    aurostd::Sleep(sleep_seconds);process_killed=(!aurostd::ProcessRunning(process,user_specific));
    //[SD20220627 - OBSOLETE]  }
    //[SD20220627 - OBSOLETE]}
  }

  //SD20220329 - overload to allow for only killing the PIDs with a specific PGID
  void ProcessKill(const string& process,const string& pgid,bool user_specific,uint signal){
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(signal<1 || signal>64){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"invalid signal specification",_VALUE_ILLEGAL_);
    }
    if(pgid.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"PGID is empty",_INPUT_ILLEGAL_);
    }
    string output_syscall="";
    vector<string> vpids=aurostd::ProcessPIDs(process,pgid,output_syscall,user_specific);
    if(vpids.empty()){return;}
    string command="kill";
    command+=" -"+aurostd::utype2string(signal)+" "+aurostd::joinWDelimiter(vpids," ")+" 2>/dev/null";
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " running command=\"" << command << "\"" << endl;}
    aurostd::execute(command);
  }

  // ***************************************************************************
  // Function ProcessRenice
  // ***************************************************************************
  //CO20210315
  bool ProcessRenice(const string& process,int nvalue,bool user_specific,const string& pgid){ //CO20210315 //CO20221029 - void->bool to catch if renice worked and added pgid to renice executable-specific processes
    vector<string> vpids;
    if(pgid.empty()){vpids=ProcessPIDs(process,user_specific);} //CO20221028
    else{string output_syscall="";vpids=ProcessPIDs(process,pgid,output_syscall,user_specific);}  //CO20221028
    if(vpids.empty()){return false;}
    string command="renice "+aurostd::utype2string(nvalue)+" "+aurostd::joinWDelimiter(vpids," ");
    string err=aurostd::execute2string(command,stderr_fsio);
    return err.empty();
  }
  
  // ***************************************************************************
  // Function ReniceAvailable
  // ***************************************************************************
  //CO20221029 - checks if renice up/down is allowed on the system
  bool ReniceAvailable(void){
    aurostd::execute("sleep 10 &"); //to background
    if(!ProcessRenice("sleep",19,true,aurostd::utype2string(getpgrp()))){return false;}  //try making more nice, this should almost always work
    if(!ProcessRenice("sleep",0,true,aurostd::utype2string(getpgrp()))){return false;}   //try making default nice, this may not always work: https://superuser.com/questions/88542/why-cant-unix-users-renice-downwards
    return true;
  }

  // ***************************************************************************
  // Function DirectoryMake
  // ***************************************************************************
  // Stefano Curtarolo
  // Make a directory by splitting each part.. returns false if something wrong
  // happens.
  bool DirectoryMake(string _Directory) {  // "" compliant SC20190401
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string Directory(CleanFileName(_Directory));
    std::deque<string> tokens;
    string dir_tmp,command;
    aurostd::string2tokens(Directory,tokens,"/");
    if(Directory.at(0)=='/') tokens[0]="/"+tokens[0]; // fix the root thing
    dir_tmp="";
    for(uint i=0;i<tokens.size();i++) {
      dir_tmp+=tokens[i]+"/";
      if(!aurostd::FileExist(dir_tmp)) {
        command="mkdir -p \""+dir_tmp+"\"";
        if(LDEBUG) cerr << "DirectoryMake creating directory=" <<  command << endl;
        aurostd::execute(command);
        if(!aurostd::FileExist(dir_tmp)) {
          // if(!QUIET) cout << "EEEEE   can not make directory: " <<  command << endl;
          return FALSE; // found some error in making directory
        }
      }
    }
    return TRUE;
  }

  // ***************************************************************************
  // Function SSH_DirectoryMake
  // ***************************************************************************
  // Stefano Curtarolo
  // Make a directory by splitting each part.. cant check much remotely
  // it starts from the assumption that you DO HAVE access to the remove account
  // which does not pretend a password
  bool SSH_DirectoryMake(string user, string machine,string _Directory) {  // "" compliant SC20190401
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string Directory(CleanFileName(_Directory));
    std::deque<string> tokens;
    string dir_tmp,command;
    aurostd::string2tokens(Directory,tokens,"/");
    if(Directory.at(0)=='/') tokens[0]="/"+tokens[0]; // fix the root thing
    dir_tmp="";
    for(uint i=0;i<tokens.size();i++) {
      dir_tmp+=tokens[i]+"/";
      command="ssh "+user+"@"+machine+" mkdir \""+dir_tmp+"\"";
      if(LDEBUG) cerr << "SSH_DirectoryMake creating directory=" <<  command << endl;
      aurostd::execute(command);
    }
    return TRUE;
  }

  // ***************************************************************************
  // Function DirectoryChmod
  // ***************************************************************************
  // Stefano Curtarolo
  // return FALSE if something got messed up
  bool DirectoryChmod(string chmod_string,string _Directory) {  // "" compliant SC20190401
    string Directory(CleanFileName(_Directory));
    ostringstream aus;
    aurostd::StringstreamClean(aus);
    aus << CHMOD_BIN << " " << chmod_string << " \"" << Directory << "\"" << endl;
    aurostd::execute(aus);
    return TRUE;
  }

  // ***************************************************************************
  // Function SubDirectoryLS
  // ***************************************************************************
  // Stefano Curtarolo
  // Returns all the subdirectories of a path
  bool SubDirectoryLS(const string& _Directory,vector<string>& vsubd){
    vector<string> vfiles;
    bool run=DirectoryLS(_Directory,vfiles);
    for(uint i=0;i<vfiles.size();i++){
      if(aurostd::IsDirectory(_Directory+"/"+vfiles[i])){
        vsubd.push_back(_Directory+"/"+vfiles[i]);
        run=(run&&SubDirectoryLS(_Directory+"/"+vfiles[i],vsubd));
      }
    }
    return run;
  }

  // ***************************************************************************
  // Function DirectoryLS
  // ***************************************************************************
  // Stefano Curtarolo
  // Returns the content of a directory without "." and ".."
  bool DirectoryLS(const string& _Directory,vector<string> &vfiles) {
    vfiles.clear();
    string Directory(CleanFileName(_Directory));
    DIR *dp;
    dp=opendir(Directory.c_str());
    if(dp!=NULL) {
      struct dirent *ep;
      string file;
      while((ep=readdir(dp))) {
        file=(ep->d_name);
        if(file!="." && file!="..")
          vfiles.push_back(file);
      }
      (void) closedir (dp);
    } else {
      cerr << "ERROR: aurostd::DirectoryLS Couldn't open the directory: " << Directory << endl;
      return FALSE;
    }
    return TRUE;
  }
  bool DirectoryLS(const string& _Directory,deque<string> &vfiles) {
    vector<string> _vfiles;
    bool run=DirectoryLS(_Directory,_vfiles);
    vfiles=aurostd::vector2deque(_vfiles);
    return run;
  }

  // ***************************************************************************
  // Function dirname
  // ***************************************************************************
  // CO20210315
  // Returns the dirname of file
  string dirname(const string& _file) {
    string file=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(CleanFileName(_file));
    if(file.empty()){return "";}
    //[no need]if(!aurostd::FileExist(file)){return "";}
    if(file[file.size()-1]=='/'){file=file.substr(0,file.size()-1);}  //remove last / for directories: matches functionality with bash dirname/basename
    if(file.find('/')==string::npos){return ".";}
    string::size_type loc=file.find_last_of('/');
    return file.substr(0,loc);
  }

  // ***************************************************************************
  // Function basename
  // ***************************************************************************
  // CO20210315
  // Returns the basename of file
  string basename(const string& _file) {  //CO20210315
    string file=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(CleanFileName(_file));
    if(file.empty()){return "";}
    //[no need]if(!aurostd::FileExist(file)){return "";}
    if(file[file.size()-1]=='/'){file=file.substr(0,file.size()-1);}  //remove last / for directories: matches functionality with bash dirname/basename
    if(file.find('/')==string::npos){return file;}
    string::size_type loc=file.find_last_of('/');
    return file.substr(loc+1);
  }

  // ***************************************************************************
  // Function DirectoryLocked
  // ***************************************************************************
  bool DirectoryLocked(string directory,string LOCK) {
    if(FileExist(directory+"/"+LOCK)) return TRUE;
    if(FileExist(directory+"/"+LOCK+_LOCK_LINK_SUFFIX_)) return TRUE;
    if(FileExist(directory+"/"+LOCK+".xz")) return TRUE;
    if(FileExist(directory+"/"+LOCK+".gz")) return TRUE;
    if(FileExist(directory+"/"+LOCK+".bz2")) return TRUE;
    return FALSE;
  }

  // ***************************************************************************
  // Function DirectorySkipped
  // ***************************************************************************
  bool DirectorySkipped(string directory) {
    if(FileExist(directory+"/SKIP")) return TRUE;
    if(FileExist(directory+"/SKIP.xz")) return TRUE;
    if(FileExist(directory+"/SKIP.gz")) return TRUE;
    if(FileExist(directory+"/SKIP.bz2")) return TRUE;
    return FALSE;
  }

  // ***************************************************************************
  // Function DirectoryWritable and DirectoryUnwritable
  // ***************************************************************************
  bool DirectoryWritable(string _Directory) {  // "" compliant SC20190401
    string Directory(CleanFileName(_Directory));
    string filename=aurostd::TmpFileCreate("DirectoryWritable",Directory,true); // SD20220223 - uses TmpFileCreate
    bool writable=aurostd::string2file("DirectoryWritable",filename);
    if(!writable || !FileExist(filename)) return FALSE;
    aurostd::RemoveFile(filename);
    return TRUE;
  }
  bool DirectoryUnwritable(string Directory) {
    return !DirectoryWritable(Directory);
  }

  // ***************************************************************************
  // Function CleanFileName
  // ***************************************************************************
  // Stefano Curtarolo
  // cleans file names from obvious things
  string CleanFileName(const string& fileIN) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(fileIN.empty()){return fileIN;}
    // ME20211001
    // Remove any control characters (below ASCII 32) while copying. This is useful
    // when fileIN is read from a file, which can have all sorts of junk and causes
    // FileExist to break. We cannot use RemoveControlCodeCharactersFromString()
    // because it keeps tabs and linebreaks and cannot use CleanStringASCII because
    // it keeps control characters.
    string fileOUT="";
    for (uint i = 0; i < fileIN.size(); i++) {
      if (fileIN[i] > 31) fileOUT += fileIN[i];
    }
    if(LDEBUG) cerr << "aurostd::CleanFileName: " << fileOUT << endl;
    // [OBSOLETE] interferes with ~/.aflow.rc   if(aurostd::substring2bool(fileOUT,"~/")) aurostd::StringSubst(fileOUT,"~/","/home/"+XHOST.user+"/");
    //ME20200922 - Cleaning // must be in a while loop or it won't clean e.g. ///
    while (fileOUT.find("//") != string::npos) aurostd::StringSubst(fileOUT,"//","/");
    aurostd::StringSubst(fileOUT,"/./","/");
    aurostd::StringSubst(fileOUT,"*","\"*\"");
    aurostd::StringSubst(fileOUT,"?","\"?\"");
    if(LDEBUG) cerr << "aurostd::CleanFileName: " << fileOUT << endl;
    return fileOUT;
  }

  // ***************************************************************************
  // Function ProperFileName
  // ***************************************************************************
  // Stefano Curtarolo
  // fix file names from obvious things
  string ProperFileName(const string& fileIN) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    // ME20211001
    // Remove any control characters (below ASCII 32) while copying. This is useful
    // when fileIN is read from a file, which can have all sorts of junk and causes
    // FileExist to break. We cannot use RemoveControlCodeCharactersFromString()
    // because it keeps tabs and linebreaks and cannot use CleanStringASCII because
    // it keeps control characters.
    string fileOUT="";
    for (uint i = 0; i < fileIN.size(); i++) {
      if (fileIN[i] > 31) fileOUT += fileIN[i];
    }
    if(LDEBUG) cerr << "aurostd::ProperFileName: " << fileOUT << endl;
    aurostd::StringSubst(fileOUT,"//",".");
    aurostd::StringSubst(fileOUT,"/",".");
    if(LDEBUG) cerr << "aurostd::ProperFileName: " << fileOUT << endl;
    return fileOUT;
  }

  // ***************************************************************************
  // Function CopyFile
  // ***************************************************************************
  // Stefano Curtarolo
  // copy the file but does not check if the directory can be made or not...
  // CO20200624 - it cannot check because we might pass in *
  bool CopyFile(const string& from,const string& to) { // "" compliant SC20190401
    stringstream command;
    command << "cp -f \"" << CleanFileName(from) << "\" \"" << CleanFileName(to) << "\" " << endl;
    aurostd::execute(command);
    return TRUE;
  }

  // ***************************************************************************
  // Function aurostd::LinkFile
  // ***************************************************************************
  // Stefano Curtarolo
  // copy the file but does not check if the directory can be made or not...
  // CO20200624 - it cannot check because we might pass in *
  bool LinkFile(const string& from,const string& to) { // "" compliant SC20190401
    stringstream command;
    command << "ln -sf \"" << CleanFileName(from) << "\" \"" << CleanFileName(to) << "\" " << endl;
    aurostd::execute(command);
    return TRUE;
  }

  // ***************************************************************************
  // Function aurostd::LinkFileAtomic
  // ***************************************************************************
  // Simon Divilov
  // Create a symbolic or hard link of a file using C++ functions
  bool LinkFileAtomic(const string& from,const string& to,bool soft) {
    int fail=0;
    string from_clean=CleanFileName(from),to_clean=CleanFileName(to);
    if(from_clean.empty() || to_clean.empty()) {return FALSE;}
    if(soft) {
      fail = symlink(from_clean.c_str(),to_clean.c_str());
    }
    else {
      fail = link(from_clean.c_str(),to_clean.c_str());
    }
    if(fail) {
      if(errno==EEXIST) {
        return FALSE;
      }
      else {
        string message = "Error linking "+from_clean+" -> "+to_clean+" | errno="+aurostd::utype2string<int>(errno);
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_FILE_ERROR_);
      }
    }
    else {
      return TRUE;
    }
  }

  // ***************************************************************************
  // Function aurostd::UnlinkFile
  // ***************************************************************************
  // Simon Divilov
  // Unlink file using C++ functions
  bool UnlinkFile(const string& link) {return (unlink(CleanFileName(link).c_str())==0);}

  //CO START
  //***************************************************************************//
  // aurostd::MatchCompressed
  //***************************************************************************//
  bool MatchCompressed(const string& CompressedFileName,const string& FileNameOUT) {
    //Corey Oses
    //Will try to mimic compression of a given file, useful for overwriting
    if(!IsCompressed(CompressedFileName)&&!IsCompressed(FileNameOUT)) {return TRUE;}
    if(IsCompressed(FileNameOUT)) {
      if(GetCompressionExtension(CompressedFileName)==GetCompressionExtension(FileNameOUT)) {return TRUE;}
      return FALSE;
    }
    if(substring2bool(CompressedFileName,".xz")) {CompressFile(FileNameOUT,"xz");return TRUE;}
    if(substring2bool(CompressedFileName,".bz2")) {CompressFile(FileNameOUT,"bzip2");return TRUE;}
    if(substring2bool(CompressedFileName,".gz")) {CompressFile(FileNameOUT,"gzip");return TRUE;}
    if(substring2bool(CompressedFileName,".zip")) {CompressFile(FileNameOUT,"zip");return TRUE;}
    return FALSE;
  }

  // [OBSOLETE] //***************************************************************************//
  // [OBSOLETE] // aurostd::DecompressFile
  // [OBSOLETE] //***************************************************************************//
  // [OBSOLETE] bool DecompressFile(const string& CompressedFileName){
  // [OBSOLETE]   //Corey Oses
  // [OBSOLETE]   //try to decompress the file the ways we know how
  // [OBSOLETE]  //try not to use unless necessary, prefer efile2stringstream, etc.
  // [OBSOLETE]   //where you might use it: open binary files
  // [OBSOLETE]   if(!IsCompressed(CompressedFileName)){return TRUE;}
  // [OBSOLETE]   if(substring2bool(CompressedFileName,".xz")) {XzipFile(CompressedFileName);return TRUE;}
  // [OBSOLETE]   if(substring2bool(CompressedFileName,".bz2")) {BzipFile(CompressedFileName);return TRUE;}
  // [OBSOLETE]   if(substring2bool(CompressedFileName,".gz")) {GzipFile(CompressedFileName);return TRUE;}
  // [OBSOLETE]   if(substring2bool(CompressedFileName,".zip")) {ZipFile(CompressedFileName);return TRUE;}
  // [OBSOLETE]   return FALSE;
  // [OBSOLETE] }

  //***************************************************************************//
  // aurostd::efile2tempfile
  //***************************************************************************//
  bool efile2tempfile(const string& _FileNameIN,string& FileNameOUT) {  //CO20210623
    bool tempfile_created=false;
    return efile2tempfile(_FileNameIN,FileNameOUT,tempfile_created);
  }
  bool efile2tempfile(const string& _FileNameIN,string& FileNameOUT,bool& tempfile_created) { //CO20210623
    //Corey Oses
    //Decompresses file to temp file, and return its path
    //SAFE: Will return FileNameIn if compression not needed
    string FileNameIN="";
    tempfile_created=false;
    if(!(aurostd::FileExist(_FileNameIN,FileNameIN) || aurostd::EFileExist(_FileNameIN,FileNameIN))) {  //CO20191110 - get FileNameIN from functions
      cerr << endl;
      cerr << "ERROR - aurostd::efile2tempfile: file=" << _FileNameIN << " not present !" << endl;  ///is empty
      cerr << endl;
      return FALSE;
    }
    //check if compressed, and decompress
    if(!aurostd::IsCompressed(FileNameIN)) {
      FileNameOUT = FileNameIN;
      return TRUE;
    }
    stringstream FileNameIN_ss;
    aurostd::efile2stringstream(FileNameIN, FileNameIN_ss);
    FileNameOUT = aurostd::TmpFileCreate();
    tempfile_created=true;
    aurostd::stringstream2file(FileNameIN_ss, FileNameOUT);
    return TRUE;
  }

  //***************************************************************************//
  // aurostd::IsCompressed
  //***************************************************************************//
  bool IsCompressed(const string& FileNameIN,string& FileNameOUT) {
    //Corey Oses
    //Given a File, it will return name of decompressed variant
    //NB: Does not actually decompress
    //FileNameOUT="";
    if(FileNameIN.find(".xz")!=string::npos) {FileNameOUT=FileNameIN;StringSubst(FileNameOUT,".xz","");return TRUE;}
    if(FileNameIN.find(".bz2")!=string::npos) {FileNameOUT=FileNameIN;StringSubst(FileNameOUT,".bz2","");return TRUE;}
    if(FileNameIN.find(".tar.gz")!=string::npos) {FileNameOUT=FileNameIN;StringSubst(FileNameOUT,".tar.gz","");return TRUE;}
    if(FileNameIN.find(".gz")!=string::npos) {FileNameOUT=FileNameIN;StringSubst(FileNameOUT,".gz","");return TRUE;}
    if(FileNameIN.find(".zip")!=string::npos) {FileNameOUT=FileNameIN;StringSubst(FileNameOUT,".zip","");return TRUE;}
    // FileNameOUT=FileNameIN; // dont touch it if not found
    return FALSE;
  }

  //***************************************************************************//
  // aurostd::IsCompressed
  //***************************************************************************//
  bool IsCompressed(const string& FileNameIN) {
    string FileNameOUT;
    return IsCompressed(FileNameIN,FileNameOUT);
  }

  //***************************************************************************//
  // aurostd::GetCompressionExtension
  //***************************************************************************//
  string GetCompressionExtension(const string& CompressedFileName) {
    //Corey Oses
    //Will determine zipped extension
    //[CO20210315 - more work]string extension="";
    //[CO20210315 - more work]if(!IsCompressed(CompressedFileName)) {return extension;}
    if(CompressedFileName.find(".xz")!=string::npos)  {return ".xz";}
    if(CompressedFileName.find(".bz2")!=string::npos) {return ".bz2";}
    if(CompressedFileName.find(".gz")!=string::npos)  {return ".gz";}
    if(CompressedFileName.find(".zip")!=string::npos) {return ".zip";}
    return "";
  }
  //CO END

  //***************************************************************************//
  // aurostd::GetCatCommand
  //***************************************************************************//
  string GetCatCommand(const string& CompressedFileName) {  //CO20210315
    //Corey Oses
    //Will determine zipped extension
    if(CompressedFileName.find(".xz")!=string::npos)  {return "xzcat";}
    if(CompressedFileName.find(".bz2")!=string::npos) {return "bzcat";}
    if(CompressedFileName.find(".gz")!=string::npos)  {return "zcat";}
    if(CompressedFileName.find(".zip")!=string::npos) {return "unzip -p";}
    return "cat";
  }

  //CO END
  // ***************************************************************************
  // Function aurostd::UncompressFile aurostd::Compress
  // ***************************************************************************
  // Bzip the file (does not check)
  bool UncompressFile(const string& _file,const string& command) {  // "" compliant SC20190401
    string file(CleanFileName(_file));
    if(command=="bunzip2" || command=="bzip2" || command=="bz2"  || command==".bz2") {
      if(!aurostd::IsCommandAvailable("bzip2")) {
        cerr << "ERROR - aurostd::UncompressFile: command \"bzip2\" is necessary !" << endl;
        return FALSE; }   
      if(file.find(".bz2")!=string::npos) aurostd::execute("bzip2 -dqf \""+file+"\"");
    }
    if(command=="xunzip" || command=="xz" || command==".xz") {
      if(!aurostd::IsCommandAvailable("xz")) {
        cerr << "ERROR - aurostd::UncompressFile: command \"xz\" is necessary !" << endl;
        return FALSE; }   
      if(file.find(".xz")!=string::npos) aurostd::execute("xz -dqf \""+file+"\"");
    }
    if(command=="gunzip" || command=="gzip" || command=="gz"  || command==".gz") {
      if(!aurostd::IsCommandAvailable("gzip")) {
        cerr << "ERROR - aurostd::UncompressFile: command \"gzip\" is necessary !" << endl;
        return FALSE; }   
      if(file.find(".gz")!=string::npos) aurostd::execute("gzip -dqf \""+file+"\"");
    }
    if(command=="unzip" || command=="zip" || command==".zip") {
      if(!aurostd::IsCommandAvailable("unzip")) {
        cerr << "ERROR - aurostd::UnzipFile: command \"unzip\" is necessary !" << endl;
        return FALSE;
      }
      if(file.find(".zip")!=string::npos) aurostd::execute("unzip -qo \""+file+"\"");
    }
    return TRUE;
  }
  bool UncompressFile(const string& _file) {  // "" compliant SC20190401
    string file(CleanFileName(_file));
    if(file.find(".xz")!=string::npos)  {return UncompressFile(file,"xz");}
    if(file.find(".bz2")!=string::npos) {return UncompressFile(file,"bzip2");}
    if(file.find(".gz")!=string::npos)  {return UncompressFile(file,"gzip");}
    if(file.find(".zip")!=string::npos) {return UncompressFile(file,"zip");}
    return FALSE;
  }

  bool CompressFile(const string& _file,const string& command) {  // "" compliant SC20190401
    string file(CleanFileName(_file));
    //  cerr << "aurostd::CompressFile FileName=[" << FileName << "]  command=[" << command << "]" << endl;
    if(aurostd::substring2bool(command,"bzip2") || aurostd::substring2bool(command,"bz2")  || aurostd::substring2bool(command,".bz2")) {
      if(!aurostd::IsCommandAvailable("bzip2")) {
        cerr << "ERROR - aurostd::CompressFile: command \"bzip2\" is necessary !" << endl;
        return FALSE;
      }   
      // [OBSOLETE]     if(FileExist(file+".bz2")) {aurostd::execute("rm -f \""+file+".bz2\"");}
      if(file.find(".bz2")==string::npos) aurostd::execute("bzip2 -9qf \""+file+"\"");
      return TRUE;
    }
    if(aurostd::substring2bool(command,"xz") || aurostd::substring2bool(command,"xzip") || aurostd::substring2bool(command,".xz")) {
      if(!aurostd::IsCommandAvailable("xz")) {
        cerr << "ERROR - aurostd::CompressFile: command \"xz\" is necessary !" << endl;
        return FALSE;
      }   
      // [OBSOLETE]     if(FileExist(file+".xz")) {aurostd::execute("rm -f \""+file+".xz\"");}
      //    cerr << "aurostd::CompressFile XZ  FileName=[" << FileName << "]  command=[" << command << "]" << endl;
      if(file.find(".xz")==string::npos) aurostd::execute("xz -9qf -q \""+file+"\""); // twice -q to avoid any verbosity
      return TRUE;
    }
    if(aurostd::substring2bool(command,"gzip") || aurostd::substring2bool(command,"gz") || aurostd::substring2bool(command,".gz")) {
      if(!aurostd::IsCommandAvailable("gzip")) {
        cerr << "ERROR - aurostd::CompressFile: command \"gzip\" is necessary !" << endl;
        return FALSE;
      }   
      // [OBSOLETE]     if(FileExist(file+".gz")) {aurostd::execute("rm -f \""+file+".gz\"");}
      if(file.find(".gz")==string::npos) aurostd::execute("gzip -9qf \""+file+"\"");
      return TRUE;
    }
    if(aurostd::substring2bool(command,"zip") || aurostd::substring2bool(command,".zip")) {
      if(!aurostd::IsCommandAvailable("zip")) {
        cerr << "ERROR - aurostd::ZipFile: command \"zip\" is necessary !" << endl;
        return FALSE;
      }
      // [OBSOLETE]     if(FileExist(file+".zip")) {cerr << file << ".zip" << endl;}
      if(FileExist(file+".zip")) {aurostd::execute("rm -f \""+file+".zip\"");}
      if(file.find(".zip")==string::npos) aurostd::execute("zip -9qm \""+file+".zip\" \""+file+"\"");
      return TRUE;
    }
    return FALSE;
  }

  // ***************************************************************************
  // aurostd::ZIP2ZIP aurostd::BZ2XZ aurostd::GZ2XZ
  // ***************************************************************************
  bool ZIP2ZIP(string _dir,string from,string to,bool VERBOSE,const string& message) {  // "" compliant SC20190401
    string from_cmd="bzip2",from_ext="bz2";
    string to_cmd="xz",to_ext="xz";
    string dir=aurostd::CleanFileName(_dir);

    if((from=="bz" || from=="bz2" || from=="bzip2") && (to=="xz")) { from_cmd="bzip2",from_ext="bz2";to_cmd="xz",to_ext="xz"; } 
    if((from=="xz") && (to=="bz" || to=="bz2" || to=="bzip2")) { from_cmd="xz",from_ext="xz";to_cmd="bzip2",to_ext="bz2"; } 
    if((from=="bz" || from=="bz2" || from=="bzip2") && (to=="gz" || to=="gzip")) { from_cmd="bzip2",from_ext="bz2";to_cmd="gzip",to_ext="gz"; } 
    if((from=="gz" || from=="gzip") && (to=="bz" || to=="bz2" || to=="bzip2")) { from_cmd="gzip",from_ext="gz";to_cmd="bzip2",to_ext="bz2"; }  
    if((from=="gz" || from=="gzip") && (to=="xz")) { from_cmd="gzip",from_ext="gz";to_cmd="xz",to_ext="xz"; } 
    if((from=="xz") && (to=="gz" || to=="gzip")) { from_cmd="xz",from_ext="xz";to_cmd="gzip",to_ext="gz"; } 
    if((from=="tbz") && (to=="xz")) { from_cmd="bzip",from_ext="tbz";to_cmd="xz",to_ext="xz"; } 
    if((from=="tgz") && (to=="xz")) { from_cmd="gzip",from_ext="tgz";to_cmd="xz",to_ext="xz"; } 

    if(VERBOSE) { cout << message << "aurostd::ZIP2ZIP: BEGIN - dir=" << dir << endl; }
    vector<string> vfile;
    //    cerr << string("ls \""+dir+"\"/* | grep "+from_ext) << endl;
    aurostd::string2vectorstring(aurostd::execute2string("ls \""+dir+"\"/* | grep "+from_ext),vfile);

    for(uint ifile=0;ifile<vfile.size();ifile++) {
      if(VERBOSE) { cout << message << "aurostd::ZIP2ZIP: vfile.at(ifile)=" << vfile.at(ifile) << endl;}
      aurostd::StringSubst(vfile.at(ifile),"."+from_ext,"");
      if(aurostd::FileExist(vfile.at(ifile)+"."+from_ext)) {
        // PATCH to be removed	if(VERBOSE)
        { cout << message << "aurostd::ZIP2ZIP: " << from_ext << "->"+to_ext+" vfile.at(" << ifile << ")=" << vfile.at(ifile) << " "; cout.flush(); }
        aurostd::UncompressFile(vfile.at(ifile)+"."+from_ext);
        // PATCH to be removedif(VERBOSE)
        { cout << "[" << from_ext << "]"; cout.flush(); }
        aurostd::CompressFile(vfile.at(ifile),to_ext);
        // PATCH to be removedif(VERBOSE)
        { cout << "["+to_ext+"]"; cout.flush(); }
        // PATCH to be removed if(VERBOSE)
        { cout << endl; cout.flush(); }
      }
    }
    if(aurostd::FileExist(dir+"/aflow.in"))
      if(aurostd::substring_present_file_FAST(dir+"/aflow.in",from_cmd)) {
        if(VERBOSE) { cout << message << "aurostd::ZIP2ZIP: " << from_ext << "->"+to_ext+" " << dir << "/aflow.in" << " " << endl; cout.flush(); }
        aurostd::execute("subst "+from_cmd+" "+to_cmd+" \""+dir+"/aflow.in\"");
        aurostd::RemoveFile("\""+dir+"/\"*~");
      }
    if(aurostd::FileExist(dir+"/agl_aflow.in"))
      if(aurostd::substring_present_file_FAST(dir+"/agl_aflow.in",from_cmd)) {
        if(VERBOSE) { cout << message << "aurostd::ZIP2ZIP: " << from_ext << "->"+to_ext+" " << dir << "/agl_aflow.in" << " " << endl; cout.flush(); }
        aurostd::execute("subst "+from_cmd+" "+to_cmd+" \""+dir+"/agl_aflow.in\"");
        aurostd::RemoveFile("\""+dir+"/\"*~");
      }
    if(aurostd::FileExist(dir+"/ael_aflow.in"))
      if(aurostd::substring_present_file_FAST(dir+"/ael_aflow.in",from_cmd)) {
        if(VERBOSE) { cout << message << "aurostd::ZIP2ZIP: " << from_ext << "->"+to_ext+" " << dir << "/ael_aflow.in" << " " << endl; cout.flush(); }
        aurostd::execute("subst "+from_cmd+" "+to_cmd+" \""+dir+"/ael_aflow.in\"");
        aurostd::RemoveFile("\""+dir+"/\"*~");
      }
    //AS20201023 BEGIN
    if(aurostd::FileExist(dir+"/aflow_qha.in"))
      if(aurostd::substring_present_file_FAST(dir+"/aflow_qha.in",from_cmd)) {
        if(VERBOSE) { cout << message << "aurostd::ZIP2ZIP: " << from_ext << "->"+to_ext+" " << dir << "/aflow_qha.in" << " " << endl; cout.flush(); }
        aurostd::execute("subst "+from_cmd+" "+to_cmd+" \""+dir+"/aflow_qha.in\"");
        aurostd::RemoveFile("\""+dir+"/\"*~");
      }
    //AS20201023 END
    if(aurostd::FileExist(dir+"/LOCK"))
      if(aurostd::substring_present_file_FAST(dir+"/LOCK",from_cmd)) {
        if(VERBOSE) { cout << message << "aurostd::ZIP2ZIP: " << from_ext << "->"+to_ext+" " << dir << "/LOCK" << " " << endl; cout.flush(); }
        aurostd::execute("subst "+from_cmd+" "+to_cmd+" \""+dir+"/LOCK\"*");
        aurostd::RemoveFile("\""+dir+"/\"*~");
      }
    if(aurostd::FileExist(dir+"/LLOCK"))
      if(aurostd::substring_present_file_FAST(dir+"/LLOCK",from_cmd)) {
        if(VERBOSE) { cout << message << "aurostd::ZIP2ZIP: " << from_ext << "->"+to_ext+" " << dir << "/LLOCK" << " " << endl; cout.flush(); }
        aurostd::execute("subst "+from_cmd+" "+to_cmd+" \""+dir+"\"/LLOCK\"*");
        aurostd::RemoveFile("\""+dir+"/\"*~");
      }

    if(VERBOSE) { cout << message << "aurostd::ZIP2ZIP: END   - dir=" << dir << endl; }
    return TRUE;
  }

  bool BZ2XZ(string dir,bool VERBOSE,const string& message) { return ZIP2ZIP(dir,"bz2","xz",VERBOSE,message); }
  bool GZ2XZ(string dir,bool VERBOSE,const string& message) { return ZIP2ZIP(dir,"gz","xz",VERBOSE,message); }

  // ***************************************************************************
  // Function FileExist
  // ***************************************************************************
  // Stefano Curtarolo
  // return a simple bool, nothing else, from a string (which is supposed to be
  // a file name)
  bool FileExist(const string& FileName) {
    if(FileName.empty()){return false;} //CO20210623
    string file(CleanFileName(FileName));
    //   cerr << file << endl;
    bool exist=FALSE;
    ifstream FileStream;
    FileStream.open(file.c_str(),std::ios::in);
    FileStream.clear();
    FileStream.close();
    if(FileStream.good()) {exist=TRUE;} else {exist=FALSE;}
    return exist;
  }

  // ***************************************************************************
  // Function FileExist
  // ***************************************************************************
  // Stefano Curtarolo
  bool FileExist(const string& FileName, string& FileNameOut) {
    if(FileExist(FileName)) {FileNameOut=FileName;return TRUE;}
    // FileNameOut=FileName;  // dont touch it if not found
    return FALSE;
  }

  // ***************************************************************************
  // Function EFileExist
  // ***************************************************************************
  // Stefano Curtarolo
  bool EFileExist(const string& FileName) {
    string FileNameOut;
    return EFileExist(FileName,FileNameOut);
  }

  // ***************************************************************************
  // Function EFileExist
  // ***************************************************************************
  // Corey Oses
  // tells you if the compressed variant exists
  bool EFileExist(const string& _FileName, string& FileNameOut){
    string FileName=aurostd::CleanFileName(_FileName); //CO20191110
    if(FileExist(FileName)) {FileNameOut=FileName;return TRUE;}
    if(FileExist(FileName+".xz")) {FileNameOut=FileName+".xz";return TRUE;}
    if(FileExist(FileName+".bz2")) {FileNameOut=FileName+".bz2";return TRUE;}
    if(FileExist(FileName+".gz")) {FileNameOut=FileName+".gz";return TRUE;}
    if(FileExist(FileName+".zip"))  {FileNameOut=FileName+".zip";return TRUE;}
    // FileNameOut=FileName;  // dont touch it if not found
    return FALSE;
  }

  // ***************************************************************************
  // Function FileSize
  // ***************************************************************************
  // Stefano Curtarolo - jan 08
  // returns in bytes the size of a file
  //ME20191001 - Changed to unsigned long long int to accommodate large files
  unsigned long long int FileSize(const string& _FileName) {
    string FileName(CleanFileName(_FileName));
    if(0){  //CO20210601 - found much faster approach below, we don't want to read the whole file if it's big
      ifstream FileStream;
      FileStream.open(FileName.c_str(),std::ios::in);
      if(!FileStream.good()){return 0;}
      unsigned long long int sizeout = 0;
      //[CO20210315 - no need to store BIG file in memory]string FileString="";
      char c; 
      while (FileStream.get(c)){
        //[CO20210315 - no need to store BIG file in memory]FileString+=c;
        sizeout+=1; //[CO20210315 - no need to store BIG file in memory]=FileString.length();
      }
      FileStream.close();
      return sizeout;
    }
    //CO20210315 - much faster approach
    //https://www.codespeedy.com/cpp-program-to-get-the-size-of-a-file/
    FILE* fp=fopen(FileName.c_str(),"r");
    if (fp==NULL){return 0;}
    fseek(fp,0L,SEEK_END);
    unsigned long long int sizeout=ftell(fp);
    fclose(fp);
    return sizeout;
  }

  bool GetMemoryUsagePercentage(double& usage_percentage_ram,double& usage_percentage_swap){ //CO20210601
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    unsigned long long int free_ram=0,total_ram=0;
    unsigned long long int free_swap=0,total_swap=0;
    usage_percentage_ram=0.0;usage_percentage_swap=0.0;
    bool memory_read=aurostd::GetMemory(free_ram,total_ram,free_swap,total_swap);
    if(memory_read){
      if(total_ram>0){usage_percentage_ram=100.0*(((double)(total_ram-free_ram))/((double)(total_ram)));}
      if(total_swap>0){usage_percentage_swap=100.0*(((double)(total_swap-free_swap))/((double)(total_swap)));}  //some qrats nodes have no swap
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " [date=" << aflow_get_time_string() << "]" << endl; //helps debugging
        cerr << __AFLOW_FUNC__ << " free_ram=" << free_ram << endl;
        cerr << __AFLOW_FUNC__ << " used_ram=" << total_ram-free_ram << endl;
        cerr << __AFLOW_FUNC__ << " total_ram=" << total_ram << endl;
        cerr << __AFLOW_FUNC__ << " usage_percentage_ram=" << usage_percentage_ram << endl;
        cerr << __AFLOW_FUNC__ << " free_swap=" << free_swap << endl;
        cerr << __AFLOW_FUNC__ << " used_swap=" << total_swap-free_swap << endl;
        cerr << __AFLOW_FUNC__ << " total_swap=" << total_swap << endl;
        cerr << __AFLOW_FUNC__ << " usage_percentage_swap=" << usage_percentage_swap << endl;
        cerr << endl; //helps debugging
      }
    }
    else{
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " unable to query memory on the node" << endl;}
    }
    return memory_read;
  }

  bool GetMemory(unsigned long long int& free_ram,unsigned long long int& total_ram,unsigned long long int& free_swap,unsigned long long int& total_swap){ //CO20210315 - only works for linux: needs `free` command
    //https://www.howtogeek.com/456943/how-to-use-the-free-command-on-linux/
    //will grab the total and the free
    //the free is the memory unused by anything
    //used column includes buff/cache, some of which the kernel can sacrifice for other applications if necessary
    //available column is an "estimate" of what could become available if needed
    //it's best to make decisions based on the free column
    //https://unix.stackexchange.com/questions/14102/real-memory-usage
    //free will follow real memory (physical RAM), using the available column will follow the actual memory (what could become available if necessary)
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    if(!aurostd::IsCommandAvailable("free")){return false;}
    string output=aurostd::execute2string("free");
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " free output:" << endl << output << endl;}
    //on most linux machines:
    //                            total        used        free      shared  buff/cache   available
    //              Mem:      395654628    41363940    30948848     4106252   323341840   349143640
    //              Swap:       2097148           0     2097148
    //on qrats:
    //                          total       used       free     shared    buffers     cached
    //             Mem:     264523076  255134588    9388488         36    1745836  232705576
    //             -/+ buffers/cache:   20683176  243839900
    //             Swap:      4194300      71436    4122864
    vector<string> vlines;
    aurostd::string2vectorstring(output,vlines);
    if(vlines.size()<3){return false;}
    vector<string> vtokens;
    uint iline=0;
    //ram
    iline=1;
    if(vlines[iline].find("Mem:")==string::npos){return false;}
    aurostd::string2tokens(vlines[iline],vtokens," ");
    if(!aurostd::isfloat(vtokens[1])){return false;}
    total_ram=aurostd::string2utype<unsigned long long int>(vtokens[1]);
    if(!aurostd::isfloat(vtokens[3])){return false;}
    free_ram=aurostd::string2utype<unsigned long long int>(vtokens[3]);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " free_ram=" << free_ram << " total_ram=" << total_ram << endl;}
    //swap
    iline=2;
    if(vlines[iline].find("Swap:")==string::npos){iline++;}  //try next line
    if(vlines[iline].find("Swap:")==string::npos){return false;}
    aurostd::string2tokens(vlines[iline],vtokens," ");
    if(!aurostd::isfloat(vtokens[1])){return false;}
    total_swap=aurostd::string2utype<unsigned long long int>(vtokens[1]);
    if(!aurostd::isfloat(vtokens[3])){return false;}
    free_swap=aurostd::string2utype<unsigned long long int>(vtokens[3]);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " free_swap=" << free_swap << " total_swap=" << total_swap << endl;}
    //
    return true;
  }

  // ***************************************************************************
  // Function FileEmpty && FileNotEmpty
  // ***************************************************************************
  // Stefano Curtarolo - jan 08
  // returns in bytes the size of a file
  bool FileEmpty(const string& _FileName) {
    string FileName(CleanFileName(_FileName));
    if(FileExist(FileName)==FALSE) return TRUE;  // does not exist hence empty
    // it exists
    if(1) {
      int i=0;
      ifstream FileStream;
      FileStream.open(FileName.c_str(),std::ios::in);
      char c; 
      while (FileStream.get(c)&&i<256) {i++;};
      // count no more that 16... it is not worth to count more
      FileStream.close();
      if(i>0) return FALSE;
      else return TRUE;
    }
    if(0) {
      if(FileSize(FileName)<=1) return TRUE;
      else return FALSE;
    }
    return FALSE;
  }
  bool FileNotEmpty(const string& FileName) {
    return !FileEmpty(FileName);
  }
  bool EFileEmpty(const string& FileName) {  //CO20190808
    string decompressed_file=""; //CO20190808
    efile2tempfile(FileName,decompressed_file);  //CO20190808
    bool fileempty=FileEmpty(decompressed_file);
    if(FileName!=decompressed_file){RemoveFile(decompressed_file);} //CO20200624 - remove tmp file IFF it is a tmp file
    return fileempty;
  }
  bool EFileNotEmpty(const string& FileName) {  //CO20190808
    return !EFileEmpty(FileName);
  }

  // ***************************************************************************
  // Function GetTimestampModified
  // ***************************************************************************
  // ME20180712
  // gets modification time and returns SECONDS since epoch (as long int)
  long int GetTimestampModified(const string& _FileName) {
    string FileName(CleanFileName(_FileName));
    if(!FileExist(FileName)){return 0;}
    time_t tm = 0;
    struct stat file_stat;
    if (stat(FileName.c_str(), &file_stat) == 0) tm = file_stat.st_mtime;
    return static_cast<long int>(tm);
  }

  // ***************************************************************************
  // Function SecondsSinceFileModified
  // ***************************************************************************
  // CO20210315
  // gets modification time and returns SECONDS since now (as long int)
  long int SecondsSinceFileModified(const string& _FileName) {
    string FileName(CleanFileName(_FileName));
    if(!FileExist(FileName)){return 0;}
    long int tmod_file=GetTimestampModified(FileName);
    if(0){  //CO20210315 - this does NOT work, current time on the machine (node) vs. NFS will cause problems
      time_t t = std::time(NULL); //DX20200319 - nullptr -> NULL
      long int tmod_curr = (long int) t;
      return max((long int)0,tmod_curr-tmod_file);  //max ensures that if something goes wrong, we return 0
    }
    //instead, write a new file and take the difference in the time stamps
    //solution inspired by ME - create a temporary file IN THE CURRENT DIRECTORY (within NFS) and check time deltas
    string dir=aurostd::dirname(FileName);
    string tmpfile=aurostd::TmpFileCreate("timestamp",dir,true); //put in current directory, make it hidden
    if(!aurostd::string2file("timestamp",tmpfile)){return 0;}  //write the file out, need a better fix here, perhaps write to /tmp?
    long int tmod_tmp=aurostd::GetTimestampModified(tmpfile);
    aurostd::RemoveFile(tmpfile);
    return max((long int)0,tmod_tmp-tmod_file); //max ensures that if something goes wrong, we return 0
  }

  // ***************************************************************************
  // Function getFileChecmSum
  // ***************************************************************************
  //ME20190219
  // Generates the checksum of a file
  // Taken from old APL/apl_hroutines.
  unsigned int getFileCheckSum(const string& filename, const string& algo) {
    ifstream infile(filename.c_str(), std::ios::in | std::ios::binary);
    if (!infile.is_open()) {
      string message = "Cannot open file " + filename + ".";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _FILE_ERROR_);
    }

    // Get file length
    infile.seekg(0, std::ios::end);
    unsigned long length = infile.tellg();
    infile.seekg(0, std::ios::beg);

    // Setup read buffer (for whole file)
    if (length % 2 != 0)
      length++;
    char* buffer = new char[length];
    buffer[length - 1] = 0x00;

    // Read it in!
    infile.read(buffer, length);
    infile.close();

    // Get checksum
    unsigned int checksum;
    if (algo == "Fletcher32") {
      checksum = getFletcher32((unsigned short*)buffer, length >> 1);
    } else {
      checksum = 0;
    }
    delete[] buffer;

    // Return value
    return checksum;
  }

  // ***************************************************************************
  // Function getFletcher32
  // ***************************************************************************
  //ME20190219
  // Generates the 32 bit checksum of a string based on Fletcher's algorithm.
  // See http://en.wikipedia.org/wiki/Fletcher%27s_checksum
  // Taken from old APL/apl_hroutines.
  unsigned int getFletcher32(unsigned short* data, size_t len) {
    unsigned int sum1 = 0xffff, sum2 = 0xffff;

    while (len) {
      unsigned tlen = len > 360 ? 360 : len;
      len -= tlen;
      do {
        sum1 += *data++;
        sum2 += sum1;
      } while (--tlen);
      sum1 = (sum1 & 0xffff) + (sum1 >> 16);
      sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    }

    // Second reduction step to reduce sums to 16 bits
    sum1 = (sum1 & 0xffff) + (sum1 >> 16);
    sum2 = (sum2 & 0xffff) + (sum2 >> 16);

    return sum2 << 16 | sum1;
  }

  // ***************************************************************************
  // Function FileToString
  // ***************************************************************************
  // Loat the content of a file into a string
  string FileToString(const string& _FileName) {
    string FileName(CleanFileName(_FileName));
    ifstream FileStream;
    stringstream strstreamout; aurostd::StringstreamClean(strstreamout);
    string strline="";
    //  cerr << FileName.c_str() << endl;  // DEBUG
    FileStream.open(FileName.c_str(),std::ios::in);
    if(FileStream.good()) // !=NULL)
    { //CO20200106 - patching for auto-indenting
      while(getline(FileStream,strline)) {
        strstreamout << strline << endl;
      }
    }
    FileStream.clear();FileStream.close();
    return strstreamout.str();
  }

  // ***************************************************************************
  // Function InFileExistCheck
  // ***************************************************************************
  // Dane Morgan
  void InFileExistCheck(const string& routine, const string& _FileName,
      ifstream& file_to_check) {
    string FileName(CleanFileName(_FileName));
    if(!file_to_check) {
      string message = "In routine " + routine + ". Cannot open file " + FileName + ".";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

  // ***************************************************************************
  // Function IsCommandAvailable
  // ***************************************************************************
  // tells you if the command is available
  bool IsCommandAvailable(const string& command, string& position) {
    // position=aurostd::execute2string("which "+command+" 2>&1 2> /dev/null");
    position=aurostd::execute2string("bash -c \"which "+command+" 2> /dev/null\"");  //CO20210315 - put stderr to /dev/null //2>&1
    // cerr << position.length() << endl;
    aurostd::StringSubst(position,"\n","");
    position=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(position); //CO20210315 - remove white spaces
    if(position.length()>0) return TRUE;
    if(aurostd::FileExist("./"+command)) {position="./"+command;return TRUE;}
    if(aurostd::FileExist("/bin/"+command)) {position="/bin/"+command;return TRUE;}
    if(aurostd::FileExist("/sbin/"+command)) {position="/sbin/"+command;return TRUE;}  // go around path
    if(aurostd::FileExist("/usr/bin/"+command)) {position="/usr/bin/"+command;return TRUE;}  // go around path
    if(aurostd::FileExist("/usr/sbin/"+command)) {position="/usr/sbin/"+command;return TRUE;}  // go around path
    if(aurostd::FileExist("/usr/local/bin/"+command)) {position="/usr/local/bin/"+command;return TRUE;}  // go around path
    if(aurostd::FileExist("/usr/local/sbin/"+command)) {position="/usr/local/sbin/"+command;return TRUE;}  // go around path
    if(aurostd::FileExist("/usr/local/maui/bin/"+command)) {position="/usr/local/maui/bin/"+command;return TRUE;}  // go around path  //CO20200526
    position="";
    return FALSE;
  }

  bool IsCommandAvailable(const string& command) {
    string position;
    return aurostd::IsCommandAvailable(command,position);
  }

  //CO20180706 - fixed this function, previously command/position trampled all over each other
  bool IsCommandAvailableModify(string& command) {
    string position;
    if(!aurostd::IsCommandAvailable(command,position)) return FALSE;
    command=position;
    return true;
  }

  // ***************************************************************************
  // Function CommandRequired
  // ***************************************************************************
  // tells you if the command is available
  bool CommandRequired(const string& command, string& position) {
    position=aurostd::execute2string("which "+command);
    aurostd::StringSubst(position,"\n","");
    if(position.length()>0) return TRUE;
    string message = "\"" + command + "\" is not available";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    return FALSE;
  }

  bool CommandRequired(const string& command) {
    string position;
    return CommandRequired(command,position);
  }

  // ***************************************************************************
  // Function IsExecutableAvailable
  // ***************************************************************************
  // tells you if the executable is available
  bool IsExecutableAvailable(const string& executable, string& position) {
    return IsCommandAvailable(executable,position);
  }

  bool IsExecutableAvailable(const string& executable) {
    string position;
    return IsCommandAvailable(executable,position);
  }

  // ***************************************************************************
  // Function ExecutableRequired
  // ***************************************************************************
  // tells you if the executable is available
  bool ExecutableRequired(const string& executable, string& position) {
    return CommandRequired(executable,position);
  }

  bool ExecutableRequired(const string& executable) {
    string position;
    return CommandRequired(executable,position);
  }

  // ***************************************************************************
  // DeleteOstringStreams
  // ***************************************************************************
  void StringstreamClean(ostringstream &aus) {
    aus.str(std::string()); //CO20200624
    aus.clear();  //CO20200624
    //[CO20200624 - fills stream with binary junk]aus.seekp(0,ios_base::beg);       // RESET
    //[CO20200624 - fills stream with binary junk]for(int i=0;i<BUFFER_MAXLEN;i++)  // RESET
    //[CO20200624 - fills stream with binary junk]  aus<<(char)0;                   // RESET
    //[CO20200624 - fills stream with binary junk]aus.seekp(0,ios_base::beg);       // RESET
    //[CO20200624 - fills stream with binary junk]// aus.str(std::string());
  }
  void StringstreamClean(stringstream &aus) {
    aus.str(std::string()); //CO20200624
    aus.clear();  //CO20200624
    //[CO20200624 - fills stream with binary junk]aus.seekp(0,ios_base::beg);       // RESET
    //[CO20200624 - fills stream with binary junk]for(int i=0;i<BUFFER_MAXLEN;i++)  // RESET
    //[CO20200624 - fills stream with binary junk]  aus<<(char)0;                   // RESET
    //[CO20200624 - fills stream with binary junk]aus.seekp(0,ios_base::beg);       // RESET
    //[CO20200624 - fills stream with binary junk]// aus.str(std::string());
  }


  // ***************************************************************************
  // Function FindIfStringInStream
  // ***************************************************************************
  //  This function returns true if string is in stream
  //  (on one line), or otherwise false. The search starts
  //  at the present file pointer location.  Note that this
  //  does alter the input string, resetting the file pointer
  //  to the input value at the end.
  // Dane Morgan style

  int FindIfStringInStream(const string& key, std::istream& instream) {
    // Get file pointer location at entry
    int loc=instream.tellg();
    int found_match=0;
    string s;
    getline(instream,s);
    int cont=0;
    if(getline(instream,s)) cont=1;
    while (cont) {
      int id=s.find(key);
      if(id!=(int) s.npos) { // Found key
        cont=0;
        found_match=1;
      }
      if(!getline(instream,s)) cont=0;
    }
    // Clear any fail bits associated with searching.
    instream.clear();
    // Set file pointer location to entry value
    instream.seekg(loc);
    return found_match;
  }

  // ***************************************************************************
  // Print Messages Errors and Warnings on and off streams.
  // ***************************************************************************
#define ErrorBarString   "EEEEE  ---------------------------------------------------------------------------------------------------------------------------- "
#define WarningBarString "WWWWW  ---------------------------------------------------------------------------------------------------------------------------- "

  //[CO20200624 - OBSOLETE]// with ostringstream
  //[CO20200624 - OBSOLETE]void PrintMessageStream(ofstream &FileMESSAGE,ostringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  FileMESSAGE << stream.str().c_str(); FileMESSAGE.flush();
  //[CO20200624 - OBSOLETE]  if(!quiet) {cout << stream.str().c_str();cout.flush();}
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  //[CO20200624 - OBSOLETE]void PrintMessageStream(std::ostream &FileMESSAGE,ostringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  FileMESSAGE << stream.str().c_str(); FileMESSAGE.flush();
  //[CO20200624 - OBSOLETE]  if(!quiet) {cout << stream.str().c_str();cout.flush();}
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  void PrintANSIEscapeSequence(const aurostd::xoption& color,FILE* fstr){
    if(color.option==FALSE){return;}
    if(color.flag("COLOR==GREEN")){cursor_fore_green(fstr);return;}
    if(color.flag("COLOR==CYAN")){cursor_fore_cyan(fstr);return;}
    if(color.flag("COLOR==YELLOW")){cursor_fore_yellow(fstr);return;}
    if(color.flag("COLOR==RED")){cursor_fore_red(fstr);return;}
  }

  void PrintMessageStream(ostringstream &stream,bool quiet,std::ostream& oss) {ofstream FileMESSAGE;return PrintMessageStream(FileMESSAGE,stream,quiet,oss);} //CO20200624
  void PrintMessageStream(ofstream &FileMESSAGE,ostringstream &stream,bool quiet,std::ostream& oss) {bool osswrite=true;return PrintMessageStream(FileMESSAGE,stream,quiet,osswrite,oss);} //CO20200624
  void PrintMessageStream(ofstream &FileMESSAGE,ostringstream &stream,bool quiet,bool osswrite,std::ostream& oss) {
    //[CO20200624 - OBSOLETE]FileMESSAGE << stream.str().c_str(); FileMESSAGE.flush();
    //[CO20200624 - OBSOLETE]if(osswrite) {if(!quiet) {oss << stream.str().c_str();oss.flush();}}
    //[CO20200624 - OBSOLETE]// cerr << stream.str().c_str(); cerr.flush();

    //CO20181226 - split by newlines and print separately
    vector<string> message_parts,_message_parts;
    string stream_str=stream.str();aurostd::StringstreamClean(stream);
    aurostd::string2vectorstring(stream_str,_message_parts);
    for(uint i=0;i<_message_parts.size();i++){
      if(!aurostd::RemoveWhiteSpacesFromTheBack(_message_parts[i]).empty()){
        message_parts.push_back(_message_parts[i]);
      }
    }
    if(message_parts.size()==0){return;}

    //CO20220129 - input quiet is always XHOST.QUIET, this is legacy: aflow.h (XHOST) was previously not made available in aurostd
    //keep redundant for now, just in case in the future we change our mind about including aflow.h
    //CO20220630 - note about osswrite, it is redundant with quiet, so it would be nice to get rid of it in the future
    //but it would require a full overhaul of many critical aflow printing functions
    //better not to touch and leave the overloads
    //bool verbose=(!XHOST.QUIET && !quiet && osswrite);  //ME20220503 - XHOST.QUIET should be part of quiet to allow for whitelisting
    bool verbose=(!quiet && osswrite);
    if(XHOST.QUIET_GLOBAL){verbose=false;}  //CO20220630
    if((&oss==&cout) && XHOST.QUIET_COUT){verbose=false;}
    if((&oss==&cerr) && XHOST.QUIET_CERR){verbose=false;}
    bool fancy_print=(!XHOST.vflag_control.flag("WWW")&&!XHOST.vflag_control.flag("NO_FANCY_PRINT"));  //CO20200404 - new web flag

    FILE* fstr=stdout;
    if(&oss==&std::cerr){fstr=stderr;}

    //COLOR CANNOT BE A STRING, this construction will cause errors for the compiler
    //string color="\033[32m";
    //printf(color.c_str());
    //the compiler needs to verify that you are not printf'ing junk
    //so it needs to be a direct injection of code that the compiler can check
    //reference aurostd.h: CO20200624 START - adding from Jahnatek
    for(uint i=0;i<message_parts.size();i++){FileMESSAGE << message_parts[i] << endl;}  //flush included in endl
    if(verbose){
      string::size_type loc;
      string str2search="";  //replicate old behavior, look for ERROR coming from logger() which has two pre spaces
      aurostd::xoption color;color.clear(); //use xoption: .option is global color flag (do we have color?), and .vxscheme tells me which color
      if(fancy_print){
        string message=stream_str;
        //COMPLETE - START
        if(color.option==FALSE){
          str2search="  COMPLETE ";  //replicate old behavior, look for ERROR coming from logger() which has two pre spaces
          if(message.find(str2search)!=string::npos){color.option=TRUE;color.flag("COLOR==GREEN",TRUE);} //green
        }
        //COMPLETE - END
        //NOTICE - START
        if(color.option==FALSE){
          str2search="  NOTICE ";  //replicate old behavior, look for ERROR coming from logger() which has two pre spaces
          if(message.find(str2search)!=string::npos){color.option=TRUE;color.flag("COLOR==CYAN",TRUE);} //cyan
        }
        //NOTICE - END
      }

      //cursor_fore_green(fstr)
      //cursor_fore_cyan(fstr)

      if(color.option==FALSE){fancy_print=false;}  //add others as needed
      if(fancy_print) PrintANSIEscapeSequence(color,fstr);
      for(uint i=0;i<message_parts.size();i++){
        loc=(!str2search.empty()?message_parts[i].find(str2search):string::npos);
        oss << message_parts[i].substr(0,loc);
        if(loc!=string::npos){
          //colors see here: https://en.m.wikipedia.org/wiki/ANSI_escape_code
          if(fancy_print) cursor_attr_none(fstr);             // turn off all cursor attributes
          if(fancy_print) PrintANSIEscapeSequence(color,fstr); // color
          if(fancy_print) {cursor_attr_blink(fstr);cursor_attr_bold(fstr);} // bold+blink
          oss << str2search;
          if(fancy_print) cursor_attr_none(fstr);             // turn off all cursor attributes
          if(fancy_print) PrintANSIEscapeSequence(color,fstr); // color
          oss << message_parts[i].substr(loc+str2search.size(),string::npos);
        }
        oss << endl;  //flush included in endl
      }
      if(fancy_print) cursor_attr_none(fstr);  // turn off all cursor attributes
    }
  }

  //[CO20200624 - OBSOLETE]void PrintMessageStream(ostringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  if(!quiet) {cout << stream.str().c_str();cout.flush();}
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  //[CO20200624 - OBSOLETE]void PrintErrorStream(ofstream &FileMESSAGE,ostringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  if(quiet) {;} // phony just to keep quiet busy
  //[CO20200624 - OBSOLETE]  FileMESSAGE << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl; FileMESSAGE.flush();
  //[CO20200624 - OBSOLETE]  cout << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl;cout.flush();
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  //[CO20200624 - OBSOLETE]void PrintErrorStream(std::ostream &FileMESSAGE,ostringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  if(quiet) {;} // phony just to keep quiet busy
  //[CO20200624 - OBSOLETE]  FileMESSAGE << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl; FileMESSAGE.flush();
  //[CO20200624 - OBSOLETE]  cout << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl;cout.flush();
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  //CO20200624 - no std::ostream& oss input: THIS MUST GO TO CERR
  void PrintErrorStream(ostringstream &stream,bool quiet) {ofstream FileMESSAGE;return PrintErrorStream(FileMESSAGE,stream,quiet);} //CO20200624
  void PrintErrorStream(ofstream &FileMESSAGE,ostringstream &stream,bool quiet) {bool osswrite=true;return PrintErrorStream(FileMESSAGE,stream,quiet,osswrite);} //CO20200624
  void PrintErrorStream(ofstream &FileMESSAGE,ostringstream &stream,bool quiet,bool osswrite) {
    //[CO20200624 - OBSOLETE]if(quiet) {;} // phony just to keep quiet busy
    //[CO20200624 - OBSOLETE]if(osswrite) {;} // phony just to keep quiet busy
    //[CO20200624 - OBSOLETE]FileMESSAGE << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl; FileMESSAGE.flush();
    //[CO20200624 - OBSOLETE]oss << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl;cout.flush();
    //[CO20200624 - OBSOLETE]// cerr << stream.str().c_str(); cerr.flush();

    //CO20181226 - split by newlines and print separately
    vector<string> message_parts,_message_parts;
    string stream_str=stream.str();aurostd::StringstreamClean(stream);
    aurostd::string2vectorstring(stream_str,_message_parts);
    for(uint i=0;i<_message_parts.size();i++){
      if(!aurostd::RemoveWhiteSpacesFromTheBack(_message_parts[i]).empty()){
        message_parts.push_back(_message_parts[i]);
      }
    }
    if(message_parts.size()==0){return;}

    //CO20220129 - input quiet is always XHOST.QUIET, this is legacy: aflow.h (XHOST) was previously not made available in aurostd
    //keep redundant for now, just in case in the future we change our mind about including aflow.h
    //CO20220630 - note about osswrite, it is redundant with quiet, so it would be nice to get rid of it in the future
    //but it would require a full overhaul of many critical aflow printing functions
    //better not to touch and leave the overloads
    //bool verbose=(!XHOST.QUIET && !quiet && osswrite);  //[CO2010315 - not always, removing for OUTCARs read during vasp runs]verbose=true; //ALWAYS! //ME20220503 - XHOST.QUIET should be part of quiet to allow for whitelisting
    bool verbose=(!quiet && osswrite);
    if(XHOST.QUIET_GLOBAL){verbose=false;}  //CO20220630
    if(XHOST.QUIET_CERR){verbose=false;}
    bool fancy_print=(!XHOST.vflag_control.flag("WWW")&&!XHOST.vflag_control.flag("NO_FANCY_PRINT"));  //CO20200404 - new web flag

    FILE* fstr=stderr;

    FileMESSAGE << ErrorBarString << endl;
    for(uint i=0;i<message_parts.size();i++){FileMESSAGE << message_parts[i] << endl;}  //flush included in endl
    FileMESSAGE << ErrorBarString << endl;
    if(verbose){
      string::size_type loc;
      string str2search="  ERROR ";  //replicate old behavior, look for ERROR coming from logger() which has two pre spaces
      std::ostream& oss=std::cerr;
      if(fancy_print) cursor_fore_red(fstr);  // red
      oss << ErrorBarString << endl;  //flush included in endl
      for(uint i=0;i<message_parts.size();i++){
        loc=message_parts[i].find(str2search);
        oss << message_parts[i].substr(0,loc);
        if(loc!=string::npos){
          if(fancy_print) cursor_attr_none(fstr);     // turn off all cursor attributes
          if(fancy_print) cursor_fore_red(fstr);      // red
          if(fancy_print) {cursor_attr_blink(fstr);cursor_attr_bold(fstr);} // bold+blink
          oss << str2search;
          if(fancy_print) cursor_attr_none(fstr);     // turn off all cursor attributes
          if(fancy_print) cursor_fore_red(fstr);      // red
          oss << message_parts[i].substr(loc+str2search.size(),string::npos);
        }
        oss << endl;  //flush included in endl
      }
      oss << ErrorBarString << endl;  //flush included in endl
      if(fancy_print) cursor_attr_none(fstr);  // turn off all cursor attributes
    }
  }

  //[CO20200624 - OBSOLETE]void PrintErrorStream(ostringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  if(quiet) {;} // phony just to keep quiet busy
  //[CO20200624 - OBSOLETE]  cout << stream.str().c_str();cout.flush();
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  //[CO20200624 - OBSOLETE]void PrintWarningStream(ofstream &FileMESSAGE,ostringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  if(quiet) {;} // phony just to keep quiet busy
  //[CO20200624 - OBSOLETE]  FileMESSAGE << stream.str().c_str(); FileMESSAGE.flush();
  //[CO20200624 - OBSOLETE]  cout << stream.str().c_str();cout.flush();
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  //[CO20200624 - OBSOLETE]void PrintWarningStream(std::ostream &FileMESSAGE,ostringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  if(quiet) {;} // phony just to keep quiet busy
  //[CO20200624 - OBSOLETE]  FileMESSAGE << stream.str().c_str(); FileMESSAGE.flush();
  //[CO20200624 - OBSOLETE]  cout << stream.str().c_str();cout.flush();
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  //CO20200624 - no std::ostream& oss input: THIS MUST GO TO CERR
  void PrintWarningStream(ostringstream &stream,bool quiet) {ofstream FileMESSAGE;return PrintWarningStream(FileMESSAGE,stream,quiet);} //CO20200624
  void PrintWarningStream(ofstream &FileMESSAGE,ostringstream &stream,bool quiet) {bool osswrite=true;return PrintWarningStream(FileMESSAGE,stream,quiet,osswrite);} //CO20200624
  void PrintWarningStream(ofstream &FileMESSAGE,ostringstream &stream,bool quiet,bool osswrite) {
    //[CO20200624 - OBSOLETE]if(quiet) {;} // phony just to keep quiet busy
    //[CO20200624 - OBSOLETE]FileMESSAGE << stream.str().c_str(); FileMESSAGE.flush();
    //[CO20200624 - OBSOLETE]if(osswrite) {oss << stream.str().c_str();oss.flush();}
    //[CO20200624 - OBSOLETE]// cerr << stream.str().c_str(); cerr.flush();

    //CO20181226 - split by newlines and print separately
    vector<string> message_parts,_message_parts;
    string stream_str=stream.str();aurostd::StringstreamClean(stream);
    aurostd::string2vectorstring(stream_str,_message_parts);
    for(uint i=0;i<_message_parts.size();i++){
      if(!aurostd::RemoveWhiteSpacesFromTheBack(_message_parts[i]).empty()){
        message_parts.push_back(_message_parts[i]);
      }
    }
    if(message_parts.size()==0){return;}

    //CO20220129 - input quiet is always XHOST.QUIET, this is legacy: aflow.h (XHOST) was previously not made available in aurostd
    //keep redundant for now, just in case in the future we change our mind about including aflow.h
    //CO20220630 - note about osswrite, it is redundant with quiet, so it would be nice to get rid of it in the future
    //but it would require a full overhaul of many critical aflow printing functions
    //better not to touch and leave the overloads
    //bool verbose=(!XHOST.QUIET && !quiet && osswrite);  //[CO2010315 - not always, removing for OUTCARs read during vasp runs]verbose=true; //ALWAYS! //ME20220503 - XHOST.QUIET should be part of quiet to allow for whitelisting
    bool verbose=(!quiet && osswrite);
    if(XHOST.QUIET_GLOBAL){verbose=false;}  //CO20220630
    if(XHOST.QUIET_CERR){verbose=false;}
    bool fancy_print=(!XHOST.vflag_control.flag("WWW")&&!XHOST.vflag_control.flag("NO_FANCY_PRINT"));  //CO20200404 - new web flag

    FILE* fstr=stderr;

    FileMESSAGE << WarningBarString << endl;
    for(uint i=0;i<message_parts.size();i++){FileMESSAGE << message_parts[i] << endl;}  //flush included in endl
    FileMESSAGE << WarningBarString << endl;
    if(verbose){
      string::size_type loc;
      string str2search="  WARNING ";  //replicate old behavior, look for WARNING coming from logger() which has two pre spaces
      std::ostream& oss=std::cerr;
      if(fancy_print) cursor_fore_yellow(fstr);   // yellow
      oss << WarningBarString << endl;  //flush included in endl
      for(uint i=0;i<message_parts.size();i++){
        loc=message_parts[i].find(str2search);
        oss << message_parts[i].substr(0,loc);
        if(loc!=string::npos){
          if(fancy_print) cursor_attr_none(fstr);     // turn off all cursor attributes
          if(fancy_print) cursor_fore_yellow(fstr);   // yellow
          if(fancy_print) {cursor_attr_blink(fstr);cursor_attr_bold(fstr);} // bold+blink
          oss << str2search;
          if(fancy_print) cursor_attr_none(fstr);     // turn off all cursor attributes
          if(fancy_print) cursor_fore_yellow(fstr);   // yellow
          oss << message_parts[i].substr(loc+str2search.size(),string::npos);
        }
        oss << endl;  //flush included in endl
      }
      oss << WarningBarString << endl;  //flush included in endl
      if(fancy_print) cursor_attr_none(fstr);  // turn off all cursor attributes
    }
  }

  //[CO20200624 - OBSOLETE]void PrintWarningStream(ostringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  if(quiet) {;} // phony just to keep quiet busy
  //[CO20200624 - OBSOLETE]  cout << stream.str().c_str();cout.flush();
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  //[CO20200624 - OBSOLETE]// with stringstream
  //[CO20200624 - OBSOLETE]void PrintMessageStream(ofstream &FileMESSAGE,stringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  FileMESSAGE << stream.str().c_str(); FileMESSAGE.flush();
  //[CO20200624 - OBSOLETE]  if(!quiet) {cout << stream.str().c_str();cout.flush();}
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  //[CO20200624 - OBSOLETE]void PrintMessageStream(std::ostream &FileMESSAGE,stringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  FileMESSAGE << stream.str().c_str(); FileMESSAGE.flush();
  //[CO20200624 - OBSOLETE]  if(!quiet) {cout << stream.str().c_str();cout.flush();}
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  void PrintMessageStream(stringstream &stream,bool quiet,std::ostream& oss) {ofstream FileMESSAGE;return PrintMessageStream(FileMESSAGE,stream,quiet,oss);} //CO20200624
  void PrintMessageStream(ofstream &FileMESSAGE,stringstream &stream,bool quiet,std::ostream& oss) {bool osswrite=true;return PrintMessageStream(FileMESSAGE,stream,quiet,osswrite,oss);} //CO20200624
  void PrintMessageStream(ofstream &FileMESSAGE,stringstream &stream,bool quiet,bool osswrite,std::ostream& oss) {ostringstream omess;omess << stream.str();aurostd::StringstreamClean(stream);return PrintMessageStream(FileMESSAGE,omess,quiet,osswrite,oss);}

  //[CO20200624 - OBSOLETE]void PrintMessageStream(stringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  if(!quiet) {cout << stream.str().c_str();cout.flush();}
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  //[CO20200624 - OBSOLETE]void PrintErrorStream(ofstream &FileMESSAGE,stringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  if(quiet) {;} // phony just to keep quiet busy
  //[CO20200624 - OBSOLETE]  FileMESSAGE << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl; FileMESSAGE.flush();
  //[CO20200624 - OBSOLETE]  cout << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl;cout.flush();
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  //[CO20200624 - OBSOLETE]void PrintErrorStream(std::ostream &FileMESSAGE,stringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  if(quiet) {;} // phony just to keep quiet busy
  //[CO20200624 - OBSOLETE]  FileMESSAGE << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl; FileMESSAGE.flush();
  //[CO20200624 - OBSOLETE]  cout << ErrorBarString << endl << stream.str().c_str() << ErrorBarString << endl;cout.flush();
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  void PrintErrorStream(stringstream &stream,bool quiet) {ofstream FileMESSAGE;return PrintErrorStream(FileMESSAGE,stream,quiet);} //CO20200624
  void PrintErrorStream(ofstream &FileMESSAGE,stringstream &stream,bool quiet) {bool osswrite=true;return PrintErrorStream(FileMESSAGE,stream,quiet,osswrite);} //CO20200624
  void PrintErrorStream(ofstream &FileMESSAGE,stringstream &stream,bool quiet,bool osswrite) {ostringstream omess;omess << stream.str();aurostd::StringstreamClean(stream);return PrintErrorStream(FileMESSAGE,omess,quiet,osswrite);}

  //[CO20200624 - OBSOLETE]void PrintErrorStream(stringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  if(quiet) {;} // phony just to keep quiet busy
  //[CO20200624 - OBSOLETE]  cout << stream.str().c_str();cout.flush();
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  //[CO20200624 - OBSOLETE]void PrintWarningStream(ofstream &FileMESSAGE,stringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  if(quiet) {;} // phony just to keep quiet busy
  //[CO20200624 - OBSOLETE]  FileMESSAGE << stream.str().c_str(); FileMESSAGE.flush();
  //[CO20200624 - OBSOLETE]  cout << stream.str().c_str();cout.flush();
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  //[CO20200624 - OBSOLETE]void PrintWarningStream(std::ostream &FileMESSAGE,stringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  if(quiet) {;} // phony just to keep quiet busy
  //[CO20200624 - OBSOLETE]  FileMESSAGE << stream.str().c_str(); FileMESSAGE.flush();
  //[CO20200624 - OBSOLETE]  cout << stream.str().c_str();cout.flush();
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  void PrintWarningStream(stringstream &stream,bool quiet) {ofstream FileMESSAGE;return PrintWarningStream(FileMESSAGE,stream,quiet);} //CO20200624
  void PrintWarningStream(ofstream &FileMESSAGE,stringstream &stream,bool quiet) {bool osswrite=true;return PrintWarningStream(FileMESSAGE,stream,quiet,osswrite);} //CO20200624
  void PrintWarningStream(ofstream &FileMESSAGE,stringstream &stream,bool quiet,bool osswrite) {ostringstream omess;omess << stream.str();aurostd::StringstreamClean(stream);return PrintWarningStream(FileMESSAGE,omess,quiet,osswrite);}

  //[CO20200624 - OBSOLETE]void PrintWarningStream(stringstream &stream,bool quiet) {
  //[CO20200624 - OBSOLETE]  if(quiet) {;} // phony just to keep quiet busy
  //[CO20200624 - OBSOLETE]  cout << stream.str().c_str();cout.flush();
  //[CO20200624 - OBSOLETE]  // cerr << stream.str().c_str(); cerr.flush();
  //[CO20200624 - OBSOLETE]  aurostd::StringstreamClean(stream);
  //[CO20200624 - OBSOLETE]}

  // ***************************************************************************
  // Execute Streams/Strings/C_strings
  // ***************************************************************************
  bool execute(ostringstream &command) {
    // cerr << "COMMAND " <<  command.str().c_str() << endl;
    //[CO20200624 - OBSOLETE]system(command.str().c_str());
    execute(command.str()); //CO20200624
    aurostd::StringstreamClean(command);
    return TRUE;
  }

  bool execute(stringstream &command) {
    // cerr << "COMMAND " <<  command.str().c_str() << endl;
    //[CO20200624 - OBSOLETE]system(command.str().c_str());
    execute(command.str()); //CO20200624
    aurostd::StringstreamClean(command);
    return TRUE;
  }

  bool execute(const string& _command) {
#ifdef AFLOW_MULTITHREADS_ENABLE  //CO+HE20221116
    std::lock_guard<std::mutex> lk(xthread_execute);  //prevents race conditions likely caused by system calls
#endif
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    // cerr << "COMMAND " <<  command.c_str() << endl;
    string command=aurostd::CleanCommand4Execute(_command); //CO20200624
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " command.c_str()=\"" << command.c_str() << "\"" << endl;}
    system(command.c_str());
    //   command="";
    return TRUE;
  }

#ifdef _stringcharstar_
  bool execute(char* _command) {
    // cerr << "COMMAND " <<  command << endl;
    //[CO20200624 - OBSOLETE]system(command);
    string command=std::string(_command); //CO20200624
    execute(command);
    return TRUE;
  }
#endif

  // ***************************************************************************
  // Execute vectors/deque of Strings
  // ***************************************************************************
  bool execute(const deque<string>& vcommand) {
    for(uint i=0;i<vcommand.size();i++)
      execute(vcommand[i]);
    return TRUE;
  }
  bool execute(const vector<string>& vcommand) {
    for(uint i=0;i<vcommand.size();i++)
      execute(vcommand[i]);
    return TRUE;
  }

  // ***************************************************************************
  // Execute & Report Streams/Strings/C_strings
  // ***************************************************************************
  string execute2string(const string& _command,FSIO fsio) { //CO20200624 - added file system IO mode
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    // bool INIT_VERBOSE=TRUE;
    // cerr << "COMMAND " <<  command << endl;

    //CO20200624 START - some command cleanup
    string command=aurostd::CleanCommand4Execute(_command);
    //if(command.find("; ")!=string::npos){command="( "+command+" )";}  //put to subshell for IO redirection; https://www.gnu.org/software/bash/manual/html_node/Command-Grouping.html#Command-Grouping
    command="( "+command+" )";  //ALWAYS put to subshell for IO redirection; https://www.gnu.org/software/bash/manual/html_node/Command-Grouping.html#Command-Grouping
    //CO20200624 END - some command cleanup

    stringstream strstream,cmdstream;
    string file=aurostd::TmpFileCreate("execute_report");
    if(fsio==stdouterr_fsio){cmdstream << "bash -c \"" << command << " &> " << file << "\"";}  //CO20200624 //SD20220311 - force bash, &> does not work in sh; be careful with quotes within quotes, althought it seems to work
    else if(fsio==stderr_fsio){cmdstream << command << " 2> " << file;} //CO20200624
    else{cmdstream << command << " > " << file;} //CO20200624
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " cmdstream=\"" << cmdstream.str() << "\"" << endl;}
    system(cmdstream.str().c_str());
    // command="";
    strstream << aurostd::file2string(file);
    aurostd::StringstreamClean(cmdstream);
#ifndef _AFLOW_TEMP_PRESERVE_
    aurostd::RemoveFile(file);
#endif
    string strout=strstream.str();
    if(strout.length()>0)
      if(strout.at(strout.length()-1)=='\n')
        strout.erase(strout.length()-1);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " strout=\"" << strout << "\"" << endl;}
    return strout;
  }

  string execute2string(ostringstream &command,FSIO fsio) { //CO20200624 - added file system IO mode
    string command_str=command.str();
    aurostd::StringstreamClean(command);
    return execute2string(command_str,fsio);  //CO20200624
  }

  string execute2string(stringstream &command,FSIO fsio) { //CO20200624 - added file system IO mode
    string command_str=command.str();
    aurostd::StringstreamClean(command);
    return execute2string(command_str,fsio);  //CO20200624
  }

  vector<string> execute2string(const vector<string>& vcommand,FSIO fsio) { //CO20200624 - added file system IO mode
    vector<string> out;
    for(uint i=0;i<vcommand.size();i++)
      out.push_back(execute2string(vcommand[i],fsio));  //CO20200624
    return out;
  }

  deque<string> execute2string(const deque<string>& vcommand,FSIO fsio) { //CO20200624 - added file system IO mode
    deque<string> out;
    for(uint i=0;i<vcommand.size();i++)
      out.push_back(execute2string(vcommand[i],fsio));  //CO20200624
    return out;
  }

#ifdef _stringcharstar_
  string execute2string(char* command,FSIO fsio) { //CO20200624 - added file system IO mode
    string command_str=string(command);
    return execute2string(command_str,fsio);  //CO20200624
  }
#endif

  string CleanCommand4Execute(const string& _command){ //CO20200624
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " command(pre )=\"" << _command << "\"" << endl;}
    //CO20200624 START - some command cleanup
    vector<string> vtokens,vtokens_new;
    aurostd::string2vectorstring(aurostd::RemoveWhiteSpacesFromTheFrontAndBack(_command),vtokens);
    string tmp="";
    uint i=0;
    for(i=0;i<vtokens.size();i++){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vtokens[i=" << i << "](pre )=\"" << vtokens[i] << "\"" << endl;}
      tmp=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vtokens[i]);
      aurostd::CleanStringASCII_InPlace(tmp);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vtokens[i=" << i << "](post)=\"" << tmp << "\"" << endl;}
      if(!tmp.empty()){vtokens_new.push_back(tmp);}
    }
    if(vtokens_new.size()==0){return "";}
    //[CO20210312 - must be smarter, could end with && or ;]string command=aurostd::joinWDelimiter(vtokens_new,"; ");
    string command="";
    uint len=0;
    bool add_semicolon=false;
    for(i=0;i<vtokens_new.size();i++){ //vtokens_new has no empty entries, so we don't need to check again
      const string& cmd=vtokens_new[i];
      command+=cmd;
      if(i<vtokens_new.size()-1){
        add_semicolon=false;
        len=cmd.size();
        if(!(cmd[len-1]=='&'||cmd[len-1]=='|'||cmd[len-1]==';')){add_semicolon=true;}
        if(add_semicolon){command+="; ";}
        else{command+=" ";} //add a space, looks good for ' this && that '
      }
    }
    //CO20200624 END - some command cleanup
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " command(post)=\"" << command << "\"" << endl;}
    return command;
  }

  // ***************************************************************************
  // Execute & Report Int Streams/Strings/C_strings
  // ***************************************************************************
  template<class utype> utype execute2utype(ostringstream &command) {
    return (utype) aurostd::string2utype<utype>(execute2string(command));
  }

  template<class utype> utype execute2utype(stringstream &command) {
    return (utype) aurostd::string2utype<utype>(execute2string(command));
  }

  template<class utype> utype execute2utype(string command) {
    return (utype) aurostd::string2utype<utype>(execute2string(command));
  }

  template<class utype> vector<utype> execute2utype(vector<utype> vcommand) {
    vector<utype> out;
    for(uint i=0;i<vcommand.size();i++)
      out.push_back((utype) execute2utype<utype>(vcommand[i]));
    return out;
  }

  template<class utype> deque<utype> execute2utype(deque<utype> vcommand) {
    deque<utype> out;
    for(uint i=0;i<vcommand.size();i++)
      out.push_back((utype) execute2utype<utype>(vcommand[i]));
    return out;
  }

#ifdef _stringcharstar_
  template<class utype> utype execute2utype(char* command) {
    return (utype) aurostd::string2utype<utype>(execute2string(command));
  }
#endif

  // ***************************************************************************
  // Sleep
  // ***************************************************************************
  unsigned int Sleep(unsigned int seconds) {
    //  ostringstream aus;
    // aus << "sleep " << (int) seconds << " " << endl;
    // aurostd::execute(aus);
    return sleep(seconds);
  }

  // *******************************************************************************************
  // *******************************************************************************************
  vector<string> GrepFile(const string& filename,const string& keyword,bool RemoveWS,bool RemoveComments){  //CO20210623 - update after integrating new substring2bool (RemoveWS,RemoveComments)
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    vector<string> vout;
    if(aurostd::FileExist(filename)==false){return vout;}

    if(RemoveWS && RemoveComments){;} //keep busy until we update substring2bool

    string strline="";
    ifstream FileStream;FileStream.open(filename.c_str(),std::ios::in);
    while(getline(FileStream,strline)){
      if(aurostd::substring2bool(strline,keyword,RemoveWS)) {
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " found matching line: \"" << strline << "\"" << endl;}
        vout.push_back(strline);
      }
    }

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}
    return vout;
  }

  // *******************************************************************************************
  // *******************************************************************************************
  //[SD20220501 - OBSOLETE]bool ExtractToFileEXPLICIT(ifstream& FileIN,const string& FileNameOUTPUT,const string& Keyword) {        // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  ofstream FileOUTPUT;
  //[SD20220501 - OBSOLETE]  FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
  //[SD20220501 - OBSOLETE]  string strline="";
  //[SD20220501 - OBSOLETE]  FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
  //[SD20220501 - OBSOLETE]  bool status=FALSE;
  //[SD20220501 - OBSOLETE]  while(getline(FileIN,strline)) {
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(strline,Keyword)) {
  //[SD20220501 - OBSOLETE]      FileOUTPUT << strline.substr(strline.find(Keyword)+Keyword.length()) << endl;
  //[SD20220501 - OBSOLETE]      status=TRUE;
  //[SD20220501 - OBSOLETE]    }
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
  //[SD20220501 - OBSOLETE]  FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
  //[SD20220501 - OBSOLETE]  return status;  // return FALSE if the keyword was never found
  //[SD20220501 - OBSOLETE]}

  //[SD20220501 - OBSOLETE]bool ExtractToFileEXPLICIT(const string& StringIN,const string& FileNameOUTPUT,const string& Keyword) {        // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  ofstream FileOUTPUT;
  //[SD20220501 - OBSOLETE]  FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
  //[SD20220501 - OBSOLETE]  string strline="";
  //[SD20220501 - OBSOLETE]  bool status=FALSE;
  //[SD20220501 - OBSOLETE]  vector<string> tokens;
  //[SD20220501 - OBSOLETE]  aurostd::string2tokens(StringIN,tokens,"\n");
  //[SD20220501 - OBSOLETE]  for(uint i=0;i<tokens.size();i++) {
  //[SD20220501 - OBSOLETE]    strline=tokens[i];
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(strline,Keyword)) {
  //[SD20220501 - OBSOLETE]      FileOUTPUT << strline.substr(strline.find(Keyword)+Keyword.length()) << endl;
  //[SD20220501 - OBSOLETE]      status=TRUE;
  //[SD20220501 - OBSOLETE]    }
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
  //[SD20220501 - OBSOLETE]  return status;  // return FALSE if the keyword was never found
  //[SD20220501 - OBSOLETE]}

  // *******************************************************************************************
  //[SD20220501 - OBSOLETE]bool ExtractToFileEXPLICIT(ifstream& FileIN,const string& FileNameOUTPUT,const string& Keyword_start,const string& Keyword_stop) {        // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  ofstream FileOUTPUT;
  //[SD20220501 - OBSOLETE]  FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
  //[SD20220501 - OBSOLETE]  string strline="";
  //[SD20220501 - OBSOLETE]  FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
  //[SD20220501 - OBSOLETE]  bool status=FALSE;
  //[SD20220501 - OBSOLETE]  while(getline(FileIN,strline)) {
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(strline,Keyword_stop))  status=FALSE;
  //[SD20220501 - OBSOLETE]    if(status) FileOUTPUT << strline << endl;
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
  //[SD20220501 - OBSOLETE]  FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
  //[SD20220501 - OBSOLETE]  return status;  // return FALSE if something got messed up
  //[SD20220501 - OBSOLETE]}

  //[SD20220501 - OBSOLETE]bool ExtractToFileEXPLICIT(const string& StringIN,const string& FileNameOUTPUT,const string& Keyword_start,const string& Keyword_stop) {        // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  ofstream FileOUTPUT;
  //[SD20220501 - OBSOLETE]  FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
  //[SD20220501 - OBSOLETE]  string strline="";
  //[SD20220501 - OBSOLETE]  bool status=FALSE;
  //[SD20220501 - OBSOLETE]  vector<string> tokens;
  //[SD20220501 - OBSOLETE]  aurostd::string2tokens(StringIN,tokens,"\n");
  //[SD20220501 - OBSOLETE]  for(uint i=0;i<tokens.size();i++) {
  //[SD20220501 - OBSOLETE]    strline=tokens[i];
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(strline,Keyword_stop))  status=FALSE;
  //[SD20220501 - OBSOLETE]    if(status) FileOUTPUT << strline << endl;
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
  //[SD20220501 - OBSOLETE]  return status;  // return FALSE if something got messed up
  //[SD20220501 - OBSOLETE]}

  // *******************************************************************************************
  // *******************************************************************************************
  //[SD20220501 - OBSOLETE]bool ExtractToStringEXPLICIT(ifstream& FileIN,string& StringOUTPUT,const string& Keyword) {        // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  string strline="";
  //[SD20220501 - OBSOLETE]  FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
  //[SD20220501 - OBSOLETE]  bool status=FALSE;
  //[SD20220501 - OBSOLETE]  while(getline(FileIN,strline)) {
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(strline,Keyword)) {
  //[SD20220501 - OBSOLETE]      StringOUTPUT=StringOUTPUT+strline.substr(strline.find(Keyword)+Keyword.length())+"\n";
  //[SD20220501 - OBSOLETE]      status=TRUE;
  //[SD20220501 - OBSOLETE]    }
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
  //[SD20220501 - OBSOLETE]  return status;  // return FALSE if the keyword was never found
  //[SD20220501 - OBSOLETE]}

  //[SD20220501 - OBSOLETE]bool ExtractToStringEXPLICIT(const string& StringIN,string& StringOUTPUT,const string& Keyword) {        // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  string strline="";
  //[SD20220501 - OBSOLETE]  bool status=FALSE;
  //[SD20220501 - OBSOLETE]  vector<string> tokens;
  //[SD20220501 - OBSOLETE]  aurostd::string2tokens(StringIN,tokens,"\n");
  //[SD20220501 - OBSOLETE]  for(uint i=0;i<tokens.size();i++) {
  //[SD20220501 - OBSOLETE]    strline=tokens[i];
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(strline,Keyword)) {
  //[SD20220501 - OBSOLETE]      StringOUTPUT=StringOUTPUT+strline.substr(strline.find(Keyword)+Keyword.length())+"\n";
  //[SD20220501 - OBSOLETE]      status=TRUE;
  //[SD20220501 - OBSOLETE]    }
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  return status;  // return FALSE if the keyword was never found
  //[SD20220501 - OBSOLETE]}

  // *******************************************************************************************
  //[SD20220501 - OBSOLETE]bool ExtractToStringEXPLICIT(ifstream& FileIN,string& StringOUTPUT,const string& Keyword_start,const string& Keyword_stop) {        // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  string strline="";
  //[SD20220501 - OBSOLETE]  FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
  //[SD20220501 - OBSOLETE]  bool status=FALSE;
  //[SD20220501 - OBSOLETE]  while(getline(FileIN,strline)) {
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(strline,Keyword_stop))  status=FALSE;
  //[SD20220501 - OBSOLETE]    if(status) StringOUTPUT=StringOUTPUT+strline+"\n";
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
  //[SD20220501 - OBSOLETE]  return status;  // return FALSE if something got messed up
  //[SD20220501 - OBSOLETE]}

  //[SD20220501 - OBSOLETE]bool ExtractToStringEXPLICIT(const string& StringIN,string& StringOUTPUT,const string& Keyword_start,const string& Keyword_stop) {        // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  string strline="";
  //[SD20220501 - OBSOLETE]  bool status=FALSE;
  //[SD20220501 - OBSOLETE]  vector<string> tokens;
  //[SD20220501 - OBSOLETE]  aurostd::string2tokens(StringIN,tokens,"\n");
  //[SD20220501 - OBSOLETE]  for(uint i=0;i<tokens.size();i++) {
  //[SD20220501 - OBSOLETE]    strline=tokens[i];
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(strline,Keyword_stop))  status=FALSE;
  //[SD20220501 - OBSOLETE]    if(status) StringOUTPUT=StringOUTPUT+strline+"\n";
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  return status;  // return FALSE if something got messed up
  //[SD20220501 - OBSOLETE]}

  // *******************************************************************************************
  // *******************************************************************************************
  //[SD20220501 - OBSOLETE]bool ExtractToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,const string& Keyword) {        // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  aurostd::StringstreamClean(StringstreamOUTPUT);
  //[SD20220501 - OBSOLETE]  string strline="";
  //[SD20220501 - OBSOLETE]  FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
  //[SD20220501 - OBSOLETE]  bool status=FALSE;
  //[SD20220501 - OBSOLETE]  while(getline(FileIN,strline)) {
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(strline,Keyword)) {
  //[SD20220501 - OBSOLETE]      StringstreamOUTPUT << strline.substr(strline.find(Keyword)+Keyword.length()) << endl;
  //[SD20220501 - OBSOLETE]      status=TRUE;
  //[SD20220501 - OBSOLETE]    }
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
  //[SD20220501 - OBSOLETE]  return status;  // return FALSE if the keyword was never found
  //[SD20220501 - OBSOLETE]}

  // *******************************************************************************************
  //[SD20220501 - OBSOLETE]bool ExtractToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,const string& Keyword_start,const string& Keyword_stop) {    // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  aurostd::StringstreamClean(StringstreamOUTPUT);
  //[SD20220501 - OBSOLETE]  string strline="";
  //[SD20220501 - OBSOLETE]  FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
  //[SD20220501 - OBSOLETE]  bool status=FALSE;
  //[SD20220501 - OBSOLETE]  while(getline(FileIN,strline)) {
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(strline,Keyword_stop))  status=FALSE;
  //[SD20220501 - OBSOLETE]    if(status) StringstreamOUTPUT << strline << endl;
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
  //[SD20220501 - OBSOLETE]  return status;  // return FALSE if something got messed up
  //[SD20220501 - OBSOLETE]}

  //[SD20220501 - OBSOLETE]bool ExtractToStringstreamEXPLICIT(stringstream& StringStreamIN,stringstream& StringstreamOUTPUT,const string& Keyword_start,const string& Keyword_stop) { // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  aurostd::StringstreamClean(StringstreamOUTPUT);
  //[SD20220501 - OBSOLETE]  string StringIN=StringStreamIN.str();
  //[SD20220501 - OBSOLETE]  return ExtractToStringstreamEXPLICIT(StringIN,StringstreamOUTPUT,Keyword_start,Keyword_stop);
  //[SD20220501 - OBSOLETE]}

  //[SD20220501 - OBSOLETE]bool ExtractToStringstreamEXPLICIT(const string& StringIN,stringstream& StringstreamOUTPUT,const string& Keyword_start,const string& Keyword_stop) {  // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  aurostd::StringstreamClean(StringstreamOUTPUT);
  //[SD20220501 - OBSOLETE]  bool status=FALSE;
  //[SD20220501 - OBSOLETE]  vector<string> tokens;
  //[SD20220501 - OBSOLETE]  aurostd::string2tokens(StringIN,tokens,"\n");
  //[SD20220501 - OBSOLETE]  for(uint i=0;i<tokens.size();i++) {
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(tokens[i],Keyword_stop))  status=FALSE;
  //[SD20220501 - OBSOLETE]    if(status) StringstreamOUTPUT << tokens[i] << endl;
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(tokens[i],Keyword_start)) status=TRUE;
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  return status;  // return FALSE if something got messed up
  //[SD20220501 - OBSOLETE]}

  //[SD20220501 - OBSOLETE]bool ExtractToStringstreamEXPLICIT(const string& StringIN,stringstream& StringstreamOUTPUT,const string& Keyword) {  // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  aurostd::StringstreamClean(StringstreamOUTPUT);
  //[SD20220501 - OBSOLETE]  bool status=FALSE;
  //[SD20220501 - OBSOLETE]  vector<string> tokens;
  //[SD20220501 - OBSOLETE]  aurostd::string2tokens(StringIN,tokens,"\n");
  //[SD20220501 - OBSOLETE]  for(uint i=0;i<tokens.size();i++) {
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(tokens[i],Keyword)) StringstreamOUTPUT << tokens[i].substr(tokens[i].find(Keyword)+Keyword.length()) << endl;
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(tokens[i],Keyword)) status=TRUE;
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  return status;  // return FALSE if something got messed up
  //[SD20220501 - OBSOLETE]}

  // *******************************************************************************************
  // *******************************************************************************************
  //[SD20220501 - OBSOLETE]bool ExtractLastToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,const string& Keyword) {   // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  return ExtractToStringstreamEXPLICIT(FileIN,StringstreamOUTPUT,Keyword);
  //[SD20220501 - OBSOLETE]}

  // *******************************************************************************************
  //[SD20220501 - OBSOLETE]bool ExtractLastToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,const string& Keyword_start,const string& Keyword_stop) { // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  bool LDEBUG=(FALSE || XHOST.DEBUG);
  //[SD20220501 - OBSOLETE]  if(LDEBUG) cerr << "LDEBUG: ExtractLastToStringstreamEXPLICIT" << endl;
  //[SD20220501 - OBSOLETE]  aurostd::StringstreamClean(StringstreamOUTPUT);
  //[SD20220501 - OBSOLETE]  vector<string> tokens;
  //[SD20220501 - OBSOLETE]  aurostd::stream2vectorstring(FileIN,tokens);
  //[SD20220501 - OBSOLETE]  int istart=-1,istop=-1;
  //[SD20220501 - OBSOLETE]  for(int i=(int) tokens.size()-1;i>0;i--) {
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(tokens[i],Keyword_stop)  && istop<1)  istop=i-1;
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(tokens[i],Keyword_start) && istart<1) istart=i+1;
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  if(LDEBUG) cerr << "LDEBUG: " << istart << " " << istop << endl;
  //[SD20220501 - OBSOLETE]  if(istart>0 && istop>0) {
  //[SD20220501 - OBSOLETE]    for(int i=istart;i<=istop;i++) StringstreamOUTPUT << tokens[i] << endl;
  //[SD20220501 - OBSOLETE]    if(LDEBUG) cerr << "StringstreamOUTPUT.str()= " << endl << StringstreamOUTPUT.str() << endl;
  //[SD20220501 - OBSOLETE]    return TRUE;
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  return FALSE;
  //[SD20220501 - OBSOLETE]}


  //[SD20220501 - OBSOLETE]bool ExtractLastToStringstreamEXPLICIT(stringstream& StringStreamIN,stringstream& StringstreamOUTPUT,const string& Keyword_start,const string& Keyword_stop) { // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  aurostd::StringstreamClean(StringstreamOUTPUT);
  //[SD20220501 - OBSOLETE]  string StringIN=StringStreamIN.str();
  //[SD20220501 - OBSOLETE]  return ExtractLastToStringstreamEXPLICIT(StringIN,StringstreamOUTPUT,Keyword_start,Keyword_stop);
  //[SD20220501 - OBSOLETE]}

  //[SD20220501 - OBSOLETE]bool ExtractLastToStringstreamEXPLICIT(const string& StringIN,stringstream& StringstreamOUTPUT,const string& Keyword_start,const string& Keyword_stop) {  // AFLOW_FUNCTION_IMPLEMENTATION
  //[SD20220501 - OBSOLETE]  bool LDEBUG=(FALSE || XHOST.DEBUG);
  //[SD20220501 - OBSOLETE]  if(LDEBUG) cerr << "LDEBUG: ExtractLastToStringstreamEXPLICIT" << endl;
  //[SD20220501 - OBSOLETE]  aurostd::StringstreamClean(StringstreamOUTPUT);
  //[SD20220501 - OBSOLETE]  vector<string> tokens;
  //[SD20220501 - OBSOLETE]  aurostd::string2vectorstring(StringIN,tokens);
  //[SD20220501 - OBSOLETE]  int istart=-1,istop=-1;
  //[SD20220501 - OBSOLETE]  for(int i=(int) tokens.size()-1;i>0;i--) {
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(tokens[i],Keyword_stop)  && istop<1)  istop=i-1;
  //[SD20220501 - OBSOLETE]    if(aurostd::substring2bool(tokens[i],Keyword_start) && istart<1) istart=i+1;
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  if(LDEBUG) cerr << "LDEBUG: " << istart << " " << istop << endl;
  //[SD20220501 - OBSOLETE]  if(istart>0 && istop>0) {
  //[SD20220501 - OBSOLETE]    for(int i=istart;i<=istop;i++) StringstreamOUTPUT << tokens[i] << endl;
  //[SD20220501 - OBSOLETE]    if(LDEBUG) cerr << "StringstreamOUTPUT.str()= " << endl << StringstreamOUTPUT.str() << endl;
  //[SD20220501 - OBSOLETE]    return TRUE;
  //[SD20220501 - OBSOLETE]  }
  //[SD20220501 - OBSOLETE]  return FALSE;
  //[SD20220501 - OBSOLETE]}

  // *******************************************************************************************
  // *******************************************************************************************

  bool ExtractJustAfterToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,const string& Keyword_start) {    // AFLOW_FUNCTION_IMPLEMENTATION
    aurostd::StringstreamClean(StringstreamOUTPUT);
    string strline="";
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    bool status=FALSE;
    while(getline(FileIN,strline)) {
      if(status) StringstreamOUTPUT << strline << endl;
      if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
    }
    FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return status;  // return FALSE if something got messed up
  }

  bool ExtractJustAfterToStringstreamEXPLICIT(stringstream& StringStreamIN,stringstream& StringstreamOUTPUT,const string& Keyword_start) { // AFLOW_FUNCTION_IMPLEMENTATION
    string StringIN=StringStreamIN.str();
    return ExtractJustAfterToStringstreamEXPLICIT(StringIN,StringstreamOUTPUT,Keyword_start);
  }

  bool ExtractJustAfterToStringstreamEXPLICIT(const string& StringIN,stringstream& StringstreamOUTPUT,const string& Keyword_start) {  // AFLOW_FUNCTION_IMPLEMENTATION
    aurostd::StringstreamClean(StringstreamOUTPUT);
    bool status=FALSE;
    vector<string> tokens;
    aurostd::string2tokens(StringIN,tokens,"\n");
    for(uint i=0;i<tokens.size();i++) {
      if(status) StringstreamOUTPUT << tokens[i] << endl;
      if(aurostd::substring2bool(tokens[i],Keyword_start)) status=TRUE;
    }
    return status;  // return FALSE if something got messed up
  }

  bool ExtractJustAfterToFileEXPLICIT(ifstream& FileIN,const string& FileNameOUTPUT,const string& Keyword_start) {        // AFLOW_FUNCTION_IMPLEMENTATION
    //[SD20220502 - OBSOLETE]ofstream FileOUTPUT;
    //[SD20220502 - OBSOLETE]FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
    //[SD20220502 - OBSOLETE]string strline="";
    //[SD20220502 - OBSOLETE]FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    //[SD20220502 - OBSOLETE]bool status=FALSE;
    //[SD20220502 - OBSOLETE]while(getline(FileIN,strline)) {
    //[SD20220502 - OBSOLETE]  if(status) FileOUTPUT << strline << endl;
    //[SD20220502 - OBSOLETE]  if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
    //[SD20220502 - OBSOLETE]}
    //[SD20220502 - OBSOLETE]FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    //[SD20220502 - OBSOLETE]FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
    //[SD20220502 - OBSOLETE]return status;  // return FALSE if something got messed up
    stringstream StringstreamOUTPUT;
    if(ExtractJustAfterToStringstreamEXPLICIT(FileIN,StringstreamOUTPUT,Keyword_start)) {
      return aurostd::stringstream2file(StringstreamOUTPUT,FileNameOUTPUT);
    }
    return FALSE;
  }

  bool ExtractJustAfterToStringEXPLICIT(ifstream& FileIN,string& StringOUTPUT,const string& Keyword_start) {        // AFLOW_FUNCTION_IMPLEMENTATION
    //[SD20220502 - OBSOLETE]string strline="";
    //[SD20220502 - OBSOLETE]FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    //[SD20220502 - OBSOLETE]bool status=FALSE;
    //[SD20220502 - OBSOLETE]while(getline(FileIN,strline)) {
    //[SD20220502 - OBSOLETE]  if(status) StringOUTPUT=StringOUTPUT+strline+"\n";
    //[SD20220502 - OBSOLETE]  if(aurostd::substring2bool(strline,Keyword_start)) status=TRUE;
    //[SD20220502 - OBSOLETE]}
    //[SD20220502 - OBSOLETE]FileIN.clear();FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    //[SD20220502 - OBSOLETE]return status;  // return FALSE if something got messed up
    stringstream StringstreamOUTPUT;
    if(ExtractJustAfterToStringstreamEXPLICIT(FileIN,StringstreamOUTPUT,Keyword_start)) {
      StringOUTPUT=StringstreamOUTPUT.str();
      return TRUE;
    }
    return FALSE;
  }

  bool ExtractJustAfterToStringEXPLICIT(const string& StringIN,string& StringOUTPUT,const string& Keyword_start) {  // AFLOW_FUNCTION_IMPLEMENTATION
    //[SD20220502 - OBSOLETE]stringstream StringstreamOUTPUT;
    //[SD20220502 - OBSOLETE]bool out=ExtractJustAfterToStringstreamEXPLICIT(StringIN,StringstreamOUTPUT,Keyword_start);
    //[SD20220502 - OBSOLETE]StringOUTPUT=StringstreamOUTPUT.str();
    //[SD20220501 - OBSOLETE]return out;  // return FALSE if something got messed up
    stringstream StringstreamOUTPUT;
    if(ExtractJustAfterToStringstreamEXPLICIT(StringIN,StringstreamOUTPUT,Keyword_start)) {
      StringOUTPUT=StringstreamOUTPUT.str();
      return TRUE;
    }
    return FALSE;
  }

  // ***************************************************************************
  // Function stream2vectorstring return UINT
  // ***************************************************************************
  // take istream into a vector strings - Stefano Curtarolo
  uint stream2vectorstring(std::istream& istreamIN,vector<string> &vstringout) {
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!istreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      istreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout.push_back(string(tmp));
    }
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.size();  // return FALSE if something got messed up
  }
  // take ifstream into a vector strings - Stefano Curtarolo
  uint stream2vectorstring(std::ifstream& ifstreamIN,vector<string> &vstringout) {
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!ifstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      ifstreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout.push_back(string(tmp));
    }
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.size();  // return FALSE if something got messed up
  } 
  uint stream2vectorstring(std::stringstream& stringstreamIN,vector<string> &vstringout) {
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!stringstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      stringstreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout.push_back(string(tmp));
    }
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.size();  // return FALSE if something got messed up
  }
  uint trimEmptyEdges(vector<string>& vstringout){  //CO20230502
    uint count=vstringout.size();
    //start with front
    while(vstringout.size()){
      if(vstringout.front().find_first_not_of("\t\n ")!=string::npos){ //CO20230502 - fast way to check for empty string //!aurostd::RemoveWhiteSpaces(vstringout.front()).empty())
        break;
      }
      vstringout.erase(vstringout.begin());
      count--;
    }
    //now back
    while(vstringout.size()){
      if(vstringout.back().find_first_not_of("\t\n ")!=string::npos){ //CO20230502 - fast way to check for empty string //aurostd::RemoveWhiteSpaces(vstringout.back()).empty())
        break;
      }
      vstringout.pop_back();
      count--;
    }
    return count;
  }
  uint string2vectorstring(const string& stringIN,vector<string> &vstringout,bool consecutive,bool trim_edges) {  //CO20170613
    //CO mods 20170613
    //we are adding functionality here, because string2tokens will treat "\n\n" same as "\n", but not "\n \n"
    //consecutive will do the following: "sssss" -> <"s","s",...>
    //trim_edges will remove delimiters from beginning and end, similar to consecutive=false behavior
    //return aurostd::string2tokens(stringIN,vstringout,"\n",true);
    uint count=aurostd::string2tokens(stringIN,vstringout,"\n",consecutive);
    if(trim_edges){count=trimEmptyEdges(vstringout);}
    return count;
  }

  // ***************************************************************************
  // Function string2vectorstring return VECTOR
  // ***************************************************************************
  // take sitring into a vector strings - Stefano Curtarolo
  vector<string> stream2vectorstring(std::istream& istreamIN) {
    vector<string> vstringout;
    aurostd::stream2vectorstring(istreamIN,vstringout);
    return vstringout;
  }
  vector<string> stream2vectorstring(std::ifstream& iftreamIN) {
    vector<string> vstringout;
    aurostd::stream2vectorstring(iftreamIN,vstringout);
    return vstringout;
  }
  vector<string> stream2vectorstring(std::stringstream& stringstreamIN) {
    vector<string> vstringout;
    aurostd::stream2vectorstring(stringstreamIN,vstringout);
    return vstringout;
  }
  vector<string> string2vectorstring(const string& stringIN,bool consecutive,bool trim_edges) { //CO20170613
    vector<string> vstringout;
    aurostd::string2vectorstring(stringIN,vstringout,consecutive,trim_edges); //CO20170613
    return vstringout;
  }

  // ***************************************************************************
  // Function liststring2string return string
  // ***************************************************************************
  string liststring2string(string s00,string s01,string s02,string s03,string s04,string s05,string s06,string s07,
      string s08,string s09,string s0A,string s0B,string s0C,string s0D,string s0E,string s0F,
      string s10,string s11,string s12,string s13,string s14,string s15,string s16,string s17,
      string s18,string s19,string s1A,string s1B,string s1C,string s1D,string s1E,string s1F,
      string s20,string s21,string s22,string s23,string s24,string s25,string s26,string s27,
      string s28,string s29,string s2A,string s2B,string s2C,string s2D,string s2E,string s2F) {
    string out="";
    if(s00!="") out+=s00+"\n";
    if(s01!="") out+=s01+"\n";
    if(s02!="") out+=s02+"\n";
    if(s03!="") out+=s03+"\n";
    if(s04!="") out+=s04+"\n";
    if(s05!="") out+=s05+"\n";
    if(s06!="") out+=s06+"\n";
    if(s07!="") out+=s07+"\n";
    if(s08!="") out+=s08+"\n";
    if(s09!="") out+=s09+"\n";
    if(s0A!="") out+=s0A+"\n";
    if(s0B!="") out+=s0B+"\n";
    if(s0C!="") out+=s0C+"\n";
    if(s0D!="") out+=s0D+"\n";
    if(s0E!="") out+=s0E+"\n";
    if(s0F!="") out+=s0F+"\n";
    if(s10!="") out+=s10+"\n";
    if(s11!="") out+=s11+"\n";
    if(s12!="") out+=s12+"\n";
    if(s13!="") out+=s13+"\n";
    if(s14!="") out+=s14+"\n";
    if(s15!="") out+=s15+"\n";
    if(s16!="") out+=s16+"\n";
    if(s17!="") out+=s17+"\n";
    if(s18!="") out+=s18+"\n";
    if(s19!="") out+=s19+"\n";
    if(s1A!="") out+=s1A+"\n";
    if(s1B!="") out+=s1B+"\n";
    if(s1C!="") out+=s1C+"\n";
    if(s1D!="") out+=s1D+"\n";
    if(s1E!="") out+=s1E+"\n";
    if(s1F!="") out+=s1F+"\n";
    if(s20!="") out+=s20+"\n";
    if(s21!="") out+=s21+"\n";
    if(s22!="") out+=s22+"\n";
    if(s23!="") out+=s23+"\n";
    if(s24!="") out+=s24+"\n";
    if(s25!="") out+=s25+"\n";
    if(s26!="") out+=s26+"\n";
    if(s27!="") out+=s27+"\n";
    if(s28!="") out+=s28+"\n";
    if(s29!="") out+=s29+"\n";
    if(s2A!="") out+=s2A+"\n";
    if(s2B!="") out+=s2B+"\n";
    if(s2C!="") out+=s2C+"\n";
    if(s2D!="") out+=s2D+"\n";
    if(s2E!="") out+=s2E+"\n";
    if(s2F!="") out+=s2F+"\n";
    return out;
  }

  // ***************************************************************************
  // Function stream2dequestring return UINT
  // ***************************************************************************
  // take istream into a deque strings - Stefano Curtarolo
  uint stream2dequestring(std::istream& istreamIN,deque<string> &vstringout) {
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!istreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      istreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout.push_back(string(tmp));
    }
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.size();  // return FALSE if something got messed up
  }
  // take ifstream into a deque strings - Stefano Curtarolo
  uint stream2dequestring(std::ifstream& ifstreamIN,deque<string> &vstringout) {
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!ifstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      ifstreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout.push_back(string(tmp));
    }
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.size();  // return FALSE if something got messed up
  }
  uint stream2dequestring(std::stringstream& stringstreamIN,deque<string> &vstringout) {
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!stringstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      stringstreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout.push_back(string(tmp));
    }
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.size();  // return FALSE if something got messed up
  }
  uint string2dequestring(const string& stringIN,deque<string> &vstringout) {
    return aurostd::string2tokens(stringIN,vstringout,"\n");
  }

  // ***************************************************************************
  // Function string2dequestring return DEQUE
  // ***************************************************************************
  // take sitring into a deque strings - Stefano Curtarolo
  deque<string> stream2dequestring(std::istream& istreamIN) {
    deque<string> vstringout;
    aurostd::stream2dequestring(istreamIN,vstringout);
    return vstringout;
  }
  deque<string> stream2dequestring(std::ifstream& iftreamIN) {
    deque<string> vstringout;
    aurostd::stream2dequestring(iftreamIN,vstringout);
    return vstringout;
  }
  deque<string> stream2dequestring(std::stringstream& stringstreamIN) {
    deque<string> vstringout;
    aurostd::stream2dequestring(stringstreamIN,vstringout);
    return vstringout;
  }
  deque<string> string2dequestring(const string& stringIN) {
    deque<string> vstringout;
    aurostd::string2dequestring(stringIN,vstringout);
    return vstringout;
  }

  // ***************************************************************************
  // Function string2file string2compressfile string2gzfile string2bz2file string2xzfile
  // ***************************************************************************
  // write string to file - Stefano Curtarolo
  bool string2file(const string& StringOUTPUT,const string& _FileNameOUTPUT,const string& mode) {
    string FileNameOUTPUT=aurostd::CleanFileName(_FileNameOUTPUT);
    bool writable=true;  //CO20190808 - captures whether we can open/write file
    if(mode=="POST" || mode=="APPEND") {
      stringstream FileINPUT;
      if(aurostd::FileExist(FileNameOUTPUT)){aurostd::file2stringstream(FileNameOUTPUT,FileINPUT);}
      ofstream FileOUTPUT;
      FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
      writable=FileOUTPUT.is_open(); //CO20190808 - captures whether we can open/write file
      FileOUTPUT << FileINPUT.str();
      FileOUTPUT << StringOUTPUT;
      FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
      return writable; //TRUE;  // return FALSE if something got messed up //CO20190808 - captures whether we can open/write file
    }
    if(mode=="PRE") {
      stringstream FileINPUT;
      if(aurostd::FileExist(FileNameOUTPUT)){aurostd::file2stringstream(FileNameOUTPUT,FileINPUT);}
      ofstream FileOUTPUT;
      FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
      writable=FileOUTPUT.is_open(); //CO20190808 - captures whether we can open/write file
      FileOUTPUT << StringOUTPUT;
      FileOUTPUT << FileINPUT.str();
      FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
      return writable; //TRUE;  // return FALSE if something got messed up //CO20190808 - captures whether we can open/write file
    }
    if(mode=="WRITE" || mode=="") {
      ofstream FileOUTPUT;
      FileOUTPUT.open(FileNameOUTPUT.c_str(),std::ios::out);
      writable=FileOUTPUT.is_open(); //CO20190808 - captures whether we can open/write file
      FileOUTPUT << StringOUTPUT;
      FileOUTPUT.flush();FileOUTPUT.clear();FileOUTPUT.close();
      return writable; //TRUE;  // return FALSE if something got messed up //CO20190808 - captures whether we can open/write file
    }
    return FALSE;
  }

  bool string2compressfile(const string& command,const string& StringOUTPUT,const string& _file,const string& mode) {
    // "" compliant SC20190401
    string file=aurostd::CleanFileName(_file);
    bool out=string2file(StringOUTPUT,file,mode);
    aurostd::execute(command+" -9fq \""+file+"\"");
    return out;
  }

  bool string2gzfile(const string& StringOUTPUT,const string& file,const string& mode) {
    return string2compressfile("gzip",StringOUTPUT,file,mode);
  }

  bool string2bz2file(const string& StringOUTPUT,const string& file,const string& mode) {
    return string2compressfile("bzip2",StringOUTPUT,file,mode);
  }

  bool string2xzfile(const string& StringOUTPUT,const string& file,const string& mode) {
    return string2compressfile("xz",StringOUTPUT,file,mode);
  }


  // ***************************************************************************
  // Function stringstream2file stringstream2compressedfile stringstream2gzfile stringstream2bz2file stringstream2xzfile
  // ***************************************************************************
  // write string to file - Stefano Curtarolo
  bool stringstream2file(const stringstream& StringstreamOUTPUT,const string& file,const string& mode) {return string2file(StringstreamOUTPUT.str(),file,mode);}  //CO20210315 - cleaned up

  bool stringstream2compressfile(const string& command,const stringstream& StringstreamOUTPUT,const string& _file,const string& mode) {
    string file=aurostd::CleanFileName(_file);
    bool out=stringstream2file(StringstreamOUTPUT,file,mode);
    aurostd::execute(command+" -9fq \""+file+"\"");
    return out;
  }

  bool stringstream2gzfile(const stringstream& StringstreamOUTPUT,const string& file,const string& mode) {
    return stringstream2compressfile("gzip",StringstreamOUTPUT,file,mode);
  }

  bool stringstream2bz2file(const stringstream& StringstreamOUTPUT,const string& file,const string& mode) {
    return stringstream2compressfile("bzip2",StringstreamOUTPUT,file,mode);
  }

  bool stringstream2xzfile(const stringstream& StringstreamOUTPUT,const string& file,const string& mode) {
    return stringstream2compressfile("xz",StringstreamOUTPUT,file,mode);
  }

  // ***************************************************************************
  // Function ostream2string  istream2string   
  // ***************************************************************************
  // convert ostream/istream to string - Stefano Curtarolo
  std::string ostream2string(std::ostream& oss) {
    std::stringstream soss;
    soss << oss.rdbuf();
    return soss.str();
  }

  // ***************************************************************************
  // Function stream2string
  // ***************************************************************************
  // take istream into a  strings - Stefano Curtarolo
  uint stream2string(std::istream& istreamIN,string &vstringout) {
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!istreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      istreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout+=string(tmp)+"\n";  // ME20210206 - fixed line break
    }
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.length();  // return FALSE if something got messed up
  }

  // take ifstream into a  strings - Stefano Curtarolo
  uint stream2string(std::ifstream& ifstreamIN,string &vstringout) {
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!ifstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      ifstreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout+=string(tmp)+"\n";  // ME20210206 - fixed line break
    }
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.length();  // return FALSE if something got messed up
  }

  uint stream2string(std::stringstream& stringstreamIN,string &vstringout) {
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while(!stringstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      stringstreamIN.getline(tmp,_CIN_LINE_BUFFER_LENGTH_-1);
      vstringout+=string(tmp)+"\n";  // ME20210206 - fixed line break
    }
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.length();  // return FALSE if something got messed up
  }

  // ***************************************************************************
  // Function getenv2string getenv2int getenv2uint getenv2double
  // ***************************************************************************
  // convert environments to string;
  string getenv2string(const string& str) {
    if(getenv(str.c_str())==NULL) return string("");
    return string(getenv(str.c_str()));
  }
  int getenv2int(const string& str) {
    if(getenv(str.c_str())==NULL) return int(0);
    return aurostd::string2utype<int>(getenv(str.c_str()));
  }
  uint getenv2uint(const string& str) {
    if(getenv(str.c_str())==NULL) return uint(0);
    return aurostd::string2utype<uint>(getenv(str.c_str()));
  }
  double getenv2double(const string& str) {
    if(getenv(str.c_str())==NULL) return double(0);
    return aurostd::string2utype<double>(getenv(str.c_str()));
  }

  // ***************************************************************************
  // Function file2string bz2file2string gzfile2string xzfile2string zipfile2string efile2string
  // ***************************************************************************
  // write file to string - Stefano Curtarolo
  uint file2string(const string& _FileNameIN,string& StringIN){
    return file2string_20220221(_FileNameIN, StringIN);
  }


  uint file2string_20220221(const string& _FileNameIN,string& StringIN){ //HE20220221
    // avoids the extra FileExist check (buffer.str() will always be empty if the file can't be open)
    // avoids reading the file char by char
    // speedup compared to file2string_20220101 10% (3000 files/ 5 warm runs)
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    std::ifstream open_file(FileNameIN);
    std::stringstream buffer;
    buffer << open_file.rdbuf();
    StringIN = buffer.str();
    return StringIN.length();
  }

  uint file2string_20220101(const string& _FileNameIN,string& StringIN) {  //CO20210624
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    if(!FileExist(FileNameIN)) {
      // cerr << "ERROR - aurostd::file2string: file=" << FileNameIN << " not present !" << endl;
      return 0;
    }
    ifstream FileIN;
    FileIN.open(FileNameIN.c_str(),std::ios::in);
    char c; while (FileIN.get(c)) StringIN+=c;
    FileIN.clear();FileIN.close();
    return StringIN.length();  // return 0 if something got messed up
  }
  uint bz2file2string(const string& _FileNameIN,string& StringIN) { //CO20210624
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    // cerr << "bz2file2string; BEGIN" << endl;
    if(!FileExist(FileNameIN)) {
      // cerr << "ERROR - aurostd::bz2file2string: file=" << FileNameIN << " not present !" << endl;
      return 0;
    }   
    if(!aurostd::IsCommandAvailable("bzcat")) {
      // cerr << "ERROR - aurostd::bz2file2string: command \"bzcat\" is necessary !" << endl;
      return 0;
    }   
    StringIN=aurostd::execute2string("bzcat \""+FileNameIN+"\"");
    return StringIN.length();  // return 0 if something got messed up
  }
  uint gzfile2string(const string& _FileNameIN,string& StringIN) {  //CO20210624
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    // cerr << "gzfile2string; BEGIN" << endl;
    if(!FileExist(FileNameIN)) {
      // cerr << "ERROR - aurostd::gzfile2string: file=" << FileNameIN << " not present !" << endl;
      return 0;
    }
    if(!aurostd::IsCommandAvailable("zcat")) {
      // cerr << "ERROR - aurostd::gzfile2string: command \"zcat\" is necessary !" << endl;
      return 0;
    }
    StringIN=aurostd::execute2string("zcat \""+FileNameIN+"\"");
    return StringIN.length();  // return 0 if something got messed up
  }
  uint xzfile2string(const string& _FileNameIN,string& StringIN) {  //CO20210624
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    // cerr << "xzfile2string; BEGIN" << endl;
    if(!FileExist(FileNameIN)) {
      // cerr << "ERROR - aurostd::xzfile2string: file=" << FileNameIN << " not present !" << endl;
      return 0;
    }
    if(!aurostd::IsCommandAvailable("xzcat")) {
      // cerr << "ERROR - aurostd::xzfile2string: command \"xzcat\" is necessary !" << endl;
      return 0;
    }
    StringIN=aurostd::execute2string("xzcat \""+FileNameIN+"\"");
    return StringIN.length();  // return 0 if something got messed up
  }
  //CO START
  uint zipfile2string(const string& _FileNameIN,string& StringIN) { //CO20210624
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    if(!FileExist(FileNameIN)) {
      return 0;
    }
    if(!aurostd::IsCommandAvailable("unzip")) {
      return 0;
    }
    StringIN=aurostd::execute2string("unzip -p \""+FileNameIN+"\"");
    return StringIN.length();  // return 0 if something got messed up
  }
  //CO END
  uint efile2string(const string& FileNameIN,string& StringIN) {  //CO20210624
    //[CO20190808 - OBSOLETE, we clean inside FileExist()]string FileNameIN=aurostd::CleanFileName(_FileNameIN),FileNameOUT;  //CO
    string FileNameOUT=""; //CO20191110
    // cerr << "efile2string; BEGIN FileNameIN=[" << FileNameIN << "]" << endl;
    if(!FileExist(FileNameIN,FileNameOUT) && !EFileExist(FileNameIN,FileNameOUT)) {
      // cerr << "ERROR - aurostd::efile2string: file=" << FileNameIN << " not present !" << endl;
      return 0;
    }   
    // cerr << aurostd::substring2bool(FileNameIN,".bz2") << endl;
    if(aurostd::substring2bool(FileNameOUT,".bz2")) {
      // cerr << "efile2string: found .bz2, using bz2file2string" << endl;
      return aurostd::bz2file2string(FileNameOUT,StringIN);
    }
    if(aurostd::substring2bool(FileNameOUT,".gz"))  {
      // cerr << "efile2string: found .gz, using gzfile2string" << endl;
      return aurostd::gzfile2string(FileNameOUT,StringIN);
    }
    if(aurostd::substring2bool(FileNameOUT,".xz"))  {
      // cerr << "efile2string: found .xz, using xzfile2string" << endl;
      return aurostd::xzfile2string(FileNameOUT,StringIN);
    }
    // cerr << "efile2string: file2string" << endl;
    return aurostd::file2string(FileNameOUT,StringIN);
  }

  // ***************************************************************************
  // Function file2vectorstring bz2file2vectorstring gzfile2vectorstring xzfile2vectorstring efile2vectorstring 
  // ***************************************************************************
  // write file to vector string - Stefano Curtarolo
  uint file2vectorstring(const string& FileNameIN,vector<string>& vline,bool consecutive,bool trim_edges) { //CO20210624
    return aurostd::string2vectorstring(file2string(aurostd::CleanFileName(FileNameIN)),vline,consecutive,trim_edges);
  }

  uint bz2file2vectorstring(const string& FileNameIN,vector<string>& vline,bool consecutive,bool trim_edges) {  //CO20210624
    return aurostd::string2vectorstring(bz2file2string(aurostd::CleanFileName(FileNameIN)),vline,consecutive,trim_edges);
  }

  uint gzfile2vectorstring(const string& FileNameIN,vector<string>& vline,bool consecutive,bool trim_edges) { //CO20210624
    return aurostd::string2vectorstring(gzfile2string(aurostd::CleanFileName(FileNameIN)),vline,consecutive,trim_edges);
  }

  uint xzfile2vectorstring(const string& FileNameIN,vector<string>& vline,bool consecutive,bool trim_edges) { //CO20210624
    return aurostd::string2vectorstring(xzfile2string(aurostd::CleanFileName(FileNameIN)),vline,consecutive,trim_edges);
  }

  uint efile2vectorstring(const string& FileNameIN,vector<string>& vline,bool consecutive,bool trim_edges) {  //CO20210624
    return aurostd::string2vectorstring(efile2string(aurostd::CleanFileName(FileNameIN)),vline,consecutive,trim_edges);
  }

  bool vectorstring2file(const vector<string>& vline,string FileNameOUT) {
    string file=aurostd::CleanFileName(FileNameOUT);
    ofstream FileOUT;
    FileOUT.open(file.c_str(),std::ios::out);
    bool writable=FileOUT.is_open(); //CO20190808 - captures whether we can open/write file
    for(uint iline=0;iline<vline.size();iline++) FileOUT << vline.at(iline) << endl;
    // FileOUT << StringstreamOUT.rdbuf();
    FileOUT.flush();FileOUT.clear();FileOUT.close();
    return writable;
  }

  // ***************************************************************************
  // Function file2dequestring bz2file2dequestring gzfile2dequestring xzfile2dequestring efile2dequestring 
  // ***************************************************************************
  // write file to deque string - Stefano Curtarolo
  uint file2dequestring(const string& FileNameIN,deque<string>& vline) {  //CO20210624
    return aurostd::string2dequestring(file2string(aurostd::CleanFileName(FileNameIN)),vline);
  }

  uint bz2file2dequestring(const string& FileNameIN,deque<string>& vline) { //CO20210624
    return aurostd::string2dequestring(bz2file2string(aurostd::CleanFileName(FileNameIN)),vline);
  }

  uint gzfile2dequestring(const string& FileNameIN,deque<string>& vline) {  //CO20210624
    return aurostd::string2dequestring(gzfile2string(aurostd::CleanFileName(FileNameIN)),vline);
  }

  uint xzfile2dequestring(const string& FileNameIN,deque<string>& vline) {  //CO20210624
    return aurostd::string2dequestring(xzfile2string(aurostd::CleanFileName(FileNameIN)),vline);
  }

  uint efile2dequestring(const string& FileNameIN,deque<string>& vline) { //CO20210624
    return aurostd::string2dequestring(efile2string(aurostd::CleanFileName(FileNameIN)),vline);
  }

  bool dequestring2file(const deque<string>& vline,string FileNameOUT) {
    string file=aurostd::CleanFileName(FileNameOUT);
    ofstream FileOUT;
    FileOUT.open(file.c_str(),std::ios::out);
    bool writable=FileOUT.is_open(); //CO20190808 - captures whether we can open/write file
    for(uint iline=0;iline<vline.size();iline++)  FileOUT << vline.at(iline) << endl;
    // FileOUT << StringstreamOUT.rdbuf();
    FileOUT.flush();FileOUT.clear();FileOUT.close();
    return writable;
  }


  // ***************************************************************************
  // Function file2vectorstring bz2file2vectorstring gzfile2vectorstring xzfile2vectorstring efile2vectorstring overloading for file2vector
  // ***************************************************************************
  // write file to deque string - Stefano Curtarolo
  uint file2vectorstring(const string& FileNameIN,deque<string>& vline) {  //CO20210624
    return aurostd::string2dequestring(file2string(aurostd::CleanFileName(FileNameIN)),vline);
  }

  uint bz2file2vectorstring(const string& FileNameIN,deque<string>& vline) {  //CO20210624
    return aurostd::string2dequestring(bz2file2string(aurostd::CleanFileName(FileNameIN)),vline);
  }

  uint gzfile2vectorstring(const string& FileNameIN,deque<string>& vline) { //CO20210624
    return aurostd::string2dequestring(gzfile2string(aurostd::CleanFileName(FileNameIN)),vline);
  }

  uint xzfile2vectorstring(const string& FileNameIN,deque<string>& vline) { //CO20210624
    return aurostd::string2dequestring(xzfile2string(aurostd::CleanFileName(FileNameIN)),vline);
  }

  uint efile2vectorstring(const string& FileNameIN,deque<string>& vline) {  //CO20210624
    return aurostd::string2dequestring(efile2string(aurostd::CleanFileName(FileNameIN)),vline);
  }

  bool vectorstring2file(const deque<string>& vline,string FileNameOUT) {
    string file=aurostd::CleanFileName(FileNameOUT);
    ofstream FileOUT;
    FileOUT.open(file.c_str(),std::ios::out);
    bool writable=FileOUT.is_open(); //CO20190808 - captures whether we can open/write file
    for(uint iline=0;iline<vline.size();iline++) FileOUT << vline.at(iline) << endl;
    // FileOUT << StringstreamOUT.rdbuf();
    FileOUT.flush();FileOUT.clear();FileOUT.close();
    return writable;
  }

  // ***************************************************************************
  // Function file2stringstream bz2file2stringstream gzfile2stringstream xzfile2stringstream zipfile2stringstream efile2stringstream
  // ***************************************************************************
  // write file to stringstream - Stefano Curtarolo
  bool file2stringstream(const string& _FileNameIN,stringstream& StringstreamIN) {  //CO20210624
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    if(!FileExist(FileNameIN)) {
      cerr << "ERROR - aurostd::file2stringstream: file=" << FileNameIN << " not present !" << endl;
      return FALSE;
    }
    aurostd::StringstreamClean(StringstreamIN);
    StringstreamIN << file2string(FileNameIN);
    return TRUE;  // return FALSE if something got messed up
  }
  bool bz2file2stringstream(const string& _FileNameIN,stringstream& StringstreamIN) { //CO20210624
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    if(!FileExist(FileNameIN)) {
      cerr << "ERROR - aurostd::bz2file2stringstream: file=" << FileNameIN << " not present !" << endl;
      return FALSE;
    }
    aurostd::StringstreamClean(StringstreamIN);
    StringstreamIN << bz2file2string(FileNameIN);
    return TRUE;  // return FALSE if something got messed up
  }
  bool gzfile2stringstream(const string& _FileNameIN,stringstream& StringstreamIN) {  //CO20210624
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    if(!FileExist(FileNameIN)) {
      cerr << "ERROR - aurostd::gzfile2stringstream: file=" << FileNameIN << " not present !" << endl;
      return FALSE;
    }
    aurostd::StringstreamClean(StringstreamIN);
    StringstreamIN << gzfile2string(FileNameIN);
    return TRUE;  // return FALSE if something got messed up
  }
  bool xzfile2stringstream(const string& _FileNameIN,stringstream& StringstreamIN) {  //CO20210624
    string FileNameIN=aurostd::CleanFileName(_FileNameIN);
    if(!FileExist(FileNameIN)) {
      cerr << "ERROR - aurostd::xzfile2stringstream: file=" << FileNameIN << " not present !" << endl;
      return FALSE;
    }
    aurostd::StringstreamClean(StringstreamIN);
    StringstreamIN << xzfile2string(FileNameIN);
    return TRUE;  // return FALSE if something got messed up
  }
  //CO START
  //zipfile 2 a string
  bool zipfile2stringstream(const string& _FileNameIN,stringstream& StringstreamIN) { //CO20210624
    string FileNameIN=CleanFileName(_FileNameIN);
    if(!FileExist(FileNameIN)) {
      cerr << "ERROR - aurostd::zipfile2stringstream: file=" << FileNameIN << " not present !" << endl;
      return FALSE;
    }
    aurostd::StringstreamClean(StringstreamIN);
    StringstreamIN << zipfile2string(FileNameIN);
    return TRUE;  // return FALSE if something got messed up
  }
  //CO END
  bool efile2stringstream(const string& _FileNameIN,stringstream& StringstreamIN) { //CO20210624
    string FileNameIN=aurostd::CleanFileName(_FileNameIN),FileNameOUT;  //CO
    if(!FileExist(FileNameIN,FileNameOUT) && !EFileExist(FileNameIN,FileNameOUT)) {
      cerr << "ERROR - aurostd::efile2stringstream: file=" << FileNameIN << " not present !" << endl;
      return FALSE;
    }
    //ME20200922 - Do not use substring2bool - it may delete paths if not properly cleaned (e.g. if they contain //)
    if(FileNameOUT.find(".bz2") != string::npos) return aurostd::bz2file2stringstream(FileNameOUT,StringstreamIN);
    if(FileNameOUT.find(".gz") != string::npos)  return aurostd::gzfile2stringstream(FileNameOUT,StringstreamIN);
    if(FileNameOUT.find(".xz") != string::npos)  return aurostd::xzfile2stringstream(FileNameOUT,StringstreamIN);
    return aurostd::file2stringstream(FileNameOUT,StringstreamIN);
  }

  // ***************************************************************************
  // Function file2string
  // ***************************************************************************
  // write file to string - Stefano Curtarolo
  string file2string(const string& FileNameIN) {  //CO20210624
    string StringIN="";
    file2string(FileNameIN,StringIN);
    return StringIN;
  }
  string bz2file2string(const string& FileNameIN) { //CO20210624
    string StringIN="";
    bz2file2string(FileNameIN,StringIN);
    return StringIN;
  }
  string gzfile2string(const string& FileNameIN) {  //CO20210624
    string StringIN="";
    gzfile2string(FileNameIN,StringIN);
    return StringIN;
  }
  string xzfile2string(const string& FileNameIN) {  //CO20210624
    string StringIN="";
    xzfile2string(FileNameIN,StringIN);
    return StringIN;
  }
  string zipfile2string(const string& FileNameIN) {  //CO20210624
    string StringIN="";
    zipfile2string(FileNameIN,StringIN);
    return StringIN;
  }
  string efile2string(const string& FileNameIN) { //CO20210624
    string StringIN="";
    efile2string(FileNameIN,StringIN);
    return StringIN;
  }

  // ***************************************************************************
  // Function url2file
  // ***************************************************************************
  // wget URL to string - Stefano Curtarolo
  bool url2file(string url,string& fileIN,bool verbose) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(!aurostd::IsCommandAvailable("wget")) {
      cerr << "ERROR - aurostd::url2file(): command \"wget\" is necessary !" << endl;
      return FALSE;
    }
    string _url=url;
    aurostd::StringSubst(_url,"http://","");
    aurostd::StringSubst(_url,"//","/");
    if(LDEBUG) cerr << "aurostd::url2file(): Loading url=" << _url << endl;
    if(verbose) cout << "aurostd::url2file(): Loading url=" << _url << endl;
#ifndef _MACOSX_
    aurostd::execute("wget --quiet --no-cache -O "+fileIN+" http://"+_url);
#else
    aurostd::execute("wget --quiet -O "+fileIN+" http://"+_url);
#endif    
    if(aurostd::FileEmpty(fileIN)) {
      aurostd::StringSubst(_url,":AFLOW","/AFLOW");
#ifndef _MACOSX_
      aurostd::execute("wget --quiet --no-cache -O "+fileIN+" http://"+_url);
#else
      aurostd::execute("wget --quiet -O "+fileIN+" http://"+_url); // _MACOSX_
#endif    
      if(aurostd::FileEmpty(fileIN)) {
        if(LDEBUG){cerr << "ERROR - aurostd::url2file(): URL not found http://" << _url << endl;} //CO20200731 - silence this, it's not an error
        return FALSE;
      }
    }
    return TRUE;
  }

  // ***************************************************************************
  // Function eurl2string
  // ***************************************************************************
  //CO20200223
  bool eurl2string(const string& url,string& stringIN,bool verbose) {
    stringIN="";
    string ext=GetCompressionExtension(url);
    if(!ext.empty()){
      string temp_file=aurostd::TmpFileCreate("eurl2string")+ext;
      url2file(url,temp_file,verbose);
      efile2string(temp_file,stringIN);
#ifndef _AFLOW_TEMP_PRESERVE_
      aurostd::RemoveFile(temp_file);
#endif
      return stringIN.length()>0;
    }
    //[CO20200404 - OBSOLETE]if(aurostd::substring2bool(url,".bz2")){
    //[CO20200404 - OBSOLETE]  string temp_file=aurostd::TmpFileCreate("eurl2string")+".bz2";
    //[CO20200404 - OBSOLETE]  url2file(url,temp_file,verbose);
    //[CO20200404 - OBSOLETE]  bz2file2string(temp_file,stringIN);
    //[CO20200404 - OBSOLETE]  return stringIN.length()>0;
    //[CO20200404 - OBSOLETE]}
    //[CO20200404 - OBSOLETE]if(aurostd::substring2bool(url,".gz")){
    //[CO20200404 - OBSOLETE]  string temp_file=aurostd::TmpFileCreate("eurl2string")+".gz";
    //[CO20200404 - OBSOLETE]  url2file(url,temp_file,verbose);
    //[CO20200404 - OBSOLETE]  gzfile2string(temp_file,stringIN);
    //[CO20200404 - OBSOLETE]  return stringIN.length()>0;
    //[CO20200404 - OBSOLETE]}
    //[CO20200404 - OBSOLETE]if(aurostd::substring2bool(url,".xz")){
    //[CO20200404 - OBSOLETE]  string temp_file=aurostd::TmpFileCreate("eurl2string")+".xz";
    //[CO20200404 - OBSOLETE]  url2file(url,temp_file,verbose);
    //[CO20200404 - OBSOLETE]  xzfile2string(temp_file,stringIN);
    //[CO20200404 - OBSOLETE]  return stringIN.length()>0;
    //[CO20200404 - OBSOLETE]}
    url2string(url,stringIN,verbose);
    return stringIN.length()>0;
  }

  // ***************************************************************************
  // Function url2string
  // ***************************************************************************
  // wget URL to string - Stefano Curtarolo
  // HE20220615 changed to aurostd http code to reduce system calls
  // Use the function in aurostd_xhttp for new code
  // Old use-cases of this function should be replaced, kept for now to maintain compatibility
  // As the original this function does not support SSL (https)
  bool url2string(const string& url,string& stringIN,bool verbose) {
    if (verbose) cerr << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    int return_code = aurostd::httpGetStatus(url, stringIN);
    if (verbose) cerr << __AFLOW_FUNC__ << " " << url << " returned " << return_code <<endl;
    if(stringIN.empty()) return false;
    else return true;
  }

  // ***************************************************************************
  // Function eurl2stringstream
  // ***************************************************************************
  //CO20200223
  bool eurl2stringstream(const string& url,stringstream& stringstreamIN,bool verbose) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    if(verbose) cout << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    string stringIN="";  //CO20200404
    bool out=eurl2string(url,stringIN,verbose);
    aurostd::StringstreamClean(stringstreamIN); stringstreamIN << stringIN;
    return out;
  }

  // ***************************************************************************
  // Function url2stringstream
  // ***************************************************************************
  // wget URL to stringstream - Stefano Curtarolo
  bool url2stringstream(const string& url,stringstream& stringstreamIN,bool verbose) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    if(verbose) cout << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    string stringIN="";  //CO20200404
    bool out=url2string(url,stringIN,verbose);
    aurostd::StringstreamClean(stringstreamIN); stringstreamIN << stringIN;
    return out;
  }

  // ***************************************************************************
  // Function eurl2vectorstring
  // ***************************************************************************
  //CO20200223
  bool eurl2vectorstring(const string& url,vector<string>& vlines,bool verbose) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    if(verbose) cout << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    string stringIN="";  //CO20200404
    bool out=eurl2string(url,stringIN,verbose);
    aurostd::string2tokens(stringIN,vlines);
    return out;
  }
  // ***************************************************************************
  // Function url2vectorstring
  // ***************************************************************************
  // wget URL to vectorstring - Stefano Curtarolo
  bool url2vectorstring(const string& url,vector<string>& vlines,bool verbose) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    if(verbose) cout << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    string stringIN="";  //CO20200404
    bool out=url2string(url,stringIN,verbose);
    aurostd::string2tokens(stringIN,vlines);
    return out;
  }

  // ***************************************************************************
  // Function eurl2dequestring
  // ***************************************************************************
  //CO20200223
  bool eurl2dequestring(const string& url,deque<string>& vlines,bool verbose) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    if(verbose) cout << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    string stringIN="";  //CO20200404
    bool out=eurl2string(url,stringIN,verbose);
    aurostd::string2tokens(stringIN,vlines);
    return out;
  }

  // ***************************************************************************
  // Function url2dequestring
  // ***************************************************************************
  // wget URL to dequestring - Stefano Curtarolo
  bool url2dequestring(const string& url,deque<string>& vlines,bool verbose) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    if(verbose) cout << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    string stringIN="";  //CO20200404
    bool out=url2string(url,stringIN,verbose);
    aurostd::string2tokens(stringIN,vlines);
    return out;
  }

  // ***************************************************************************
  // Function eurl2tokens
  // ***************************************************************************
  //CO20200223
  template<typename utype> uint eurl2tokens(const string& url,vector<utype>& tokens,const string& delimiters) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    if(!aurostd::IsCommandAvailable("wget")) {
      cerr << "ERROR - " << __AFLOW_FUNC__ << ": command \"wget\" is necessary !" << endl;
      return 0;}	
    tokens.clear(); 
    string content = "";
    aurostd::eurl2string(url,content);
    if(LDEBUG) { //CO20180627
      cerr << __AFLOW_FUNC__ << " content=" << endl;
      cerr << content << endl;
    }
    if(content.empty()) {
      if(LDEBUG){cerr << "ERROR - " << __AFLOW_FUNC__ << ": URL empty http://" << url << endl;} //CO20200731 - silence this, it's not an error
      return 0;
    }

    vector<string> stokens;
    aurostd::string2tokens(content,stokens,delimiters);
    for(uint i=0;i<stokens.size();i++)
      if(!stokens[i].empty()) 
        tokens.push_back(aurostd::string2utype<utype>(stokens[i]));
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " [5] tokens.size()=" << tokens.size() << endl;
    return tokens.size();
  }

  // ***************************************************************************
  // Function url2tokens
  // ***************************************************************************
  // wget URL to vector of tokens - Stefano Curtarolo
  template<typename utype> uint url2tokens(const string& url,vector<utype>& tokens,const string& delimiters) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    if(!aurostd::IsCommandAvailable("wget")) {
      cerr << "ERROR - " << __AFLOW_FUNC__ << ": command \"wget\" is necessary !" << endl;
      return 0;}	
    tokens.clear(); 
    string content = "";
    aurostd::url2string(url,content);
    if(LDEBUG) { //CO20180627
      cerr << __AFLOW_FUNC__ << " content=" << endl;
      cerr << content << endl;
    }
    if(content.empty()) {
      if(LDEBUG){cerr << "ERROR - " << __AFLOW_FUNC__ << ": URL empty http://" << url << endl;} //CO20200731 - silence this, it's not an error
      return 0;
    }

    vector<string> stokens;
    aurostd::string2tokens(content,stokens,delimiters);
    for(uint i=0;i<stokens.size();i++)
      if(stokens[i]!="") 
        tokens.push_back(aurostd::string2utype<utype>(stokens[i]));
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " [5] tokens.size()=" << tokens.size() << endl;
    return tokens.size();
  }

  // ***************************************************************************
  // Function eurl2tokens
  // ***************************************************************************
  //CO20200223
  template<typename utype> uint eurl2tokens(const string& url,deque<utype>& tokens,const string& delimiters) {
    vector<utype> vtokens;
    aurostd::eurl2tokens(url,vtokens,delimiters);
    for(uint i=0;i<vtokens.size();i++) tokens.push_back(vtokens[i]);
    return tokens.size();
  }

  // ***************************************************************************
  // Function url2tokens
  // ***************************************************************************
  // wget URL to deque of tokens - Stefano Curtarolo
  template<typename utype> uint url2tokens(const string& url,deque<utype>& tokens,const string& delimiters) {
    vector<utype> vtokens;
    aurostd::url2tokens(url,vtokens,delimiters);
    for(uint i=0;i<vtokens.size();i++) tokens.push_back(vtokens[i]);
    return tokens.size();
  }

  // ***************************************************************************
  // Function eurl2string
  // ***************************************************************************
  //CO20200223
  string eurl2string(const string& url) {
    string stringIN="";    //CO20200404
    eurl2string(url,stringIN);
    return stringIN;
  }

  // ***************************************************************************
  // Function url2string
  // ***************************************************************************
  // wget URL to stringstream - Stefano Curtarolo
  string url2string(const string& url) {
    string stringIN="";  //CO20200404
    url2string(url,stringIN);
    return stringIN;
  }

  // ***************************************************************************
  // Function ChmodFile
  // ***************************************************************************
  // change mod of a file - Stefano Curtarolo
  bool ChmodFile(string chmod_string,string _file) { // "" compliant SC20190401
    string file=aurostd::CleanFileName(_file);
    if(chmod_string.empty()){return false;} //CO20190321
    if(file.empty()){return false;} //CO20190321
    ostringstream aus;
    aurostd::StringstreamClean(aus);
    aus << CHMOD_BIN << " " << chmod_string << " \"" << file << "\"" << endl;
    aurostd::execute(aus);
    return TRUE;  // return FALSE if something got messed up
  }

  //CO START
  //***************************************************************************//
  // aurostd::file2directory
  //***************************************************************************//
  // Corey Oses
  // Move file to destination
  bool file2directory(const string& _file,const string& _destination) { // "" compliant SC20190401
    string file=CleanFileName(_file);
    string destination=CleanFileName(_destination);
    if(!aurostd::IsCommandAvailable("mv")) {
      cerr << "ERROR - aurostd::file2directory: command \"mv\" is necessary !" << endl;
      return FALSE;
    }
    if(!aurostd::FileExist(file)) {
      cerr << "ERROR - aurostd::file2directory: " << file << " cannot be found !" << endl;
      //      cerr << "ERROR - aurostd::file2directory: _file=" << _file << " " << endl;
      //      cerr << "ERROR - aurostd::file2directory:  file=" << file << " " << endl;
      //      cerr << "ERROR - aurostd::file2directory: _destination=" << _destination << " " << endl;
      //     cerr << "ERROR - aurostd::file2directory:  destination=" << destination << " " << endl;
      return FALSE;
    }
    if(!aurostd::IsDirectory(destination)) //CO20180220 //FileExist(destination))
    { //CO20200106 - patching for auto-indenting
      cerr << "ERROR - aurostd::file2directory: " << destination << " cannot be found !" << endl;
      return FALSE;
    }
    aurostd::execute("mv \""+file+"\" \""+destination+"/\"");
    return TRUE;
  }

  //***************************************************************************//
  // aurostd::file2directory
  //***************************************************************************//
  // Corey Oses
  // Move files to destination, do one at a time to check if files exist
  bool file2directory(const vector<string>& files,const string& destination) {
    for(uint i=0;i<files.size();i++) {
      if(!file2directory(files[i],destination)) {return FALSE;}
    }
    return TRUE;
  }

  //***************************************************************************//
  // aurostd::file2file
  //***************************************************************************//
  // Corey Oses
  bool file2file(const string& _file,const string& _destination) {
    string file=CleanFileName(_file);
    string destination=CleanFileName(_destination);
    if(!aurostd::IsCommandAvailable("mv")) {
      cerr << "ERROR - aurostd::file2file: command \"mv\" is necessary !" << endl;
      return FALSE;
    }
    if(!aurostd::FileExist(file)) {
      cerr << "ERROR - aurostd::file2file: " << file << " cannot be found !" << endl;
      return FALSE;
    }
    aurostd::execute("mv \""+file+"\" \""+destination+"\"");
    return TRUE;
  }

  //***************************************************************************//
  // aurostd::file2md5sum
  //***************************************************************************//
  // Stefano Curtarolo
  string file2md5sum(const string& file) { //SC20200326
    vector<string> vtokens;
    if(aurostd::FileExist(file)) {
      aurostd::string2tokens(aurostd::execute2string("md5sum "+file),vtokens," ");
      if(vtokens.size()>0) return vtokens.at(0);
    }
    return "";
  }

  //***************************************************************************//
  // aurostd::file2auid
  //***************************************************************************//
  // Stefano Curtarolo
  string file2auid(const string& file) { //SC20200326
    vector<string> vtokens;
    if(aurostd::FileExist(file)) {
      uint64_t crc=0;
      crc=aurostd::crc64(crc,aurostd::efile2string(file)); // DONT TOUCH THIS
      return aurostd::crc2string(crc);
    }
    return "";
  }

  // ***************************************************************************
  // Function IsDirectory
  // ***************************************************************************
  // true if path is directory
  // http://stackoverflow.com/questions/146924/how-can-i-tell-if-a-given-path-is-a-directory-or-a-file-c-c
  bool IsDirectory(const string& _path){
    string path=CleanFileName(_path);
    struct stat s;
    if( stat(path.c_str(),&s) != 0 ){return FALSE;} //error
    if( s.st_mode & S_IFDIR ){      //CO20200531 - still trips if directory exists but is not accessible
      DIR *dp;                      //CO20200531 - see if directory is accessible by user
      dp=opendir(path.c_str());     //CO20200531 - see if directory is accessible by user
      if(dp==NULL) {return FALSE;}  //CO20200531 - see if directory is accessible by user
      (void) closedir (dp);         //CO20200531 - see if directory is accessible by user
      return TRUE;
    }
    return FALSE;
  }

  // ***************************************************************************
  // Function IsFile
  // ***************************************************************************
  // true if path is file
  // http://stackoverflow.com/questions/146924/how-can-i-tell-if-a-given-path-is-a-directory-or-a-file-c-c
  bool IsFile(const string& _path){
    string path=CleanFileName(_path);
    struct stat s;
    if( stat(path.c_str(),&s) != 0 ){return FALSE;} //error
    if( s.st_mode & S_IFREG ){      //CO20200531 - still trips if file exists but is not accessible
      ifstream f(path.c_str());     //CO20200531 - see if file is accessible by user
      if(!f.good()){return FALSE;}  //CO20200531 - see if file is accessible by user
      return TRUE;
    }
    return FALSE;
  }
  //CO END

  // ***************************************************************************
  // Function RemoveFile
  // ***************************************************************************
  // change remove a file - Stefano Curtarolo
  bool RemoveFile(string _file) {
    string file=CleanFileName(_file);
    // [OBSOLETE]   ostringstream aus;
    // [OBSOLETE] aurostd::StringstreamClean(aus);
    // [OBSOLETE] //    aus << "rm -f \"" << file << "\"" << endl; //   cerr << aus.str() << endl;
    // [OBSOLETE] //    aurostd::execute(aus);
    // [OBSOLETE]    remove(file.c_str()); // this works only with well defined names and not with * which are build by sh/bash
    execute("rm -f \""+file+"\""); // the simplest possible call th sh/bash
    return TRUE;  // return FALSE if something got messed up
  }

  //CO START
  //***************************************************************************//
  // aurostd::RemoveFile
  //***************************************************************************//
  bool RemoveFile(const vector<string>& files) {
    for(uint i=0;i<files.size();i++) {
      if(!RemoveFile(files[i])) {return FALSE;}
    }
    return TRUE;
  }

  //***************************************************************************//
  // aurostd::RemoveDirectory
  //***************************************************************************//
  bool RemoveDirectory(const string& _path) {
    string path=CleanFileName(_path);
    if(!IsCommandAvailable("rm")) {
      cerr << "ERROR - aurostd::RemoveDirectory: command \"rm\" is necessary !" << endl;
      return FALSE;
    }
    aurostd::execute("rm -r \""+path+"\"");
    return TRUE;
  }
  //CO END

//#define DEBUG_STRING2TOKENS
  // ***************************************************************************
  // Function string2tokens string2tokens<utype>
  // ***************************************************************************
  // Finds string2tokens to split strings in tokens
  // Stefano Curtarolo
  // void string2tokens(const string& str,vector<string>& tokens,const string& delimiters=" ") {  //[CO20200106 - close bracket for indenting]}
  uint string2tokens(const string& str,std::vector<string>& tokens,const string& delimiters,bool consecutive) { //CO20170613
    //CO mods 20170613
    //we are adding functionality here, because string2tokens will treat "\n\n" same as "\n", but not "\n \n"
    //consecutive will do the following: "sssss" -> <"s","s",...>
    //consecutive ALSO starts at 0, not first_not_of
    //return aurostd::string2tokens(stringIN,vstringout,"\n",true);
    tokens.clear(); // clear in the case there was something already in!
    string::size_type lastPos=(consecutive ? 0 : str.find_first_not_of(delimiters,0));   // Skip delimiters at beginning.
    string::size_type pos=str.find_first_of(delimiters,lastPos);     // Find first "non-delimiter".
#ifdef DEBUG_STRING2TOKENS
    cerr << __AFLOW_FUNC__ << " delimiters=" << delimiters << endl;
    cerr << __AFLOW_FUNC__ << " consecutive=" << consecutive << endl;
    cerr << __AFLOW_FUNC__ << " lastPos=" << lastPos << endl;
    cerr << __AFLOW_FUNC__ << " pos=" << pos << endl;
#endif
    while (pos!=string::npos || lastPos!=string::npos) {
      tokens.push_back(str.substr(lastPos,pos-lastPos));             // Found a token, add it to the vector.
#ifdef DEBUG_STRING2TOKENS
      cerr << __AFLOW_FUNC__ << " tokens.back()=\"" << tokens.back() << "\"" << endl;
#endif
      if(consecutive){lastPos=(pos!=string::npos ? pos+1 : string::npos);}
      else{lastPos=str.find_first_not_of(delimiters,pos);}           // Skip delimiters.  Note the "not_of"
      pos=str.find_first_of(delimiters,lastPos);                     // Find next "non-delimiter"
#ifdef DEBUG_STRING2TOKENS
      cerr << __AFLOW_FUNC__ << " lastPos=" << lastPos << endl;
      cerr << __AFLOW_FUNC__ << " pos=" << pos << endl;
#endif
    }
    return tokens.size();
  }

  // void string2tokens(const string& str,deque<string>& tokens,const string& delimiters=" ") { //[CO20200106 - close bracket for indenting]}
  uint string2tokens(const string& str,std::deque<string>& tokens,const string& delimiters,bool consecutive) { //CO20170613
    vector<string> vtokens;
    uint i=aurostd::string2tokens(str,vtokens,delimiters,consecutive);  //CO20170613
    tokens.clear();
    for(i=0;i<vtokens.size();i++) tokens.push_back(vtokens[i]);
    return tokens.size();
  }
  template<class utype> uint string2tokens(const string& str,std::vector<utype>& tokens,const string& delimiters,bool consecutive) {  //CO20170613
    vector<string> stokens;
    uint out=aurostd::string2tokens(str,stokens, delimiters,consecutive); //CO20170613
    tokens.clear();
    for(uint i=0;i<stokens.size();i++)
      tokens.push_back(aurostd::string2utype<utype>(stokens[i]));
    return out;
  }
  template<class utype> uint string2tokens(const string& str,std::deque<utype>& tokens,const string& delimiters,bool consecutive) { //CO20170613
    deque<string> stokens;
    uint out=aurostd::string2tokens(str,stokens, delimiters,consecutive); //CO20170613
    tokens.clear();
    for(uint i=0;i<stokens.size();i++)
      tokens.push_back(aurostd::string2utype<utype>(stokens[i]));
    return out;
  }

  //SD20220504 - string2tokens, but using a single delimiter that can have more than one char. Makes use of "find" rather than "find_first_of"
  uint string2tokensByDelimiter(const string& str, vector<string>& tokens, const string& delimiter) { //SD20220504
    tokens.clear();
    uint dlen = delimiter.length();
    if (dlen == 1) {return aurostd::string2tokens(str, tokens, delimiter);}
    string::size_type lpos = 0, cpos = str.find(delimiter, lpos);
    while (cpos != string::npos) {
      tokens.push_back(str.substr(lpos, cpos - lpos));
      lpos = cpos + dlen;
      cpos = str.find(delimiter, lpos);
    }
    tokens.push_back(str.substr(lpos, str.length() - lpos));
    return tokens.size();
  }
  uint string2tokensByDelimiter(const string& str, deque<string>& tokens, const string& delimiter) { //SD20220504
    tokens.clear();
    vector<string> vtokens;
    aurostd::string2tokensByDelimiter(str, vtokens, delimiter);
    tokens = aurostd::vector2deque(vtokens);
    return tokens.size();
  }

  // ***************************************************************************
  // Function string2tokensAdd string2tokensAdd<utype>
  // ***************************************************************************

  uint string2tokensAdd(const string& str,std::vector<string>& tokens,const string& delimiters) {
    vector<string> vtokens;
    uint i=aurostd::string2tokens(str,vtokens,delimiters);
    for(i=0;i<vtokens.size();i++) tokens.push_back(vtokens[i]);
    return tokens.size();
  }
  uint string2tokensAdd(const string& str,std::deque<string>& tokens,const string& delimiters) {
    vector<string> vtokens;
    uint i=aurostd::string2tokens(str,vtokens,delimiters);
    for(i=0;i<vtokens.size();i++) tokens.push_back(vtokens[i]);
    return tokens.size();
  }
  template<class utype> uint string2tokensAdd(const string& str,std::vector<utype>& tokens,const string& delimiters) {
    vector<string> vtokens;
    uint i=aurostd::string2tokens(str,vtokens,delimiters);
    for(i=0;i<vtokens.size();i++) tokens.push_back(aurostd::string2utype<utype>(vtokens[i]));
    return tokens.size();
  }
  template<class utype> uint string2tokensAdd(const string& str,std::deque<utype>& tokens,const string& delimiters) {
    deque<string> vtokens;
    uint i=aurostd::string2tokens(str,vtokens,delimiters);
    for(i=0;i<vtokens.size();i++) tokens.push_back(aurostd::string2utype<utype>(vtokens[i]));
    return tokens.size();
  }

  //[CO20210315 - OBSOLETE use stream2stream()]// ***************************************************************************
  //[CO20210315 - OBSOLETE use stream2stream()]// Function StringStreamConvert
  //[CO20210315 - OBSOLETE use stream2stream()]// ***************************************************************************
  //[CO20210315 - OBSOLETE use stream2stream()]// convert whatever into a string !
  //[CO20210315 - OBSOLETE use stream2stream()]template<typename typeTo, typename typeFrom> typeTo StringStreamConvert(const typeFrom& from) {  //CO20210315 - cleaned up
  //[CO20210315 - OBSOLETE use stream2stream()]  std::stringstream temp;
  //[CO20210315 - OBSOLETE use stream2stream()]  temp << from;
  //[CO20210315 - OBSOLETE use stream2stream()]  typeTo to=typeTo();
  //[CO20210315 - OBSOLETE use stream2stream()]  temp >> to;
  //[CO20210315 - OBSOLETE use stream2stream()]  return to;
  //[CO20210315 - OBSOLETE use stream2stream()]}
  //[CO20210315 - OBSOLETE use stream2stream()]
  //[CO20210315 - OBSOLETE use stream2stream()]template<typename typeFrom> std::string StringConvert(const typeFrom& from) { //CO20210315 - cleaned up
  //[CO20210315 - OBSOLETE use stream2stream()]  return StringStreamConvert<std::string>((typeFrom) from);
  //[CO20210315 - OBSOLETE use stream2stream()]}
  //[CO20210315 - OBSOLETE use stream2stream()]// to initialize the templates....
  //[CO20210315 - OBSOLETE use stream2stream()]void _StringConvert(void) {
  //[CO20210315 - OBSOLETE use stream2stream()]  cerr << StringStreamConvert<std::string>((int) 1);
  //[CO20210315 - OBSOLETE use stream2stream()]  cerr << StringStreamConvert<std::string>((float) 1);
  //[CO20210315 - OBSOLETE use stream2stream()]  cerr << StringStreamConvert<std::string>((double) 1);
  //[CO20210315 - OBSOLETE use stream2stream()]  cerr << StringConvert<int>((int) 1);
  //[CO20210315 - OBSOLETE use stream2stream()]  cerr << StringConvert<float>((float) 1);
  //[CO20210315 - OBSOLETE use stream2stream()]  cerr << StringConvert<double>((double) 1);
  //[CO20210315 - OBSOLETE use stream2stream()]}

  // ***************************************************************************
  // Function stream2stream
  // ***************************************************************************
  // convert whatever into a string !
  template<typename typeTo, typename typeFrom> typeTo stream2stream(const typeFrom& from,int precision,char FORMAT) { //CO20210315 - cleaned up
    std::stringstream temp;
    if(FORMAT==DEFAULT_STREAM){;} //default
    if(FORMAT==FIXED_STREAM){temp << std::fixed;}
    if(FORMAT==SCIENTIFIC_STREAM){temp << std::scientific;}
    temp.precision(precision);
    temp << from;
    typeTo to=typeTo();
    temp >> to;
    return to;
  }
  template<typename typeTo, typename typeFrom> typeTo stream2stream(const typeFrom& from,int precision) {  //CO20210315 - cleaned up
    return (typeTo) stream2stream<typeTo>(from,precision,DEFAULT_STREAM);
  }
  template<typename typeTo, typename typeFrom> typeTo stream2stream(const typeFrom& from) { //CO20210315 - cleaned up
    return (typeTo) stream2stream<typeTo>(from,AUROSTD_DEFAULT_PRECISION,DEFAULT_STREAM);
  }
  template<typename utype> utype string2utype(const string& from, const uint base) {
    if(from.empty()){return (utype) stream2stream<utype>("0",AUROSTD_DEFAULT_PRECISION,DEFAULT_STREAM);} //CO20210315 - stream2stream behavior is not defined for empty string input: https://stackoverflow.com/questions/4999650/c-how-do-i-check-if-the-cin-buffer-is-empty
    if(from.find("T")!=string::npos||from.find("t")!=string::npos||from.find("F")!=string::npos||from.find("f")!=string::npos){ //CO20221112
      //CO20221112 - gdb died here when running chull in parallel
      //seems toupper might have some issues with threaded processes
      //adding this guard to mitigate the issue
      string FROM=aurostd::toupper(from); //CO20210315
      if(FROM=="TRUE"||FROM=="T"||FROM==".TRUE."){return (utype) stream2stream<utype>("1",AUROSTD_DEFAULT_PRECISION,DEFAULT_STREAM);;}  //CO20210315 - safe because inputs are generally digits
      if(FROM=="FALSE"||FROM=="F"||FROM==".FALSE."){return (utype) stream2stream<utype>("0",AUROSTD_DEFAULT_PRECISION,DEFAULT_STREAM);;}  //CO20210315 - safe because inputs are generally digits
    }
    if (base != 10) { //HE20220324 add non-decimal bases (will ignore positions behind a point)
      std::stringstream temp;
      temp << std::stoll(from, nullptr, base); // stoll -> string to long long
      return (utype) stream2stream<utype>(temp.str(),AUROSTD_DEFAULT_PRECISION,DEFAULT_STREAM);
    }
    //[CO20210315 - doesn't work]if(!aurostd::isfloat(from)){return (utype) 0;} //CO20210315 - stream2stream undefined behavior
    return (utype) stream2stream<utype>(from,AUROSTD_DEFAULT_PRECISION,DEFAULT_STREAM);
  }

  string string2string(const string& from) {return from;}

  template<typename utype> vector<utype> vectorstring2vectorutype(const vector<string>& from) { //SD20220520
    vector<utype> vout;for(uint i=0;i<from.size();i++) vout.push_back(aurostd::string2utype<utype>(from[i]));return vout;
  }
  template<typename utype> vector<utype> vectorstring2vectorutype(const deque<string>& from) { //SD20220520
    vector<utype> vout;for(uint i=0;i<from.size();i++) vout.push_back(aurostd::string2utype<utype>(from[i]));return vout;
  }

  //[SD20220525 - OBSOLETE]vector<double> vectorstring2vectordouble(const vector<string>& from) {  //CO20210315 - cleaned up
  //[SD20220525 - OBSOLETE]  vector<double> vout;for(uint i=0;i<from.size();i++) vout.push_back(aurostd::string2utype<double>(from[i]));return vout;
  //[SD20220525 - OBSOLETE]}

  //[SD20220525 - OBSOLETE]vector<int> vectorstring2vectorint(const vector<string>& from) {  //CO20210315 - cleaned up
  //[SD20220525 - OBSOLETE]  vector<int> vout;for(uint i=0;i<from.size();i++) vout.push_back(aurostd::string2utype<int>(from[i]));return vout;
  //[SD20220525 - OBSOLETE]}

  //[SD20220525 - OBSOLETE]vector<uint> vectorstring2vectoruint(const vector<string>& from) { //CO20210315 - cleaned up
  //[SD20220525 - OBSOLETE]  vector<uint> vout;for(uint i=0;i<from.size();i++) vout.push_back(aurostd::string2utype<uint>(from[i]));return vout;
  //[SD20220525 - OBSOLETE]}

  //[SD20220525 - OBSOLETE]vector<float> vectorstring2vectorfloat(const vector<string>& from) { //CO20210315 - cleaned up
  //[SD20220525 - OBSOLETE]  vector<float> vout;for(uint i=0;i<from.size();i++) vout.push_back(aurostd::string2utype<float>(from[i]));return vout;
  //[SD20220525 - OBSOLETE]}

  string vectorstring2string(const vector<string>& vstrings) {
    string out="";
    for(uint istr=0;istr<vstrings.size();istr++) out+=vstrings.at(istr);
    return out;
  }
  string vectorstring2string(const deque<string>& vstrings) {
    string out="";
    for(uint istr=0;istr<vstrings.size();istr++) out+=vstrings.at(istr);
    return out;
  }

  // ***************************************************************************
  // Function utype2string
  // ***************************************************************************

  //  template<typename string> string utype2string(const string& from) {
  //   return (string) stream2stream<string>(from);
  // }

  //DX20210128 [OBSOLETE - use default arguments] template<typename utype> string utype2string(const utype& from) {
  //DX20210128 [OBSOLETE - use default arguments]   return (string) stream2stream<string>(from);
  //DX20210128 [OBSOLETE - use default arguments] }
  //DX20210128 [OBSOLETE - use default arguments] template<typename utype> string utype2string(const utype& from,int precision) {
  //DX20210128 [OBSOLETE - use default arguments]   return (string) stream2stream<string>(from,precision);
  //DX20210128 [OBSOLETE - use default arguments]}
  template<typename utype> string utype2string(const utype& from,int precision,char FORMAT) { //see DEFAULT_STREAM, FIXED_STREAM, SCIENTIFIC_STREAM
    return (string) stream2stream<string>(from,precision,FORMAT);
  }
  //  string utype2string(const string& from) {
  //    return (string) from;
  //  }
  //  string utype2string(const std::basic_string<char, std::char_traits<char>, std::allocator<char> >& from) {    return (string) from;  }
  //  string utype2string(std::basic_string<char, std::char_traits<char>, std::allocator<char> > from) {    return (string) from;  }

  //cannot template this the same as others, char's don't make sense with roff
  string utype2string(double from,bool roff) {return utype2string(from,AUROSTD_DEFAULT_PRECISION,roff,DEFAULT_STREAM);}
  string utype2string(double from,int precision,bool roff) {return utype2string(from,precision,roff,AUROSTD_ROUNDOFF_TOL,DEFAULT_STREAM);}
  string utype2string(double from,bool roff,double tol) {return utype2string(from,AUROSTD_DEFAULT_PRECISION,roff,tol,DEFAULT_STREAM);}
  string utype2string(double from,int precision,bool roff,double tol) {return utype2string(from,precision,roff,tol,DEFAULT_STREAM);}
  string utype2string(double from,bool roff,char FORMAT) {return utype2string(from,AUROSTD_DEFAULT_PRECISION,roff,FORMAT);}
  string utype2string(double from,int precision,char FORMAT,bool roff) {return utype2string(from,precision,roff,FORMAT);}  //CO20200624
  string utype2string(double from,int precision,bool roff,char FORMAT) {return utype2string(from,precision,roff,AUROSTD_ROUNDOFF_TOL,FORMAT);}
  string utype2string(double from,bool roff,double tol,char FORMAT) {return utype2string(from,AUROSTD_DEFAULT_PRECISION,roff,tol,FORMAT);}
  string utype2string(double from,int precision,bool roff,double tol,char FORMAT) {
    double tmp=from;
    if(roff){tmp=roundoff(from,tol);}
    return (string) stream2stream<string>(tmp,precision,FORMAT);
  }
  string bool2string(bool from) {
    if(from) return "TRUE";
    return "FALSE";
  }

  // ***************************************************************************
  // Function utypes2deque
  // ***************************************************************************
  template<class utype> deque<utype> utypes2deque(utype u1) {
    deque<utype> out;
    out.push_back(u1);
    return out;}
  template<class utype> deque<utype> utypes2deque(utype u1,utype u2) {
    deque<utype> out;
    out.push_back(u1);out.push_back(u2);
    return out;}
  template<class utype> deque<utype> utypes2deque(utype u1,utype u2,utype u3) {
    deque<utype> out;
    out.push_back(u1);out.push_back(u2);out.push_back(u3);
    return out;}
  template<class utype> deque<utype> utypes2deque(utype u1,utype u2,utype u3,utype u4) {
    deque<utype> out;
    out.push_back(u1);out.push_back(u2);out.push_back(u3);out.push_back(u4);
    return out;}

  // ***************************************************************************
  // Function StringCommasColumsVectorInt
  // ***************************************************************************
  void StringCommasColumsVectorInt(string vstring,vector<int> &vint) {
    vector<string> tokens_commas,tokens_colums;
    vint.clear();
    string2tokens(vstring,tokens_commas,",");
    for(uint i=0;i<tokens_commas.size();i++) {
      tokens_colums.clear();
      if(aurostd::substring2bool(tokens_commas[i],":")) {
        string2tokens(tokens_commas[i],tokens_colums,":");
        for(int j=aurostd::string2utype<int>(tokens_colums.at(0));
            j<=aurostd::string2utype<int>(tokens_colums.at(tokens_colums.size()-1));
            j++)
          vint.push_back(j);
      } else {
        vint.push_back(aurostd::string2utype<int>(tokens_commas[i]));
      }	
    }
  }

  // ***************************************************************************
  // Function StringCommasColumsVectorUnsignedInt
  // ***************************************************************************
  void StringCommasColumsVectorUnsignedInt(string vstring,vector<uint> &vuint) {
    vector<string> tokens_commas,tokens_colums;
    vuint.clear();
    string2tokens(vstring,tokens_commas,",");
    for(uint i=0;i<tokens_commas.size();i++) {
      tokens_colums.clear();
      if(aurostd::substring2bool(tokens_commas[i],":")) {
        string2tokens(tokens_commas[i],tokens_colums,":");
        for(uint j=aurostd::string2utype<uint>(tokens_colums.at(0));
            j<=(uint) aurostd::string2utype<uint>(tokens_colums.at(tokens_colums.size()-1));
            j++)
          vuint.push_back(j);
      } else {
        vuint.push_back(aurostd::string2utype<uint>(tokens_commas[i]));
      }	
    }
  }

  // ***************************************************************************
  // Function StringCommasColumsVectorFloat
  // ***************************************************************************
  void StringCommasColumsVectorFloat(string vstring,vector<float> &vfloat) {
    vector<string> tokens_commas,tokens_colums;
    vfloat.clear();
    string2tokens(vstring,tokens_commas,",");
    for(uint i=0;i<tokens_commas.size();i++) {
      tokens_colums.clear();
      if(aurostd::substring2bool(tokens_commas[i],":")) {
        string2tokens(tokens_commas[i],tokens_colums,":");
        for(float j=aurostd::string2utype<float>(tokens_colums.at(0));
            j<=aurostd::string2utype<float>(tokens_colums.at(tokens_colums.size()-1));
            j++)
          vfloat.push_back(j);
      } else {
        vfloat.push_back(aurostd::string2utype<float>(tokens_commas[i]));
      }	
    }
  }

  // ***************************************************************************
  // Function StringCommasColumsVectorDouble
  // ***************************************************************************
  void StringCommasColumsVectorDouble(string vstring,vector<double> &vdouble) {
    vector<string> tokens_commas,tokens_colums;
    vdouble.clear();
    string2tokens(vstring,tokens_commas,",");
    for(uint i=0;i<tokens_commas.size();i++) {
      tokens_colums.clear();
      if(aurostd::substring2bool(tokens_commas[i],":")) {
        string2tokens(tokens_commas[i],tokens_colums,":");
        for(double j=aurostd::string2utype<double>(tokens_colums.at(0));
            j<=aurostd::string2utype<double>(tokens_colums.at(tokens_colums.size()-1));
            j++)
          vdouble.push_back(j);
      } else {
        vdouble.push_back(aurostd::string2utype<double>(tokens_commas[i]));
      }	
    }
  }

  // ***************************************************************************
  // Function StringsAlphabetic
  // ***************************************************************************
  // says if two strings are in alphabetical order
  bool StringsAlphabetic(const string& A,const string& B,bool allow_identical) { //CO20181019
    if(A<B) return TRUE; // cerr << "A<B" << "  " << A << "<" << B << endl;
    if(A>B) return FALSE; //cerr << "A>B" << "  " << A << ">" << B << endl;
    if(A==B) return allow_identical; //cerr << "A==B" << "  " << A << "==" << B << endl; //CO20181019
    return TRUE;
  }

  // ***************************************************************************
  // Function StringsAlphabetic
  // ***************************************************************************
  // says if two strings are in alphabetical order  //CO20180801
  bool StringsAlphabetic(const vector<string>& input,bool allow_identical) {
    for(uint i=1;i<input.size();i++){ //CO20190218
      if(!StringsAlphabetic(input[i-1],input[i],allow_identical)){return false;}
    }
    return true;
  }
  bool StringsAlphabetic(const deque<string>& input,bool allow_identical) { //CO20190218
    for(uint i=1;i<input.size();i++){
      if(!StringsAlphabetic(input[i-1],input[i],allow_identical)){return false;}
    }
    return true;
  }

  // ***************************************************************************
  // Function StringSubst
  // ***************************************************************************
  // Substitute strings here and there
  // Stefano Curtarolo
  string StringSubst(string &strstring, const string &strfind, const string &strreplace) {
    if(strfind.empty()) return strstring;
    string::size_type pos=0;
    while((pos=strstring.find(strfind, pos))!=string::npos) {
      strstring.erase(pos, strfind.length());
      strstring.insert(pos, strreplace);
      pos+=strreplace.length();
    }
    return strstring;
  }

  //HE20220321 overload for const strings
  string StringSubst(const string &strstring, const string &strfind, const string &strreplace) {
    std::string work_copy = strstring;
    return StringSubst(work_copy, strfind, strreplace);
  }

  string StringSubst(string &strstring, const char &charfind, const char &charreplace) {
    string stroutput;
    for (uint i=0;i<strstring.size();i++)
      if(strstring[i]==charfind)
        stroutput+=charreplace;
      else
        stroutput+=strstring[i];
    strstring=stroutput;
    return strstring;
  }

  //HE20220321 overload for const strings
  string StringSubst(const string &strstring, const char &charfind, const char &charreplace) {
    std::string work_copy = strstring;
    return StringSubst(work_copy, charfind, charreplace);
  }

  void StringStreamSubst(stringstream &strstringstream, const string &strfind, const string &strreplace) {
    string strstring=strstringstream.str();
    StringSubst(strstring,strfind,strreplace);
    aurostd::StringstreamClean(strstringstream);
    strstringstream << strstring;
  }

  //  string StringSubst(string &strstring, const string &strfind0, const string &strfind1, const string &strfind2, const string &strfind3, const string &strreplace) {
  //    StringSubst(strstring,strfind0,strreplace);
  //    StringSubst(strstring,strfind1,strreplace);
  //    StringSubst(strstring,strfind2,strreplace);
  //    StringSubst(strstring,strfind3,strreplace);
  //    return strstring;
  //  }


  // ***************************************************************************
  // Function SubStrings
  // ***************************************************************************
  // Finds strings here and there.
  // Stefano Curtarolo

  int GetNLinesString(const string& str) {
    // SLOW
    //     string _str(str);
    //     int N=1;
    //     while(_str.find("\n")!=string::npos) {
    //       N++;_str=_str.substr(_str.find("\n")+1);
    //     }
    //     return N;
    stringstream strstream(str);
    return GetNLinesString(strstream);
  }
  int GetNLinesString(const stringstream& strstream) {
    // VERY SLOW    return aurostd::GetNLinesString(strstream.str());
    // FAST
    stringstream strstream_new(strstream.str());  //copy stringstream
    string line;
    int count_line = 0;
    while(getline(strstream_new, line)) {count_line++;}
    return count_line;
  }

  int GetNLinesFile(const string& file_name) {
    stringstream streamFILE;
    aurostd::file2stringstream(file_name, streamFILE);
    return GetNLinesString(streamFILE);
  }

  string GetLineString(const string& strstream,int line) {
    string _strstream(strstream),_strline;
    //  if(line>aurostd::GetNLinesString(_strstream)) return (string) "";   // TOO SLOW IF THE STRING IS LONG !
    for(int i=0;i<line;i++) {
      _strline=_strstream.substr(0,_strstream.find("\n"));
      _strstream=_strstream.substr(_strstream.find("\n")+1);
    }
    return _strline;
  }
  string GetLineString(const stringstream& strstream,int line) {
    return aurostd::GetLineString(strstream.str(),line);
  }

  // ***************************************************************************
  // Function SubStringsPresent
  // ***************************************************************************
  bool substring2bool(const string& strstream,const string& strsub1,bool RemoveWS,bool RemoveComments) { //CO20210315 - cleaned up
    //substring2bool and kvpair2bool are similar but distinct
    //substring2bool will match any strsub1 and return true
    //kvpair2bool will assume the line is written KEY+DELIMITER+VALUE, if it matches KEY and DELIMITER exactly, it will return true
    //matching KEY exactly is useful, e.g.:
    //_FILE_START_
    //IALGO=48
    //_FILE_END_
    //strsub1="ALGO": substring2bool will return true
    //keyword="ALGO,delim="=": kvpair2bool will return false
    //kvpair2bool must match KEY exactly! skips the rest
    //substring2bool is good for aflow.in's which has no set delimiter style: [AFLOW_BIN_XZ] vs. [AFLOW_BIN=XZ] vs. [AFLOW_BIN]XZ vs. [AFLOW]BIN=XZ
    bool LDEBUG=FALSE;//TRUE;
    if(LDEBUG) cerr << XPID << "aurostd::substring2bool(): BEGIN [substring=\"" << strsub1 << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
    string _strstream(strstream);
    if(RemoveWS==TRUE) _strstream=aurostd::RemoveWhiteSpaces(_strstream,'"');
    if(LDEBUG) cerr << XPID << "aurostd::substring2bool(): [input=\"" << strstream << "\"], [substring=\"" << strsub1 << "\"]" << endl;
    if(_strstream.find(strsub1)==string::npos) return false;
    if(RemoveComments){ //SD20220403 - substring exists, but now check if it exists outside of comments 
      vector<string> tokens;
      aurostd::string2tokens(_strstream,tokens,"\n");
      string strline="";
      for(uint i=0;i<tokens.size();i++) {
        strline=aurostd::RemoveComments(tokens[i]);  //CO20210315
        if(strline.find(strsub1)!=string::npos) {
          if(LDEBUG) cerr << XPID << "aurostd::substring2bool(): END [substring=\"" << strsub1 << "\" found] [RemoveWS=" << RemoveWS << "]" << endl;
          return true;
        }
      }
      if(LDEBUG) cerr << XPID << "aurostd::substring2bool(): END [substring=" << strsub1 << " NOT found] [RemoveWS=" << RemoveWS << "]" << endl;
      return false;
    }
    return true; //SD20220403 - since substring exists, return true
  }

  bool substring2bool(const vector<string>& vstrstream,const string& strsub1,bool RemoveWS,bool RemoveComments) {
    for(uint i=0;i<vstrstream.size();i++)
      if(aurostd::substring2bool(vstrstream[i],strsub1,RemoveWS,RemoveComments)) return TRUE;
    return FALSE;
  }
  bool substring2bool(const deque<string>& vstrstream,const string& strsub1,bool RemoveWS,bool RemoveComments) {
    for(uint i=0;i<vstrstream.size();i++)
      if(aurostd::substring2bool(vstrstream[i],strsub1,RemoveWS,RemoveComments)) return TRUE;
    return FALSE;
  }

  // ME20220505
  // Matches a list of substrings to a string
  // match_all: only substring must be inside the string
  bool substringlist2bool(const string& strin, const vector<string>& substrings, bool match_all) {
    for (uint i = 0; i < substrings.size(); i++) {
      if (strin.find(substrings[i]) == string::npos) {
        // Didn't find substring, but need all
        if (match_all) return false;
      } else if (!match_all) {
        // Found something and only need one
        return true;
      }
    }
    // Code only gets here when all substrings are
    // in the string (match_all) or none have (!match_all)
    return match_all;
  }
  
  bool substring2bool(const stringstream& strstream,const string& strsub1,bool RemoveWS,bool RemoveComments) {
    return aurostd::substring2bool(strstream.str(),strsub1,RemoveWS,RemoveComments);
  }

  bool substringlist2bool(const string& strin, const deque<string>& substrings, bool match_all) {
    for (uint i = 0; i < substrings.size(); i++) {
      if (strin.find(substrings[i]) == string::npos) {
        // Didn't find substring, but need all
        if (match_all) return false;
      } else if (!match_all) {
        // Found something and only need one
        return true;
      }
    }
    // Code only gets here when all substrings are
    // in the string (match_all) or none have (!match_all)
    return match_all;
  }

  /// @brief Determines whether a specific value exists inside a list
  ///
  /// @param list List of values
  /// @param input Specific value to match the values within the list
  /// @param index Index of the element that matched the input
  /// @param sorted Whether the list is sorted
  /// @return Boolean if the input was found in the list or not
  ///
  /// @note index is assigned the value of 0 if the input is not found in the list
  template<class utype> bool WithinList(const vector<utype>& list, const utype& input, size_t& index, bool sorted) { //SD20220705
    for (size_t i = 0; i < list.size(); i++) {
      if(sorted && list[i]>input){break;}
      if(aurostd::isequal(list[i],input)) {
        index = i;
        return true;
      }
    }
    index = 0;
    return false;
  }
  bool WithinList(const vector<string>& list, const string& input, size_t& index, bool sorted) { //SD20220705
    for (size_t i = 0; i < list.size(); i++) {
      if(sorted && list[i]>input){break;}
      if(aurostd::isequal(list[i],input)) {
        index = i;
        return true;
      }
    }
    index = 0;
    return false;
  }
  template<class utype> bool WithinList(const deque<utype>& list, const utype& input, size_t& index, bool sorted) { //SD20220705
    return WithinList(aurostd::deque2vector(list), input, index, sorted);
  }
  bool WithinList(const deque<string>& list, const string& input, size_t& index, bool sorted) { //SD20220705
    return WithinList(aurostd::deque2vector(list), input, index, sorted);
  }
  template<class utype> bool WithinList(const vector<utype>& list, const utype& input, bool sorted) { //SD20220705
    size_t index=0;
    return WithinList(list, input, index, sorted);
  }
  bool WithinList(const vector<string>& list, const string& input, bool sorted) { //SD20220705
    size_t index=0;
    return WithinList(list, input, index, sorted);
  }
  template<class utype> bool WithinList(const deque<utype>& list, const utype& input, bool sorted) { //SD20220705
    return WithinList(aurostd::deque2vector(list), input, sorted);
  }
  bool WithinList(const deque<string>& list, const string& input, bool sorted) { //SD20220705
    return WithinList(aurostd::deque2vector(list), input, sorted);
  }

  template<class utype> bool WithinList(const vector<utype>& list, const utype& input, vector<size_t>& index, bool sorted) { //SD20220705
    index.clear();
    for (size_t i = 0; i < list.size(); i++) {
      if(sorted && list[i]>input){break;}
      if(aurostd::isequal(list[i],input)) {
        index.push_back(i);
      }
    }
    return !index.empty();
  }
  bool WithinList(const vector<string>& list, const string& input, vector<size_t>& index, bool sorted) { //SD20220705
    index.clear();
    for (size_t i = 0; i < list.size(); i++) {
      if(sorted && list[i]>input){break;}
      if(aurostd::isequal(list[i],input)) {
        index.push_back(i);
      }
    }
    return !index.empty();
  }
  template<class utype> bool WithinList(const deque<utype>& list, const utype& input, vector<size_t>& index, bool sorted) { //SD20220705
    return WithinList(aurostd::deque2vector(list), input, index, sorted);
  }
  bool WithinList(const deque<string>& list, const string& input, vector<size_t>& index, bool sorted) { //SD20220705
    return WithinList(aurostd::deque2vector(list), input, index, sorted);
  }

  //[SD20220705 - OBSOLETE]bool WithinList(const vector<string>& list,const string& input,bool sorted) { //CO20181010
  //[SD20220705 - OBSOLETE]  //for(uint i=0;i<list.size();i++){if(list[i]==input){return true;}}  OBSOLETE ME20190905
  //[SD20220705 - OBSOLETE]  //return false;  OBSOLETE ME20190905
  //[SD20220705 - OBSOLETE]  int index=-1;
  //[SD20220705 - OBSOLETE]  return WithinList(list, input, index, sorted);
  //[SD20220705 - OBSOLETE]}
  //[SD20220705 - OBSOLETE]bool WithinList(const deque<string>& list,const string& input,bool sorted) { //CO20181010
  //[SD20220705 - OBSOLETE]  int index=-1;
  //[SD20220705 - OBSOLETE]  return WithinList(aurostd::deque2vector(list), input, index, sorted);
  //[SD20220705 - OBSOLETE]}
  //[SD20220705 - OBSOLETE]bool WithinList(const vector<int>& list,int input,bool sorted) {  //CO20181010
  //[SD20220705 - OBSOLETE]  //for(uint i=0;i<list.size();i++){if(list[i]==input){return true;}}  OBSOLETE ME20190905
  //[SD20220705 - OBSOLETE]  //return false;  OBSOLETE ME20190905
  //[SD20220705 - OBSOLETE]  int index=-1;
  //[SD20220705 - OBSOLETE]  return WithinList(list, input, index, sorted);
  //[SD20220705 - OBSOLETE]}
  //[SD20220705 - OBSOLETE]bool WithinList(const vector<uint>& list,uint input,bool sorted) {  //CO20181010
  //[SD20220705 - OBSOLETE]  //for(uint i=0;i<list.size();i++){if(list[i]==input){return true;}}  OBSOLETE ME20190905
  //[SD20220705 - OBSOLETE]  //return false;  OBSOLETE ME20190905
  //[SD20220705 - OBSOLETE]  int index=-1;
  //[SD20220705 - OBSOLETE]  return WithinList(list, input, index, sorted);
  //[SD20220705 - OBSOLETE]}

  //[SD20220705 - OBSOLETE]//ME20190813 - added versions that also determine the index of the item in the list
  //[SD20220705 - OBSOLETE]bool WithinList(const vector<string>& list, const string& input, int& index, bool sorted) {
  //[SD20220705 - OBSOLETE]  for (int i = 0, nlist = (int) list.size(); i < nlist; i++) {
  //[SD20220705 - OBSOLETE]    if(sorted && list[i]>input){break;} //CO20201111
  //[SD20220705 - OBSOLETE]    if(list[i]==input) {
  //[SD20220705 - OBSOLETE]      index = i;
  //[SD20220705 - OBSOLETE]      return true;
  //[SD20220705 - OBSOLETE]    }
  //[SD20220705 - OBSOLETE]  }
  //[SD20220705 - OBSOLETE]  index = -1;
  //[SD20220705 - OBSOLETE]  return false;
  //[SD20220705 - OBSOLETE]}

  //[SD20220705 - OBSOLETE]bool WithinList(const vector<int>& list, int input, int& index, bool sorted) {
  //[SD20220705 - OBSOLETE]  for (int i = 0, nlist = (int) list.size(); i < nlist; i++) {
  //[SD20220705 - OBSOLETE]    if(sorted && list[i]>input){break;} //CO20201111
  //[SD20220705 - OBSOLETE]    if(list[i]==input) {
  //[SD20220705 - OBSOLETE]      index = i;
  //[SD20220705 - OBSOLETE]      return true;
  //[SD20220705 - OBSOLETE]    }
  //[SD20220705 - OBSOLETE]  }
  //[SD20220705 - OBSOLETE]  index = -1;
  //[SD20220705 - OBSOLETE]  return false;
  //[SD20220705 - OBSOLETE]}

  //[SD20220705 - OBSOLETE]bool WithinList(const vector<uint>& list, uint input, int& index, bool sorted) {
  //[SD20220705 - OBSOLETE]  for (int i = 0, nlist = (int) list.size(); i < nlist; i++) {
  //[SD20220705 - OBSOLETE]    if(sorted && list[i]>input){break;} //CO20201111
  //[SD20220705 - OBSOLETE]    if(list[i]==input) {
  //[SD20220705 - OBSOLETE]      index = i;
  //[SD20220705 - OBSOLETE]      return true;
  //[SD20220705 - OBSOLETE]    }
  //[SD20220705 - OBSOLETE]  }
  //[SD20220705 - OBSOLETE]  index = -1;
  //[SD20220705 - OBSOLETE]  return false;
  //[SD20220705 - OBSOLETE]}

  //ME20220503
  bool SubstringWithinList(const deque<string>& list, const string& input) {
    int index = -1;
    return SubstringWithinList(list, input, index);
  }

  bool SubstringWithinList(const deque<string>& list, const string& input, int& index) {
    for (deque<string>::const_iterator it = list.begin(); it != list.end(); ++it) {
      if ((*it).find(input) != string::npos) {
        index = std::distance(list.begin(), it);
        return true;
      }
    }
    index = -1;
    return false;
  }

  bool SubstringWithinList(const vector<string>& list, const string& input) {
    int index = -1;
    return SubstringWithinList(list, input, index);
  }

  bool SubstringWithinList(const vector<string>& list, const string& input, int& index) {
    for (vector<string>::const_iterator it = list.begin(); it != list.end(); ++it) {
      if ((*it).find(input) != string::npos) {
        index = std::distance(list.begin(), it);
        return true;
      }
    }
    index = -1;
    return false;
  }

  bool EWithinList(const vector<string>& list,const string& input) { //CO20200223
    string output="";
    return EWithinList(list, input, output);
  }
  bool EWithinList(const vector<string>& list, const string& input, string& output) { //CO20200223
    output="";
    for (uint i = 0, nlist = list.size(); i < nlist; i++) {
      if(list[i]==input){output=input;return true;}
      if(list[i]==input+".xz"){output=input+".xz";return true;}
      if(list[i]==input+".gz"){output=input+".gz";return true;}
      if(list[i]==input+".bz2"){output=input+".bz2";return true;}
    }
    return false;
  }

  // ***************************************************************************
  bool substring_present_file(const string& FileName,const string& strsub1,bool RemoveWS,bool RemoveComments) {
    string StringFile;
    ifstream FileFile;
    FileFile.open(FileName.c_str(),std::ios::in);
    FileFile.clear();FileFile.seekg(0);
    if(!FileFile) {
      cerr << "ERROR  FileName=" << FileName << "   not found" << endl;
      return FALSE;
    }
    StringFile="";char c; while (FileFile.get(c)) StringFile+=c;
    FileFile.close();
    return aurostd::substring2bool(StringFile,strsub1,RemoveWS,RemoveComments);
  }

  bool substring_present_file_FAST(const string& FileName, const string& _strsub1, bool RemoveWS, bool case_insensitive,bool expect_near_end,unsigned long long int size_max) {
    //be careful, this does not filter-out # comments
    //CO20210315 - this function is not only fast, but it enables grepping through VERY large files,
    //whereas reading the whole file into memory can lead to out-of-memory issues for aflow
    //leverages regex power of grep (better than reading in cpp line by line and processing the string)
    //CO20210315 - adding case_insensitive
    //https://stackoverflow.com/questions/13913014/grepping-a-huge-file-80gb-any-way-to-speed-it-up
    //use LC_ALL=C to speed up grep
    //use fgrep if possible
    //use -m 1 to stop at the first match
    //use cat/tac if you know it comes at the end of the file
    //adding size_max: if the file is bigger than size_max, then do not search and return FALSE
    //files that are too big will freeze-up the grep command
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string message = "";

    if(!aurostd::FileExist(FileName)) {
      message = "file input not found =" + FileName;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    if(size_max!=AUROSTD_MAX_ULLINT){
      unsigned long long int fsize=aurostd::FileSize(FileName);
      if(fsize>=size_max){return false;}  //CO20210315 - in the future, consider throwing instead, and create a special function for vasp.out greps in kvasp
    }

    bool use_regex=true;  //do we want to use regex?
    bool is_regex_used=false; //do we actually use regex?

    string strsub1=_strsub1;
    if(case_insensitive==true){strsub1=aurostd::tolower(strsub1);}

    if(RemoveWS){
      if(use_regex){  //use regex
        aurostd::StringSubst(strsub1,"\t"," ");
        string strsub1_keep=strsub1;
        vector<string> vtokens;
        aurostd::string2tokens(strsub1,vtokens," ");
        strsub1=aurostd::joinWDelimiter(vtokens,"\\s*");
        if(!strsub1_keep.empty() && strsub1_keep[0]==' '){strsub1="\\s*"+strsub1;}  //insert at the front
        if(!strsub1_keep.empty() && strsub1_keep[strsub1_keep.size()-1]==' '){strsub1+="\\s*";} //append at the back
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " strsub1(regex)=\"" << strsub1 << "\"" << endl;}
        if(strsub1.find("\\s*")!=string::npos){is_regex_used=true;}
      }else{
        strsub1=aurostd::RemoveWhiteSpaces(strsub1);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " strsub1(!regex)=\"" << strsub1 << "\"" << endl;}
      }
    }

    string temp_file=aurostd::TmpFileCreate("substring");
    ostringstream aus;
    ifstream FileFile;
    int found=0;
    aurostd::StringstreamClean(aus);
    aus << "rm -f " << temp_file  << endl;
    if(1){aus << "LC_ALL=C ";}  //speed up grep

    string cat_command="cat";
    if(expect_near_end && aurostd::IsCommandAvailable("tac")){cat_command="tac";}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " cat_command=\"" << cat_command << "\"" << endl;}

    aus << cat_command << " \"" << FileName << "\"";
    if(use_regex==false && RemoveWS==true){aus << " | sed \"s/ //g\" | sed \"s/\\t//g\"";}

    //decide is_regex_used above here

    //if using regex, protect against incoming ()
    //e.g., "Total CPU time used (sec)"
    if(is_regex_used){
      aurostd::StringSubst(strsub1,"(","\\(");
      aurostd::StringSubst(strsub1,")","\\)");
    }

    string grep_command="grep";
    if(aurostd::IsCommandAvailable("fgrep") && is_regex_used==false){grep_command="fgrep";}  //fixed string search
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " grep_command=\"" << grep_command << "\"" << endl;}

    string grep_flags="-c"; //return count instead of strings matching
    if(case_insensitive==true){grep_flags="-ic";} //case-insensitive, should work with both grep and fgrep (tested by CO20210601)
    if(is_regex_used){grep_flags="-E "+grep_flags;} //regex flag
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " grep_flags=\"" << grep_flags << "\"" << endl;}

    aus << " | " << grep_command << " " << grep_flags << " -m 1 \"" << strsub1 << "\" > " << temp_file  << endl;  //-m 1: stop when you find a match
    aus << "echo >> " << temp_file << endl; // to give EOL
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " command=\"" << aus.str() << "\"" << endl;}
    aurostd::execute(aus);
    if(!aurostd::FileExist(temp_file)) {
      message = "file output not found =" + FileName;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    FileFile.open(temp_file.c_str(),std::ios::in);
    FileFile.clear();FileFile.seekg(0);
    FileFile >> found;
    FileFile.close();
    aurostd::StringstreamClean(aus);
#ifndef _AFLOW_TEMP_PRESERVE_
    aurostd::RemoveFile(temp_file);
#endif
    if(found>0) return true;
    return false;
  }

  // ***************************************************************************
  // Function SubStringsPresent and EXTRACT
  // ***************************************************************************
  //SD20220520
  //Rewritten substring2string to be more understandable, accept more input types, and incorporate extracting Nth entries
  //n==0 returns all entries, starts counting from 1, negative numbers go backwards
  //Original substring2string always returned just the first entry
  //When two substrings are present strsub1 is the start keyword and strsub2 is the stop keyword
  //CO20230502 - another rewrite for strsub1 and strsub2
  //substring2string() with strsub1 grabs everything on the line
  //substring2string() with strsub1 and strsub2 grabs everything in between
  //substring2strings() with strsub1 and strsub2 now handles newlines within entries of the return vector, as well as start/stop tags within the same line
  string substring2string(ifstream& input,const string& strsub1,const int instance,bool RemoveWS,bool RemoveComments) {
    //substring2string and kvpair2string are similar but distinct
    //substring2string will match any strsub1 and return everything AFTER strsub1
    //kvpair2string will assume the line is written KEY+DELIMITER+VALUE, if it matches KEY and DELIMITER exactly, it will return VALUE
    //matching KEY exactly is useful, e.g.:
    //_FILE_START_
    //IALGO==48
    //ALGO==FAST
    //_FILE_END_
    //strsub1="ALGO",instance=1: substring2string will return "==48"
    //keyword="ALGO,delim="==": kvpair2string will return "FAST"
    //kvpair2string must match KEY exactly! skips the rest
    //substring2string is good for aflow.in's which has no set delimiter style: [AFLOW_BIN_XZ] vs. [AFLOW_BIN=XZ] vs. [AFLOW_BIN]XZ vs. [AFLOW]BIN=XZ
    bool LDEBUG=FALSE;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "BEGIN [substring=\"" << strsub1 << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;}
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "[input=\"" << input.rdbuf() << "\"] [substring=\"" << strsub1 << "\"]" << endl;}
    string strline="";
    input.clear();input.seekg(0);
    vector<string> tokens;
    int iter=0;
    while(getline(input,strline) && (instance==0 || iter!=instance)) {
      if(RemoveWS) {strline=aurostd::RemoveWhiteSpaces(strline,'"');}
      if(RemoveComments) {strline=aurostd::RemoveComments(strline);}
      if(strline.find(strsub1)!=string::npos) {
        tokens.push_back(strline.substr(strline.find(strsub1)+strsub1.length()));
        iter++;
      }
    }
    input.clear();input.seekg(0);
    if(tokens.size()==0 || (uint)aurostd::abs(instance)>tokens.size()) {return "";}
    stringstream strstream;
    if(instance==0) {
      for(uint i=0;i<tokens.size();i++) {strstream << tokens[i] << endl;}
    }
    else if(instance>0) {
      strstream << tokens[tokens.size()-1];
    }
    else if(instance<0) {
      uint i=(uint)aurostd::boundary_conditions_periodic(0,tokens.size()-1,instance);
      strstream << tokens[i];
    }
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "END [substring=\"" << strsub1 << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;}
    return strstream.str();
  }

  string substring2string(const string& _input,const string& strsub1,const int instance,bool RemoveWS,bool RemoveComments) {
    bool LDEBUG=FALSE;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " BEGIN [substring=\"" << strsub1 << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;}
    string input=_input;
    if(RemoveWS) {input=aurostd::RemoveWhiteSpaces(_input,'"');}
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " [input=\"" << input << "\"] [substring=\"" << strsub1 << "\"]" << endl;}
    if(input.find(strsub1)==string::npos) {return "";}
    stringstream strstream;
    vector<string> tokens;
    aurostd::string2vectorstring(input,tokens);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " [tokens.size()=" << tokens.size() << "]" << endl;}
    int iter=0;
    if(instance>0) {
      for(uint i=0;i<tokens.size();i++) {
        if(RemoveComments) {tokens[i]=aurostd::RemoveComments(tokens[i]);}
        if(tokens[i].find(strsub1)!=string::npos) {
          iter++;
          if(instance==iter) {strstream << tokens[i].substr(tokens[i].find(strsub1)+strsub1.length());break;}
        }
      }
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " END [substring=\"" << strsub1 << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;}
      return strstream.str();
    }
    else if(instance<0) {
      for(int i=tokens.size()-1;i>=0;i--) {
        if(RemoveComments) {tokens[i]=aurostd::RemoveComments(tokens[i]);}
        if(tokens[i].find(strsub1)!=string::npos) {
          iter--;
          if(instance==iter) {strstream << tokens[i].substr(tokens[i].find(strsub1)+strsub1.length());break;}
        }
      }
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " END [substring=\"" << strsub1 << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;}
      return strstream.str();
    }
    else { //instance==0
      for(uint i=0;i<tokens.size();i++) {
        if(RemoveComments) {tokens[i]=aurostd::RemoveComments(tokens[i]);}
        if(tokens[i].find(strsub1)!=string::npos) {
          strstream << tokens[i].substr(tokens[i].find(strsub1)+strsub1.length()) << endl;
        }
      }
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " END [substring=\"" << strsub1 << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;}
      return strstream.str();
    }
    return "";
  }

  string substring2string(const stringstream& input,const string& strsub1,const int instance,bool RemoveWS,bool RemoveComments) {
    return substring2string(input.str(),strsub1,instance,RemoveWS,RemoveComments);
  }

  uint substring2stringInnerLoop(const string& _strline,bool& read,vector<string>& vlines,vector<string>& tokens,const string& strsub1,const string& strsub2,bool RemoveComments,bool trim_edges) { //CO20230502 - helper function for substring2strings()
    bool LDEBUG=false;
    string strline=_strline,strline2="",strline3="";
    bool found_strsub1=false,found_strsub2=false;
    if(RemoveComments) {strline=aurostd::RemoveComments(strline);}
    if(read==FALSE && strline.find(strsub1)!=string::npos) {
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " strline.find(strsub1)!=string::npos" << endl;}
      //CO20230501 - fixing in case strsub1 and strsub2 are in the same line
      while(strline.find(strsub1)!=string::npos){
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " (WHILE LOOP START)" << endl;
          cerr << __AFLOW_FUNC__ << " strline=" << strline << endl;}
        strline2=strline.substr(strline.find(strsub1)+strsub1.length());
        strline="";
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " strline2=" << strline2 << endl;
          cerr << __AFLOW_FUNC__ << " strline=" << strline << endl;
        }
        if(strline2.find(strsub2)==string::npos){
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " strline2.find(strsub2)==string::npos" << endl;}
          read=TRUE;
          vlines.push_back(strline2);
        }
        else{
          strline3=strline2.substr(0,strline2.find(strsub2));
          strline =strline2.substr(strline2.find(strsub2)+strsub2.length());
          strline2=strline3;
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " strline2.find(strsub2)!=string::npos" << endl;
            cerr << __AFLOW_FUNC__ << " strline2=" << strline2 << endl;
            cerr << __AFLOW_FUNC__ << " strline=" << strline << endl;
          }
          vlines.push_back(strline2);
          if(trim_edges){trimEmptyEdges(vlines);} //CO20230502
          tokens.push_back(aurostd::joinWDelimiter(vlines,"\n"));
          vlines.clear();
        }
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " (WHILE LOOP END)" << endl;}
      }
    }
    else if(read==TRUE && strline.find(strsub2)!=string::npos) {
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " strline.find(strsub2)!=string::npos)" << endl;}
      //CO20230501 - fixing in case strsub1 and strsub2 are in the same line
      while(!strline.empty() || strline.find(strsub2)!=string::npos){
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " (WHILE LOOP START)" << endl;
          cerr << __AFLOW_FUNC__ << " strline=" << strline << endl;
        }
        found_strsub2=(strline.find(strsub2)!=string::npos);
        strline2=strline;
        strline="";
        if(found_strsub2){
          strline3=strline2.substr(0,strline2.find(strsub2));
          strline =strline2.substr(strline2.find(strsub2)+strsub2.length());
          strline2=strline3;
        }
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " strline2=" << strline2 << endl;
          cerr << __AFLOW_FUNC__ << " strline=" << strline << endl;
        }
        vlines.push_back(strline2);
        if(found_strsub2){
          if(trim_edges){trimEmptyEdges(vlines);} //CO20230502
          tokens.push_back(aurostd::joinWDelimiter(vlines,"\n"));
          if(LDEBUG){
            for(uint i=0;i<vlines.size();i++){cerr << __AFLOW_FUNC__ << " vlines[i=" << i << "]=" << vlines[i] << endl;}
            for(uint i=0;i<tokens.size();i++){cerr << __AFLOW_FUNC__ << " tokens[i=" << i << "]=" << tokens[i] << endl;}
          }
          vlines.clear();
          if(strline.empty()){read=FALSE;}
        }
        if(!strline.empty()){
          found_strsub1=(strline.find(strsub1)!=string::npos);
          if(found_strsub1){
            strline=strline.substr(strline.find(strsub1)+strsub1.length());
          }
        }
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " read=" << read << endl;
          cerr << __AFLOW_FUNC__ << " (WHILE LOOP END)" << endl;
        }
      }
    }
    else if(read) {
      vlines.push_back(strline);
    }
    return tokens.size();
  }

  string substring2stringConsolidate(const vector<string>& tokens,const int instance){  //CO20230502 - helper function for substring2strings()
    if(tokens.size()==0 || (uint)aurostd::abs(instance)>tokens.size()) {return "";}
    stringstream strstream;
    if(instance==0) {
      for(uint i=0;i<tokens.size();i++) {strstream << tokens[i] << endl;}
    }
    else if(instance>0) {
      strstream << tokens[instance-1]; //[CO20230501 - OBSOLETE]tokens[tokens.size()-1];
    }
    else if(instance<0) { //CO20230501 - if instance<0, then we need to get full tokens to see how many to subtract
      uint i=(uint)aurostd::boundary_conditions_periodic(0,tokens.size()-1,instance);
      strstream << tokens[i];
    }
    return strstream.str();
  }

  string substring2string(ifstream& input,const string& strsub1,const string& strsub2,const int instance,bool RemoveWS,bool RemoveComments) {
    //CO20230502 - another rewrite to handle newlines within start/stop tags, and start/stop tags within the same line
    bool LDEBUG=FALSE;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "BEGIN [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;}
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "[input=\"" << input.rdbuf() << "\"] [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"]" << endl;}
    vector<string> tokens;
    aurostd::substring2strings(input,tokens,strsub1,strsub2,instance,RemoveWS,RemoveComments);
    string output=substring2stringConsolidate(tokens,instance);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "END [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [instance=" << instance << "] [output=\"" << output << "\"] [RemoveWS=" << RemoveWS << "]" << endl;}
    return output;
  }

  string substring2string(const string& input,const string& strsub1,const string& strsub2,const int instance,bool RemoveWS,bool RemoveComments) {
    //CO20230502 - another rewrite to handle newlines within start/stop tags, and start/stop tags within the same line
    bool LDEBUG=FALSE;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "BEGIN [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;}
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "[input=\"" << input << "\"] [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"]" << endl;}
    vector<string> tokens;
    aurostd::substring2strings(input,tokens,strsub1,strsub2,instance,RemoveWS,RemoveComments);
    string output=substring2stringConsolidate(tokens,instance);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "END [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [instance=" << instance << "] [output=\"" << output << "\"] [RemoveWS=" << RemoveWS << "]" << endl;}
    return output;
  }
  
  //[CO20230501 - OBSOLETE]string substring2string(const string& _input,const string& strsub1,const string& strsub2,const int instance,bool RemoveWS,bool RemoveComments) {
  //[CO20230501 - OBSOLETE]  bool LDEBUG=FALSE;
  //[CO20230501 - OBSOLETE]  if(LDEBUG) {cerr << __AFLOW_FUNC__ << "BEGIN [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;}
  //[CO20230501 - OBSOLETE]  string input=_input;
  //[CO20230501 - OBSOLETE]  if(RemoveWS) {input=aurostd::RemoveWhiteSpaces(_input,'"');}
  //[CO20230501 - OBSOLETE]  if(LDEBUG) {cerr << __AFLOW_FUNC__ << "[input=\"" << input << "\"] [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"]" << endl;}
  //[CO20230501 - OBSOLETE]  if(input.find(strsub1)==string::npos || input.find(strsub2)==string::npos) {return "";}
  //[CO20230501 - OBSOLETE]  stringstream strstream;
  //[CO20230501 - OBSOLETE]  vector<string> tokens;
  //[CO20230501 - OBSOLETE]  vector<uint> vstart, vstop;
  //[CO20230501 - OBSOLETE]  aurostd::string2vectorstring(input,tokens);
  //[CO20230501 - OBSOLETE]  uint istart=0,istop=0;
  //[CO20230501 - OBSOLETE]  int iter=0;
  //[CO20230501 - OBSOLETE]  if(instance>0) {
  //[CO20230501 - OBSOLETE]    for(uint i=0;i<tokens.size();i++) {
  //[CO20230501 - OBSOLETE]      if(RemoveComments) {tokens[i]=aurostd::RemoveComments(tokens[i]);}
  //[CO20230501 - OBSOLETE]      if(tokens[i].find(strsub2)!=string::npos) {
  //[CO20230501 - OBSOLETE]        istop=i-1;
  //[CO20230501 - OBSOLETE]        if(instance==iter) {vstop.push_back(istop);break;}
  //[CO20230501 - OBSOLETE]      }
  //[CO20230501 - OBSOLETE]      if(tokens[i].find(strsub1)!=string::npos) {
  //[CO20230501 - OBSOLETE]        iter++;
  //[CO20230501 - OBSOLETE]        istart=i+1;
  //[CO20230501 - OBSOLETE]        if(instance==iter) {vstart.push_back(istart);}
  //[CO20230501 - OBSOLETE]      }
  //[CO20230501 - OBSOLETE]    }
  //[CO20230501 - OBSOLETE]  }
  //[CO20230501 - OBSOLETE]  else if(instance<0) {
  //[CO20230501 - OBSOLETE]    for(int i=tokens.size()-1;i>=0;i--) {
  //[CO20230501 - OBSOLETE]      if(RemoveComments) {tokens[i]=aurostd::RemoveComments(tokens[i]);}
  //[CO20230501 - OBSOLETE]      if(tokens[i].find(strsub2)!=string::npos) {
  //[CO20230501 - OBSOLETE]        iter--;
  //[CO20230501 - OBSOLETE]        istop=i-1;
  //[CO20230501 - OBSOLETE]        if(instance==iter) {vstop.push_back(istop);}
  //[CO20230501 - OBSOLETE]      }
  //[CO20230501 - OBSOLETE]      if(tokens[i].find(strsub1)!=string::npos) {
  //[CO20230501 - OBSOLETE]        istart=i+1;
  //[CO20230501 - OBSOLETE]        if(instance==iter) {vstart.push_back(istart);break;}
  //[CO20230501 - OBSOLETE]      }
  //[CO20230501 - OBSOLETE]    }
  //[CO20230501 - OBSOLETE]  }
  //[CO20230501 - OBSOLETE]  else { //instance==0
  //[CO20230501 - OBSOLETE]    for(uint i=0;i<tokens.size();i++) {
  //[CO20230501 - OBSOLETE]      if(RemoveComments) {tokens[i]=aurostd::RemoveComments(tokens[i]);}
  //[CO20230501 - OBSOLETE]      if(tokens[i].find(strsub2)!=string::npos) {
  //[CO20230501 - OBSOLETE]        istop=i-1;
  //[CO20230501 - OBSOLETE]        vstop.push_back(istop);
  //[CO20230501 - OBSOLETE]      }
  //[CO20230501 - OBSOLETE]      if(tokens[i].find(strsub1)!=string::npos) {
  //[CO20230501 - OBSOLETE]        istart=i+1;
  //[CO20230501 - OBSOLETE]        vstart.push_back(istart);
  //[CO20230501 - OBSOLETE]      }
  //[CO20230501 - OBSOLETE]    }
  //[CO20230501 - OBSOLETE]  }
  //[CO20230501 - OBSOLETE]  if(vstop.empty()) {return "";}
  //[CO20230501 - OBSOLETE]  if(vstop.size()!=vstart.size()) {vstart.pop_back();} // remove unfinished start
  //[CO20230501 - OBSOLETE]  for(uint i=0;i<vstart.size();i++) {
  //[CO20230501 - OBSOLETE]    for(uint j=vstart[i];j<=vstop[i];j++) {strstream << tokens[j] << endl;}
  //[CO20230501 - OBSOLETE]  }
  //[CO20230501 - OBSOLETE]  if(LDEBUG) {cerr << __AFLOW_FUNC__ << "END [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;}
  //[CO20230501 - OBSOLETE]  return strstream.str();
  //[CO20230501 - OBSOLETE]}

  string substring2string(const stringstream& input,const string& strsub1,const string& strsub2,const int instance,bool RemoveWS,bool RemoveComments) { 
    return substring2string(input.str(),strsub1,strsub2,instance,RemoveWS,RemoveComments);
  }

  //[SD20220520 - OBSOLETE]string substring2string(const string& strstream,const string& strsub1,bool RemoveWS,bool RemoveComments) { //CO20210315 - cleaned up
  //[SD20220520 - OBSOLETE]  //substring2string and kvpair2string are similar but distinct
  //[SD20220520 - OBSOLETE]  //substring2string will match any strsub1 and return everything AFTER strsub1
  //[SD20220520 - OBSOLETE]  //kvpair2string will assume the line is written KEY+DELIMITER+VALUE, if it matches KEY and DELIMITER exactly, it will return VALUE
  //[SD20220520 - OBSOLETE]  //matching KEY exactly is useful, e.g.:
  //[SD20220520 - OBSOLETE]  //_FILE_START_
  //[SD20220520 - OBSOLETE]  //IALGO=48
  //[SD20220520 - OBSOLETE]  //ALGO=FAST
  //[SD20220520 - OBSOLETE]  //_FILE_END_
  //[SD20220520 - OBSOLETE]  //strsub1="ALGO": substring2string will return "=48"
  //[SD20220520 - OBSOLETE]  //keyword="ALGO,delim="=": kvpair2string will return "FAST"
  //[SD20220520 - OBSOLETE]  //kvpair2string must match KEY exactly! skips the rest
  //[SD20220520 - OBSOLETE]  //substring2string is good for aflow.in's which has no set delimiter style: [AFLOW_BIN_XZ] vs. [AFLOW_BIN=XZ] vs. [AFLOW_BIN]XZ vs. [AFLOW]BIN=XZ
  //[SD20220520 - OBSOLETE]  bool LDEBUG=FALSE;//TRUE;
  //[SD20220520 - OBSOLETE]  if(LDEBUG) cerr << XPID << "aurostd::substring2string(): BEGIN [substring=\"" << strsub1 << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
  //[SD20220520 - OBSOLETE]  string _strstream(strstream);
  //[SD20220520 - OBSOLETE]  if(RemoveWS==TRUE) _strstream=aurostd::RemoveWhiteSpaces(_strstream,'"');
  //[SD20220520 - OBSOLETE]  if(LDEBUG) cerr << XPID << "aurostd::substring2string(): [input=\"" << strstream << "\"], [substring=\"" << strsub1 << "\"]" << endl;
  //[SD20220520 - OBSOLETE]  if(_strstream.find(strsub1)==string::npos) return "";
  //[SD20220520 - OBSOLETE]
  //[SD20220520 - OBSOLETE]  vector<string> tokens;
  //[SD20220520 - OBSOLETE]  aurostd::string2tokens(_strstream,tokens,"\n");
  //[SD20220520 - OBSOLETE]  string strline="";
  //[SD20220520 - OBSOLETE]  string::size_type idxS1;
  //[SD20220520 - OBSOLETE]  for(uint i=0;i<tokens.size();i++) {
  //[SD20220520 - OBSOLETE]    if(RemoveComments){strline=aurostd::RemoveComments(tokens[i]);}  //CO20210315
  //[SD20220520 - OBSOLETE]    idxS1=strline.find(strsub1);
  //[SD20220520 - OBSOLETE]    if(idxS1!=string::npos) {
  //[SD20220520 - OBSOLETE]      strline=strline.substr(idxS1+strsub1.length());
  //[SD20220520 - OBSOLETE]      strline=aurostd::RemoveWhiteSpacesFromTheBack(strline);
  //[SD20220520 - OBSOLETE]      if(LDEBUG) cerr << XPID << "aurostd::substring2string(): END [substring=\"" << strsub1 << "\" found] [output=\"" << strline << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
  //[SD20220520 - OBSOLETE]      return strline;
  //[SD20220520 - OBSOLETE]    }
  //[SD20220520 - OBSOLETE]  }
  //[SD20220520 - OBSOLETE]  if(LDEBUG) cerr << XPID << "aurostd::substring2string(): END [substring=" << strsub1 << " NOT found] [RemoveWS=" << RemoveWS << "]" << endl;
  //[SD20220520 - OBSOLETE]  return "";
  //[SD20220520 - OBSOLETE]}

  //[CO20210315 - not used, not sure the purpose of strsub2]string substring2string(const string& strstream, const string& strsub1, const string& strsub2, bool RemoveWS) {
  //[CO20210315 - not used, not sure the purpose of strsub2]  bool LDEBUG=(FALSE || XHOST.DEBUG);
  //[CO20210315 - not used, not sure the purpose of strsub2]  if(LDEBUG) cerr << "DEBUG substring2string5: (BEGIN) " << strsub1 << " " << RemoveWS << endl;
  //[CO20210315 - not used, not sure the purpose of strsub2]  string _strstream(strstream),_strline,_strsub1(strsub1),_strsub2(strsub2);
  //[CO20210315 - not used, not sure the purpose of strsub2]  string strout="";
  //[CO20210315 - not used, not sure the purpose of strsub2]  string::size_type idxS1,idxS2;
  //[CO20210315 - not used, not sure the purpose of strsub2]  if(RemoveWS==TRUE) _strstream=aurostd::RemoveWhiteSpaces(_strstream,'"');
  //[CO20210315 - not used, not sure the purpose of strsub2]  if(_strstream.find(_strsub1)==string::npos) return (string) strout;
  //[CO20210315 - not used, not sure the purpose of strsub2]  if(_strstream.find(_strsub2)==string::npos) return (string) strout;
  //[CO20210315 - not used, not sure the purpose of strsub2]  //  transform(_strstream.begin(),_strstream.end(),_strstream.begin(),toupper); // pout everything UPPER
  //[CO20210315 - not used, not sure the purpose of strsub2]  vector<string> tokens;
  //[CO20210315 - not used, not sure the purpose of strsub2]  aurostd::string2tokens(_strstream,tokens,"\n");
  //[CO20210315 - not used, not sure the purpose of strsub2]  for(uint i=0;i<tokens.size();i++) {
  //[CO20210315 - not used, not sure the purpose of strsub2]    _strline=tokens[i];
  //[CO20210315 - not used, not sure the purpose of strsub2]    if(_strline.find(COMMENT_NEGLECT_1)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_1));
  //[CO20210315 - not used, not sure the purpose of strsub2]    if(_strline.find(COMMENT_NEGLECT_2)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_2));
  //[CO20210315 - not used, not sure the purpose of strsub2]    if(_strline.find(COMMENT_NEGLECT_3)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_3));
  //[CO20210315 - not used, not sure the purpose of strsub2]    idxS1=_strline.find(_strsub1);
  //[CO20210315 - not used, not sure the purpose of strsub2]    idxS2=_strline.find(_strsub2);
  //[CO20210315 - not used, not sure the purpose of strsub2]    if(idxS1!=string::npos && idxS2!=string::npos && idxS1<=idxS2) {
  //[CO20210315 - not used, not sure the purpose of strsub2]      strout=_strline.substr(std::max(_strline.find(_strsub2)+_strsub2.length(),_strline.find(_strsub1)+_strsub1.length()));
  //[CO20210315 - not used, not sure the purpose of strsub2]      strout=aurostd::RemoveWhiteSpacesFromTheBack(strout);
  //[CO20210315 - not used, not sure the purpose of strsub2]      if(LDEBUG) cerr << "DEBUG substring2string5: (END) " << strsub1 << " " << RemoveWS << endl;
  //[CO20210315 - not used, not sure the purpose of strsub2]      return (string) strout;
  //[CO20210315 - not used, not sure the purpose of strsub2]    }
  //[CO20210315 - not used, not sure the purpose of strsub2]  }
  //[CO20210315 - not used, not sure the purpose of strsub2]  return (string) strout;
  //[CO20210315 - not used, not sure the purpose of strsub2]}

  template<typename utype> utype substring2utype(ifstream& input,const string& strsub1,const int instance,bool RemoveWS,bool RemoveComments) {
    return string2utype<utype>(substring2string(input,strsub1,instance,RemoveWS,RemoveComments));
  }
  template<typename utype> utype substring2utype(const string& input,const string& strsub1,const int instance,bool RemoveWS,bool RemoveComments) {
    return string2utype<utype>(substring2string(input,strsub1,instance,RemoveWS,RemoveComments));
  }
  template<typename utype> utype substring2utype(const stringstream& input,const string& strsub1,const int instance,bool RemoveWS,bool RemoveComments) {
    return substring2utype<utype>(input.str(),strsub1,instance,RemoveWS,RemoveComments);
  }
  template<typename utype> utype substring2utype(ifstream& input,const string& strsub1,const string& strsub2,const int instance,bool RemoveWS,bool RemoveComments) {
    return string2utype<utype>(substring2string(input,strsub1,strsub2,instance,RemoveWS,RemoveComments));
  }
  template<typename utype> utype substring2utype(const string& input,const string& strsub1,const string& strsub2,const int instance,bool RemoveWS,bool RemoveComments) {
    return string2utype<utype>(substring2string(input,strsub1,strsub2,instance,RemoveWS,RemoveComments));
  }
  template<typename utype> utype substring2utype(const stringstream& input,const string& strsub1,const string& strsub2,const int instance,bool RemoveWS,bool RemoveComments) {
    return substring2utype<utype>(input.str(),strsub1,strsub2,instance,RemoveWS,RemoveComments);
  }

  bool kvpair2bool(ifstream& input,const string& keyword,const string& delim,bool RemoveWS,bool RemoveComments) { //SD20220520
    //substring2bool and kvpair2bool are similar but distinct
    //substring2bool will match any strsub1 and return true
    //kvpair2bool will assume the line is written KEY+DELIMITER+VALUE, if it matches KEY and DELIMITER exactly, it will return true
    //matching KEY exactly is useful, e.g.:
    //_FILE_START_
    //IALGO==48
    //_FILE_END_
    //strsub1="ALGO": substring2bool will return true
    //keyword="ALGO,delim="==": kvpair2bool will return false
    //kvpair2bool must match KEY exactly! skips the rest
    //substring2bool is good for aflow.in's which has no set delimiter style: [AFLOW_BIN_XZ] vs. [AFLOW_BIN=XZ] vs. [AFLOW_BIN]XZ vs. [AFLOW]BIN=XZ
    bool LDEBUG=FALSE;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [RemoveWS=" << RemoveWS << "]" << endl;}
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "[input=\"" << input.rdbuf() << "\"] [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"]" << endl;}
    string strline="",_keyword="";
    string::size_type idx=0;
    while(getline(input,strline)) {
      if(RemoveWS) {strline=aurostd::RemoveWhiteSpaces(strline,'"');}
      if(RemoveComments) {strline=aurostd::RemoveComments(strline);}
      idx=strline.find(delim);
      if(idx!=string::npos){
        _keyword=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(strline.substr(0,idx));
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << "_keyword=\"" << _keyword << "\"" << endl;}
        if(_keyword==keyword){
          if(LDEBUG) {cerr << __AFLOW_FUNC__ << "END [keyword=\"" << keyword << "\" found] [RemoveWS=" << RemoveWS << "]" << endl;}
          return true;
        }
      }
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << "END [keyword=\"" << keyword << "\" NOT found] [RemoveWS=" << RemoveWS << "]" << endl;
    return false;
  }

  bool kvpair2bool(const string& input,const string& keyword,const string& delim,bool RemoveWS,bool RemoveComments) { //CO20210315
    bool LDEBUG=FALSE;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [RemoveWS=" << RemoveWS << "]" << endl;}
    string _input(input);
    if(RemoveWS==TRUE) _input=aurostd::RemoveWhiteSpaces(_input,'"');
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "[input=\"" << input << "\"], [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"]" << endl;}
    if(_input.find(keyword)==string::npos) return false;

    if(RemoveComments){ //SD20220403 - substring exists, but now check if it exists outside of comments
      vector<string> tokens;
      aurostd::string2tokens(_input,tokens,"\n");
      string strline="",_keyword="";
      string::size_type idx=0;
      for(uint i=0;i<tokens.size();i++) {
        strline=aurostd::RemoveComments(tokens[i]);
        idx=strline.find(delim);
        if(idx!=string::npos){
          _keyword=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(strline.substr(0,idx));
          if(LDEBUG) {cerr << __AFLOW_FUNC__ << "_keyword=\"" << _keyword << "\"" << endl;}
          if(_keyword==keyword){
            if(LDEBUG) {cerr << __AFLOW_FUNC__ << "END [keyword=\"" << keyword << "\" found] [RemoveWS=" << RemoveWS << "]" << endl;}
            return true;
          }
        }
      }
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << "END [keyword=\"" << keyword << "\" NOT found] [RemoveWS=" << RemoveWS << "]" << endl;}
      return false;
    }
    return true; //SD20220403 - since substring exists, return true
  }

  bool kvpair2bool(const stringstream& input,const string& keyword,const string& delim,bool RemoveWS,bool RemoveComments) { //CO20210315 - cleaned up
    return kvpair2bool(input.str(),keyword,delim,RemoveWS,RemoveComments);
  }

  //SD20220520
  //Rewritten kvpair2string to be more understandable, accept more input types, allow for multi-char delimiters and incorporate extracting Nth entries
  //n==0 returns all entries, starts counting from 1, negative numbers go backwards
  //Original kvpair2string always returned just the first entry
  string kvpair2string(ifstream& input,const string& keyword,const string& delim,const int instance,bool RemoveWS,bool RemoveComments) {
    //substring2string and kvpair2string are similar but distinct
    //substring2string will match any strsub1 and return everything AFTER strsub1
    //kvpair2string will assume the line is written KEY+DELIMITER+VALUE, if it matches KEY and DELIMITER exactly, it will return VALUE
    //matching KEY exactly is useful, e.g.:
    //_FILE_START_
    //IALGO==48
    //ALGO==FAST
    //_FILE_END_
    //strsub1="ALGO",instance=1: substring2string will return "==48"
    //keyword="ALGO,delim="==": kvpair2string will return "FAST"
    //kvpair2string must match KEY exactly! skips the rest
    //substring2string is good for aflow.in's which has no set delimiter style: [AFLOW_BIN_XZ] vs. [AFLOW_BIN=XZ] vs. [AFLOW_BIN]XZ vs. [AFLOW]BIN=XZ
    bool LDEBUG=FALSE;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;}
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "[input=\"" << input.rdbuf() << "\"] [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"]" << endl;}
    string strline="",_keyword="";
    input.clear();input.seekg(0);
    vector<string> tokens;
    string::size_type idx=0;
    int iter=0;
    while(getline(input,strline) && (instance==0 || iter!=instance)) {
      if(RemoveWS) {strline=aurostd::RemoveWhiteSpaces(strline,'"');}
      if(RemoveComments) {strline=aurostd::RemoveComments(strline);}
      idx=strline.find(delim);
      if(idx!=string::npos) {
        _keyword=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(strline.substr(0,idx));
        if(_keyword==keyword) {
          tokens.push_back(aurostd::RemoveWhiteSpacesFromTheFrontAndBack(strline.substr(idx+delim.length())));
          iter++;
        }
      }
    }
    input.clear();input.seekg(0);
    if(tokens.size()==0 || (uint)aurostd::abs(instance)>tokens.size()) {return "";}
    stringstream strstream;
    if(instance==0) {
      for(uint i=0;i<tokens.size();i++) {strstream << tokens[i] << endl;}
    }
    else if(instance>0) {
      strstream << tokens[tokens.size()-1];
    }
    else if(instance<0) {
      uint i=(uint)aurostd::boundary_conditions_periodic(0,tokens.size()-1,instance);
      strstream << tokens[i];
    }
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;}
    return strstream.str();
  }

  string kvpair2string(const string& _input,const string& keyword,const string& delim,const int instance,bool RemoveWS,bool RemoveComments) {
    bool LDEBUG=FALSE;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;}
    string input=_input;
    if(RemoveWS) {input=aurostd::RemoveWhiteSpaces(_input,'"');}
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "[input=\"" << input << "\"] [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"]" << endl;}
    if(input.find(delim)==string::npos || input.find(keyword)==string::npos) {return "";}
    stringstream strstream;
    vector<string> tokens;
    aurostd::string2vectorstring(input,tokens);
    string _keyword="";
    string::size_type idx=0;
    int iter=0;
    if(instance>0) {
      for(uint i=0;i<tokens.size();i++) {
        if(RemoveComments) {tokens[i]=aurostd::RemoveComments(tokens[i]);}
        idx=tokens[i].find(delim);
        if(idx!=string::npos) {
          _keyword=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[i].substr(0,idx));
          if(_keyword==keyword) {
            iter++;
            if(instance==iter) {strstream << aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[i].substr(idx+delim.length()));break;}
          }
        }
      }
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;}
      return strstream.str();
    }
    else if(instance<0) {
      for(int i=tokens.size()-1;i>=0;i--) {
        if(RemoveComments) {tokens[i]=aurostd::RemoveComments(tokens[i]);}
        idx=tokens[i].find(delim);
        if(idx!=string::npos) {
          _keyword=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[i].substr(0,idx));
          if(_keyword==keyword) {
            iter--;
            if(instance==iter) {strstream << aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[i].substr(idx+delim.length()));break;}
          }
        }
      }
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;}
      return strstream.str();
    }
    else { //instance==0
      for(uint i=0;i<tokens.size();i++) {
        if(RemoveComments) {tokens[i]=aurostd::RemoveComments(tokens[i]);}
        idx=tokens[i].find(delim);
        if(idx!=string::npos) {
          _keyword=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[i].substr(0,idx));
          if(_keyword==keyword) {
            strstream << aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[i].substr(idx+delim.length())) << endl;
          }
        }
      }
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;}
      return strstream.str();
    }
    return "";
  }

  string kvpair2string(const stringstream& input,const string& keyword,const string& delim,const int instance,bool RemoveWS,bool RemoveComments) {
    return kvpair2string(input.str(),keyword,delim,instance,RemoveWS,RemoveComments);
  }

  //[SD20220520 - OBSOLETE]string kvpair2string(const string& strstream,const string& keyword,const string& delim,bool RemoveWS,bool RemoveComments) { //CO20210315
  //[SD20220520 - OBSOLETE]  //substring2string and kvpair2string are similar but distinct
  //[SD20220520 - OBSOLETE]  //substring2string will match any strsub1 and return everything AFTER strsub1
  //[SD20220520 - OBSOLETE]  //kvpair2string will assume the line is written KEY+DELIMITER+VALUE, if it matches KEY and DELIMITER exactly, it will return VALUE
  //[SD20220520 - OBSOLETE]  //matching KEY exactly is useful, e.g.:
  //[SD20220520 - OBSOLETE]  //_FILE_START_
  //[SD20220520 - OBSOLETE]  //IALGO=48
  //[SD20220520 - OBSOLETE]  //ALGO=FAST
  //[SD20220520 - OBSOLETE]  //_FILE_END_
  //[SD20220520 - OBSOLETE]  //strsub1="ALGO": substring2string will return "=48"
  //[SD20220520 - OBSOLETE]  //keyword="ALGO,delim="=": kvpair2string will return "FAST"
  //[SD20220520 - OBSOLETE]  //kvpair2string must match KEY exactly! skips the rest
  //[SD20220520 - OBSOLETE]  //substring2string is good for aflow.in's which has no set delimiter style: [AFLOW_BIN_XZ] vs. [AFLOW_BIN=XZ] vs. [AFLOW_BIN]XZ vs. [AFLOW]BIN=XZ
  //[SD20220520 - OBSOLETE]  bool LDEBUG=FALSE;//TRUE;
  //[SD20220520 - OBSOLETE]  if(LDEBUG) cerr << XPID << "aurostd::kvpair2string(): BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
  //[SD20220520 - OBSOLETE]  string _strstream(strstream);
  //[SD20220520 - OBSOLETE]  if(RemoveWS==TRUE) _strstream=aurostd::RemoveWhiteSpaces(_strstream,'"');
  //[SD20220520 - OBSOLETE]  if(LDEBUG) cerr << XPID << "aurostd::kvpair2string(): [input=\"" << strstream << "\"], [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"]" << endl;
  //[SD20220520 - OBSOLETE]  if(_strstream.find(keyword)==string::npos) return "";
  //[SD20220520 - OBSOLETE]
  //[SD20220520 - OBSOLETE]  vector<string> tokens;
  //[SD20220520 - OBSOLETE]  aurostd::string2tokens(_strstream,tokens,"\n");
  //[SD20220520 - OBSOLETE]  string strline="",_keyword="",value="";
  //[SD20220520 - OBSOLETE]  string::size_type idxS1;
  //[SD20220520 - OBSOLETE]  for(uint i=0;i<tokens.size();i++) {
  //[SD20220520 - OBSOLETE]    if(RemoveComments){strline=aurostd::RemoveComments(tokens[i]);}
  //[SD20220520 - OBSOLETE]    idxS1=strline.find(delim);
  //[SD20220520 - OBSOLETE]    if(idxS1!=string::npos){
  //[SD20220520 - OBSOLETE]      _keyword=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(strline.substr(0,idxS1));
  //[SD20220520 - OBSOLETE]      if(LDEBUG) cerr << XPID << "aurostd::kvpair2string(): _keyword=\"" << _keyword << "\"" << endl;
  //[SD20220520 - OBSOLETE]      if(_keyword==keyword){
  //[SD20220520 - OBSOLETE]        value=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(strline.substr(idxS1+1));
  //[SD20220520 - OBSOLETE]        if(LDEBUG) cerr << XPID << "aurostd::kvpair2string(): END [keyword=\"" << keyword << "\" found] [output=\"" << value << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
  //[SD20220520 - OBSOLETE]        return value;
  //[SD20220520 - OBSOLETE]      }
  //[SD20220520 - OBSOLETE]    }
  //[SD20220520 - OBSOLETE]  }
  //[SD20220520 - OBSOLETE]  if(LDEBUG) cerr << XPID << "aurostd::kvpair2string(): END [keyword=" << keyword << " NOT found] [RemoveWS=" << RemoveWS << "]" << endl;
  //[SD20220520 - OBSOLETE]  return "";
  //[SD20220520 - OBSOLETE]}

  template<typename utype> utype kvpair2utype(ifstream& input,const string& keyword,const string& delim,const int instance,bool RemoveWS,bool RemoveComments) {
    return string2utype<utype>(kvpair2string(input,keyword,delim,instance,RemoveWS,RemoveComments));
  }
  template<typename utype> utype kvpair2utype(const string& input,const string& keyword,const string& delim,const int instance,bool RemoveWS,bool RemoveComments) {
    return string2utype<utype>(kvpair2string(input,keyword,delim,instance,RemoveWS,RemoveComments));
  }
  template<typename utype> utype kvpair2utype(const stringstream& input,const string& keyword,const string& delim,const int instance,bool RemoveWS,bool RemoveComments) {
    return kvpair2utype<utype>(input.str(),keyword,delim,instance,RemoveWS,RemoveComments);
  }

  // ***************************************************************************
  // Function SubStringsPresentExtractString and other
  // ***************************************************************************
  uint substring2strings(ifstream& input,vector<string> &vstringout,const string& strsub1,bool RemoveWS,bool RemoveComments) { //SD20220520
    vstringout=aurostd::string2vectorstring(aurostd::substring2string(input,strsub1,0,RemoveWS,RemoveComments));
    return vstringout.size();
  }
  uint substring2strings(const string& input,vector<string> &vstringout,const string& strsub1,bool RemoveWS,bool RemoveComments) {
    //[SD20220520 - OBSOLETE]//an AFLOW-specific substring matcher: will remove comments before searching
    //[SD20220520 - OBSOLETE]bool LDEBUG=FALSE;//TRUE;
    //[SD20220520 - OBSOLETE]if(LDEBUG) cerr << XPID << "aurostd::substring2strings(): BEGIN [substring=\"" << strsub1 << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
    //[SD20220520 - OBSOLETE]vstringout.clear();
    //[SD20220520 - OBSOLETE]string _strstream(strstream);
    //[SD20220520 - OBSOLETE]if(RemoveWS==TRUE) _strstream=aurostd::RemoveWhiteSpaces(_strstream,'"');
    //[SD20220520 - OBSOLETE]if(LDEBUG) cerr << XPID << "aurostd::substring2strings(): [input=\"" << strstream << "\"], [substring=\"" << strsub1 << "\"]" << endl;
    //[SD20220520 - OBSOLETE]if(_strstream.find(strsub1)==string::npos) return 0;
    //[SD20220520 - OBSOLETE]
    //[SD20220520 - OBSOLETE]vector<string> tokens;
    //[SD20220520 - OBSOLETE]aurostd::string2tokens(_strstream,tokens,"\n");
    //[SD20220520 - OBSOLETE]string strline="";
    //[SD20220520 - OBSOLETE]string::size_type idxS1;
    //[SD20220520 - OBSOLETE]for(uint i=0;i<tokens.size();i++) {
    //[SD20220520 - OBSOLETE]  if(RemoveComments){strline=aurostd::RemoveComments(tokens[i]);}  //CO20210315
    //[SD20220520 - OBSOLETE]  idxS1=strline.find(strsub1);
    //[SD20220520 - OBSOLETE]  if(idxS1!=string::npos) {
    //[SD20220520 - OBSOLETE]    strline=strline.substr(idxS1+strsub1.length());
    //[SD20220520 - OBSOLETE]    strline=aurostd::RemoveWhiteSpacesFromTheBack(strline);
    //[SD20220520 - OBSOLETE]    if(LDEBUG) cerr << XPID << "aurostd::substring2strings(): [substring=\"" << strsub1 << "\" found] [output=\"" << strline << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
    //[SD20220520 - OBSOLETE]    vstringout.push_back(strline);
    //[SD20220520 - OBSOLETE]  }
    //[SD20220520 - OBSOLETE]}
    //[SD20220520 - OBSOLETE]if(LDEBUG) cerr << XPID << "aurostd::substring2strings(): END [substring=\"" << strsub1 << "\"] [hits=" << vstringout.size() << "] [RemoveWS=" << RemoveWS << "]" << endl;
    //[SD20220520 - OBSOLETE]return vstringout.size();
    vstringout=aurostd::string2vectorstring(aurostd::substring2string(input,strsub1,0,RemoveWS,RemoveComments));
    return vstringout.size();
  }
  uint substring2strings(const stringstream& input,vector<string> &vstringout,const string& strsub1,bool RemoveWS,bool RemoveComments) {
    vstringout=aurostd::string2vectorstring(aurostd::substring2string(input,strsub1,0,RemoveWS,RemoveComments));
    return vstringout.size();
  }

  uint substring2strings(ifstream& input,vector<string> &vstringout,const string& strsub1,const string& strsub2,const int instance,bool RemoveWS,bool RemoveComments,bool trim_edges) { //CO20230502
    //[CO20230502 - OBSOLETE]vstringout=aurostd::string2vectorstring(aurostd::substring2string(input,strsub1,strsub2,0,RemoveWS,RemoveComments));
    //[CO20230502 - OBSOLETE]return vstringout.size();
    //CO20230502 - another rewrite to handle newlines within start/stop tags, and start/stop tags within the same line
    bool LDEBUG=false;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "BEGIN [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;}
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "[input=\"" << input.rdbuf() << "\"] [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"]" << endl;}
    vstringout.clear();
    string strline="";
    input.clear();input.seekg(0);
    vector<string> vlines;
    int iter=0;
    bool read=FALSE;
    while(getline(input,strline) && (instance<=0 || iter<instance)) {
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " strline=" << strline << endl;}
      if(RemoveWS) {strline=aurostd::RemoveWhiteSpaces(strline,'"');}
      //
      iter=substring2stringInnerLoop(strline,read,vlines,vstringout,strsub1,strsub2,RemoveComments,trim_edges); //CO20230502
      //
    }
    input.clear();input.seekg(0);
    uint count=vstringout.size();
    if(LDEBUG){
      for(uint i=0;i<count;i++) {
        cerr << __AFLOW_FUNC__ << " vstringout[i=" << i << "]=\"" <<  vstringout[i] << "\"" << endl;
      }
    }
    return count;
  }
  uint substring2strings(const string& _input,vector<string> &vstringout,const string& strsub1,const string& strsub2,const int instance,bool RemoveWS,bool RemoveComments,bool trim_edges) {  //CO20230502
    //[SD20220520 - OBSOLETE]bool LDEBUG=(FALSE || XHOST.DEBUG);
    //[SD20220520 - OBSOLETE]if(LDEBUG) cerr << "DEBUG substring2strings5: (BEGIN) " << strsub1 << " " << strsub2 << " " << RemoveWS << endl;
    //[SD20220520 - OBSOLETE]string _strstream(strstream),_strline,_strsub1(strsub1),_strsub2(strsub2);
    //[SD20220520 - OBSOLETE]string::size_type idxS1,idxS2;
    //[SD20220520 - OBSOLETE]vstringout.clear(); // clear so it is empty
    //[SD20220520 - OBSOLETE]if(RemoveWS==TRUE) _strstream=aurostd::RemoveWhiteSpaces(_strstream,'"');
    //[SD20220520 - OBSOLETE]if(_strstream.find(_strsub1)==string::npos) return 0; // there is not
    //[SD20220520 - OBSOLETE]if(_strstream.find(_strsub2)==string::npos) return 0; // there is not
    //[SD20220520 - OBSOLETE]//  transform(_strstream.begin(),_strstream.end(),_strstream.begin(),toupper); // pout everything UPPER
    //[SD20220520 - OBSOLETE]vector<string> tokens;
    //[SD20220520 - OBSOLETE]aurostd::string2tokens(_strstream,tokens,"\n");
    //[SD20220520 - OBSOLETE]for(uint i=0;i<tokens.size();i++) {
    //[SD20220520 - OBSOLETE]  _strline=tokens[i];
    //[SD20220520 - OBSOLETE]  if(_strline.find(COMMENT_NEGLECT_1)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_1));
    //[SD20220520 - OBSOLETE]  if(_strline.find(COMMENT_NEGLECT_2)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_2));
    //[SD20220520 - OBSOLETE]  if(_strline.find(COMMENT_NEGLECT_3)!=string::npos) _strline=_strline.substr(0,_strline.find(COMMENT_NEGLECT_3));
    //[SD20220520 - OBSOLETE]  idxS1=_strline.find(_strsub1);
    //[SD20220520 - OBSOLETE]  idxS2=_strline.find(_strsub2);
    //[SD20220520 - OBSOLETE]  if(idxS1!=string::npos && idxS2!=string::npos) {
    //[SD20220520 - OBSOLETE]    vstringout.push_back(aurostd::RemoveWhiteSpacesFromTheBack(_strline.substr(std::max(_strline.find(_strsub2)+_strsub2.length(),_strline.find(_strsub1)+_strsub1.length()))));
    //[SD20220520 - OBSOLETE]  }
    //[SD20220520 - OBSOLETE]}
    //[SD20220520 - OBSOLETE]if(LDEBUG) cerr << "DEBUG substring2string5: (END) " << strsub1 << " " << strsub2 << " " << vstringout.size() << " " << RemoveWS << endl;
    //[SD20220520 - OBSOLETE]return vstringout.size();
    //[CO20230502 - OBSOLETE]vstringout=aurostd::string2vectorstring(aurostd::substring2string(input,strsub1,strsub2,0,RemoveWS,RemoveComments));
    //[CO20230502 - OBSOLETE]return vstringout.size();
    //CO20230502 - another rewrite to handle newlines within start/stop tags, and start/stop tags within the same line
    bool LDEBUG=FALSE;
    string input=_input;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "BEGIN [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;}
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << "[input=\"" << input << "\"] [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"]" << endl;}
    vstringout.clear();
    if(RemoveWS) {input=aurostd::RemoveWhiteSpaces(input,'"');}
    if(input.find(strsub1)==string::npos && input.find(strsub2)==string::npos) {return vstringout.size();}
    string strline="";
    vector<string> vlines,vlines2;
    aurostd::string2vectorstring(input,vlines2);
    int iter=0;
    bool read=false;
    for(uint iline=0;iline<vlines2.size() && (instance<=0 || iter<instance);iline++){
      strline=vlines2[iline];
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " strline=" << strline << endl;}
      //
      iter=substring2stringInnerLoop(strline,read,vlines,vstringout,strsub1,strsub2,RemoveComments,trim_edges); //CO20230502
      //
    }
    uint count=vstringout.size();
    if(LDEBUG){
      for(uint i=0;i<count;i++) {
        cerr << __AFLOW_FUNC__ << " vstringout[i=" << i << "]=\"" <<  vstringout[i] << "\"" << endl;
      }
    }
    return count;
  }
  uint substring2strings(const stringstream& input,vector<string> &vstringout,const string& strsub1,const string& strsub2,const int instance,bool RemoveWS,bool RemoveComments,bool trim_edges) { //CO20230502
    //[CO20230502 - OBSOLETE]vstringout=aurostd::string2vectorstring(aurostd::substring2string(input,strsub1,strsub2,0,RemoveWS,RemoveComments));
    //[CO20230502 - OBSOLETE]return vstringout.size();
    return aurostd::substring2strings(input.str(),vstringout,strsub1,strsub2,instance,RemoveWS,RemoveComments,trim_edges);  //CO20230502
  }

  template<typename utype> uint substring2utypes(ifstream& input,vector<utype> &vutypeout,const string& strsub1,bool RemoveWS,bool RemoveComments) {  //SD202205020
    vector<string> vstringout;
    aurostd::substring2strings(input,vstringout,strsub1,RemoveWS,RemoveComments);
    vutypeout=aurostd::vectorstring2vectorutype<utype>(vstringout);
    return vutypeout.size();
  }
  template<typename utype> uint substring2utypes(const string& input,vector<utype> &vutypeout,const string& strsub1,bool RemoveWS,bool RemoveComments) {  //CO20210315 - cleaned up //SD202205020 - added utype
    vector<string> vstringout;
    aurostd::substring2strings(input,vstringout,strsub1,RemoveWS,RemoveComments);
    vutypeout=aurostd::vectorstring2vectorutype<utype>(vstringout);
    return vutypeout.size();
  }
  template<typename utype> uint substring2utypes(const stringstream& input,vector<utype> &vutypeout,const string& strsub1,bool RemoveWS,bool RemoveComments) {  //CO20210315 - cleaned up //SD202205020 - added utype
    return substring2utypes<utype>(input.str(),vutypeout,strsub1,RemoveWS,RemoveComments);
  }

  template<typename utype> uint substring2utypes(ifstream& input,vector<utype> &vutypeout,const string& strsub1,const string& strsub2,bool RemoveWS,bool RemoveComments) { //SD202205020
    vector<string> vstringout;
    aurostd::substring2strings(input,vstringout,strsub1,strsub2,RemoveWS,RemoveComments);
    vutypeout=aurostd::vectorstring2vectorutype<utype>(vstringout);
    return vutypeout.size();
  }
  template<typename utype> uint substring2utypes(const string& input,vector<utype> &vutypeout,const string& strsub1,const string& strsub2,bool RemoveWS,bool RemoveComments) { //SD202205020 - added utype
    vector<string> vstringout;
    aurostd::substring2strings(input,vstringout,strsub1,strsub2,RemoveWS,RemoveComments);
    vutypeout=aurostd::vectorstring2vectorutype<utype>(vstringout);
    return vutypeout.size();
  }
  template<typename utype> uint substring2utypes(const stringstream& input,vector<utype> &vutypeout,const string& strsub1,const string& strsub2,bool RemoveWS,bool RemoveComments) { //SD202205020
    return substring2utypes<utype>(input.str(),vutypeout,strsub1,strsub2,RemoveWS,RemoveComments);
  }
}

// ***************************************************************************
// FUNCTION HTML LATEX TXT

namespace aurostd {

  // http://www.w3schools.com/html/html_entities.asp
  // http://en.wikibooks.org/wiki/LaTeX/Accents

  //ME20200921 - Replaces HTML special characters with the correct entity name
  string text2html(const string& str) {
    string out = str;
    // Ampersand must come first since it is in the entity name!
    aurostd::StringSubst(out, "&", "&amp;");
    aurostd::StringSubst(out, "<", "&lt;");
    aurostd::StringSubst(out, ">", "&gt;");
    aurostd::StringSubst(out, "\"", "&quot;");
    aurostd::StringSubst(out, "'", "&apos;");
    return out;
  }

  string html2latex(const string& str) {
    string out=str;
    aurostd::StringSubst(out,"_","\\_");
    aurostd::StringSubst(out,"<sub>","$_{");aurostd::StringSubst(out,"</sub>","}$");
    aurostd::StringSubst(out,"<i>","\\textit{");aurostd::StringSubst(out,"</i>","}");
    aurostd::StringSubst(out,"<b>","\\textbf{"); aurostd::StringSubst(out,"</b>","}");
    aurostd::StringSubst(out,"<blink>","\\textbf{"); aurostd::StringSubst(out,"</blink>","}");
    aurostd::StringSubst(out,"MgB2","MgB$_2$");
    aurostd::StringSubst(out,"Schuttler","Sch\\\"uttler");
    aurostd::StringSubst(out,"Csányi","Cs\\'anyi");aurostd::StringSubst(out,"Csanyi","Cs\\'anyi");
    aurostd::StringSubst(out,"Pólya","P\\'{o}lya");
    if(!aurostd::substring2bool(out,"Rosenbrock")) aurostd::StringSubst(out,"Rosen","Ros\\'en");
    // string bar="";//;bar.at(0)=92;

    // http://en.wikibooks.org/wiki/LaTeX/Accents
    // umlaut
    aurostd::StringSubst(out,"&auml;","\\\"{a}");aurostd::StringSubst(out,"&Auml;","\\\"{A}");
    aurostd::StringSubst(out,"&euml;","\\\"{e}");aurostd::StringSubst(out,"&Euml;","\\\"{E}");
    aurostd::StringSubst(out,"&iuml;","\\\"{i}");aurostd::StringSubst(out,"&Iuml;","\\\"{I}");
    aurostd::StringSubst(out,"&ouml;","\\\"{o}");aurostd::StringSubst(out,"&Ouml;","\\\"{O}");
    aurostd::StringSubst(out,"&uuml;","\\\"{u}");aurostd::StringSubst(out,"&Uuml;","\\\"{U}");
    // grave accent
    aurostd::StringSubst(out,"&agrave;","\\`{a}");aurostd::StringSubst(out,"&Agrave;","\\`{A}");
    aurostd::StringSubst(out,"&egrave;","\\`{e}");aurostd::StringSubst(out,"&Egrave;","\\`{E}");
    aurostd::StringSubst(out,"&igrave;","\\`{i}");aurostd::StringSubst(out,"&Igrave;","\\`{I}");
    aurostd::StringSubst(out,"&ograve;","\\`{o}");aurostd::StringSubst(out,"&Ograve;","\\`{O}");
    aurostd::StringSubst(out,"&ugrave;","\\`{u}");aurostd::StringSubst(out,"&Ugrave;","\\`{U}");
    // acute accent
    aurostd::StringSubst(out,"&aacute;","\\'{a}");aurostd::StringSubst(out,"&Aacute;","\\'{A}");
    aurostd::StringSubst(out,"&eacute;","\\'{e}");aurostd::StringSubst(out,"&Eacute;","\\'{E}");
    aurostd::StringSubst(out,"&iacute;","\\'{i}");aurostd::StringSubst(out,"&Iacute;","\\'{I}");
    aurostd::StringSubst(out,"&oacute;","\\'{o}");aurostd::StringSubst(out,"&Oacute;","\\'{O}");
    aurostd::StringSubst(out,"&uacute;","\\'{u}");aurostd::StringSubst(out,"&Uacute;","\\'{U}");
    // tilde
    aurostd::StringSubst(out,"&atilde;","\\~{a}");aurostd::StringSubst(out,"&Atilde;","\\~{A}");
    aurostd::StringSubst(out,"&etilde;","\\~{e}");aurostd::StringSubst(out,"&Etilde;","\\~{E}");
    aurostd::StringSubst(out,"&itilde;","\\~{i}");aurostd::StringSubst(out,"&Itilde;","\\~{I}");
    aurostd::StringSubst(out,"&otilde;","\\~{o}");aurostd::StringSubst(out,"&Otilde;","\\~{O}");
    aurostd::StringSubst(out,"&utilde;","\\~{u}");aurostd::StringSubst(out,"&Utilde;","\\~{U}");
    // circ
    aurostd::StringSubst(out,"&acirc;","\\^{a}");aurostd::StringSubst(out,"&Acirc;","\\^{A}");
    aurostd::StringSubst(out,"&ecirc;","\\^{e}");aurostd::StringSubst(out,"&Ecirc;","\\^{E}");
    aurostd::StringSubst(out,"&icirc;","\\^{i}");aurostd::StringSubst(out,"&Icirc;","\\^{I}");
    aurostd::StringSubst(out,"&ocirc;","\\^{o}");aurostd::StringSubst(out,"&Ocirc;","\\^{O}");
    aurostd::StringSubst(out,"&ucirc;","\\^{u}");aurostd::StringSubst(out,"&Ucirc;","\\^{U}");
    // ring
    aurostd::StringSubst(out,"&aring;","\\r{a}");aurostd::StringSubst(out,"&Aring;","\\r{A}");
    aurostd::StringSubst(out,"&ering;","\\r{e}");aurostd::StringSubst(out,"&Ering;","\\r{E}");
    aurostd::StringSubst(out,"&iring;","\\r{i}");aurostd::StringSubst(out,"&Iring;","\\r{I}");
    aurostd::StringSubst(out,"&oring;","\\r{o}");aurostd::StringSubst(out,"&Oring;","\\r{O}");
    aurostd::StringSubst(out,"&uring;","\\r{u}");aurostd::StringSubst(out,"&Uring;","\\r{U}");
    // cedil
    aurostd::StringSubst(out,"&acedil;","\\c{a}");aurostd::StringSubst(out,"&Acedil;","\\c{A}");
    aurostd::StringSubst(out,"&ecedil;","\\c{e}");aurostd::StringSubst(out,"&Ecedil;","\\c{E}");
    aurostd::StringSubst(out,"&icedil;","\\c{i}");aurostd::StringSubst(out,"&Icedil;","\\c{I}");
    aurostd::StringSubst(out,"&ocedil;","\\c{o}");aurostd::StringSubst(out,"&Ocedil;","\\c{O}");
    aurostd::StringSubst(out,"&ucedil;","\\c{u}");aurostd::StringSubst(out,"&Ucedil;","\\c{U}");
    // caron
    aurostd::StringSubst(out,"&zcaron;","{\\v{z}}");aurostd::StringSubst(out,"&Zcaron;","{\\v{Z}}");
    // slash
    aurostd::StringSubst(out,"&oslash;","{\\o}");aurostd::StringSubst(out,"&Oslash;","{\\O}");
    // math
    aurostd::StringSubst(out,"&Alpha;","$\\Alpha$");aurostd::StringSubst(out,"&alpha;","$\\alpha$");
    aurostd::StringSubst(out,"&Beta;","$\\Βeta$");aurostd::StringSubst(out,"&beta;","$\\beta$");
    aurostd::StringSubst(out,"&Gamma;","$\\Gamma$");aurostd::StringSubst(out,"&gamma;","$\\gamma$");
    aurostd::StringSubst(out,"&Delta;","$\\Delta$");aurostd::StringSubst(out,"&delta;","$\\delta$");
    aurostd::StringSubst(out,"&Epsilon;","$\\Εpsilon$");aurostd::StringSubst(out,"&epsilon;","$\\epsilon$");
    aurostd::StringSubst(out,"&Zeta;","$\\Ζeta$");aurostd::StringSubst(out,"&zeta;","$\\zeta$");
    aurostd::StringSubst(out,"&Eta;","$\\Eta$");aurostd::StringSubst(out,"&eta;","$\\eta$");
    aurostd::StringSubst(out,"&Theta;","$\\Theta$");aurostd::StringSubst(out,"&theta;","$\\theta$");
    aurostd::StringSubst(out,"&Iota;","$\\Ιiota$");aurostd::StringSubst(out,"&iota;","$\\iota$");
    aurostd::StringSubst(out,"&Kappa;","$\\Kappa$");aurostd::StringSubst(out,"&kappa;","$\\kappa$");
    aurostd::StringSubst(out,"&Lambda;","$\\Lambda$");aurostd::StringSubst(out,"&lambda;","$\\lambda$");
    aurostd::StringSubst(out,"&Mu;","$\\Mu$");aurostd::StringSubst(out,"&mu;","$\\mu$");
    aurostd::StringSubst(out,"&Nu;","$\\Νu$");aurostd::StringSubst(out,"&nu;","$\\nu$");
    aurostd::StringSubst(out,"&Xi;","$\\Xi$");aurostd::StringSubst(out,"&xi;","$\\xi$");
    aurostd::StringSubst(out,"&Omicron;","$\\Omicron$");aurostd::StringSubst(out,"&omicron;","$\\omicron$");
    aurostd::StringSubst(out,"&Pi;","$\\Pi$");aurostd::StringSubst(out,"&pi;","$\\pi$");
    aurostd::StringSubst(out,"&Rho;","$\\Rho$");aurostd::StringSubst(out,"&rho;","$\\rho$");
    aurostd::StringSubst(out,"&Sigma;","$\\Sigma$");aurostd::StringSubst(out,"&sigma;","$\\sigma$");
    aurostd::StringSubst(out,"&Tau;","$\\Tau$");aurostd::StringSubst(out,"&tau;","$\\tau$");
    aurostd::StringSubst(out,"&Upsilon;","$\\Upsilon$");aurostd::StringSubst(out,"&upsilon;","$\\upsilon$");
    aurostd::StringSubst(out,"&Phi;","$\\Phi$");aurostd::StringSubst(out,"&phi;","$\\phi$");
    aurostd::StringSubst(out,"&Chi;","$\\Chi$");aurostd::StringSubst(out,"&chi;","$\\chi$");
    aurostd::StringSubst(out,"&Psi;","$\\Psi$");aurostd::StringSubst(out,"&psi;","$\\psi$");
    aurostd::StringSubst(out,"&Omega;","$\\Omega$");aurostd::StringSubst(out,"&omega;","$\\omega$");
    aurostd::StringSubst(out,"&thetasym","$\\thetasym$");
    // FINAL
    aurostd::StringSubst(out,"&","\\&");

    return out;
  }

  string html2txt(const string& str) {
    string out=str;
    aurostd::StringSubst(out,"<sub>","");aurostd::StringSubst(out,"</sub>","");
    aurostd::StringSubst(out,"<i>","");aurostd::StringSubst(out,"</i>","");
    aurostd::StringSubst(out,"<b>",""); aurostd::StringSubst(out,"</b>","");
    aurostd::StringSubst(out,"MgB2","MgB2");
    aurostd::StringSubst(out,"&","&");
    aurostd::StringSubst(out,"_","");aurostd::StringSubst(out,"\\","");
    return out;
  }


  // ***************************************************************************
  // Function aurostd::string2latex
  // ***************************************************************************
  string string2latex(const string& str) {
    string out=str;
    aurostd::StringSubst(out,"_pv","_{pv}");aurostd::StringSubst(out,"_sv","_{sv}");aurostd::StringSubst(out,"_h","_{h}");
    aurostd::StringSubst(out,"_d","_{d}");aurostd::StringSubst(out,"_s","_{s}");
    aurostd::StringSubst(out,"_1","_{1}");aurostd::StringSubst(out,"_2","_{2}");aurostd::StringSubst(out,"_3","_{3}");
    return out;
  }

  // ***************************************************************************
  // Function aurostd::latex2html
  // ***************************************************************************
  string latex2html(const string& str) {
    string out=str;
    aurostd::StringSubst(out,"\\alpha","&alpha;");aurostd::StringSubst(out,"\\Alpha","&Alpha;");
    aurostd::StringSubst(out,"\\beta","&beta;");aurostd::StringSubst(out,"\\Beta","&Beta;");
    aurostd::StringSubst(out,"\\epsilon","&epsilon;");aurostd::StringSubst(out,"\\Epsilon","&Epsilon;");
    aurostd::StringSubst(out,"\\eta","&eta;");aurostd::StringSubst(out,"\\Eta","&Eta;");
    aurostd::StringSubst(out,"\\gamma","&gamma;");aurostd::StringSubst(out,"\\Gamma","&Gamma;");
    aurostd::StringSubst(out,"\\delta","&delta;");aurostd::StringSubst(out,"\\Delta","&Delta;");
    aurostd::StringSubst(out,"\\omega","&omega;");aurostd::StringSubst(out,"\\Omega","&Omega;");
    aurostd::StringSubst(out,"\\sigma","&sigma;");aurostd::StringSubst(out,"\\Sigma","&Sigma;");
    aurostd::StringSubst(out,"_{a}","<sub>a</sub>");aurostd::StringSubst(out,"_a","<sub>a</sub>"); 
    aurostd::StringSubst(out,"_{b}","<sub>b</sub>");aurostd::StringSubst(out,"_b","<sub>b</sub>");
    aurostd::StringSubst(out,"_{c}","<sub>d</sub>");aurostd::StringSubst(out,"_c","<sub>d</sub>");
    aurostd::StringSubst(out,"_{d}","<sub>d</sub>");aurostd::StringSubst(out,"_d","<sub>d</sub>");
    aurostd::StringSubst(out,"_{h}","<sub>h</sub>");aurostd::StringSubst(out,"_h","<sub>h</sub>");
    aurostd::StringSubst(out,"_{s}","<sub>s</sub>");aurostd::StringSubst(out,"_s","<sub>s</sub>");
    aurostd::StringSubst(out,"_{v}","<sub>v</sub>");aurostd::StringSubst(out,"_v","<sub>v</sub>");
    aurostd::StringSubst(out,"_{AB}","<sub>AB</sub>");
    aurostd::StringSubst(out,"_{AB2}","<sub>AB2</sub>");
    aurostd::StringSubst(out,"_{A2B2}","<sub>A2B2</sub>");
    aurostd::StringSubst(out,"_{AB3}","<sub>AB3</sub>"); 
    for(uint i=0;i<100;i++) aurostd::StringSubst(out,"_{"+aurostd::utype2string(i)+"}","<sub>"+aurostd::utype2string(i)+"</sub>");
    for(uint i=0;i<10;i++) aurostd::StringSubst(out,"_"+aurostd::utype2string(i)+"","<sub>"+aurostd::utype2string(i)+"</sub>"); // patch
    for(uint i1=0;i1<=3;i1++)
      for(uint i2=0;i2<=3;i2++)
        for(uint i3=0;i3<=3;i3++)
          aurostd::StringSubst(out,
              "^{["+aurostd::utype2string(i1)+aurostd::utype2string(i2)+aurostd::utype2string(i3)+"]}",
              "<sup>"+aurostd::utype2string(i1)+aurostd::utype2string(i2)+aurostd::utype2string(i3)+"</sup>");
    string s="AB";
    stringstream ss;
    for(uint i1=0;i1<=1;i1++)
      for(uint i2=0;i2<=1;i2++)
        for(uint i3=0;i3<=1;i3++)
          for(uint i4=0;i4<=1;i4++)
            for(uint i5=0;i5<=1;i5++) {
              aurostd::StringstreamClean(ss);
              ss << s.at(i1) << s.at(i2) << s.at(i3) << s.at(i4) << s.at(i5); 
              aurostd::StringSubst(out,"_{"+ss.str()+"}","<sub>"+ss.str()+"</sub>");
            }
    //    return out;
    //  string latex2html(const string& str) {  //[CO20200106 - close bracket for indenting]}
    // string out=str;
    aurostd::StringSubst(out,"\\&","&");
    aurostd::StringSubst(out,"MgB$_2$","MgB<sub>2</sub>");
    //  aurostd::StringSubst(out,"<sub>","$_{");aurostd::StringSubst(out,"</sub>","}$");
    //  aurostd::StringSubst(out,"<i>","\\textit{");aurostd::StringSubst(out,"</i>","}");
    // aurostd::StringSubst(out,"<b>","\\textbf{"); aurostd::StringSubst(out,"</b>","}");
    // aurostd::StringSubst(out,"&","\\&");
    //  aurostd::StringSubst(out,"Schuttler","Sch\\\"uttler");
    //  aurostd::StringSubst(out,"Csányi","Cs\\'anyi");aurostd::StringSubst(out,"Csanyi","Cs\\'anyi");
    // aurostd::StringSubst(out,"Rosen","Ros\\'en");
    // umlaut
    aurostd::StringSubst(out,"\\:a","&auml;");aurostd::StringSubst(out,"\\:A","&Auml;");
    aurostd::StringSubst(out,"\\:e","&euml;");aurostd::StringSubst(out,"\\:E","&Euml;");
    aurostd::StringSubst(out,"\\:i","&iuml;");aurostd::StringSubst(out,"\\:I","&Iuml;");
    aurostd::StringSubst(out,"\\:o","&ouml;");aurostd::StringSubst(out,"\\:O","&Ouml;");
    aurostd::StringSubst(out,"\\:u","&uuml;");aurostd::StringSubst(out,"\\:U","&Uuml;");
    // grave accent
    aurostd::StringSubst(out,"\\`a","&agrave;");aurostd::StringSubst(out,"\\`A","&Agrave;");
    aurostd::StringSubst(out,"\\`e","&egrave;");aurostd::StringSubst(out,"\\`E","&Egrave;");
    aurostd::StringSubst(out,"\\`i","&igrave;");aurostd::StringSubst(out,"\\`I","&Igrave;");
    aurostd::StringSubst(out,"\\`o","&ograve;");aurostd::StringSubst(out,"\\`O","&Ograve;");
    aurostd::StringSubst(out,"\\`u","&ugrave;");aurostd::StringSubst(out,"\\`U","&Ugrave;");
    // acute accent
    aurostd::StringSubst(out,"\\'a","&aacute;");aurostd::StringSubst(out,"\\'A","&Aacute;");
    aurostd::StringSubst(out,"\\'e","&eacute;");aurostd::StringSubst(out,"\\'E","&Eacute;");
    aurostd::StringSubst(out,"\\'i","&iacute;");aurostd::StringSubst(out,"\\'I","&Iacute;");
    aurostd::StringSubst(out,"\\'o","&oacute;");aurostd::StringSubst(out,"\\'O","&Oacute;");
    aurostd::StringSubst(out,"\\'u","&uacute;");aurostd::StringSubst(out,"\\'U","&Uacute;");
    // tilde
    aurostd::StringSubst(out,"\\~a","&atilde;");aurostd::StringSubst(out,"\\~A","&Atilde;");
    aurostd::StringSubst(out,"\\~e","&etilde;");aurostd::StringSubst(out,"\\~E","&Etilde;");
    aurostd::StringSubst(out,"\\~i","&itilde;");aurostd::StringSubst(out,"\\~I","&Itilde;");
    aurostd::StringSubst(out,"\\~o","&otilde;");aurostd::StringSubst(out,"\\~O","&Otilde;");
    aurostd::StringSubst(out,"\\~u","&utilde;");aurostd::StringSubst(out,"\\~U","&Utilde;");

    // caron
    aurostd::StringSubst(out,"\\v{z}","&zcaron;"); aurostd::StringSubst(out,"\\v{Z}","&Zcaron;");
    // slash
    aurostd::StringSubst(out,"\\o","&oslash;"); aurostd::StringSubst(out,"\\O","&Oslash;");

    return out;
  }

  string latex2txt(const string& str) {
    string out=str;
    aurostd::StringSubst(out,"\\&","&");
    aurostd::StringSubst(out,"MgB$_2$","MgB2");
    aurostd::StringSubst(out,"<sub>","");aurostd::StringSubst(out,"</sub>","");
    aurostd::StringSubst(out,"<i>","");aurostd::StringSubst(out,"</i>","");
    aurostd::StringSubst(out,"<b>","");aurostd::StringSubst(out,"</b>","");
    return out;
  }

  //CO20190419 - moved from chull
  string fixStringLatex(const string& input, bool double_back_slash,bool symmetry_string) {
    // deals with special characters for LaTeX, like some characters in prototype
    // see http://tex.stackexchange.com/questions/34580/escape-character-in-latex
    // double_back_slash was needed SOMETIMES for gnuplot output, as one backslash
    // went away when writing to file, and  -- OBSOLETE NOW
    string output;
    vector<char> problem_characters;
    problem_characters.push_back('&');
    problem_characters.push_back('%');
    problem_characters.push_back('$');
    problem_characters.push_back('#');
    if(!symmetry_string) {
      problem_characters.push_back('_');
      problem_characters.push_back('{');
      problem_characters.push_back('}');
    }
    problem_characters.push_back('~');  // different fix
    problem_characters.push_back('^');  // different fix
    string solution_string;
    solution_string = "\\\\";  // has to be string, \\ char does not work
    bool found_escaped_char;
    bool found_hyphen_symmetry = false;
    bool solved_hyphen_symmetry = false;
    for(uint i=0;i<input.length();i++) {
      // we first enter this loop because symmetry_string and input[i]=='-'
      // second enter loop because symmetry_string and found_hyphen_symmetry
      if(symmetry_string && (input[i] == '-' || found_hyphen_symmetry)) {
        if(!found_hyphen_symmetry) {
          // first enter loop, come here
          found_hyphen_symmetry = true;
          output.append((double_back_slash?string("\\"):string(""))+string("\\overline{"));
          // very important, we don't want to add hyphen, just replace
          // with overline, so continue
          continue;
        } else {
          // second enter loop, do nothing but turn this flag on
          // allow us to add input[i]
          found_hyphen_symmetry = false;
          solved_hyphen_symmetry = true;
        }
      } else {
        if(symmetry_string && solved_hyphen_symmetry) {
          // last step of symmetry_string fix, but we have to do this in part of
          // the loop to allow for next character to be identified as problem
          // character as well
          output.append(1, '}');
          solved_hyphen_symmetry = false;
        }
        // go through all problem characters
        for(uint j=0,fl_size_j=problem_characters.size();j<fl_size_j;j++) {
          if(input[i] == problem_characters[j]) {
            if(double_back_slash) {
              // if we find one, but it has double backslash, leave alone
              // doesn't matter what it is, if it has double backslash it's good
              // if we find one, but it only has single backslash, add one
              if(i && i - 1 && input[i - 1] == '\\' && input[i - 2] == '\\') {break;}
              else if(i && input[i - 1] == '\\') {
                output.append(1, '\\');  // just add one
                break;
              }
              // if we find one, give two backslashes
              output.append("\\\\");
              break;
            } else {
              // if we find one, but it has single backslash, leave alone
              // doesn't matter what it is, if it has single backslash it's good
              // if we find one, give single backslash
              if(i && input[i - 1] == '\\') {break;}  
              output.append(1, '\\');
              break;
            }
          }
        }
        // we also have to add {} for these characters
        if(input[i] == '~' || input[i] == '^') {output.append("{}");}
        found_escaped_char = false;
        if(input[i] == '\\') {
          for(uint j=0,fl_size_j=problem_characters.size();j<fl_size_j;j++) {
            // the only way this works if it's serving as an escape for a character
            // don't worry about double backslash here, we get to that when we find
            // the actual character
            if(i != (input.length() - 1) && input[i+1] == problem_characters[j]) {
              found_escaped_char = true;
              break;  // doesn't matter what it is, if it has backslash it's good
            }
          }
          // this is a problem, no way around it--we cannot output single backslash
          if(!found_escaped_char) {
            stringstream message;
            message << "Extraneous backslash found in \"" << input << "\" which may cause problems for LaTeX/gnuplot";
            //[moved from chull]pflow::logger(__AFLOW_FUNC__,message,m_aflags, *p_FileMESSAGE,*p_oss,_LOGGER_WARNING_);
            cerr << __AFLOW_FUNC__ << " ERROR - " << message.str() << endl;
            return input;
          }
        }
      }
      // add in character from input
      output.append(1, input[i]);
    }
    return output;
  }
}

// ***************************************************************************
// SORT WORLD
// ----------------------------------------------------------------------------
// sort for vector (starting from xvector)

namespace aurostd {
  template<class utype1> // function quicksort
    void sort(vector<utype1>& arr) {
      xvector<utype1> xarr = aurostd::vector2xvector(arr);
      aurostd::sort(xarr);
      arr = aurostd::xvector2vector(xarr);
    }

  template<class utype1,class utype2> // function quicksort
    void sort(vector<utype1>& arr, vector<utype2>& brr) {
    xvector<utype1> xarr = aurostd::vector2xvector(arr);
    xvector<utype2> xbrr = aurostd::vector2xvector(brr);
    aurostd::sort2(xarr.rows, xarr, xbrr);
    arr = aurostd::xvector2vector(xarr);
    brr = aurostd::xvector2vector(xbrr);
  }


  template<class utype1,class utype2> // function quicksort //CO20200915
    void sort(deque<utype1>& arr, deque<utype2>& brr) {
      xvector<utype1> xarr(arr.size());
      xvector<utype2> xbrr(brr.size());
      for(uint i=0;i<arr.size();i++) xarr[i+1]=arr[i];
      for(uint i=0;i<brr.size();i++) xbrr[i+1]=brr[i];
      aurostd::sort2(xarr.rows,xarr,xbrr);
      // aurostd::sort2(xarr,xbrr);
      arr.clear();brr.clear();
      for(int i=0;i<xarr.rows;i++) {
        arr.push_back(xarr[i+1]);
        brr.push_back(xbrr[i+1]);
      }
    }

  template<class utype1,class utype2,class utype3> // function quicksort
    void sort(vector<utype1>& arr, vector<utype2>& brr, vector<utype3>& crr) {
    xvector<utype1> xarr = aurostd::vector2xvector(arr);
    xvector<utype2> xbrr = aurostd::vector2xvector(brr);
    xvector<utype3> xcrr = aurostd::vector2xvector(crr);
    aurostd::sort3(xarr.rows,xarr,xbrr, xcrr);
    arr = aurostd::xvector2vector(xarr);
    brr = aurostd::xvector2vector(xbrr);
    crr = aurostd::xvector2vector(xcrr);
    }

  template<class utype1,class utype2,class utype3,class utype4> // function quicksort
    void sort(vector<utype1>& arr, vector<utype2>& brr, vector<utype3>& crr, vector<utype4>& drr) {
    xvector<utype1> xarr = aurostd::vector2xvector(arr);
    xvector<utype2> xbrr = aurostd::vector2xvector(brr);
    xvector<utype3> xcrr = aurostd::vector2xvector(crr);
    xvector<utype4> xdrr = aurostd::vector2xvector(drr);
    aurostd::sort4(xarr.rows,xarr,xbrr,xcrr, xdrr);
    arr = aurostd::xvector2vector(xarr);
    brr = aurostd::xvector2vector(xbrr);
    crr = aurostd::xvector2vector(xcrr);
    drr = aurostd::xvector2vector(xdrr);
    }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// sort for vector of strings
namespace aurostd {
  void sort(vector<string>& arg) {
    sort(arg.begin(),arg.end(),aurostd::_sort_string_());
  }
  void sort(deque<string>& arg) {
    std::sort(arg.begin(),arg.end(),aurostd::_sort_string_());
  }
  void rsort(vector<string>& arg) {
    std::reverse(arg.begin(),arg.end());//,aurostd::_sort_string_());
  }
  void rsort(deque<string>& arg) {
    std::reverse(arg.begin(),arg.end());//,aurostd::_sort_string_());
  }
}

// sort_remove_duplicates for vector of strings
namespace aurostd {
  void sort_remove_duplicates(vector<string>& arg) {
    sort(arg.begin(),arg.end(),aurostd::_sort_string_());
    arg.erase(std::unique(arg.begin(),arg.end()),arg.end());
  }
  void sort_remove_duplicates(deque<string>& arg) {
    std::sort(arg.begin(),arg.end(),aurostd::_sort_string_());
    arg.erase(std::unique(arg.begin(),arg.end()),arg.end());
  }
  void rsort_remove_duplicates(vector<string>& arg) {
    std::reverse(arg.begin(),arg.end());//,aurostd::_sort_string_());
    arg.erase(std::unique(arg.begin(),arg.end()),arg.end());
  }
  void rsort_remove_duplicates(deque<string>& arg) {
    std::reverse(arg.begin(),arg.end());//,aurostd::_sort_string_());
    arg.erase(std::unique(arg.begin(),arg.end()),arg.end());
  }
}


// ----------------------------------------------------------------------------
// sort for vector/deque of string_int
namespace aurostd {
  void sort(vector<string>& varg1,vector<int>& varg2) {
    vector<aurostd::_string_int_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];}
    sort(vv.begin(),vv.end(),_sort_string_int_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;}
  }
  void sort(deque<string>& varg1,deque<int>& varg2) {
    deque<aurostd::_string_int_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];}
    sort(vv.begin(),vv.end(),_sort_string_int_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of string_double
namespace aurostd {
  void sort(vector<string>& varg1,vector<double>& varg2) {
    vector<aurostd::_string_double_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];}
    sort(vv.begin(),vv.end(),_sort_string_double_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;}
  }
  void sort(deque<string>& varg1,deque<double>& varg2) {
    deque<aurostd::_string_double_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];}
    sort(vv.begin(),vv.end(),_sort_string_double_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of string_string
namespace aurostd {
  void sort(vector<string>& varg1,vector<string>& varg2) {
    vector<aurostd::_string_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];}
    sort(vv.begin(),vv.end(),_sort_string_string_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;}
  }
  void sort(deque<string>& varg1,deque<string>& varg2) {
    deque<aurostd::_string_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];}
    sort(vv.begin(),vv.end(),_sort_string_string_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;}
  }
}


// ----------------------------------------------------------------------------
// sort for vector/deque of double_int
// HERE THEY ARE

namespace aurostd {
  void sort(vector<double>& varg1,vector<int>& varg2) {
    vector<aurostd::_double_int_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];}
    sort(vv.begin(),vv.end(),_sort_double_int_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;}
  }
  void sort(deque<double>& varg1,deque<int>& varg2) {
    deque<aurostd::_double_int_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];}
    sort(vv.begin(),vv.end(),_sort_double_int_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of double_double
namespace aurostd {
  void sort(vector<double>& varg1,vector<double>& varg2) {
    vector<aurostd::_double_double_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];}
    sort(vv.begin(),vv.end(),_sort_double_double_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;}
  }
  void sort(deque<double>& varg1,deque<double>& varg2) {
    deque<aurostd::_double_double_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];}
    sort(vv.begin(),vv.end(),_sort_double_double_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of double_string
namespace aurostd {
  void sort(vector<double>& varg1,vector<string>& varg2) {
    vector<aurostd::_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];}
    sort(vv.begin(),vv.end(),_sort_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;}
  }
  void sort(deque<double>& varg1,deque<string>& varg2) {
    deque<aurostd::_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];}
    sort(vv.begin(),vv.end(),_sort_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of string_int_string
namespace aurostd {
  void sort(vector<string>& varg1,vector<int>& varg2,vector<string>& varg3) {
    vector<aurostd::_string_int_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];vv[i].arg3=varg3[i];}
    sort(vv.begin(),vv.end(),_sort_string_int_string_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;varg3[i]=vv[i].arg3;}
  }
  void sort(deque<string>& varg1,deque<int>& varg2,deque<string>& varg3) {
    deque<aurostd::_string_int_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];vv[i].arg3=varg3[i];}
    sort(vv.begin(),vv.end(),_sort_string_int_string_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;varg3[i]=vv[i].arg3;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of string_double_string
namespace aurostd {
  void sort(vector<string>& varg1,vector<double>& varg2,vector<string>& varg3) {
    vector<aurostd::_string_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];vv[i].arg3=varg3[i];}
    sort(vv.begin(),vv.end(),_sort_string_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;varg3[i]=vv[i].arg3;}
  }
  void sort(deque<string>& varg1,deque<double>& varg2,deque<string>& varg3) {
    deque<aurostd::_string_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];vv[i].arg3=varg3[i];}
    sort(vv.begin(),vv.end(),_sort_string_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;varg3[i]=vv[i].arg3;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of string_string_string
namespace aurostd {
  void sort(vector<string>& varg1,vector<string>& varg2,vector<string>& varg3) {
    vector<aurostd::_string_string_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];vv[i].arg3=varg3[i];}
    sort(vv.begin(),vv.end(),_sort_string_string_string_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;varg3[i]=vv[i].arg3;}
  }
  void sort(deque<string>& varg1,deque<string>& varg2,deque<string>& varg3) {
    deque<aurostd::_string_string_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];vv[i].arg3=varg3[i];}
    sort(vv.begin(),vv.end(),_sort_string_string_string_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;varg3[i]=vv[i].arg3;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of string_string_double_string
namespace aurostd {
  void sort(vector<string>& varg1,vector<string>& varg2,vector<double>& varg3,vector<string>& varg4) {
    vector<aurostd::_string_string_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];vv[i].arg3=varg3[i];vv[i].arg4=varg4[i];}
    sort(vv.begin(),vv.end(),_sort_string_string_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;varg3[i]=vv[i].arg3;varg4[i]=vv[i].arg4;}
  }
  void sort(deque<string>& varg1,deque<string>& varg2,deque<double>& varg3,deque<string>& varg4) {
    deque<aurostd::_string_string_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];vv[i].arg3=varg3[i];vv[i].arg4=varg4[i];}
    sort(vv.begin(),vv.end(),_sort_string_string_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;varg3[i]=vv[i].arg3;varg4[i]=vv[i].arg4;}
  }
}

// ----------------------------------------------------------------------------
// sort for vector/deque of string_string_double_double_string
namespace aurostd {
  void sort(vector<string>& varg1,vector<string>& varg2,vector<double>& varg3,vector<double>& varg4,vector<string>& varg5) {
    vector<aurostd::_string_string_double_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];vv[i].arg3=varg3[i];vv[i].arg4=varg4[i];vv[i].arg5=varg5[i];}
    sort(vv.begin(),vv.end(),_sort_string_string_double_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;varg3[i]=vv[i].arg3;varg4[i]=vv[i].arg4;varg5[i]=vv[i].arg5;}
  }
  void sort(deque<string>& varg1,deque<string>& varg2,deque<double>& varg3,deque<double>& varg4,deque<string>& varg5) {
    deque<aurostd::_string_string_double_double_string_> vv(varg1.size());
    for(uint i=0;i<varg1.size();i++) {vv[i].arg1=varg1[i];vv[i].arg2=varg2[i];vv[i].arg3=varg3[i];vv[i].arg4=varg4[i];vv[i].arg5=varg5[i];}
    sort(vv.begin(),vv.end(),_sort_string_string_double_double_string_());
    for(uint i=0;i<varg1.size();i++) {varg1[i]=vv[i].arg1;varg2[i]=vv[i].arg2;varg3[i]=vv[i].arg3;varg4[i]=vv[i].arg4;varg5[i]=vv[i].arg5;}
  }
}

// ----------------------------------------------------------------------------
// reorder vector //CO20221111
namespace aurostd {
  template<class utype> // function quicksort
    void reorder(vector<utype>& vec,vector<uint>& vorder,uint mode){	//CO20221111
			//algorithms and discussion from here: https://stackoverflow.com/questions/838384/reorder-vector-using-a-vector-of-indices (very good!)
			//solution by chmike
      //reorder a vector given input indices
      //there are two ways this can be done depending on what is inside vorder
      //input: vec={7,5,9,6}; vorder={1,3,0,2}
      //
      //mode 1: ``draw the elements of vector from the position of the indices''
      //result: {5,6,7,9}
      //NOTE: this is the default mode
      //
      //mode 2: ``move elements of vector to the position of the indices''
      //result: {9,7,6,5}
      //NOTE: this can also be accomplished with aurostd::sort(vorder,vec) but it requires vec to be C++ type or string
      //this function seems to run faster than aurostd::sort() as well
			uint i=0,j=0;
			if(mode==1){
				for(i=0;i<vec.size()-1;i++){
					if(vorder[i]==i){continue;}
					for(j=i+1;j<vorder.size();j++){
						if(vorder[j]==i){break;}
					}
					std::iter_swap(vec.begin()+i,vec.begin()+vorder[i]);
					std::iter_swap(vorder.begin()+i,vorder.begin()+j);
				}
        return;
      }
      else if(mode==2) {
				uint alt=0;
				// for all elements to put in place
				for(i=0;i<vec.size()-1;++i){
					// while the element i is not yet in place 
					while(i!=vorder[i]){
						// swap it with the element at its final place
						alt=vorder[i];
            std::iter_swap(vec.begin()+i,vec.begin()+alt);
            std::iter_swap(vorder.begin()+i,vorder.begin()+alt);
					}
				}
			}
			else{throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Unknown mode", _INPUT_ILLEGAL_);}
    }
}

// ***************************************************************************
// Function some statistical stuff
// combinations
// ***************************************************************************
namespace aurostd  {
  template<class utype> utype combinations(utype n,utype k) { // http://en.wikipedia.org/wiki/Combination // C^n_k=n!/k!(n-k)!   hard to calculate
    double cnk=1.0;
    for(utype i=0;i<=k-1;i++) cnk=cnk*(n-i)/(k-i);
    return (utype) cnk;
  }
  template<class utype> utype Cnk(utype n,utype k) { return combinations(n,k);}  // http://en.wikipedia.org/wiki/Combination
}

// ***************************************************************************
// aurostd::ShiftFirstColumn(const vector<vector<double> >& a, const double& value)
// ***************************************************************************
namespace aurostd  {
  vector<vector<double> > ShiftFirstColumn(const vector<vector<double> >& vva, const double& value) {
    //change value in the first column (usually menas energy in DOS)
    vector<vector<double> > vvb=vva;
    for (uint i=0; i<vvb.size(); i++) {
      vvb[i].at(0)=vvb[i].at(0) - value;
    }
    return vvb;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::ShrinkValuesExceptFirstColumn(const vector<vector<double> >& vva, const double& value)
// ***************************************************************************
namespace aurostd  {
  vector<vector<double> > ShrinkValuesExceptFirstColumn(const vector<vector<double> >& vva, const double& Fi) {
    //shrink Fis (usually means DOS Fis in DOS); Fi means probability
    vector<vector<double> > vvb=vva;
    for (uint i=0; i<vvb.size(); i++) {
      for (uint j=1; j<vvb[i].size();j++) {
        vvb[i][j]*=Fi;
      }
    }
    return vvb;
  }
} // namespace aurostd

// ***************************************************************************
// vector<vector<double> > aurostd::NormalizeAndSum3DVector(const vector<vector<vector<double> > >& vvva, const vector<vector<double> >& vFi)
// ***************************************************************************
namespace aurostd  {
  vector<vector<double> > NormalizeAndSum3DVector(const vector<vector<vector<double> > >& vvva, const vector<double>& vFi) {
    //normalize DOS and sum
    if(vvva.size()!=vFi.size()) {
      string message = "Vector sizes are not equal.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    vector<vector<double> > vvb, vv_tmp, vv_tmp_shrinked;
    vector<vector<vector<double> > > vvvc;
    double Fi;
    for (uint i=0; i<vvva.size();i++) {
      vv_tmp=vvva[i];
      Fi=vFi[i];
      vv_tmp_shrinked=aurostd::ShrinkValuesExceptFirstColumn(vv_tmp, Fi);
      vvvc.push_back(vv_tmp_shrinked);
    }
    vvb=aurostd::Sum3DVectorAndReduce2D(vvvc);
    return vvb;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::Sum3DVectorAndReduce2D(const vector<vector<vector<double> > >& vvva)
// ***************************************************************************
namespace aurostd  {
  vector<vector<double> > Sum3DVectorAndReduce2D(const vector<vector<vector<double> > >& vvva) {
    //The first column will not change! (For example, PDOS into TOTALPDOS)
    vector<vector<double> > vvtmp, vv_sum; 
    vv_sum=vvva.at(0);
    for (uint i=1; i<vvva.size();i++) {
      vvtmp=vvva[i];
      vv_sum=aurostd::Sum2DVectorExceptFirstColumn(vv_sum, vvtmp);
    }
    return vv_sum;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::Sum3DVectorAndReduce2D(const vector<vector<vector<double> > >& vvva)
// ***************************************************************************
namespace aurostd  {
  vector<vector<double> > Sum2DVectorExceptFirstColumn(const vector<vector<double> >& vva, const vector<vector<double> >& vvb) {
    if((vva.size()!=vvb.size()) && (vva.at(0).size() != vvb.at(0).size())) {
      string message = "Vector sizes are not equal.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }

    vector<vector<double> > vv_sum; vv_sum.resize(vva.size());
    for (uint i=0; i<vva.size(); i++) {
      int N=vva[i].size();
      vv_sum[i].resize(N);
    }

    for (uint i=0; i<vva.size();i++) {
      vv_sum[i][0]=vva[i].at(0);
      for (uint j=1; j<vva[i].size();j++) {
        vv_sum[i][j]=vva[i][j] + vvb[i][j];
      }
    }
    return vv_sum;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::ReduceVector(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd  {
  vector<vector<double> > ReduceVector(const vector<vector<double> >& vva, const int& n) {
    //Pick up the first (begin from 0) and the nth column of 2D vector
    vector<vector<double> > vvb; vvb.clear();
    vector<double> vtmp;
    for (uint i=0; i<vva.size();i++) {
      vtmp.clear();
      vtmp.push_back(vva[i].at(0));
      vtmp.push_back(vva[i].at(n));
      vvb.push_back(vtmp);
    }
    return vvb;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::CalculateIntegrate(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd  {
  double CalculateIntegrate(const vector<vector<double> >& vva, const int& n) {
    //Calculate integration of vva, the 0st column is x0, x1..., the n column is y1, y2 ...
    //begin from 0
    vector<vector<double> > vvb=aurostd::ReduceVector(vva, n);
    return aurostd::CalculateIntegrate(vvb);
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::CalculateIntegrate(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd  {
  double CalculateIntegrate(const vector<vector<double> >& vva, const int& n, const double& Emin, const double& Emax) {
    //Calculate integration of vva, the 0st column is x0, x1..., the n column is y1, y2 ...
    //begin from 0
    vector<vector<double> > vvb=aurostd::ReduceVector(vva, n);
    return aurostd::CalculateIntegrate(vvb, Emin, Emax);
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::CalculateIntegrate(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd  {
  double CalculateIntegrate(const vector<vector<double> >& vva) {
    double Emin=-100; double Emax=0.0; //default setting
    return aurostd::CalculateIntegrate(vva, Emin, Emax);
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::CalculateIntegrate(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd  {
  double CalculateIntegrate(const vector<vector<double> >& vva, const double& Emin, const double& Emax) {
    //Integral function
    //format of vva: x0, y0; x1, y1; x2, y2
    double integral_result=0.0;
    double area_tmp =0.0;
    double xbeg, xend, ybeg, yend;
    for (uint i=0; i<vva.size()-1;i++) {
      xbeg=vva[i].at(0); xend=vva.at(i+1).at(0);
      ybeg=vva[i].at(1); yend=vva.at(i+1).at(1);
      if(xbeg >= Emin && xend <= Emax) {
        area_tmp=0.5*(ybeg + yend)*(xend - xbeg);
        integral_result += area_tmp;
      }
    }
    return integral_result;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::vector2string(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd  {
  string vector2string(const vector<vector<double> >& vva) {
    stringstream ss_vva; aurostd::StringstreamClean(ss_vva);
    ss_vva << std::scientific;
    for (uint i=0; i<vva.size();i++) {
      for (uint j=0; j<vva[i].size();j++) {
        ss_vva << vva[i][j] << "   ";
      }
      ss_vva << endl;
    }
    return ss_vva.str();
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::vector2deque(const vector<utype>& vin)
// ***************************************************************************
//CO20181226
namespace aurostd  {
  template<class utype> deque<utype> vector2deque(const vector<utype>& vin){
    deque<utype> dout;
    for(uint i=0;i<vin.size();i++){dout.push_back(vin[i]);}
    return dout;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::vector2deque(const vector<utype>& vin)
// ***************************************************************************
//CO20181226
namespace aurostd  {
  template<class utype> vector<utype> deque2vector(const deque<utype>& din){
    vector<utype> vout;
    for(uint i=0;i<din.size();i++){vout.push_back(din[i]);}
    return vout;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::FindMaxIn2DvectorExcept1stColumn(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd  {
  double FindMaxIn2DvectorExcept1stColumn(const vector<vector<double> >& vva) {
    double min=-10;  //default
    double max=10;
    return aurostd::FindMaxIn2DvectorExcept1stColumn(vva, min, max);
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::FindMaxIn2DvectorExcept1stColumn(const vector<vector<double>& vva, const double& min, const double& max)
// ***************************************************************************
namespace aurostd  {
  double FindMaxIn2DvectorExcept1stColumn(const vector<vector<double> >& vva, const double& min, const double& max) {
    double max_value=0.0;
    for (uint i=0; i<vva.size();i++) {
      double E_tmp=vva[i].at(0);
      if(E_tmp >= min && E_tmp <= max) {
        for (uint j=1; j<vva[i].size();j++) {
          double db_tmp=vva[i][j];
          if(abs(db_tmp) > max_value) max_value=abs(db_tmp);
        }
      }
    }
    return max_value;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::FindMaxInTDOS(const vector<vector<double> >& vva, const double& min, const double& max)
// ***************************************************************************
namespace aurostd  {
  double FindMaxInTDOS(const vector<vector<double> >& vva, const double& min, const double& max) {
    double max_value=0.0;
    for (uint i=0; i<vva.size();i++) {
      double E_tmp=vva[i].at(0);
      if(E_tmp >= min && E_tmp <= max) {
        int column_max=0; // some default
        if(vva.at(0).size()==3) column_max=2; //get rid of the sum of TDOS
        if(vva.at(0).size()==5) column_max=3;
        for (int j=1; j<column_max;j++) {
          double db_tmp=vva[i][j];
          if(abs(db_tmp) > max_value) max_value=db_tmp;
        }
      }
    }
    return max_value;
  }
} // namespace aurostd


namespace aurostd {
  //***************************************************************************//
  // aurostd::joinWDelimiter(vector<uint>& uientries,const stringstream&
  // delimiter,const stringstream& m_delimiter,const stringstream& l_delimiter)
  //***************************************************************************//
  // joinWDelimiters int/uint type of objects together by a delimiter
  // no point for double objects, faster to just do it on the spot with
  // setprecision,fixed, etc.
  // m_delimiter is used if input is exactly length 2
  // l_delimiter otherwise
  // ME20220324 - added missing uint variant for xvector
  string joinWDelimiter(const xvector<int>& ientries, const char& _delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const xvector<int>& ientries, const char& _delimiter,
      const char& _l_delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const xvector<int>& ientries, const char& _delimiter,
      const char& _m_delimiter, const char& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(ientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const xvector<int>& ientries, const string& _delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const xvector<int>& ientries, const string& _delimiter,
      const string& _l_delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const xvector<int>& ientries, const string& _delimiter,
      const string& _m_delimiter, const string& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(ientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const xvector<int>& ientries, const stringstream& delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const xvector<int>& ientries, const stringstream& delimiter,
      const stringstream& l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const xvector<int>& ientries, const stringstream& delimiter,
      const stringstream& m_delimiter,
      const stringstream& l_delimiter) {
    stringstream output;
    string delim = delimiter.str();
    string mDelim = m_delimiter.str();
    string lDelim = l_delimiter.str();

    if (ientries.rows > 2) {
      for (int i =ientries.lrows; i <= ientries.urows; i++) {
        output << ientries[i];
        if (i == ientries.urows - 1) {  //CO20180216 - added -1
          output << lDelim;
        } else if (i !=ientries.urows) {
          output << delim;
        }
      }
    } else {
      for (int i = ientries.lrows; i <= ientries.urows; i++) {
        output << ientries[i];
        if (i == ientries.urows - 1) {  //CO20180216 - added -1
          output << mDelim;
        } else if (i != ientries.urows) {
          output << delim;
        }
      }
    }
    return output.str();
  }
  string joinWDelimiter(const xvector<uint>& ientries, const char& _delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const xvector<uint>& ientries, const char& _delimiter,
      const char& _l_delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const xvector<uint>& ientries, const char& _delimiter,
      const char& _m_delimiter, const char& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(ientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const xvector<uint>& ientries, const string& _delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const xvector<uint>& ientries, const string& _delimiter,
      const string& _l_delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const xvector<uint>& ientries, const string& _delimiter,
      const string& _m_delimiter, const string& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(ientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const xvector<uint>& ientries, const stringstream& delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const xvector<uint>& ientries, const stringstream& delimiter,
      const stringstream& l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const xvector<uint>& ientries, const stringstream& delimiter,
      const stringstream& m_delimiter,
      const stringstream& l_delimiter) {
    stringstream output;
    string delim = delimiter.str();
    string mDelim = m_delimiter.str();
    string lDelim = l_delimiter.str();

    if (ientries.rows > 2) {
      for (int i =ientries.lrows; i <= ientries.urows; i++) {
        output << ientries[i];
        if (i == ientries.urows - 1) {  //CO20180216 - added -1
          output << lDelim;
        } else if (i !=ientries.urows) {
          output << delim;
        }
      }
    } else {
      for (int i = ientries.lrows; i <= ientries.urows; i++) {
        output << ientries[i];
        if (i == ientries.urows - 1) {  //CO20180216 - added -1
          output << mDelim;
        } else if (i != ientries.urows) {
          output << delim;
        }
      }
    }
    return output.str();
  }
  string joinWDelimiter(const vector<int>& ientries, const char& _delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const vector<int>& ientries, const char& _delimiter,
      const char& _l_delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const vector<int>& ientries, const char& _delimiter,
      const char& _m_delimiter, const char& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(ientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<int>& ientries, const string& _delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const vector<int>& ientries, const string& _delimiter,
      const string& _l_delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const vector<int>& ientries, const string& _delimiter,
      const string& _m_delimiter, const string& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(ientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<int>& ientries, const stringstream& delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const vector<int>& ientries, const stringstream& delimiter,
      const stringstream& l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<int>& ientries, const stringstream& delimiter,
      const stringstream& m_delimiter,
      const stringstream& l_delimiter) {
    stringstream output;
    string delim = delimiter.str();
    string mDelim = m_delimiter.str();
    string lDelim = l_delimiter.str();

    if (ientries.size() > 2) {
      for (uint i = 0; i < ientries.size(); i++) {
        output << ientries[i];
        if (i == ientries.size() - 2) {
          output << lDelim;
        } else if (i != ientries.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (uint i = 0; i < ientries.size(); i++) {
        output << ientries[i];
        if (i == ientries.size() - 2) {
          output << mDelim;
        } else if (i != ientries.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
  string joinWDelimiter(const vector<uint>& uientries, const char& _delimiter) {
    return joinWDelimiter(uientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const vector<uint>& uientries, const char& _delimiter,
      const char& _l_delimiter) {
    return joinWDelimiter(uientries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const vector<uint>& uientries, const char& _delimiter,
      const char& _m_delimiter, const char& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(uientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<uint>& uientries, const string& _delimiter) {
    return joinWDelimiter(uientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const vector<uint>& uientries, const string& _delimiter,
      const string& _l_delimiter) {
    return joinWDelimiter(uientries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const vector<uint>& uientries, const string& _delimiter,
      const string& _m_delimiter, const string& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(uientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<uint>& uientries, const stringstream& delimiter) {
    return joinWDelimiter(uientries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const vector<uint>& uientries, const stringstream& delimiter,
      const stringstream& l_delimiter) {
    return joinWDelimiter(uientries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<uint>& uientries, const stringstream& delimiter,
      const stringstream& m_delimiter,
      const stringstream& l_delimiter) {
    stringstream output;
    string delim = delimiter.str();
    string mDelim = m_delimiter.str();
    string lDelim = l_delimiter.str();
    if (uientries.size() > 2) {
      for (uint i = 0; i < uientries.size(); i++) {
        output << uientries[i];
        if (i == uientries.size() - 2) {
          output << lDelim;
        } else if (i != uientries.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (uint i = 0; i < uientries.size(); i++) {
        output << uientries[i];
        if (i == uientries.size() - 2) {
          output << mDelim;
        } else if (i != uientries.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
} // namespace aurostd

namespace aurostd {
  //***************************************************************************//
  // aurostd::joinWDelimiter(vector<string>& _sentries,const stringstream&
  // delimiter,const stringstream& m_delimiter,const stringstream& l_delimiter)
  //***************************************************************************//
  // joinWDelimiters string type of objects together by a delimiter
  // m_delimiter is used if input is exactly length 2
  // l_delimiter otherwise
  string joinWDelimiter(const vector<string>& _sentries, const char& _delimiter) {
    return joinWDelimiter(_sentries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const vector<string>& _sentries, const char& _delimiter,
      const char& _l_delimiter) {
    return joinWDelimiter(_sentries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const vector<string>& _sentries, const char& _delimiter,
      const char& _m_delimiter, const char& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(_sentries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<string>& _sentries, const string& _delimiter) {
    return joinWDelimiter(_sentries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const vector<string>& _sentries, const string& _delimiter,
      const string& _l_delimiter) {
    return joinWDelimiter(_sentries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const vector<string>& _sentries, const string& _delimiter,
      const string& _m_delimiter, const string& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(_sentries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<string>& _sentries, const stringstream& delimiter) {
    return joinWDelimiter(_sentries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const vector<string>& _sentries, const stringstream& delimiter,
      const stringstream& l_delimiter) {
    return joinWDelimiter(_sentries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<string>& _sentries, const stringstream& delimiter,
      const stringstream& m_delimiter,
      const stringstream& l_delimiter) {
    stringstream output;
    vector<string> sentries;
    string delim = delimiter.str();
    string mDelim = m_delimiter.str();
    string lDelim = l_delimiter.str();
    // go through once to eliminate empty strings
    for (uint i = 0; i < _sentries.size(); i++) {
      if (_sentries[i].length()) {
        sentries.push_back(_sentries[i]);
      }
    }
    if (sentries.size() > 2) {
      for (uint i = 0; i < sentries.size(); i++) {
        output << sentries[i];
        if (i == sentries.size() - 2) {
          output << lDelim;
        } else if (i != sentries.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (uint i = 0; i < sentries.size(); i++) {
        output << sentries[i];
        if (i == sentries.size() - 2) {
          output << mDelim;
        } else if (i != sentries.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
} // namespace aurostd

namespace aurostd {
  //***************************************************************************//
  // aurostd::joinWDelimiter(deque<uint>& uientries,const stringstream&
  // delimiter,const stringstream& m_delimiter,const stringstream& l_delimiter)
  //***************************************************************************//
  // joinWDelimiters int/uint type of objects together by a delimiter
  // no point for double objects, faster to just do it on the spot with
  // setprecision,fixed, etc.
  // m_delimiter is used if input is exactly length 2
  // l_delimiter otherwise
  string joinWDelimiter(const deque<int>& ientries, const char& _delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const deque<int>& ientries, const char& _delimiter,
      const char& _l_delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const deque<int>& ientries, const char& _delimiter,
      const char& _m_delimiter, const char& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(ientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<int>& ientries, const string& _delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const deque<int>& ientries, const string& _delimiter,
      const string& _l_delimiter) {
    return joinWDelimiter(ientries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const deque<int>& ientries, const string& _delimiter,
      const string& _m_delimiter, const string& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(ientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<int>& ientries, const stringstream& delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const deque<int>& ientries, const stringstream& delimiter,
      const stringstream& l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<int>& ientries, const stringstream& delimiter,
      const stringstream& m_delimiter,
      const stringstream& l_delimiter) {
    stringstream output;
    string delim = delimiter.str();
    string mDelim = m_delimiter.str();
    string lDelim = l_delimiter.str();
    if (ientries.size() > 2) {
      for (uint i = 0; i < ientries.size(); i++) {
        output << ientries[i];
        if (i == ientries.size() - 2) {
          output << lDelim;
        } else if (i != ientries.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (uint i = 0; i < ientries.size(); i++) {
        output << ientries[i];
        if (i == ientries.size() - 2) {
          output << mDelim;
        } else if (i != ientries.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
  string joinWDelimiter(const deque<uint>& uientries, const char& _delimiter) {
    return joinWDelimiter(uientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const deque<uint>& uientries, const char& _delimiter,
      const char& _l_delimiter) {
    return joinWDelimiter(uientries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const deque<uint>& uientries, const char& _delimiter,
      const char& _m_delimiter, const char& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(uientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<uint>& uientries, const string& _delimiter) {
    return joinWDelimiter(uientries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const deque<uint>& uientries, const string& _delimiter,
      const string& _l_delimiter) {
    return joinWDelimiter(uientries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const deque<uint>& uientries, const string& _delimiter,
      const string& _m_delimiter, const string& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(uientries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<uint>& uientries, const stringstream& delimiter) {
    return joinWDelimiter(uientries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const deque<uint>& uientries, const stringstream& delimiter,
      const stringstream& l_delimiter) {
    return joinWDelimiter(uientries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<uint>& uientries, const stringstream& delimiter,
      const stringstream& m_delimiter,
      const stringstream& l_delimiter) {
    stringstream output;
    string delim = delimiter.str();
    string mDelim = m_delimiter.str();
    string lDelim = l_delimiter.str();
    if (uientries.size() > 2) {
      for (uint i = 0; i < uientries.size(); i++) {
        output << uientries[i];
        if (i == uientries.size() - 2) {
          output << lDelim;
        } else if (i != uientries.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (uint i = 0; i < uientries.size(); i++) {
        output << uientries[i];
        if (i == uientries.size() - 2) {
          output << mDelim;
        } else if (i != uientries.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
}

namespace aurostd {
  //***************************************************************************//
  // aurostd::joinWDelimiter(deque<string>& _sentries,const stringstream&
  // delimiter,const stringstream& m_delimiter,const stringstream& l_delimiter)
  //***************************************************************************//
  // joinWDelimiters string type of objects together by a delimiter
  // m_delimiter is used if input is exactly length 2
  // l_delimiter otherwise
  string joinWDelimiter(const deque<string>& _sentries, const char& _delimiter) {
    return joinWDelimiter(_sentries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const deque<string>& _sentries, const char& _delimiter,
      const char& _l_delimiter) {
    return joinWDelimiter(_sentries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const deque<string>& _sentries, const char& _delimiter,
      const char& _m_delimiter, const char& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(_sentries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<string>& _sentries, const string& _delimiter) {
    return joinWDelimiter(_sentries, _delimiter, _delimiter, _delimiter);
  }
  string joinWDelimiter(const deque<string>& _sentries, const string& _delimiter,
      const string& _l_delimiter) {
    return joinWDelimiter(_sentries, _delimiter, _delimiter, _l_delimiter);
  }
  string joinWDelimiter(const deque<string>& _sentries, const string& _delimiter,
      const string& _m_delimiter, const string& _l_delimiter) {
    stringstream delimiter, m_delimiter, l_delimiter;
    delimiter << _delimiter;
    m_delimiter << _m_delimiter;
    l_delimiter << _l_delimiter;
    return joinWDelimiter(_sentries, delimiter, m_delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<string>& _sentries, const stringstream& delimiter) {
    return joinWDelimiter(_sentries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const deque<string>& _sentries, const stringstream& delimiter,
      const stringstream& l_delimiter) {
    return joinWDelimiter(_sentries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<string>& _sentries, const stringstream& delimiter,
      const stringstream& m_delimiter,
      const stringstream& l_delimiter) {
    stringstream output;
    vector<string> sentries;  // no point working with deque
    string delim = delimiter.str();
    string mDelim = m_delimiter.str();
    string lDelim = l_delimiter.str();
    // go through once to eliminate empty strings
    for (uint i = 0; i < _sentries.size(); i++)
      //DX - Should not be an "!"; we want it to push back if it has a length [OBSOLETE] if (!_sentries[i].length())
    { //CO20200106 - patching for auto-indenting
      if (_sentries[i].length()) {
        sentries.push_back(_sentries[i]);
      }
    }
    if (sentries.size() > 2) {
      for (uint i = 0; i < sentries.size(); i++) {
        output << sentries[i];
        if (i == sentries.size() - 2) {
          output << lDelim;
        } else if (i != sentries.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (uint i = 0; i < sentries.size(); i++) {
        output << sentries[i];
        if (i == sentries.size() - 2) {
          output << mDelim;
        } else if (i != sentries.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
}

namespace aurostd {
  string wrapString(const string& input,const string& wrapper){return wrapString(input,wrapper,wrapper);}
  string wrapString(const string& input,const string& wrapper_start,const string& wrapper_end){
    if(input.empty()){return input;}
    return wrapper_start+input+wrapper_end;
  }
}

//DX20180118 START: XCOMPLEX TO JSON
namespace aurostd {
  //***************************************************************************//
  // aurostd::xcomplex2json
  //***************************************************************************//
  template<typename utype> string _xcomplex2json(xcomplex<utype>& number){
    string eendl="";
    bool roff=true; //round off
    stringstream sss;
    stringstream sscontent_json;
    vector<string> vcontent_json;
    // real
    sscontent_json << "\"real\":\"" << aurostd::utype2string(number.re,5,roff) << "\"" << eendl;
    vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
    // imaginary
    sscontent_json << "\"imag\":\"" << aurostd::utype2string(number.im,5,roff) << "\"" << eendl;
    vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

    sss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}" << eendl;
    return sss.str();    
  }
}

//Need to initalize 
namespace aurostd {
  string xcomplex2json(xcomplex<double>& number){ return _xcomplex2json(number); }
}

//DX20180118 END: XCOMPLEX TO JSON

//DX20170803 START: Matrix to JSON
namespace aurostd {
  //***************************************************************************//
  // aurostd::xmatDouble2String(xmatrix<double>& xmat_in)
  //***************************************************************************//
  // converts xmatrix<double> to json string
  // [OBSOLETE] string xmatDouble2String(const xmatrix<double>& xmat_in, bool roff){
  // [OBSOLETE]   stringstream output;
  // [OBSOLETE]   vector<string> rows;
  // [OBSOLETE]   for(uint i=1;i<(uint)xmat_in.rows+1;i++){
  // [OBSOLETE]     stringstream row;
  // [OBSOLETE]     xvector<double> xvec = xmat_in(i); //DX20170822 - added roundoff
  // [OBSOLETE]     if(roff){ xvec = roundoff(xvec,1e-8);} //DX20170822 - added roundoff
  // [OBSOLETE]     row << "[" << joinWDelimiter(xvecDouble2vecString(xvec),",") << "]";
  // [OBSOLETE]     rows.push_back(row.str());
  // [OBSOLETE]   }
  // [OBSOLETE]   output << joinWDelimiter(rows,",");
  // [OBSOLETE]   return output.str();
  // [OBSOLETE] }
  string xmatDouble2String(const xmatrix<double>& xmat_in, int precision, bool roff, double tol, char FORMAT){
    stringstream output;
    vector<string> rows;
    for(int i=xmat_in.lrows;i<=xmat_in.urows;i++){ //DX20180323 - fixed typo for initial index "int i=1" not "int i=xmat_in.urows" //ME20220324 - changed to lrows
      stringstream row;
      xvector<double> xvec = xmat_in(i); //DX20170822 - added roundoff
      //if(roff){ xvec = roundoff(xvec,tol);} //DX20170822 - added roundoff
      row << "[" << joinWDelimiter(xvecDouble2vecString(xvec,precision,roff,tol,FORMAT),",") << "]";
      rows.push_back(row.str());
      //cerr << i << "row.str(): " << row.str() << endl;
    }
    output << joinWDelimiter(rows,",");
    return output.str();
  }

  //ME20220324
  template <typename utype>
    string xmat2String(const xmatrix<utype>& xmat_in) {
      vector<string> rows;
      for (int i = xmat_in.lrows; i <= xmat_in.urows; i++) {
        rows.push_back("[" + joinWDelimiter(xmat_in(i), ",") + "]");
      }
      return joinWDelimiter(rows, ",");
    }
  template string xmat2String(const xmatrix<int>&);
  template string xmat2String(const xmatrix<uint>&);
}
//DX20170803 START: Matrix to END

namespace aurostd {
  //***************************************************************************//
  // aurostd::vecDouble2vecString(vector<double>& vin,int precision)
  //***************************************************************************//
  // converts vector<double> to vector<string> with precision
  // also works for xvectors and deques
  // [OBSOLETE] vector<string> vecDouble2vecString(const vector<double>& vin, bool roff) {
  // [OBSOLETE]   vector<string> vout;
  // [OBSOLETE]   for(uint i=0;i<vin.size();i++){
  // [OBSOLETE]     double tmp = vin[i]; //DX20170822 - add roundoff
  // [OBSOLETE]     if(roff){ tmp=scalar_roundoff(tmp,1e-8); } //DX20170822 - add roundoff
  // [OBSOLETE]     vout.push_back(aurostd::utype2string(tmp)); //DX20170822 - add roundoff
  // [OBSOLETE]   }
  // [OBSOLETE]   return vout;
  // [OBSOLETE] }
  vector<string> vecDouble2vecString(const vector<double>& vin,int precision, bool roff, double tol, char FORMAT) {
    vector<string> vout;
    for(uint i=0;i<vin.size();i++){
      //double tmp = vin[i]; //DX20170822 - add roundoff
      //if(roff){ tmp=roundoff(tmp,tol); } //DX20170822 - add roundoff
      vout.push_back(aurostd::utype2string(vin[i],precision,roff,tol,FORMAT)); //DX20170822 - add roundoff
    }
    return vout;
  }
  string vecDouble2String(const vector<double>& vin,int precision, bool roff, double tol, char FORMAT) {
    return aurostd::joinWDelimiter(vecDouble2vecString(vin, precision, roff, tol, FORMAT), ",");
  }
  // [OBSOLETE] vector<string> xvecDouble2vecString(const xvector<double>& vin, bool roff) {
  // [OBSOLETE]   vector<string> vout;
  // [OBSOLETE]   for(uint i=1;i<(uint)vin.rows+1;i++){
  // [OBSOLETE]     double tmp = vin(i); //DX20170822 - add roundoff
  // [OBSOLETE]    if(roff){ tmp=scalar_roundoff(tmp,1e-8); } //DX20170822 - add roundoff
  // [OBSOLETE]     vout.push_back(aurostd::utype2string(tmp)); //DX20170822 - add roundoff
  // [OBSOLETE]   }
  // [OBSOLETE]   return vout;
  // [OBSOLETE] }
  vector<string> xvecDouble2vecString(const xvector<double>& vin,int precision, bool roff, double tol, char FORMAT) {
    vector<string> vout;
    for(int i=vin.lrows;i<=vin.urows;i++){
      //double tmp = vin(i); //DX20170822 - add roundoff
      //if(roff){ tmp=roundoff(tmp,tol); } //DX20170822 - add roundoff
      vout.push_back(aurostd::utype2string(vin[i],precision,roff,tol,FORMAT)); //DX20170822 - add roundoff
    }
    return vout;
  }
  string xvecDouble2String(const xvector<double>& vin,int precision, bool roff, double tol, char FORMAT) {
    return aurostd::joinWDelimiter(xvecDouble2vecString(vin, precision, roff, tol, FORMAT), ",");
  }
  // [OBSOLETE] deque<string> deqDouble2deqString(const deque<double>& vin, bool roff) {
  // [OBSOLETE]   deque<string> vout;
  // [OBSOLETE]   for(uint i=0;i<vin.size();i++){
  // [OBSOLETE]     double tmp = vin[i]; //DX20170822 - add roundoff
  // [OBSOLETE]     if(roff){ tmp=scalar_roundoff(tmp,1e-8); } //DX20170822 - add roundoff
  // [OBSOLETE]     vout.push_back(aurostd::utype2string(tmp)); //DX20170822 - add roundoff
  // [OBSOLETE]   }
  // [OBSOLETE]   return vout;
  // [OBSOLETE] }
  // [OBSOLETE]  deque<string> deqDouble2deqString(const deque<double>& vin,int precision, bool roff, double tol, char FORMAT)  // USE OVERLOADING
  deque<string> vecDouble2vecString(const deque<double>& vin,int precision, bool roff, double tol, char FORMAT) { //SC20200330
    deque<string> vout;
    for(uint i=0;i<vin.size();i++){
      //double tmp = vin[i]; //DX20170822 - add roundoff
      //if(roff){ tmp=roundoff(tmp,tol); } //DX20170822 - add roundoff
      vout.push_back(aurostd::utype2string(vin[i],precision,roff,tol,FORMAT)); //DX20170822 - add roundoff
    }
    return vout;
  }
  string vecDouble2String(const deque<double>& vin,int precision, bool roff, double tol, char FORMAT) {
    return aurostd::joinWDelimiter(vecDouble2vecString(vin, precision, roff, tol, FORMAT), ",");
  }
}

namespace aurostd {
  //***************************************************************************//
  // aurostd::wrapVecEntries(vector<string>& vin,string wrap)
  //***************************************************************************//
  // individually wraps entries of vector with specified string
  // converts <a,b,c> to <'a','b','c'>
  // also works for deques
  vector<string> wrapVecEntries(const vector<string>& vin,const string& wrap){
    return wrapVecEntries(vin,wrap,wrap);
  }
  vector<string> wrapVecEntries(const vector<string>& vin,const string& wrap_start,const string& wrap_end){
    vector<string> vout;
    for(uint i=0;i<vin.size();i++){
      if(vin[i].length()){
        vout.push_back(wrap_start+vin[i]+wrap_end);
      }
    }
    return vout;
  }
  deque<string> wrapVecEntries(const deque<string>& vin,const string& wrap){
    return wrapVecEntries(vin,wrap,wrap);
  }
  deque<string> wrapVecEntries(const deque<string>& vin,const string& wrap_start,const string& wrap_end){
    deque<string> vout;
    for(uint i=0;i<vin.size();i++){
      if(vin[i].length()){
        vout.push_back(wrap_start+vin[i]+wrap_end);
      }
    }
    return vout;
  }
}

//base64 stuff
//CO START
namespace aurostd {
  // ***************************************************************************
  // aurostd::isBase64(unsigned char c)
  // ***************************************************************************
  // determines if char is base64
  // http://www.adp-gmbh.ch/cpp/common/base64.html
  //static inline bool isBase64(unsigned char c)
  inline bool isBase64(unsigned char c)
  { //CO20200106 - patching for auto-indenting
    return (isalnum(c) || (c == '+') || (c == '/'));
  }

  // ***************************************************************************
  // aurostd::base64Encoder(unsigned char const* bytes_to_encode, unsigned int in_len)
  // ***************************************************************************
  // encodes bytes to base64
  // http://www.adp-gmbh.ch/cpp/common/base64.html
  std::string base64Encoder(unsigned char const* bytes_to_encode, unsigned int in_len) {
    std::string ret;
    int i = 0;
    int j = 0;
    unsigned char char_array_3[3];
    unsigned char char_array_4[4];

    while (in_len--) {
      char_array_3[i++] = *(bytes_to_encode++);
      if (i == 3) {
        char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
        char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
        char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
        char_array_4[3] = char_array_3[2] & 0x3f;

        for(i = 0; (i <4) ; i++) {
          ret += base64_chars[char_array_4[i]];
        }
        i = 0;
      }
    }

    if (i) {
      for(j = i; j < 3; j++) {
        char_array_3[j] = '\0';
      }

      char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
      char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
      char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
      char_array_4[3] = char_array_3[2] & 0x3f;

      for (j = 0; (j < i + 1); j++) {
        ret += base64_chars[char_array_4[j]];
      }

      while((i++ < 3)) {
        ret += '=';
      }
    }

    return ret;
  }

  // ***************************************************************************
  // aurostd::base64Decoder(std::string const& encoded_string)
  // ***************************************************************************
  // decodes base64 to bytes
  // http://www.adp-gmbh.ch/cpp/common/base64.html
  std::string base64Decoder(std::string const& encoded_string) {
    int in_len = encoded_string.size();
    int i = 0;
    int j = 0;
    int in_ = 0;
    unsigned char char_array_4[4], char_array_3[3];
    std::string ret;

    while (in_len-- && ( encoded_string[in_] != '=') && isBase64(encoded_string[in_])) {
      char_array_4[i++] = encoded_string[in_]; in_++;
      if (i ==4) {
        for (i = 0; i <4; i++) {
          char_array_4[i] = base64_chars.find(char_array_4[i]);
        }

        char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
        char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
        char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

        for (i = 0; (i < 3); i++) {
          ret += char_array_3[i];
        }
        i = 0;
      }
    }

    if (i) {
      for (j = i; j <4; j++) {
        char_array_4[j] = 0;
      }

      for (j = 0; j <4; j++) {
        char_array_4[j] = base64_chars.find(char_array_4[j]);
      }

      char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
      char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
      char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

      for (j = 0; (j < i - 1); j++) {
        ret += char_array_3[j];
      }
    }

    return ret;
  }

  // ***************************************************************************
  // aurostd::bin2base64(const std::string& b_file, std::string& b64String)
  // ***************************************************************************
  // converts binary file to base64 string
  bool bin2base64(const std::string& b_file, std::string& b64String) {
    stringstream output;
    if (!aurostd::FileExist(b_file)) {
      cerr << "ERROR - aurostd::bin2base64: Binary file " << b_file << " does not exist!";
      return FALSE;
    }
    ifstream file(b_file.c_str(), std::ios::in | std::ios::binary );
    output << b64_encoder << file;
    b64String=output.str();
    return TRUE;
  }

  // ***************************************************************************
  // aurostd::base642bin(const std::string& b64String, const std::string& b_file)
  // ***************************************************************************
  // converts base64 string to binary file
  bool base642bin(const std::string& b64String, const std::string& b_file) {
    ofstream output;
    output.open(b_file.c_str(),std::ios::out | std::ios::binary);
    output << b64_decoder << b64String;
    output.flush();output.clear();output.close();
    return TRUE;
  }

  b64_encoder_proxy operator<<(std::ostream & os, b64_encoder_creator) {
    return b64_encoder_proxy(os);
  }

  b64_decoder_proxy operator<<(std::ostream & os, b64_decoder_creator) {
    return b64_decoder_proxy(os);
  }

}  // namespace aurostd
//CO END

#endif  // _AURO_IMPLEMENTATIONS_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
