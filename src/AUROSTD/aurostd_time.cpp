
#include "aurostd_time.h"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

#include <sys/time.h>
#include <time.h>

#include "aurostd.h"

#include "aflow_xhost.h" // for XHOST.DEBUG

using std::cerr;
using std::endl;
using std::string;
using std::vector;

// ***************************************************************************
// TIME evolution stuff
namespace aurostd {
  int get_day() {
    const time_t t = time(nullptr);
    struct tm* ptr_now = localtime(&t);
    return get_day(*ptr_now);
  } // CO20200624
  int get_day(const tm& tstruct) {
    return tstruct.tm_mday;
  } // CO20200624
  int get_month() {
    const time_t t = time(nullptr);
    struct tm* ptr_now = localtime(&t);
    return get_month(*ptr_now);
  } // CO20200624
  int get_month(const tm& tstruct) {
    return tstruct.tm_mon + 1;
  }  // CO20200624
  int get_year() {
    const time_t t = time(nullptr);
    struct tm* ptr_now = localtime(&t);
    return get_year(*ptr_now);
  }
  int get_year(const tm& tstruct) {
    return tstruct.tm_year + 1900;
  } // CO20200624
  void get_offset_utc(int& offset_hours, int& offset_mins) {
    const time_t t = time(nullptr);
    struct tm* ptr_now = localtime(&t);
    return get_offset_utc(*ptr_now, offset_hours, offset_mins);
  }  // CO20210601: https://codereview.stackexchange.com/questions/175353/getting-current-timezone
  void get_offset_utc(const tm& _tstruct_inp, int& offset_hours, int& offset_mins) {  // CO20210601
    // https://codereview.stackexchange.com/questions/175353/getting-current-timezone
    const bool LDEBUG = (false || XHOST.DEBUG);
    char buffer[30];
    tm tstruct_inp = _tstruct_inp; // mktime modifies tstruct, make copy
    const time_t t_inp = std::mktime(&tstruct_inp);
    if (LDEBUG) {
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
      strftime(buffer, 30, "%F %T %Z", &tstruct_inp);
      cerr << __AFLOW_FUNC__ << " tstruct_inp=" << buffer << endl;  //%Y:%m:%d %H:%M:%S
      cerr << __AFLOW_FUNC__ << " ///////////////////////////////////////////////" << endl;
    }
    //
    time_t t_local = t_inp;
    const bool fix_utc_2_now = false; // this is good for debugging different time zones
    if (fix_utc_2_now) {
      t_local = time(nullptr);
    } // struct tm *tstruct_now=localtime(&t_now);
    struct tm* ptr_tstruct_gmt = std::gmtime(&t_local); // get gmt wrt to local (now vs. input)
    if (LDEBUG) {
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
      tm tstruct_tmp = *ptr_tstruct_gmt;  // mktime modifies tstruct
      cerr << __AFLOW_FUNC__ << " mktime(ptr_tstruct_gmt)=" << std::mktime(&tstruct_tmp) << endl;
      strftime(buffer, 30, "%F %T %Z", ptr_tstruct_gmt);
      cerr << __AFLOW_FUNC__ << " tstruct_gmt=" << buffer << endl;  //%Y:%m:%d %H:%M:%S
      cerr << __AFLOW_FUNC__ << " ///////////////////////////////////////////////" << endl;
    }
    // NB: before the following DST change, the %Z of ptr_struct_gmt is GMT, after it is EST (or EDT)
    ptr_tstruct_gmt->tm_isdst = -1; // VERY IMPORTANT, forces mktime to figure out dst
    const time_t t_gmt = std::mktime(ptr_tstruct_gmt);
    if (LDEBUG) {
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
      strftime(buffer, 30, "%F %T %Z", ptr_tstruct_gmt);
      cerr << __AFLOW_FUNC__ << " tstruct_gmt=" << buffer << endl;  //%Y:%m:%d %H:%M:%S
      cerr << __AFLOW_FUNC__ << " ///////////////////////////////////////////////" << endl;
    }
    //
    const long int t_diff = static_cast<long int>(t_inp - t_gmt); // flip to get right sign
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " t_diff=" << t_diff << endl;
    }
    const double offset = (double) t_diff / 3600.0;
    offset_hours = (int) std::floor(offset);
    offset_mins = (int) (std::round((offset - (double) offset_hours) * 60.0));
  }
  long int get_date() {
    const time_t t = time(nullptr);
    struct tm* ptr_now = localtime(&t);
    return get_date(*ptr_now);
  }  // CO20200624
  long int get_date(const tm& tstruct) {
    return aurostd::get_year(tstruct) * 10000 + aurostd::get_month(tstruct) * 100 + aurostd::get_day(tstruct);
  }  // CO20200624
  int get_hour() {
    const time_t t = time(nullptr);
    struct tm* ptr_now = localtime(&t);
    return get_hour(*ptr_now);
  }  // CO20200624
  int get_hour(const tm& tstruct) {
    return tstruct.tm_hour;
  }  // CO20200624
  int get_min() {
    const time_t t = time(nullptr);
    struct tm* ptr_now = localtime(&t);
    return get_min(*ptr_now);
  }  // CO20200624
  int get_min(const tm& tstruct) {
    return tstruct.tm_min;
  }  // CO20200624
  int get_sec() {
    const time_t t = time(nullptr);
    struct tm* ptr_now = localtime(&t);
    return get_sec(*ptr_now);
  }  // CO20200624
  int get_sec(const tm& tstruct) {
    return tstruct.tm_sec;
  }  // CO20200624
  long double get_seconds() {
    timeval tim;
    gettimeofday(&tim, nullptr);
    return tim.tv_sec + tim.tv_usec / 1e6;
  }
  long double get_seconds(long double reference_seconds) {
    return get_seconds() - reference_seconds;
  }
  long double get_delta_seconds(long double& seconds_begin) {
    const long double out = get_seconds() - seconds_begin;
    seconds_begin = get_seconds();
    return out;
  }
  long double get_mseconds() {
    timeval tim;
    gettimeofday(&tim, nullptr);
    return tim.tv_usec / 1000.0;
  }
  long double get_mseconds(long double reference_useconds) {
    return (aurostd::get_useconds() - reference_useconds) / 1000.0;
  }
  long double get_delta_mseconds(long double& useconds_begin) {
    const long double out = (aurostd::get_useconds() - useconds_begin) / 1000.0;
    useconds_begin = aurostd::get_useconds() / 1000.0;
    return out;
  }
  long double get_useconds() {
    timeval tim;
    gettimeofday(&tim, nullptr);
    return tim.tv_usec;
  }
  long double get_useconds(long double reference_useconds) {
    return aurostd::get_useconds() - reference_useconds;
  }
  long double get_delta_useconds(long double& useconds_begin) {
    const long double out = aurostd::get_useconds() - useconds_begin;
    useconds_begin = aurostd::get_useconds();
    return out;
  }
  string get_time() {
    const time_t t = time(nullptr);
    struct tm* ptr_now = localtime(&t);
    return get_time(*ptr_now);
  } // CO20200624
  string get_time(const tm& tstruct) {
    const int h = get_hour(tstruct);
    const int m = get_min(tstruct);
    const int s = get_sec(tstruct);
    return (h < 10 ? "0" : "") + aurostd::utype2string(h) + ":" + (m < 10 ? "0" : "") + aurostd::utype2string(m) + ":" + (s < 10 ? "0" : "") + aurostd::utype2string(s);
  } // CO20200624
  string get_datetime(bool include_utc_offset) {
    const time_t t = time(nullptr);
    struct tm* ptr_now = localtime(&t);
    return get_datetime(*ptr_now, include_utc_offset);
  } // CO20200624
  string get_datetime(const tm& tstruct, bool include_utc_offset) {  // CO20200624
    string datetime = utype2string(get_date(tstruct)) + "_" + get_time(tstruct);
    if (include_utc_offset) { // CO20210624
      int offset_hours = 0;
      int offset_mins = 0;
      get_offset_utc(tstruct, offset_hours, offset_mins);
      datetime += "_GMT";
      datetime += (std::signbit(offset_hours) ? string("-") : string("+")); // sign +/-
      const bool pad_hours = false; // default AFLOW behavior does not pad hours
      datetime += aurostd::PaddedNumString(std::abs(offset_hours), (pad_hours ? 2 : 1)); // CO20210624 - PaddedNumString() struggles with negative numbers
      const bool print_mins = false;  // default AFLOW behavior does not print mins
      if (print_mins || offset_mins != 0) {
        datetime += ":" + aurostd::PaddedNumString(offset_mins, 2);
      }
    }
    return datetime;
  }
  string get_datetime_formatted(const string& date_delim, bool include_time, const string& date_time_sep, const string& time_delim) {
    const time_t t = time(nullptr);
    struct tm* ptr_now = localtime(&t);
    return get_datetime_formatted(*ptr_now, date_delim, include_time, date_time_sep, time_delim);
  }  // CO20171215  //CO20200624
  string get_datetime_formatted(const tm& tstruct, const string& date_delim, bool include_time, const string& date_time_sep, const string& time_delim) {  // CO20171215 //CO20200624
    stringstream misc_ss;
    const int y = aurostd::get_year(tstruct);
    const int b = aurostd::get_month(tstruct);
    const int d = aurostd::get_day(tstruct); // CO20200624
    misc_ss << y << date_delim << (b < 10 ? "0" : "") << b << date_delim << (d < 10 ? "0" : "") << d;
    if (include_time) {
      const int h = get_hour(tstruct);
      const int m = get_min(tstruct);
      const int s = get_sec(tstruct);
      misc_ss << date_time_sep << (h < 10 ? "0" : "") << h << time_delim << (m < 10 ? "0" : "") << m << time_delim << (s < 10 ? "0" : "") << s; // CO20200624
    }
    return misc_ss.str();
  }
  bool beep(uint freq, uint duration) {
    return aurostd::execute("beep -f " + aurostd::utype2string<uint>(freq) + " -l " + aurostd::utype2string<uint>(duration));
  }
} // namespace aurostd

// ***************************************************************************
// aflow_get_time_string
// ***************************************************************************
string aflow_get_time_string() {
  // OUTPUT: http://www.cplusplus.com/reference/ctime/ctime/
  // Www Mmm dd hh:mm:ss yyyy
  // Where Www is the weekday, Mmm the month (in letters), dd the day of the month, hh:mm:ss the time, and yyyy the year.
  // The string is followed by a new-line character ('\n') and terminated with a null-character.
  const long ltime = time(nullptr);

  string date = string(ctime(&ltime));
  if (!date.empty()) {
    if (date.at(date.length() - 1) == '\n') {
      date.erase(date.length() - 1);
    }
  }
  return date;
}
string aflow_convert_time_ctime2aurostd(const string& time_LOCK) { // CO20200624
  // refer to aflow_get_time_string()
  // convert Www Mmm dd hh:mm:ss yyyy style to aurostd::get_datetime() one
  const bool LDEBUG = (false || XHOST.DEBUG);

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
  }

  vector<string> tokens;
  aurostd::string2tokens(time_LOCK, tokens);
  if (tokens.size() != 5) {
    return "";
  }

  // https://en.cppreference.com/w/c/chrono/strftime
  //'Www Mmm dd hh:mm:ss yyyy' === '%a %b %d %H:%M:%S %Y'
  // https://stackoverflow.com/questions/19524720/using-strptime-converting-string-to-time-but-getting-garbage
  tm tstruct;
  if (!strptime(time_LOCK.c_str(), "%a %b %d %H:%M:%S %Y", &tstruct)) {
    return "";
  }
  tstruct.tm_isdst = -1;  // let computer figure it out
  std::mktime(&tstruct);  // get is_dst

  if (LDEBUG) {
    char buffer[30];
    strftime(buffer, 30, "%F %T %Z", &tstruct);
    cerr << __AFLOW_FUNC__ << " tstruct=" << buffer << endl;
  }

  if (LDEBUG) {
    cerr << __AFLOW_FUNC__ << " END" << endl;
  }

  const bool include_utc_offset = true;
  return aurostd::get_datetime(tstruct, include_utc_offset);
}

// ***************************************************************************
// aflow_get_time_string_short
// ***************************************************************************
string aflow_get_time_string_short() {
  string date;
  vector<string> tokens;
  aurostd::string2tokens(aflow_get_time_string(), tokens);
  date = "na";
  if (tokens.size() > 4) {
    if (tokens.at(2).size() > 1) {
      date = tokens.at(4).substr(2, 2) + tokens.at(1) + tokens.at(2);
    } else {
      date = tokens.at(4).substr(2, 2) + tokens.at(1) + "0" + tokens.at(2);
    }
    aurostd::StringSubstInPlace(date, "Jan", "01");
    aurostd::StringSubstInPlace(date, "Feb", "02");
    aurostd::StringSubstInPlace(date, "Mar", "03");
    aurostd::StringSubstInPlace(date, "Apr", "04");
    aurostd::StringSubstInPlace(date, "May", "05");
    aurostd::StringSubstInPlace(date, "Jun", "06");
    aurostd::StringSubstInPlace(date, "Jul", "07");
    aurostd::StringSubstInPlace(date, "Aug", "08");
    aurostd::StringSubstInPlace(date, "Sep", "09");
    aurostd::StringSubstInPlace(date, "Oct", "10");
    aurostd::StringSubstInPlace(date, "Nov", "11");
    aurostd::StringSubstInPlace(date, "Dec", "12");
    date = date.substr(0, 6);
  }
  if (!date.empty()) {
    if (date.at(date.length() - 1) == '\n') {
      date.erase(date.length() - 1);
    }
  }
  return date;
}
