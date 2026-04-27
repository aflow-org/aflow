
#ifndef AUROSTD_TIME_H
#define AUROSTD_TIME_H

#include <ctime>
#include <string>

#ifndef uint
typedef unsigned uint;
#endif

// ----------------------------------------------------------------------------
// TIME stuff
namespace aurostd {
  int get_day();
  int get_day(const tm& tstruct); // CO20200624
  int get_month();
  int get_month(const tm& tstruct);  // CO20200624
  int get_year();
  int get_year(const tm& tsruct);  // CO20200624
  void get_offset_utc(int& offset_hours, int& offset_mins); // CO20210601
  void get_offset_utc(const tm& tstruct, int& offset_hours, int& offset_mins);  // CO20210601
  long int get_date();
  long int get_date(const tm& tsruct);  // CO20200624
  int get_hour();
  int get_hour(const tm& tsruct);  // CO20200624
  int get_min();
  int get_min(const tm& tsruct); // CO20200624
  int get_sec();
  int get_sec(const tm& tsruct); // CO20200624
  long double get_seconds();
  long double get_seconds(long double reference_seconds);
  long double get_delta_seconds(long double& seconds_begin);
  long double get_mseconds();
  long double get_mseconds(long double reference_useconds);
  long double get_delta_mseconds(long double& useconds_begin);
  long double get_useconds();
  long double get_useconds(long double reference_useconds);
  long double get_delta_useconds(long double& useconds_begin);
  std::string get_time();
  std::string get_time(const tm& tsruct); // CO20200624
  std::string get_datetime(bool include_utc_offset = false);
  std::string get_datetime(const tm& tsruct, bool include_utc_offset = false); // CO20200624
  std::string get_datetime_formatted(const std::string& date_delim = "/", bool include_time = true, const std::string& date_time_sep = " ", const std::string& time_delim = ":");  // CO20171215
  std::string get_datetime_utc_formatted(const std::string& date_delim="/",bool include_time=true,const std::string& date_time_sep=" ",const std::string& time_delim=":");  //CO20171215
  std::string get_datetime_formatted(const tm& tsruct, const std::string& date_delim = "/", bool include_time = true, const std::string& date_time_sep = " ", const std::string& time_delim = ":"); // CO20171215 //CO20200624
  bool beep(uint = 2000, uint = 100); // standard values
} // namespace aurostd

std::string aflow_get_time_string();
std::string aflow_convert_time_ctime2aurostd(const std::string& time_LOCK); // CO20200624
std::string aflow_get_time_string_short();

#endif // AUROSTD_TIME_H
