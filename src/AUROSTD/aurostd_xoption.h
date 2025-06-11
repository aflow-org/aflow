// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 2013-2014
// added template<class utype> bool xoption::args2addattachedscheme SC 2017

#ifndef _AUROSTD_XOPTION_H_
#define _AUROSTD_XOPTION_H_

#include <ostream>
#include <string>
#include <vector>

#include <sys/types.h>

// --------------------------------------------------------------------------
// general flag for xoption to take/manipulate options
#define aurostd_xoptionONOFF int(-1)
#define aurostd_xoptionMULTI int(-2)

namespace aurostd {
  // namespace aurostd
  class xoption {
  public:
      // trivial constructurs/destuctors/operators
    xoption();                                         // default, just allocate
    ~xoption();                                        // kill everything
    xoption(const xoption& b);                         // constructor copy
    const xoption& operator=(const xoption& b);        // copy
    friend std::ostream& operator<<(std::ostream&, const xoption&);       // ostream
    void clear();                                  // clear
      // CONTENT
    std::string keyword;            // the keyword found (we can provide a bunch with |) //CO20180404
    bool isentry;              // is the entry available
    std::string content_string;     // the content
    double content_double;     // the content
    int content_int;           // the content
    uint content_uint;         // the content
    bool option;               // the output
    bool option_default;       // the default
    std::string xscheme;            // the content
    std::vector<std::string> vxscheme;   // tokenized "," content
    std::vector<std::string> vxsghost;   // tokenized "," content
    bool preserved;            // the output
      // LOAD BOOLS FUNCTIONS
    void options2entry(const std::string&, const std::string&, int = aurostd_xoptionONOFF, const std::string& xscheme_DEFAULT = "");  // CO20210805 - const&
    void scheme2scheme(char, const std::string&); // CO20210805 - const&
    void scheme2scheme(const std::string&, const std::string&);  // CO20210805 - const&
    [[nodiscard]] bool isscheme(const std::string&) const; // check if available //CO20180101 //SC20191227 //CO20210805 - const&
    uint opscheme(const std::string&, bool);  // add/remove scheme then returns vscheme.size()  //CO20210805 - const&
    uint push(const std::string&);           // add scheme then returns vscheme.size() //CO20210805 - const&
    uint pop(const std::string&);            // remove scheme then returns vscheme.size()  //CO20210805 - const&
      // for plain flags
    bool flag(const std::string&, bool);      // if bool=true/false => add/remove "string"  //CO20210805 - const&
    [[nodiscard]] bool flag(const std::string&) const;     // interrogate=true/false, same as ischeme //CO20180101  //SC20191227 //CO20210805 - const&
    [[nodiscard]] bool flag() const;       // return if there is any scheme inside //CO20180101 //SC20191227
      // attached stuff..
    [[nodiscard]] bool isdefined(const std::string&) const;                                            // SC20200114  //CO20210805 - const&
    uint opattachedscheme(const std::string&, const std::string&, bool);                        // add/remove attached_scheme if flag=true, then returns vghost.size() //CO20210805 - const&
    uint addattachedscheme(const std::string& scheme, const std::string& attached, bool flag);  // add attached_scheme if flag=true, then returns vghost.size()  //CO20210805 - const&
    uint push_attached(const std::string& scheme, const std::string& attached);        // add attached_scheme, then returns vghost.size() - like addattachedscheme with flag=true //CO20210805 - const&
    uint pop_attached(const std::string& check);                                 // remove attached_scheme, then returns vghost.size() //CO20210805 - const&
    [[nodiscard]] std::string getattachedscheme(const std::string& scheme) const; // CO20180101
    template <class utype> [[nodiscard]] utype getattachedutype(const std::string& scheme) const; // CO20200731
    bool args2addattachedscheme(std::vector<std::string>& argv, const std::string scheme, const std::string& _s_search, std::string string_default);
    bool args2addattachedscheme(std::vector<std::string>& argv, std::vector<std::string>& cmds, const std::string scheme, const std::string& _s_search, std::string string_default);
    bool args2addattachedscheme(std::vector<std::string>& argv, const std::string scheme, const std::string& _s_search, const char* string_default);
    bool args2addattachedscheme(std::vector<std::string>& argv, std::vector<std::string>& cmds, const std::string scheme, const std::string& _s_search, const char* string_default);
    template <class utype> bool args2addattachedscheme(std::vector<std::string>& argv, const std::string scheme, const std::string& _s_search, utype utype_default);
    template <class utype> bool args2addattachedscheme(std::vector<std::string>& argv, std::vector<std::string>& cmds, const std::string scheme, const std::string& _s_search, utype utype_default);
    bool refresh();

  private:                                              //
    void free();                                        // free space
    void copy(const xoption& b);                        //
  };
} // namespace aurostd

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

#endif  // _AUROSTD_XOPTION_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
