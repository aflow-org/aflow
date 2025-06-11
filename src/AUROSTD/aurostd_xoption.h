// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 2013-2014
// added template<class utype> bool xoption::args2addattachedscheme SC 2017

#ifndef _AUROSTD_XOPTION_H_
#define _AUROSTD_XOPTION_H_

// --------------------------------------------------------------------------
// general flag for xoption to take/manipulate options
#define aurostd_xoptionONOFF int(-1)
#define aurostd_xoptionMULTI int(-2)

using std::ostream;
using std::vector;
using std::string;

namespace aurostd {
  // namespace aurostd
  class xoption {
    public:
      // trivial constructurs/destuctors/operators
      xoption();                                         // default, just allocate
      ~xoption();                                        // kill everything
      xoption(const xoption& b);                         // constructor copy
      const xoption& operator=(const xoption &b);        // copy
      friend std::ostream& operator<<(std::ostream&,const xoption&);       // ostream
      void clear(void);                                  // clear
      // CONTENT
      string keyword;            // the keyword found (we can provide a bunch with |) //CO20180404
      bool isentry;              // is the entry available
      string content_string;     // the content
      double content_double;     // the content
      int content_int;           // the content
      uint content_uint;         // the content
      bool option;               // the output
      bool option_default;       // the default 
      string xscheme;            // the content
      vector<string> vxscheme;   // tokenized "," content
      vector<string> vxsghost;   // tokenized "," content
      bool preserved;            // the output
      // LOAD BOOLS FUNCTIONS
      void options2entry(const string&,const string&,int=aurostd_xoptionONOFF,const string& xscheme_DEFAULT="");  //CO20210805 - const&
      void scheme2scheme(char,const string&); //CO20210805 - const&
      void scheme2scheme(const string&,const string&);  //CO20210805 - const&
      bool isscheme(const string&) const; // check if available //CO20180101 //SC20191227 //CO20210805 - const&
      // [OBSOLETE] uint addscheme(string);      // add scheme then returns vscheme.size()
      // [OBSOLETE] uint purgescheme(string);    // remove scheme then returns vscheme.size()
      uint opscheme(const string&,bool);  // add/remove scheme then returns vscheme.size()  //CO20210805 - const&
      uint push(const string&);           // add scheme then returns vscheme.size() //CO20210805 - const&
      uint pop(const string&);            // remove scheme then returns vscheme.size()  //CO20210805 - const&
      // for plain flags
      bool flag(const string&,bool);      // if bool=TRUE/FALSE => add/remove "string"  //CO20210805 - const&
      bool flag(const string&) const;     // interrogate=TRUE/FALSE, same as ischeme //CO20180101  //SC20191227 //CO20210805 - const&
      bool flag(void) const;       // return if there is any scheme inside //CO20180101 //SC20191227
      // attached stuff..
      bool isdefined(const string&) const;                                            //SC20200114  //CO20210805 - const&
      uint opattachedscheme(const string&,const string&,bool);                        // add/remove attached_scheme if flag=TRUE, then returns vghost.size() //CO20210805 - const&
      uint addattachedscheme(const string& scheme,const string& attached,bool flag);  // add attached_scheme if flag=TRUE, then returns vghost.size()  //CO20210805 - const&
      // [OBSOLETE] uint purgeattachedscheme(string check);            // remove attached_scheme, then returns vghost.size() - same as pop_attached
      uint push_attached(const string& scheme,const string& attached);        // add attached_scheme, then returns vghost.size() - like addattachedscheme with flag=TRUE //CO20210805 - const&
      uint pop_attached(const string& check);                                 // remove attached_scheme, then returns vghost.size() //CO20210805 - const&
      string getattachedscheme(const string& scheme) const; //CO20180101
      template<class utype> utype getattachedutype(const string& scheme) const;  //CO20200731
      bool args2addattachedscheme(vector<string> &argv,const string scheme,const string& _s_search,string string_default); 
      bool args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string scheme,const string& _s_search,string string_default);
      bool args2addattachedscheme(vector<string> &argv,const string scheme,const string& _s_search,char const* string_default); 
      bool args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string scheme,const string& _s_search,char const* string_default);
      template<class utype> bool args2addattachedscheme(vector<string> &argv,const string scheme,const string& _s_search,utype utype_default); 
      template<class utype> bool args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string scheme,const string& _s_search,utype utype_default);
      bool refresh(void);
    private:                                              //
      void free();                                        // free space
      void copy(const xoption& b);                        //
  };
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

#endif  // _AUROSTD_XOPTION_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

