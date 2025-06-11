// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo - Corey Oses

#ifndef _AFLOW_OAIMS_CPP_
#define _AFLOW_OAIMS_CPP_

#include "aflow.h"

//---------------------------------------------------------------------------------
// class xAIMSOUT
//---------------------------------------------------------------------------------
xAIMSOUT::xAIMSOUT() {
  content="";
  filename="";
  natoms=0;
  //[CO20200404 - OBSOLETE]ERROR="";
  free();
}

xAIMSOUT::~xAIMSOUT() {
  free();
}

void xAIMSOUT::free() {
  vcontent.clear();
  vforces.clear();
}

void xAIMSOUT::copy(const xAIMSOUT& b){
  content=b.content;
  vcontent.clear(); for(uint i=0;i<b.vcontent.size();i++){vcontent.push_back(b.vcontent[i]);}
  filename=b.filename;
  vforces.clear(); for(uint i=0;i<b.vforces.size();i++){vforces.push_back(b.vforces[i]);}
  natoms=b.natoms;
}

const xAIMSOUT& xAIMSOUT::operator=(const xAIMSOUT& b) {
  if(this!=&b) {free();copy(b);}
  return *this;
}

xAIMSOUT::xAIMSOUT(const string& fileIN,bool QUIET) {
  clear(); // so it does not mess up vector/deque
  filename=fileIN;
  GetPropertiesFile(fileIN,QUIET);
}

xAIMSOUT::xAIMSOUT(const xAIMSOUT& b) { // copy PUBLIC
  //  free(); *this=b;
  copy(b);
}

void xAIMSOUT::clear() {  // clear PRIVATE
  xAIMSOUT _temp;
  string filename_aus=filename;
  copy(_temp);
  filename=filename_aus;
}

bool xAIMSOUT::GetProperties(const string& stringIN,bool QUIET) {
  stringstream sss; sss.str(stringIN);
  if(filename=="") filename="string";
  return xAIMSOUT::GetProperties(sss,QUIET);
}

bool xAIMSOUT::GetPropertiesFile(const string& fileIN,bool QUIET) {
  stringstream sss;
  aurostd::efile2stringstream(fileIN,sss);
  if(filename=="") filename=fileIN;
  return xAIMSOUT::GetProperties(sss,QUIET);
}

bool xAIMSOUT::GetPropertiesFile(const string& fileIN,uint natoms_check,bool QUIET) {
  bool flag=GetPropertiesFile(fileIN,QUIET);
  if(aurostd::abs(natoms_check-(double) natoms)>0.1) {
    stringstream message;
    message << "natoms_check(" << natoms_check << ")!= (int) natoms(" << natoms << ")";
    throw aurostd::xerror(__AFLOW_FILE__,"xAIMSOUT::GetPropertiesFile():",message,_FILE_CORRUPT_); //CO20200624
  }
  return flag;
}

bool xAIMSOUT::GetPropertiesUrlFile(const string& url,const string& file,bool VERBOSE) {
  //[CO20200502 - OBSOLETE]string tmpfile=XHOST.tmpfs+"/_aflow_"+XHOST.user+".pid"+XHOST.ostrPID.str()+".tid"+XHOST.ostrTID.str()+".a"+string(AFLOW_VERSION)+".rnd"+aurostd::utype2string(uint((double) std::floor((double)100000*aurostd::ran0())))+".u"+aurostd::utype2string(uint((double) aurostd::get_useconds()))+"_"+file;  //CO20200502 - threadID
  string tmpfile=aurostd::TmpFileCreate("xAIMSOUT_GetProperties"); //CO20200502 - threadID
  aurostd::url2file(url+"/"+file,tmpfile,VERBOSE);
  bool out=GetPropertiesFile(tmpfile);
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xAIMSOUT::GetProperties(const stringstream& stringstreamIN,bool QUIET) {
  bool LVERBOSE=(FALSE || XHOST.DEBUG || !QUIET);
  bool ERROR_flag=FALSE;
  clear();
  stringstream sss; sss.str(stringstreamIN.str());
  content=stringstreamIN.str();
  vcontent.clear();
  vector<string> vline,tokens;
  aurostd::string2vectorstring(content,vcontent);
  string line;
  if(filename=="") filename="stringstream";
  if(LVERBOSE) cerr << "xAIMSOUT::GetProperties: ---------------------------------" << endl;
  if(LVERBOSE) cerr << "xAIMSOUT::GetProperties: BEGIN" << endl;
  if(LVERBOSE) cerr << "xAIMSOUT::GetProperties: filename=[" << filename << "]" << endl;
  if(LVERBOSE) cerr.precision(12);
  // ----------------------------------------------------------------------
  //grab natoms first
  if(LVERBOSE) cerr << "xAIMSOUT::GetProperties: ---------------------------------" << endl;
  if(LVERBOSE) cerr << "xAIMSOUT::GetProperties: LOAD NATOMS" << endl;
  line="";
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--){  // NEW - FROM THE BACK 
    if(aurostd::substring2bool(vcontent.at(iline),"Number of atoms")){ // VASP
      aurostd::string2tokens(vcontent.at(iline),tokens,":");
      if(tokens.size()==2){natoms=aurostd::string2utype<double>(tokens[1]);}
    }
  }
  if(natoms==0){
    if(LVERBOSE) cerr << "WARNING - xAIMSOUT::GetProperties:" << " no natoms tag found";
    ERROR_flag=TRUE;
  }
  if(LVERBOSE){cerr << "xAIMSOUT::GetProperties: natoms=" << natoms << endl;}
  // ----------------------------------------------------------------------
  //now grab vforces
  if(LVERBOSE) cerr << "xAIMSOUT::GetProperties: ---------------------------------" << endl;
  if(LVERBOSE) cerr << "xAIMSOUT::GetProperties: LOAD VFORCES" << endl;
  int iforce;
  for(int iline=(int)vcontent.size()-1;iline>=0;iline--){  // NEW - FROM THE BACK 
    if(vforces.size()==0 && aurostd::substring2bool(vcontent.at(iline),"Total atomic forces")){ // VASP
      for(uint i=1;i<natoms+1&&iline+i<vcontent.size()-1;i++){
        line=aurostd::RemoveCharacter(vcontent[iline+i],'|');  //remove pesky | in the beginning
        aurostd::string2tokens(line,tokens," ");
        if(tokens.size()!=4){
          if(LVERBOSE) {
            cerr << "WARNING - xAIMSOUT::GetProperties: ";
            cerr << "ill-written forces line, tokens.size()=" << tokens.size() << ";";
            cerr <<"should be 4" << endl;
          }
          ERROR_flag=TRUE;
        } else {
          iforce=aurostd::string2utype<int>(tokens[0]);
          if((iforce!=(int)i)){
            if(LVERBOSE) {
              cerr << "WARNING - xAIMSOUT::GetProperties: ";
              cerr << "missing force " << i << "; ";
              cerr << "found force " << iforce << " instead"  << endl;
            }
            ERROR_flag=TRUE;
          } else {
            xvector<double> force(3);
            force[1]=aurostd::string2utype<double>(tokens[1]);
            force[2]=aurostd::string2utype<double>(tokens[2]);
            force[3]=aurostd::string2utype<double>(tokens[3]);
            vforces.push_back(force);
            if(LVERBOSE) {cerr << "xAIMSOUT::GetProperties: found new force[" << vforces.size() << "]=" << force << endl;}
          }
        }
      }
    }
  }
  if(vforces.size()==0 || vforces.size()!=natoms){
    if(LVERBOSE){cerr << "WARNING - xAIMSOUT::GetProperties: not all forces found, vforces.size()=" << vforces.size() << endl;}
  }
  // DONE NOW RETURN
  if(ERROR_flag && LVERBOSE) cerr << "WARNING - xAIMSOUT::GetProperties: ERROR_flag set in xAIMSOUT" << endl;
  if(ERROR_flag) return FALSE;
  return TRUE;
}

#endif //  _AFLOW_OAIMS_CPP_
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
