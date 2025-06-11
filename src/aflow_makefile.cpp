// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow COREY OSES - Duke University 2013-2021                  *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses
// corey.oses@duke.edu

#ifndef _AFLOW_MAKEFILE_CPP_
#define _AFLOW_MAKEFILE_CPP_

#include "aflow.h"

#define _DEBUG_MAKEFILE_ true
#define EXCLUDE_DATA true

//DX20200831 - revert back to hard-coded prototypes and do not compile symbolic math
#define COMPILE_HARDCODED_PROTOTYPES false
#define COMPILE_SYMBOLIC !(COMPILE_HARDCODED_PROTOTYPES)

namespace makefile {
  void trimPath(string& filename){
    //get rid of ".", "..", etc.
    filename=aurostd::CleanFileName(filename);  //get rid of /./
    aurostd::CleanStringASCII_InPlace(filename);
    if(filename.size()>1 && filename[0]=='.' && filename[1]=='/'){filename=filename.substr(2,string::npos);} //removes ./ from beginning so j can start at 1 below
    if(filename.find("..")!=string::npos){ //might go up and down in tree structure
      vector<string> vtokens;
      aurostd::string2tokens(filename,vtokens,"/");

      uint j=0;
      for(j=1;j<vtokens.size();j++){  //start with 1
        if(vtokens[j]==".."){vtokens.erase(vtokens.begin()+j-1,vtokens.begin()+j+1);} //get rid of AUROSTD/../
      }
      filename=aurostd::joinWDelimiter(vtokens,"/");
    }
  }
  void getDependencies(const string& _filename,vector<string>& files_already_explored,vector<string>& dfiles){
    bool mt_required=false;
    return getDependencies(_filename,files_already_explored,dfiles,mt_required);
  }
  void getDependencies(const string& _filename,vector<string>& files_already_explored,vector<string>& dfiles,bool& mt_required){
    //find #include and get files
    bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);

    string filename=aurostd::CleanFileName(_filename);
    aurostd::CleanStringASCII_InPlace(filename);
    trimPath(filename);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " fetching dependencies for " << filename << endl;};
    if(!aurostd::FileExist(filename)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"File \""+filename+"\" does not exist",_VALUE_ILLEGAL_);}

    if(aurostd::WithinList(files_already_explored,filename)){return;}
    files_already_explored.push_back(filename);

    string dirname=".";
    if(filename.find("/")!=string::npos){ //this is a directory
      vector<string> vtokens;
      aurostd::string2tokens(filename,vtokens,"/");
      vtokens.pop_back(); //remove filename
      dirname=aurostd::joinWDelimiter(vtokens,"/");
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " dirname=" << dirname << endl;}

    vector<string> vlines;
    aurostd::file2vectorstring(filename,vlines);

    uint i=0;
    string::size_type loc;
    string line="",dfile="";
    mt_required=false;
    //only check for line comments. does not work for block comments
    for(i=0;i<vlines.size();i++){
      line = vlines[i];
      loc=line.find("//");line=line.substr(0,loc);line=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(line); //remove comments and left/right whitespaces
      if(line.find("#include")==0 && (line.find("<")==string::npos && line.find(">")==string::npos)){ //needs to be at the beginning of the line, files with <> are standard libraries
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " line with #include = \"" << line << "\"" << endl;}
        dfile=line;
        aurostd::StringSubst(dfile,"#include","");
        aurostd::StringSubst(dfile,"\"","");
        dfile=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(dfile);
        dfile=dirname+"/"+dfile;
        trimPath(dfile);  //resolve relative path as trim as possible
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " dfile=" << dfile << endl;}
        dfiles.push_back(dfile);
        getDependencies(dfile,files_already_explored,dfiles); //do not propagate mt_required from sub-dependencies
      }
      if(line.find("#define")==0 && line.find("AFLOW_")!=string::npos && line.find("_MULTITHREADS_ENABLE")!=string::npos){mt_required=true;}
      if(line.find("AFLOW_MULTITHREADS_ENABLE")!=string::npos){mt_required=true;}  //ME20220204
    }
    std::sort(dfiles.begin(),dfiles.end());dfiles.erase( std::unique( dfiles.begin(), dfiles.end() ), dfiles.end() );  //get unique set of dependent files
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " dfiles=" << aurostd::joinWDelimiter(dfiles,",") << endl;}
  }

  //[CO20200508 - OBSOLETE]void replaceMakefileDefinitions(const vector<string>& vvariables,vector<string>& vdefinitions){
  //[CO20200508 - OBSOLETE]  bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
  //[CO20200508 - OBSOLETE]
  //[CO20200508 - OBSOLETE]  bool replacement_made=false;
  //[CO20200508 - OBSOLETE]  string::size_type loc_var_first,loc_var_last;
  //[CO20200508 - OBSOLETE]  string variable="";
  //[CO20200508 - OBSOLETE]  uint i=0,j=0;
  //[CO20200508 - OBSOLETE]  for(j=0;j<vdefinitions.size();j++){
  //[CO20200508 - OBSOLETE]    string& definitions=vdefinitions[j];
  //[CO20200508 - OBSOLETE]    loc_var_first=definitions.find("$(");
  //[CO20200508 - OBSOLETE]    if(loc_var_first!=string::npos){ //variable found in definition
  //[CO20200508 - OBSOLETE]      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vdefinitions[j=" << j << "]=\"" << definitions << "\"" << endl;}
  //[CO20200508 - OBSOLETE]      loc_var_last=definitions.find(')',loc_var_first);
  //[CO20200508 - OBSOLETE]      variable=definitions.substr(loc_var_first,loc_var_last-loc_var_first);
  //[CO20200508 - OBSOLETE]      aurostd::StringSubst(variable,"$(","");
  //[CO20200508 - OBSOLETE]      aurostd::StringSubst(variable,")","");
  //[CO20200508 - OBSOLETE]      //aurostd::StringSubst(variable,"\\\"","");  //weird functionality: $(TIME)\" needs to remove \"
  //[CO20200508 - OBSOLETE]      if(LDEBUG){cerr << __AFLOW_FUNC__ << " variable inside vdefinitions[j=" << j << "]: \"" << variable << "\"" << endl;}
  //[CO20200508 - OBSOLETE]      if(variable.find("shell ")!=string::npos){continue;} //skip calls to shell
  //[CO20200508 - OBSOLETE]      if(variable.find(":")!=string::npos && variable.find("=")!=string::npos){continue;} //skip subst_ref
  //[CO20200508 - OBSOLETE]      replacement_made=false;
  //[CO20200508 - OBSOLETE]      for(i=0;i<vvariables.size();i++){
  //[CO20200508 - OBSOLETE]        if(variable==vvariables[i]){
  //[CO20200508 - OBSOLETE]          if(LDEBUG){cerr << __AFLOW_FUNC__ << " vdefinitions[j=" << j << "](pre )=\"" << definitions << "\"" << endl;}
  //[CO20200508 - OBSOLETE]          string& definitions_sub=vdefinitions[i];
  //[CO20200508 - OBSOLETE]          if(LDEBUG){cerr << __AFLOW_FUNC__ << " vdefinitions[i=" << i << "]=\"" << definitions_sub << "\"" << endl;}
  //[CO20200508 - OBSOLETE]          if(definitions_sub.find("shell ")!=string::npos){break;} //skip replacing calls to shell
  //[CO20200508 - OBSOLETE]          if(definitions_sub.find(":")!=string::npos && definitions_sub.find("=")!=string::npos){break;} //skip replacing subst_ref
  //[CO20200508 - OBSOLETE]          aurostd::StringSubst(definitions,"$("+variable+")",definitions_sub);
  //[CO20200508 - OBSOLETE]          replacement_made=true;
  //[CO20200508 - OBSOLETE]          if(LDEBUG){cerr << __AFLOW_FUNC__ << " vdefinitions[j=" << j << "](post)=\"" << definitions << "\"" << endl;}
  //[CO20200508 - OBSOLETE]          break;
  //[CO20200508 - OBSOLETE]        }
  //[CO20200508 - OBSOLETE]      }
  //[CO20200508 - OBSOLETE]      if(replacement_made){j--;continue;}
  //[CO20200508 - OBSOLETE]    }
  //[CO20200508 - OBSOLETE]  }
  //[CO20200508 - OBSOLETE]}
  //[CO20200508 - OBSOLETE]void replaceMakefileDefinitions(const vector<string>& vvariables,vector<vector<string> >& vvdefinitions){
  //[CO20200508 - OBSOLETE]  bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
  //[CO20200508 - OBSOLETE]  if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}
  //[CO20200508 - OBSOLETE]
  //[CO20200508 - OBSOLETE]  vector<string> vdefinitions;
  //[CO20200508 - OBSOLETE]  for(uint j=0;j<vvdefinitions.size();j++){vdefinitions.push_back(aurostd::joinWDelimiter(vvdefinitions[j]," "));}
  //[CO20200508 - OBSOLETE]  replaceMakefileDefinitions(vvariables,vdefinitions);
  //[CO20200508 - OBSOLETE]  splitMakefileDefinitions(vdefinitions,vvdefinitions);
  //[CO20200508 - OBSOLETE]}
  //[CO20200508 - OBSOLETE]void splitMakefileDefinitions(const string& definitions,vector<string>& vdefinitions){
  //[CO20200508 - OBSOLETE]  bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
  //[CO20200508 - OBSOLETE]
  //[CO20200508 - OBSOLETE]  //split by spaces, except for $(shell date ...)
  //[CO20200508 - OBSOLETE]  vdefinitions.clear();
  //[CO20200508 - OBSOLETE]  if(LDEBUG){cerr << __AFLOW_FUNC__ << " definitions=" << definitions << endl;}
  //[CO20200508 - OBSOLETE]  string definition="";
  //[CO20200508 - OBSOLETE]  uint count_start=0,count_end=0;
  //[CO20200508 - OBSOLETE]  //aurostd::string2tokens(definitions,vdefinitions," ");
  //[CO20200508 - OBSOLETE]  string::size_type loc_var_first,loc_var_middle,loc_var_last;
  //[CO20200508 - OBSOLETE]  loc_var_first=loc_var_last=0;
  //[CO20200508 - OBSOLETE]  loc_var_last=definitions.find(' ',loc_var_last+1);
  //[CO20200508 - OBSOLETE]  uint i=0;
  //[CO20200508 - OBSOLETE]  while(loc_var_last!=string::npos){
  //[CO20200508 - OBSOLETE]    if(LDEBUG){
  //[CO20200508 - OBSOLETE]      cerr << __AFLOW_FUNC__ << " loc_var_first=" << loc_var_first << endl;
  //[CO20200508 - OBSOLETE]      cerr << __AFLOW_FUNC__ << " loc_var_last=" << loc_var_last << endl;
  //[CO20200508 - OBSOLETE]    }
  //[CO20200508 - OBSOLETE]    definition=definitions.substr(loc_var_first,loc_var_last-loc_var_first);
  //[CO20200508 - OBSOLETE]    loc_var_middle=loc_var_last+1;
  //[CO20200508 - OBSOLETE]    loc_var_last=definitions.find(' ',loc_var_middle);
  //[CO20200508 - OBSOLETE]    //check that you have enclosing $()
  //[CO20200508 - OBSOLETE]    count_start=count_end=0;
  //[CO20200508 - OBSOLETE]    for(i=0;i<definition.size();i++){
  //[CO20200508 - OBSOLETE]      if(definition[i]=='$' && i+1<definition.size() && definition[i+1]=='('){count_start++;}
  //[CO20200508 - OBSOLETE]      if(definition[i]==')'){count_end++;}
  //[CO20200508 - OBSOLETE]    }
  //[CO20200508 - OBSOLETE]    if(count_start!=count_end){continue;}
  //[CO20200508 - OBSOLETE]    //if(definition.find("$(")!=string::npos && definition.find(")")==string::npos){continue;}
  //[CO20200508 - OBSOLETE]    //if(definition.find("$(")==string::npos && definition.find(")")!=string::npos){continue;}
  //[CO20200508 - OBSOLETE]    definition=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(definition); //repetita iuvant
  //[CO20200508 - OBSOLETE]    if(definition.empty()){continue;}
  //[CO20200508 - OBSOLETE]    if(LDEBUG){cerr << __AFLOW_FUNC__ << " splitting: " << definition << endl;}
  //[CO20200508 - OBSOLETE]    vdefinitions.push_back(definition);
  //[CO20200508 - OBSOLETE]    loc_var_first=loc_var_middle;
  //[CO20200508 - OBSOLETE]  }
  //[CO20200508 - OBSOLETE]  //get end
  //[CO20200508 - OBSOLETE]  definition=definitions.substr(loc_var_first,string::npos);
  //[CO20200508 - OBSOLETE]  definition=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(definition); //repetita iuvant
  //[CO20200508 - OBSOLETE]  if(!definition.empty()){
  //[CO20200508 - OBSOLETE]    if(LDEBUG){cerr << __AFLOW_FUNC__ << " splitting: " << definition << endl;}
  //[CO20200508 - OBSOLETE]    vdefinitions.push_back(definition);
  //[CO20200508 - OBSOLETE]  }
  //[CO20200508 - OBSOLETE]}
  //[CO20200508 - OBSOLETE]void splitMakefileDefinitions(const vector<string>& vdefinitions,vector<vector<string> >& vvdefinitions){
  //[CO20200508 - OBSOLETE]  bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
  //[CO20200508 - OBSOLETE]  
  //[CO20200508 - OBSOLETE]  //split by spaces, except for $(shell date ...)
  //[CO20200508 - OBSOLETE]  for(uint j=0;j<vdefinitions.size();j++){
  //[CO20200508 - OBSOLETE]    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vdefinitions[j=" << j << "]=\"" << vdefinitions[j] << "\"" << endl;}
  //[CO20200508 - OBSOLETE]    vvdefinitions.push_back(vector<string>(0));
  //[CO20200508 - OBSOLETE]    splitMakefileDefinitions(vdefinitions[j],vvdefinitions.back());
  //[CO20200508 - OBSOLETE]  }
  //[CO20200508 - OBSOLETE]}
  //[CO20200508 - OBSOLETE]void readMakefileVariables(const string& directory,vector<string>& vvariables,vector<vector<string> >& vvdefinitions){
  //[CO20200508 - OBSOLETE]  bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
  //[CO20200508 - OBSOLETE]
  //[CO20200508 - OBSOLETE]  if(LDEBUG){cerr << __AFLOW_FUNC__ << " directory=" << directory << endl;}
  //[CO20200508 - OBSOLETE]  
  //[CO20200508 - OBSOLETE]  string filename=aurostd::CleanFileName(directory+"/Makefile");
  //[CO20200508 - OBSOLETE]  aurostd::CleanStringASCII_InPlace(filename);
  //[CO20200508 - OBSOLETE]  if(!aurostd::FileExist(filename)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"File \""+filename+"\" does not exist",_VALUE_ILLEGAL_);}
  //[CO20200508 - OBSOLETE]  vector<string> vlines;
  //[CO20200508 - OBSOLETE]  aurostd::file2vectorstring(filename,vlines);
  //[CO20200508 - OBSOLETE]  return readMakefileVariables(vlines,vvariables,vvdefinitions);
  //[CO20200508 - OBSOLETE]}
  //[CO20200508 - OBSOLETE]void readMakefileVariables(const vector<string>& vlines,vector<string>& vvariables,vector<vector<string> >& vvdefinitions){
  //[CO20200508 - OBSOLETE]  bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
  //[CO20200508 - OBSOLETE]
  //[CO20200508 - OBSOLETE]  //find raw variables and definitions
  //[CO20200508 - OBSOLETE]  vector<string> vtokens,vdefinitions;
  //[CO20200508 - OBSOLETE]  string::size_type loc,loc_var_first,loc_var_last;
  //[CO20200508 - OBSOLETE]  string line="",variable="",definition="",subst_ref_string="";
  //[CO20200508 - OBSOLETE]  bool found_subst_ref=false;
  //[CO20200508 - OBSOLETE]  uint i=0,j=0;
  //[CO20200508 - OBSOLETE]  for(i=0;i<vlines.size();i++){
  //[CO20200508 - OBSOLETE]    line = vlines[i];
  //[CO20200508 - OBSOLETE]    loc=line.find("#");line=line.substr(0,loc);line=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(line);  //remove comments
  //[CO20200508 - OBSOLETE]    if(LDEBUG){cerr << __AFLOW_FUNC__ << " line=\"" << line << "\"" << endl;}
  //[CO20200508 - OBSOLETE]    if(line.size()>0 && line[0]=='\t'){continue;} //a recipe
  //[CO20200508 - OBSOLETE]    if(line.find("=")!=string::npos){
  //[CO20200508 - OBSOLETE]      //better to define definition first because you could have: TODAY=-DTODAY=\"$(TIME)\"
  //[CO20200508 - OBSOLETE]      //definition
  //[CO20200508 - OBSOLETE]      found_subst_ref=true;
  //[CO20200508 - OBSOLETE]      loc=string::npos;
  //[CO20200508 - OBSOLETE]      while(found_subst_ref){
  //[CO20200508 - OBSOLETE]        loc=line.find_last_of("=",loc);
  //[CO20200508 - OBSOLETE]        if(loc==string::npos){break;}
  //[CO20200508 - OBSOLETE]        //check that it is not a substitution ref: $(AFLOW_CPP:.cpp=.o)
  //[CO20200508 - OBSOLETE]        loc_var_first=line.find_last_of("$(",loc);
  //[CO20200508 - OBSOLETE]        if(loc_var_first==string::npos){break;}
  //[CO20200508 - OBSOLETE]        loc_var_last=line.find(")",loc);
  //[CO20200508 - OBSOLETE]        if(LDEBUG){
  //[CO20200508 - OBSOLETE]          cerr << __AFLOW_FUNC__ << " loc=" << loc << endl;
  //[CO20200508 - OBSOLETE]          cerr << __AFLOW_FUNC__ << " loc_var_first=" << loc_var_first << endl;
  //[CO20200508 - OBSOLETE]          cerr << __AFLOW_FUNC__ << " loc_var_last=" << loc_var_last << endl;
  //[CO20200508 - OBSOLETE]        }
  //[CO20200508 - OBSOLETE]        subst_ref_string=line.substr(loc_var_first-1,loc_var_last-loc_var_first+2);  //-1 because "$(" is length 2
  //[CO20200508 - OBSOLETE]        if(LDEBUG){cerr << __AFLOW_FUNC__ << " subst_ref_string=" << subst_ref_string << endl;}
  //[CO20200508 - OBSOLETE]        if(subst_ref_string.find(":")==string::npos && subst_ref_string.find("=")==string::npos){break;} //not a subst_ref
  //[CO20200508 - OBSOLETE]        loc-=1;
  //[CO20200508 - OBSOLETE]      }
  //[CO20200508 - OBSOLETE]      if(loc==string::npos){continue;}
  //[CO20200508 - OBSOLETE]      definition=line.substr(loc+1,string::npos); //+1 for equals sign
  //[CO20200508 - OBSOLETE]      definition=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(definition);
  //[CO20200508 - OBSOLETE]      if(0&&LDEBUG){cerr << __AFLOW_FUNC__ << " definition=\"" << definition << "\"" << endl;}
  //[CO20200508 - OBSOLETE]      //variable
  //[CO20200508 - OBSOLETE]      variable=line.substr(0,loc);  //no +1 because of equals sign
  //[CO20200508 - OBSOLETE]      variable=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(variable);
  //[CO20200508 - OBSOLETE]      if(variable.back()==':'){variable.pop_back();}  //in case the definition uses := instead of =
  //[CO20200508 - OBSOLETE]      aurostd::string2tokens(variable,vtokens,"="); //double variable assignment: TODAY=-DTODAY=\"$(TIME)\"
  //[CO20200508 - OBSOLETE]      for(j=0;j<vtokens.size();j++){
  //[CO20200508 - OBSOLETE]        variable=vtokens[j];
  //[CO20200508 - OBSOLETE]        variable=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(variable);
  //[CO20200508 - OBSOLETE]        if(variable.empty()){continue;}
  //[CO20200508 - OBSOLETE]        vvariables.push_back(variable);
  //[CO20200508 - OBSOLETE]        vdefinitions.push_back(definition);
  //[CO20200508 - OBSOLETE]        if(LDEBUG){
  //[CO20200508 - OBSOLETE]          cerr << __AFLOW_FUNC__ << " variable=\"" << vvariables.back() << "\"" << endl;
  //[CO20200508 - OBSOLETE]          cerr << __AFLOW_FUNC__ << " definition=\"" << vdefinitions.back() << "\"" << endl;
  //[CO20200508 - OBSOLETE]        }
  //[CO20200508 - OBSOLETE]      }
  //[CO20200508 - OBSOLETE]    }
  //[CO20200508 - OBSOLETE]  }
  //[CO20200508 - OBSOLETE]
  //[CO20200508 - OBSOLETE]  //perform replacements
  //[CO20200508 - OBSOLETE]  //NB: some variables appear twice, e.g., ARCH
  //[CO20200508 - OBSOLETE]  //so this replaces with the first instance (not necessarily the correct one)
  //[CO20200508 - OBSOLETE]  //this does not affect current functionality, as we only care about prerequisites
  //[CO20200508 - OBSOLETE]  //you are warned
  //[CO20200508 - OBSOLETE]  replaceMakefileDefinitions(vvariables,vdefinitions);
  //[CO20200508 - OBSOLETE]
  //[CO20200508 - OBSOLETE]  splitMakefileDefinitions(vdefinitions,vvdefinitions);
  //[CO20200508 - OBSOLETE]}

  void updateDependenciesAUROSTD(vector<string>& vdeps){
    bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vdeps(pre )=" << aurostd::joinWDelimiter(vdeps,",") << endl;}
    //first look for AUROSTD dependencies and replace with single AUROSTD/aurostd.o
    uint i=0;
    bool found_aurostd=false;
    for(i=vdeps.size()-1;i<vdeps.size();i--){  //go backwards so you can erase
      if(vdeps[i].find("AUROSTD/")!=string::npos){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " found AUROSTD dependencies" << endl;}
        vdeps.erase(vdeps.begin()+i);
        found_aurostd=true;
      }
    }
    if(found_aurostd){vdeps.insert(vdeps.begin(),"AUROSTD/aurostd.o");}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vdeps(post)=" << aurostd::joinWDelimiter(vdeps,",") << endl;}
  }
  void updateDependenciesVariable(const vector<string>& vdeps_var,const string& var,vector<string>& vdeps_replace){
    //ASSUMES vdeps_replace[0] IS $<
    bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vdeps_replace(pre )=" << aurostd::joinWDelimiter(vdeps_replace,",") << endl;}
    //first look for if ALL vdeps_var in vdeps_replace
    uint i=0,index_first=AUROSTD_MAX_UINT;
    bool found_all=true;
    for(i=0;i<vdeps_var.size()&&found_all;i++){
      if(!aurostd::WithinList(vdeps_replace,vdeps_var[i])){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " " << vdeps_var[i] << " NOT found in vdeps_replace" << endl;}
        found_all=false;
      }
    }
    if(!found_all){return;}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " found ALL $(" << var << ") in vdeps_replace" << endl;}
    for(i=vdeps_replace.size()-1;i<vdeps_replace.size();i--){  //go backwards so you can erase
      const string& dep=vdeps_replace[i];
      if(aurostd::WithinList(vdeps_var,dep)){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " found " << var << " dependency:" << dep << endl;}
        if(!(i==0 && dep!=vdeps_var[0])){ //SUPER PROTECTION FOR $< : protect $< at all costs! even if it means having duplicates in dependencies
          vdeps_replace.erase(vdeps_replace.begin()+i);
          if(i<index_first){index_first=i;}
        }
      }
    }
    //replace the variable at the first entry it was found, this ensures AUROSTD/aurostd.cpp stays at position $<
    vdeps_replace.insert(vdeps_replace.begin()+index_first,"$("+var+")");
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vdeps_replace(post)=" << aurostd::joinWDelimiter(vdeps_replace,",") << endl;}
  }
  void createMakefileAFLOW(const string& directory){
    bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " directory=" << directory << endl;}

    if(!aurostd::IsDirectory(directory+"/AUROSTD")){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,directory+"/AUROSTD not found",_VALUE_ILLEGAL_);}  //check that AUROSTD exists (fast check that ANY object files exist)
    if(!aurostd::FileExist(directory+"/aflow.h")){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,directory+"/aflow.h not found",_VALUE_ILLEGAL_);} //check we are in an aflow directory
    //[not a good idea, circular flow of information]if(!aurostd::FileExist(directory+"/aflow.o")){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,directory+"/aflow.o not found",_VALUE_ILLEGAL_);}  //check that aflow.o exists (fast check that ANY object files exist)

    //[CO20200508 - OBSOLETE]vector<string> vvariables;
    //[CO20200508 - OBSOLETE]vector<vector<string> > vvdefinitions;
    //[CO20200508 - OBSOLETE]readMakefileVariables(directory,vvariables,vvdefinitions);
    //[CO20200508 - OBSOLETE]if(LDEBUG){
    //[CO20200508 - OBSOLETE]  for(uint i=0;i<vvariables.size();i++){
    //[CO20200508 - OBSOLETE]    cerr << __AFLOW_FUNC__ << " variable  =\"" << vvariables[i] << "\"" << endl;
    //[CO20200508 - OBSOLETE]    cerr << __AFLOW_FUNC__ << " definition=\"" << aurostd::joinWDelimiter(vvdefinitions[i],",") << "\"" << endl;
    //[CO20200508 - OBSOLETE]  }
    //[CO20200508 - OBSOLETE]}

    vector<string> vsubdirectories;
    vsubdirectories.push_back(directory);
    vsubdirectories.push_back(directory+"/APL");
    if(COMPILE_HARDCODED_PROTOTYPES){vsubdirectories.push_back(directory+"/ANRL");} //DX20200623 - add if-statement
    vsubdirectories.push_back(directory+"/SQLITE");
    std::sort(vsubdirectories.begin(),vsubdirectories.end()); //sort before adding AUROSTD, which must come first

    //[do AUROSTD separately]vsubdirectories.insert(vsubdirectories.begin(),directory+"/AUROSTD");  //do first

    uint i=0,j=0;
    string dir="",file="",_file="";
    vector<string> vfs,vfiles,files_already_explored;
    vector<vector<string> > vvdependencies;
    vector<bool> vmt_required;
    bool mt_required=false;

    string Makefile_aflow="Makefile.aflow",Makefile_aflow_OLD="Makefile.aflow.OLD";

    //do AUROSTD first - it gets compiled separately
    file=directory+"/AUROSTD/aurostd.cpp";
    trimPath(file);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " building dependency for " << file << endl;}
    vfiles.push_back(file);
    vvdependencies.push_back(vector<string>(1, Makefile_aflow));  // ME20220210 - Add Makefile.aflow as dependency or files with changed flags/dependencies won't compile
    files_already_explored.clear();
    mt_required=false;
    getDependencies(vfiles.back(),files_already_explored,vvdependencies.back(),mt_required); //[didn't compile for some reason]vmt_required.back());
    vvdependencies.back().insert(vvdependencies.back().begin(),file);  //put file at the BEGINNING for $<
    vmt_required.push_back(mt_required);

    //SC variables - hack
    //AUROSTD
    string var_vcpp_aurostd="AUROSTD_CPPS",var_vhpp_aurostd="AUROSTD_HPPS";
    vector<string> vcpp_aurostd,vhpp_aurostd;
    //there are some AUROSTD files that are OBSOLETE, e.g., aurostd_xtensor_template.cpp
    //[GET vcpp_aurostd and vhpp_aurostd FROM aurostd.cpp dependencies, close the loop]dir=directory+"/AUROSTD";
    //[GET vcpp_aurostd and vhpp_aurostd FROM aurostd.cpp dependencies, close the loop]trimPath(dir);
    //[GET vcpp_aurostd and vhpp_aurostd FROM aurostd.cpp dependencies, close the loop]aurostd::DirectoryLS(dir,vfs);
    //[GET vcpp_aurostd and vhpp_aurostd FROM aurostd.cpp dependencies, close the loop]std::sort(vfs.begin(),vfs.end());
    //[GET vcpp_aurostd and vhpp_aurostd FROM aurostd.cpp dependencies, close the loop]for(j=0;j<vfs.size();j++){
    //[GET vcpp_aurostd and vhpp_aurostd FROM aurostd.cpp dependencies, close the loop]  file=dir+"/"+vfs[j];
    //[GET vcpp_aurostd and vhpp_aurostd FROM aurostd.cpp dependencies, close the loop]  if(aurostd::IsFile(file)){
    //[GET vcpp_aurostd and vhpp_aurostd FROM aurostd.cpp dependencies, close the loop]    if(vfs[j].size()>4 && vfs[j].find(".cpp")==vfs[j].size()-4){vcpp_aurostd.push_back(file);}
    //[GET vcpp_aurostd and vhpp_aurostd FROM aurostd.cpp dependencies, close the loop]    else if(vfs[j].size()>2 && vfs[j].find(".h")==vfs[j].size()-2){vhpp_aurostd.push_back(file);}
    //[GET vcpp_aurostd and vhpp_aurostd FROM aurostd.cpp dependencies, close the loop]  }
    //[GET vcpp_aurostd and vhpp_aurostd FROM aurostd.cpp dependencies, close the loop]}
    for(i=0;i<vvdependencies.back().size();i++){
      const string& dep=vvdependencies.back()[i];
      if(dep.find("AUROSTD")!=string::npos && dep.size()>4 && dep.find(".cpp")==dep.size()-4){vcpp_aurostd.push_back(dep);}
      else if(dep.find("AUROSTD")!=string::npos && dep.size()>2 && dep.find(".h")==dep.size()-2){vhpp_aurostd.push_back(dep);}
    }
    //stupidity test: make sure vcpp_aurostd[0]==directory+"/AUROSTD/aurostd.cpp" for $<
    if(vcpp_aurostd.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"vcpp_aurostd.size()==0",_RUNTIME_ERROR_);}
    if(vcpp_aurostd[0]!=file){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,file+" not set to first entry of vcpp_aurostd",_RUNTIME_ERROR_);}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vcpp_aurostd=" << aurostd::joinWDelimiter(vcpp_aurostd,",") << endl;}

    //DX20200801 - SYMBOLIC MATH - START
#if COMPILE_SYMBOLIC
    //do SYMBOLICCPLUSPLUS second - it gets compiled separately
    file=directory+"/SYMBOLICCPLUSPLUS/symbolic_main.cpp";
    trimPath(file);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " building dependency for " << file << endl;}
    vfiles.push_back(file);
    vvdependencies.push_back(vector<string>(1, Makefile_aflow));  // ME20220210 - Add Makefile.aflow as dependency or files with changed flags/dependencies won't compile
    files_already_explored.clear();
    mt_required=false;
    getDependencies(vfiles.back(),files_already_explored,vvdependencies.back(),mt_required); //[didn't compile for some reason]vmt_required.back());
    vvdependencies.back().insert(vvdependencies.back().begin(),file);  //put file at the BEGINNING for $<
    vmt_required.push_back(mt_required);

    //SYMBOLIC
    string var_vcpp_symbolic="SYMBOLIC_CPPS",var_vhpp_symbolic="SYMBOLIC_HPPS";
    //string var_vhpp_symbolic="SYMBOLIC_HPPS";
    vector<string> vcpp_symbolic,vhpp_symbolic;
    for(i=0;i<vvdependencies.back().size();i++){
      const string& dep=vvdependencies.back()[i];
      if(dep.find("SYMBOLICCPLUSPLUS")!=string::npos && dep.size()>4 && dep.find(".cpp")==dep.size()-4){vcpp_symbolic.push_back(dep);}
      else if(dep.find("SYMBOLICCPLUSPLUS")!=string::npos && dep.size()>2 && dep.find(".h")==dep.size()-2){vhpp_symbolic.push_back(dep);}
    }
    if(vcpp_symbolic.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"vcpp_symbolic.size()==0",_RUNTIME_ERROR_);}
    if(vcpp_symbolic[0]!=file){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,file+" not set to first entry of vcpp_symbolic",_RUNTIME_ERROR_);}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vcpp_symbolic=" << aurostd::joinWDelimiter(vcpp_symbolic,",") << endl;}
#endif
    //DX20200801 - SYMBOLIC MATH - END

    //ANRL
    vector<string> vcpp_anrl;
    if(COMPILE_HARDCODED_PROTOTYPES){
      dir=directory+"/ANRL";
      trimPath(dir);
      aurostd::DirectoryLS(dir,vfs);
      std::sort(vfs.begin(),vfs.end());
      for(j=0;j<vfs.size();j++){
        file=dir+"/"+vfs[j];
        if(aurostd::IsFile(file)){
          if(vfs[j].size()>4 && vfs[j].find(".cpp")==vfs[j].size()-4 && 
              vfs[j].find("aflow_anrl_A")!=string::npos && 
              vfs[j].find("aflow_anrl_sigma")!=string::npos &&  //DX - you may need to add here eventually
              true){vcpp_anrl.push_back(file);}
        }
      }
    }

    //do aflow.h dependencies next
    string var_vdep_aflowh="AFLOW_H_DEPS";
    vector<string> vdep_aflowh;
    file=directory+"/aflow.h";
    trimPath(file);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " building dependency for " << file << endl;}
    vdep_aflowh.push_back(file);  //add aflow.h to dependency
    files_already_explored.clear();
    getDependencies(file,files_already_explored,vdep_aflowh);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vdep_aflowh=" << aurostd::joinWDelimiter(vdep_aflowh,",") << endl;}

    vector<string> vhpp_aflow; //SC variables - hack
    vector<string> vdep_aflow; //populate me with missing skipped files

    bool skip_file=false; //these files do NOT get .o
    for(i=0;i<vsubdirectories.size();i++){
      dir=vsubdirectories[i];
      trimPath(dir);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " dir=" << dir << endl;}
      if(!aurostd::IsDirectory(dir)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,dir+" directory not found",_VALUE_ILLEGAL_);}
      aurostd::DirectoryLS(dir,vfs);
      std::sort(vfs.begin(),vfs.end());
      for(j=0;j<vfs.size();j++){
        file=dir+"/"+vfs[j];
        trimPath(file);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " file=" << file << endl;}
        if(aurostd::IsFile(file)){
          //if(file.find(".cpp")!=string::npos) //NO - ignore .cpp.orig
          if(file.size()>2 && file.find(".h")==file.size()-2){vhpp_aflow.push_back(file);} //SC variables - hack
          else if(file.size()>4 && file.find(".cpp")==file.size()-4){
            skip_file=false;
            //BEGIN skipping
            ////check that it has a header file
            //_file=file;aurostd::StringSubst(_file,".cpp",".h");if(!aurostd::FileExist(_file)){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING non-headed file=" << _file << endl;}skip_file=true;}
            //check that it's not a .cpp generated by another file
            _file=file;aurostd::StringSubst(_file,".cpp",".js");if(aurostd::FileExist(_file)||aurostd::FileExist("DATA/"+_file)){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING generated file=" << _file << endl;}skip_file=true;}  //CO20230507 - checking for DATA directory
            _file=file;aurostd::StringSubst(_file,".cpp",".json");if(aurostd::FileExist(_file)||aurostd::FileExist("DATA/"+_file)){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING generated file=" << _file << endl;}skip_file=true;}  //CO20230507 - checking for DATA directory
            _file=file;aurostd::StringSubst(_file,".cpp",".py");if(aurostd::FileExist(_file)||aurostd::FileExist("DATA/"+_file)){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING generated file=" << _file << endl;}skip_file=true;}  //CO20230507 - checking for DATA directory
            _file=file;aurostd::StringSubst(_file,".cpp",".txt");if(aurostd::FileExist(_file)||aurostd::FileExist("DATA/"+_file)){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING generated file=" << _file << endl;}skip_file=true;} //CO20230507 - checking for DATA directory
            //check that it's not a git file
            if(file.find("_BASE_")!=string::npos){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING git file=" << file << endl;}continue;} //continue, we don't want these in vdep_aflow
            if(file.find("_REMOTE_")!=string::npos){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING git file=" << file << endl;}continue;} //continue, we don't want these in vdep_aflow
            if(file.find("_LOCAL_")!=string::npos){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING git file=" << file << endl;}continue;}  //continue, we don't want these in vdep_aflow
            if(file.find("_BACKUP_")!=string::npos){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING git file=" << file << endl;}continue;} //continue, we don't want these in vdep_aflow
            //remove anything with aflow_data
            if(EXCLUDE_DATA){if(file.find("aflow_data")!=string::npos){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING " << file << endl;}continue;}} //continue, we don't want these in vdep_aflow
            //skip misc cpp files
            if(file.find("aflow_nomix.2014-01-15.cpp")!=string::npos){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING " << file << endl;}continue;} //continue, we don't want these in vdep_aflow  //no reason why we have this file in the build
            if(file.find("aflow_xproto_gus_lib.cpp")!=string::npos){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING " << file << endl;}skip_file=true;}
            if(file.find("aflow_test.cpp")!=string::npos){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING " << file << endl;}skip_file=true;}
            if(file.find("aflow_matlab_funcs.cpp")!=string::npos){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING " << file << endl;}skip_file=true;}
            if(file.find("aflow_gnuplot_funcs.cpp")!=string::npos){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING " << file << endl;}skip_file=true;}
            if(file.find("aflow_xpseudopotentials_data.cpp")!=string::npos){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING " << file << endl;}skip_file=true;}
            //[CO20200521 - OBSOLETE]if(file.find("aflow_aflowrc.cpp")!=string::npos){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING " << file << endl;}skip_file=true;}
            if(file.find("aflow_xproto_library_default.cpp")!=string::npos){if(LDEBUG){cerr << __AFLOW_FUNC__ << " SKIPPING " << file << endl;}skip_file=true;}
            if(skip_file){
              vdep_aflow.push_back(file);  //add to aflow dependencies
              continue;
            }
            //END skipping
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " building dependency for " << file << endl;}
            vfiles.push_back(file);
            vvdependencies.push_back(vector<string>(1, Makefile_aflow));  // ME20220210 - Add Makefile.aflow as dependency or files with changed flags/dependencies won't compile
            files_already_explored.clear();
            //[didn't compile for some reason]vmt_required.push_back(false);
            mt_required=false;
            getDependencies(vfiles.back(),files_already_explored,vvdependencies.back(),mt_required); //[didn't compile for some reason]vmt_required.back());
            vvdependencies.back().insert(vvdependencies.back().begin(),file);  //put file at the BEGINNING for $<
            //[CO20200508 - OBSOLETE]//unfortunate hacks that are needed
            //[CO20200508 - OBSOLETE]if(file=="aflowlib_libraries.cpp"){
            //[CO20200508 - OBSOLETE]  if(!aurostd::WithinList(vvdependencies.back(),"aflow_matlab_funcs.cpp")){vvdependencies.back().push_back("aflow_matlab_funcs.cpp");}  //safety
            //[CO20200508 - OBSOLETE]  if(!aurostd::WithinList(vvdependencies.back(),"aflow_gnuplot_funcs.cpp")){vvdependencies.back().push_back("aflow_gnuplot_funcs.cpp");}  //safety
            //[CO20200508 - OBSOLETE]}
            vmt_required.push_back(mt_required);
            //get rid of any dependencies with AUROSTD and replace with $(AUROSTD_OBJ)
            //[CO20200521 - do NOT replace AUROSTD .cpp and .h with .o, this will slow down make unnecessarily]updateDependenciesAUROSTD(vvdependencies.back());
            updateDependenciesVariable(vdep_aflowh,var_vdep_aflowh,vvdependencies.back());
          }
        }
      }
    }
    //[CO20200508 - OBSOLETE]//unfortunate hacks that are needed
    //[CO20200508 - OBSOLETE]if(!aurostd::WithinList(vdep_aflow,"aflow_matlab_funcs.cpp")){vdep_aflow.push_back("aflow_matlab_funcs.cpp");}  //safety
    //[CO20200508 - OBSOLETE]if(!aurostd::WithinList(vdep_aflow,"aflow_gnuplot_funcs.cpp")){vdep_aflow.push_back("aflow_gnuplot_funcs.cpp");}  //safety
    stringstream makefile_rules_ss;
    vector<string> vfile_obj;
    string obj_file="";
    for(i=0;i<vfiles.size();i++){
      if(vvdependencies[i].empty()){continue;} //we will define generic one at the end
      obj_file=vfiles[i];aurostd::StringSubst(obj_file,".cpp",".o");
      if(vfiles[i].find("aflow_data")==string::npos){ //CO20200508 - do not include anything for aflow_data in aflow dependencies
        vfile_obj.push_back(obj_file);
        vdep_aflow.push_back(obj_file);
      }
      //[not a good idea, circular flow of information]if(!aurostd::FileExist(obj_file)){continue;}  //since we can only run this code once aflow is compiled, we can check that the obj_file is a real target or not
      //[CO20200521 - OBSOLETE, we already injected vfiles[i] above]vvdependencies[i].insert(vvdependencies[i].begin(),vfiles[i]); //put .cpp at the front of dependencies for $<
      if(obj_file.find("aurostd.o")!=string::npos){ //exception for aurostd.o
        updateDependenciesVariable(vcpp_aurostd,var_vcpp_aurostd,vvdependencies[i]);
        updateDependenciesVariable(vhpp_aurostd,var_vhpp_aurostd,vvdependencies[i]);
      }
      //DX20200801 - SYMBOLIC MATH - START
#if COMPILE_SYMBOLIC
      if(obj_file.find("symbolic_main.o")!=string::npos){ //exception for symbolic_main.o
        updateDependenciesVariable(vcpp_symbolic,var_vcpp_symbolic,vvdependencies[i]);
        updateDependenciesVariable(vhpp_symbolic,var_vhpp_symbolic,vvdependencies[i]);
        updateDependenciesVariable(vhpp_aurostd,var_vhpp_aurostd,vvdependencies[i]);
      }
      if(obj_file.find("aflow_symbolic.o")!=string::npos || obj_file.find("aflow_anrl.o")!=string::npos){ //exception for aflow_symbolic.o and aflow_anrl.o
        updateDependenciesVariable(vhpp_symbolic,var_vhpp_symbolic,vvdependencies[i]);
      }
#endif
      //DX20200801 - SYMBOLIC MATH - END
      //makefile_rules_ss << obj_file << ": " << vfiles[i] << " " << aurostd::joinWDelimiter(vvdependencies[i]," ") << endl;
      makefile_rules_ss << obj_file << ": " << aurostd::joinWDelimiter(vvdependencies[i]," ") << endl;
      makefile_rules_ss << "\t" << "$(CPP) $(VERS) $(INCLUDE) $(CCFLAGS" << (vmt_required[i] ? "_MT" : "") << ") $(OPTS" << (vmt_required[i] ? "_MT" : "") << ") $(ARCH) $< -c -o $@" << endl;  //(obj_file=="AUROSTD/aurostd.o"?"$^":"$<")  //ME20200514 - Added CCFLAGS_MT
    }
    //[CO20200521 - OBSOLETE]//unfortunate hacks that are needed
    //[CO20200508 - OBSOLETE]if(!aurostd::WithinList(vfile_obj,"aflow_matlab_funcs.cpp")){vfile_obj.push_back("aflow_matlab_funcs.cpp");}  //safety
    //[CO20200508 - OBSOLETE]if(!aurostd::WithinList(vfile_obj,"aflow_gnuplot_funcs.cpp")){vfile_obj.push_back("aflow_gnuplot_funcs.cpp");}  //safety
    //[CO20200521 - OBSOLETE]if(!aurostd::WithinList(vfile_obj,"aflow_aflowrc.cpp")){vfile_obj.push_back("aflow_aflowrc.cpp");}  //safety

    //make some variable replacements before printing
    updateDependenciesVariable(vcpp_aurostd,var_vcpp_aurostd,vdep_aflowh);
    updateDependenciesVariable(vhpp_aurostd,var_vhpp_aurostd,vdep_aflowh);

    //this approach is not sustainable - do you think SC won't demand more hacks?
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]//unfortunate hacks that are needed
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]//vdep_aflowh needs aflow_aflowrc.cpp first for $<
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]//aflow_aflowrc.cpp is half .cpp half .h
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]//since it has the same dependencies as aflow.h, we need to make sure the AFLOW_H_DEPS variables is defined with aflow_aflowrc.cpp first
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]//if another file 
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]file=directory+"/aflow_aflowrc.cpp";
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]trimPath(file);
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]if(vdep_aflowh.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"vdep_aflowh.size()==0",_RUNTIME_ERROR_);}
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]bool found=false;
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]for(i=0;i<vdep_aflowh.size()&&found==false;i++){
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]  const string& dep=vdep_aflowh[i];
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]  if(dep.find(file)!=string::npos){
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]    if(i!=0){
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]      vdep_aflowh.erase(vdep_aflowh.begin()+i);
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]      vdep_aflowh.insert(vdep_aflowh.begin(),file);
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]    }
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]    found=true;
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]  }
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]}
    //[CO20200521 - OBSOLETE WITH SUPER PROECTION FOR $< ABOVE]if(vdep_aflowh[0]!=file){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,file+" not set to first entry of vdep_aflowh",_RUNTIME_ERROR_);}

    stringstream makefile_definitions_ss;
    //MAIN DEPENDENCIES START
    if(vcpp_aurostd.size()){makefile_definitions_ss << var_vcpp_aurostd << "=" << aurostd::joinWDelimiter(vcpp_aurostd," ") << endl;} //SC variables - hack
    if(vhpp_aurostd.size()){makefile_definitions_ss << var_vhpp_aurostd << "=" << aurostd::joinWDelimiter(vhpp_aurostd," ") << endl;} //SC variables - hack
#if COMPILE_SYMBOLIC
    if(vcpp_symbolic.size()){makefile_definitions_ss << var_vcpp_symbolic << "=" << aurostd::joinWDelimiter(vcpp_symbolic," ") << endl;} //DX20200831 - symbolic math
    if(vhpp_symbolic.size()){makefile_definitions_ss << var_vhpp_symbolic << "=" << aurostd::joinWDelimiter(vhpp_symbolic," ") << endl;} //DX20200831 - symbolic math
#endif
    if(vdep_aflowh.size()){makefile_definitions_ss << var_vdep_aflowh << "=" << aurostd::joinWDelimiter(vdep_aflowh," ") << endl;} //SC variables - hack
    //MAIN DEPENDENCIES END
    if(vdep_aflow.size()){makefile_definitions_ss << "AFLOW_DEPS=" << aurostd::joinWDelimiter(vdep_aflow," ") << endl;}  //SC variables - hack
    if(vfile_obj.size()){makefile_definitions_ss << "AFLOW_OBJS=" << aurostd::joinWDelimiter(vfile_obj," ") << endl;}  //SC variables - hack
    if(vhpp_aflow.size()){makefile_definitions_ss << "AFLOW_HPPS=" << aurostd::joinWDelimiter(vhpp_aflow," ") << endl;} //SC variables - hack
    if(COMPILE_HARDCODED_PROTOTYPES&&vcpp_anrl.size()){makefile_definitions_ss << "ANRL_CPPS=" << aurostd::joinWDelimiter(vcpp_anrl," ") << endl;} //SC variables - hack //DX20200623 - added USE_HARDCODED_PROTOTYPES

    //join together
    stringstream makefile_aflow;
    makefile_aflow << "# AUTOMATICALLY GENERATED AFLOW - BEGIN" << endl << endl;
    makefile_aflow << makefile_definitions_ss.str() << endl;
    makefile_aflow << makefile_rules_ss.str() << endl;
    makefile_aflow << "# AUTOMATICALLY GENERATED AFLOW - END" << endl;
    if(aurostd::FileExist(directory+"/"+Makefile_aflow)){aurostd::file2file(directory+"/"+Makefile_aflow,directory+"/"+Makefile_aflow_OLD);}  //saving
    aurostd::stringstream2file(makefile_aflow,directory+"/"+Makefile_aflow);
  }
}

#endif  // _AFLOW_MAKEFILE_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow COREY OSES - Duke University 2013-2021                  *
// *                                                                         *
// ***************************************************************************
