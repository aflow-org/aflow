// ***************************************************************************
// *                                                                         *
// *               AFlow SHIDONG WANG - Duke University 2010-2011            *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_shidong_funs.cpp
// functions written by
// 2010-2011: shidong.wang@duke.edu

#ifndef _AFLOW_CONTRIB_SHIDONG_AUXILIARY_CPP_
#define _AFLOW_CONTRIB_SHIDONG_AUXILIARY_CPP_

#include "aflow.h"
#include "aflow_contrib_shidong_auxiliary.h"
#include <fcntl.h>

//using namespace std;

//[CO20200624 - OBSOLETE]//*************************************
//[CO20200624 - OBSOLETE]// Output the error message to cerr
//[CO20200624 - OBSOLETE]//*************************************
//[CO20200624 - OBSOLETE]void ErrorMessage(const string & errmsg, int errnr) {
//[CO20200624 - OBSOLETE]  cerr << errmsg << endl;
//[CO20200624 - OBSOLETE]  cerr << "Abort the program ... \n";
//[CO20200624 - OBSOLETE]}

//*************************************
// Wrapper for generate new training item
//*************************************

void GenerateNewTrainingItem(string & SL_name) {
#ifdef TESTCE
  CalculateNewStateTest(SL_name);
# else
  CalculateNewStateAFLOW(SL_name);
#endif
}

//*************************************
// Generate new training item
//*************************************
void CalculateNewStateAFLOW(string & SL_name) {

  pid_t pID;
  int child_exit_status;

  // creat a directory for aflow calculation

  aurostd::execute("rm -rf " + SL_name);
  aurostd::execute("mkdir " + SL_name);

  // setup the input file
  aurostd::execute("cp "+_AFLOWIN_+" " + SL_name);

  // ***********************************************
  // for Unix/Linux
  // ***********************************************
  // first get the current working directory
  //char *pwd_old, *pwd;  //CO20181019 - avoid strncat warning of size()
  //int pwd_old_size;

  string pwd,pwd_old;
  pwd=pwd_old=aurostd::getPWD(); //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd");
  //pwd_old = getcwd(NULL, 0);  //CO20181019 - avoid strncat warning of size()
  //pwd_old_size = strlen(pwd_old);
  //pwd = new char[strlen(pwd_old) + SL_name.size()+2]; //CO20181019 - avoid strncat warning of size()
  //strcpy(pwd, pwd_old); //CO20181019 - avoid strncat warning of size()
  //strncat(pwd, "/", 1); //CO20181019 - avoid strncat warning of size()
  //strncat(pwd, SL_name.c_str(), SL_name.size());  //CO20181019 - avoid strncat warning of size()
  pwd+='/';
  pwd+=SL_name;

  pID= fork();
  if( pID < 0 ) {
    // fail to fork
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"fail to fork",_RUNTIME_ERROR_); //CO20200624
  } else if( pID == 0 ) {
    // child process

    // change to the new working directory
    chdir(pwd.c_str()); //CO20181019 - avoid strncat warning of size()

    // run the program generate new state
    char cmd_name[] = "aflow";
    char *args_list[] = {NULL, NULL};

    args_list[0] = cmd_name;

    // redirect stdout and stderr in child process
    int fout = open("out.dat", O_WRONLY | O_CREAT, 0600); //FK20180515
    int ferr = open("err.dat", O_WRONLY | O_CREAT, 0600); //FK20180515

    dup2(fout, STDOUT_FILENO);
    dup2(ferr, STDOUT_FILENO);

    close(fout);
    close(ferr);

    execvp(cmd_name, args_list);

    // clean pointers
    //args_list[0] = NULL;

  } else {
    waitpid(pID, &child_exit_status, 0);
    //delete [] pwd;  //CO20181019 - avoid strncat warning of size()
    //delete [] pwd_old;  //CO20181019 - avoid strncat warning of size()

    aurostd::execute("cd "+SL_name+" &&  chmod "+ _MODE + " *");
    aurostd::execute("cp " + SL_name + "/" + AFLOW_RESULT_FILE + ".bz2 ./");
    cerr << "parent finishes preparing new fitting structure \n";
  }
}

//void CalculateNewStateTest(string & SL_name)
//{
//
//    pid_t pID;
//    int child_exit_status;
//    string command;
//
//    // creat a directory for aflow calculation
//
//    aurostd::execute("rm -rf " + SL_name);
//
//    aurostd::execute("mkdir " + SL_name);
//
//    // setup the input file
//    aurostd::execute("cp "+_AFLOWIN_+" " + SL_name);
//
//    // ***********************************************
//    // for Unix/Linux
//    // ***********************************************
//    // first get the current working directory
//    char *pwd_old, *pwd;
//    int pwd_old_size;
//
//    pwd_old = getcwd(NULL, NULL);
//    pwd_old_size = strlen(pwd_old);
//    pwd = new char[strlen(pwd_old) + SL_name.size()+2];
//    strcpy(pwd, pwd_old);
//    strncat(pwd, "/", 1);
//    strncat(pwd, SL_name.c_str(), SL_name.size());
//
//
//    pID= fork();
//    if( pID < 0 ) {
//        // fail to fork
//        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"fail to fork",_RUNTIME_ERROR_); //CO20200624
//    } else if( pID == 0 ) {
//        // child process
//
//        // change to the new working directory
//        chdir(pwd);
//
//        // run the program generate new state
//        char *cmd;
//        char cmd_name[] = "generateState.py";
//        int cmd_name_size = strlen(cmd_name);
//        char *args_list[] = {NULL, NULL, NULL};
//
//        cmd = new char[pwd_old_size + cmd_name_size + 2];
//        strcpy(cmd, pwd_old);
//        strncat(cmd, "/", 1);
//        strncat(cmd, cmd_name, cmd_name_size);
//
//        args_list[0] = cmd;
//
//        args_list[1] = new char[SL_name.size()+1];
//        strcpy(args_list[1], SL_name.c_str());
//
//        execvp(cmd, args_list);
//
//        // clean pointers
//        args_list[0] = NULL;
//        delete [] cmd;
//        delete [] args_list[1];
//
//    } else {
//        waitpid(pID, &child_exit_status, 0);
//        delete [] pwd;
//        delete [] pwd_old;
//
//        aurostd::execute("cd " + SL_name + " && chmod "+ _MODE + " *");
//        aurostd::execute("cp " + SL_name + "/" + AFLOW_RESULT_FILE + ".bz2 ./");
//        cerr << "parent finishes preparing new fitting structure \n";
//    }
//}

//*************************************
// Manipulate generated files
//*************************************

double GetResultFromAFLOW() {
  // get the calculated formation_energy_per_atom from aflow
  // the results are read from DEFAULT_AFLOW_QMVASP_OUT.bz2
  // commands used here depend on the output format of aflow

  // change everything back to the mode "655"

  // unzip the result file
  aurostd::execute("bzip -dq " + AFLOW_RESULT_FILE + ".bz2");

  ifstream fin_result;
  fin_result.open(AFLOW_RESULT_FILE.c_str());

  string line;

  string _anchor = "E0/N", _anchor2 = "=";

  double energy = 0;

  while ( !fin_result.eof() && getline(fin_result, line) ) {
    if( line.find(_anchor) != string::npos ) {
      int pos = line.find(_anchor2);
      // remove the leading spaces
      line.erase(line.begin(), line.begin()+pos+1);
      for (uint i=0; i<line.size(); i++) {
        if( line.at(i) == ' ' ) {
          line.erase(line.begin()+i);
        } else {
          break;
        }
      }
      energy = aurostd::string2utype<long double>(line);
      break;
    }
  }
  aurostd::execute("bzip2 " + AFLOW_RESULT_FILE);

  return energy;

}

void GenerateGNUplotScript() {

  ofstream fout;
  fout.open(_GNUPLOTPTFILE.c_str());

  fout << "set terminal postscript eps enhanced color solid" << endl;
  fout << "set output \"" + _GNUPLOTPTOUTFILE + "\" " << endl;
  fout << endl;
  fout << "set xlabel \"x_b\"" << endl;
  fout << "set ylabel \"enthalpy\"" << endl;
  fout << "plot \"" + _TOTALPTFILE + "\" u 2:3 w p pt 6 lc 1, \""
    + _HULLPTFILE +"\"  u 2:3 w l \\" << endl;
  fout << ", \"" + _FITCOMPARISONFILE + "\" u 2:3 w p pt 2 lc 3, \""
    + _FITCOMPARISONFILE+ "\" u 2:4 w p pt 3 lc 5,  \\" << endl;
  fout << "\"" + _RALLOYRESULT + "\" u 2:4 w l lc 3 lw 2" << endl;
  fout << endl;

  fout.close();
}

//*************************************
// deal with files
//*************************************

void RenameFiles(int count) {
  // rename all output file by attaching the calculation number
  vector<string> backup_file_names;

  string filename;

  backup_file_names.push_back(_TOTALPTFILE);
  backup_file_names.push_back(_HULLPTFILE);
  backup_file_names.push_back(_GNDPTFILE);
  backup_file_names.push_back(_GNUPLOTPTOUTFILE);
  backup_file_names.push_back(_FITSTRUCTUREFILE);
  backup_file_names.push_back(_RALLOYRESULT);
  backup_file_names.push_back(_FITCOMPARISONFILE);
  backup_file_names.push_back(_ECIFILERESULT);
  backup_file_names.push_back(_SLFILERESULT);

  //string affix = aurostd::utype2string(count);

  string com;
  vector<string>::iterator name_itr;
  for (name_itr = backup_file_names.begin();
      name_itr < backup_file_names.end(); name_itr++) {
    string old_name = *name_itr;

    if( old_name == _FITSTRUCTUREFILE ) {
      com = "cp";
    } else {
      com = "mv";
    }

    BackUpFile(old_name, com, count);
  }

}

void BackUpFile(string old_name, string com, int count) {
  // change file name from xxx.xx to xxx_count.xx

  string postfix = aurostd::utype2string(count);

  int pos = old_name.find(".");

  string new_name;
  new_name = old_name.substr(0,pos) + "_" + postfix + old_name.substr(pos);

  string command;
  aurostd::execute(com + " " + old_name + " " + new_name);

}

void MoveFiles(int count) {
  // create a directory with name "count" and move
  // all files end with "_count.dat" to it
  string postfix = aurostd::utype2string(count);
  string file_wildcard = "*_" + postfix + ".*";
  // remove the old directory if existed
  aurostd::execute("rm -rf " + postfix);
  aurostd::execute("mkdir " + postfix);
  aurostd::execute("mv " + file_wildcard + " " + postfix);
  aurostd::execute("cp "+_AFLOWIN_+" " + postfix);
  aurostd::execute("gzip -r " + postfix);
}

void CleanUp() {
  // Do not care about the success of removing files
  aurostd::RemoveFile(_LOCVFILE);
  aurostd::RemoveFile(_SCOREFILE);
  aurostd::RemoveFile(_AFLOWIN_);
}

#endif

