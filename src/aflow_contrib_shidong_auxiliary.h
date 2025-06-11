// ***************************************************************************
// *                                                                         *
// *               AFlow SHIDONG WANG - Duke University 2010-2011            *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_shidong_funs.cpp
// functions written by
// 2010-2011: shidong.wang@duke.edu
//

#ifndef _AFLOW_CONTRIB_SHIDONG_AUXILIARY_H_
#define _AFLOW_CONTRIB_SHIDONG_AUXILIARY_H_

#include <iostream>
#include <iomanip>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <cstdlib>
using namespace std;

const string _SLFILE="SL.dat";
const string _SLFILENAME="SLname.dat";
const string _SLFILERESULT="SLresult.dat";
const string _SLFILECOR="SLcor.dat";
const string _SLFILESQS="SLSQS.dat";

const string _GNUPLOTPTFILE="ce-enthalpy.gp";
const string _FITSTRUCTUREFILE = "fitstructure.dat";

const string _TOTALPTFILE="total_pt.dat";
const string _HULLPTFILE="hull_pt.dat";
const string _GNDPTFILE="groundstate_pt.dat";
const string _GNUPLOTPTOUTFILE="enthalpy.eps";

const string _RALLOYRESULT="ralloy.dat";
const string _FITCOMPARISONFILE="fit_comparison.dat";
const string _ECIFILERESULT="eci_result.dat";

const string _LOCVFILE="LOCV.dat";
const string _SCOREFILE="ga_score.dat";

void ErrorMessage(const string & errmsg, int errnr);
const int _EXIT_FAIL_FORK = 100;
const string AFLOW_RESULT_FILE = "aflow.qmvasp.out";
const string _MODE = "644"; // mode rw-r--r--

void GenerateNewTrainingItem(string & SL_name);
void CalculateNewStateTest(string & SL_name);
void CalculateNewStateAFLOW(string & SL_name);
void BackUpFile(string old_name, string com, int count);
void MoveFiles(int count);
double GetResultFromAFLOW();
void RenameFiles(int count);
void CleanUp();

void GenerateGNUplotScript();
#endif
