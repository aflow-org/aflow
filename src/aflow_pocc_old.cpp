// ***************************************************************************
// *                                                                         *
// *              AFlow KESONG YANG  Duke University 2010-2011               *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_kesong_ipocc.cpp
// functions written by KESONG YANG
// 2010-2011: kesong.yang@gmail.com

// This file contains the routines to prepare partial occupation input files.

#ifndef _AFLOW_CONTRIB_KESONG_IPOCC_CPP
#define _AFLOW_CONTRIB_KESONG_IPOCC_CPP

//#include "aflow_contrib_kesong.h"
#include "aflow_pocc_old.h"
#include "aflow_pocc.h" //CO20180409

//ipocc.cpp
using aurostd::StringSubst;
using std::setfill;
//const int MaxNumberPOSCAR=999;
const unsigned long long int MaxNumberPOSCAR=999;

//pocc_basic.cpp
const double Tol_DecimalPart=1E-6;
const double Tol_FracDeciDiff=1E-6;
//This is the optimized tolerance
const int MaxDenominator=999999;

//hnfcell.cpp
bool CALCULATE_ENERGY_UFF = true;

//poccupation_forcefield.cpp
//const double unit_kcaltoev=4.336443203200000E-002; // 1(kcal/mol) = 4.33644E-2 eV
const double KCAL_TO_EV =4.336443203200000E-002; // 1(kcal/mol) = 4.33644E-2 eV
//const double KJ_TO_EV =1.0357416650425146E-002 ; // 1 (KJ/mol) = 0.010357 eV
//const double KCAL_TO_KJ = 4.1868; // 1 Kilocalories (Kcal) = 4.1868 Kilojoules (Kj)

bool PRINT_AVG_IDOS=false; //CO20190808 - quick patch for IDOS, keep this off in general

// ***************************************************************************
// pflow::POCC_INPUT(void)
// ***************************************************************************
namespace pflow {
  void POCC_INPUT(void) {
    _aflags aflags;aflags.Directory="./";
    ofstream oss;
    string aflowin, MESSAGE="pflow::POCC_INPUT ERROR";
    aflowin=string(aflags.Directory+_AFLOWIN_);
    if(!aurostd::FileExist(aflowin)) {
      throw aurostd::xerror(__AFLOW_FILE__,"pflow::POCC_INPUT():","file not found: "+aflowin,_FILE_CORRUPT_); //CO20200624
    }
    POCC_GENERATE_INPUT(oss,aflags);
  }
}

// ***************************************************************************
// bool POCC_GENERATE_INPUT(ofstream &FileMESSAGE,_aflags &aflags)
// ***************************************************************************
bool POCC_GENERATE_INPUT(ofstream &FileMESSAGE,_aflags &aflags) {
  ostringstream aus;
  aus << "00000  MESSAGE running POCC_GENERATE_INPUT files " << Message(__AFLOW_FILE__,aflags) << endl;
  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

  string aflowin=aflags.Directory+"/"+_AFLOWIN_;

  if(!aurostd::FileExist(aflowin)) {
    FileMESSAGE << "ERROR" << ": file not found " << aflowin << endl;
    return FALSE;
  }

  string AflowIn;aurostd::file2string(aflowin,AflowIn);
  stringstream sspoccSTRUCTURE;

  // CHECK FOR INSIDE STUFF
  if(!pocc::POCC_Already_Calculated_Input(AflowIn,aflags.Directory)) {
    stringstream input_file; input_file.clear();input_file.str("");
    stringstream input_file_aus; input_file_aus.clear();input_file_aus.str("");

    sspoccSTRUCTURE.str(aurostd::substring2string(AflowIn,"[POCC_MODE_EXPLICIT]START.POCC_STRUCTURE","[POCC_MODE_EXPLICIT]STOP.POCC_STRUCTURE"));
    //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn,sspoccSTRUCTURE, "[POCC_MODE_EXPLICIT]START.POCC_STRUCTURE", "[POCC_MODE_EXPLICIT]STOP.POCC_STRUCTURE");

    xstructure  xstr_pocc(sspoccSTRUCTURE, IOVASP_POSCAR);
    ofstream FileMESSAGE;
    ofstream file_aflowin;
    file_aflowin.open(aflowin.c_str(),ios_base::app);

    file_aflowin << AFLOWIN_SEPARATION_LINE << endl; 
    if(!aurostd::substring2bool(AflowIn,"[AFLOW_POCC]CALC")) {
      file_aflowin << "[AFLOW_POCC]CALC "<< endl;
    }

    vector<xstructure> vecgroupxstr_sorted;
    stringstream ssxstr_sorted;
    stringstream ss;
    ssxstr_sorted.str("");

    unsigned long long int Num_calculated;
    unsigned long long int Num_xstr;

    if(1){  //OLD KY
      vecgroupxstr_sorted = Partial2Supercell(xstr_pocc);  //CO20170721 - KY OLD code
      Num_xstr=vecgroupxstr_sorted.size();

      if(Num_xstr<MaxNumberPOSCAR) {
        Num_calculated=Num_xstr;
      }
      else {
        Num_calculated=MaxNumberPOSCAR;
      }

      for(unsigned long long int i=0;i<Num_calculated;i++) {
        ss.str("");
        ss << "POCC_" << setfill('0') << setw(2) <<(i+1);
        ssxstr_sorted << AFLOWIN_SEPARATION_LINE<< endl;
        ssxstr_sorted << _VASP_POSCAR_MODE_EXPLICIT_START_ << ss.str() << endl;  //CO20200624
        ssxstr_sorted << vecgroupxstr_sorted.at(i);
        ssxstr_sorted << _VASP_POSCAR_MODE_EXPLICIT_STOP_ << ss.str() << endl; //CO20200624
        ssxstr_sorted << AFLOWIN_SEPARATION_LINE<< endl;
      }
    } else {  //START CO
      pocc::POccCalculator pcalc(xstr_pocc,aflags,FileMESSAGE,cout);
      if(!pcalc.m_initialized){
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POCC algorithm failed to initialize",_RUNTIME_ERROR_); //CO20200624
      }
      pcalc.calculate();

      //vector<xstructure> vecgroupxstr_sorted = pcalc.getUniqueDerivativeStructures(); //too much memory

      Num_xstr=pcalc.getUniqueSuperCellsCount();//vecgroupxstr_sorted.size();

      //NO NEED TO LIMIT HERE, just print everything out, who cares?
      //if (Num_xstr<MaxNumberPOSCAR) {
      Num_calculated=Num_xstr;
      //}
      //else {
      //    Num_calculated=MaxNumberPOSCAR;
      //}

      stringstream message;
      ostream& oss=cout;
      //ofstream FileMESSAGE;
      message << "Creating list of (primitivized) unique derivative supercells. Please be patient as primitivization can be slow";
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,FileMESSAGE,oss,_LOGGER_MESSAGE_);

      for(unsigned long long int i=0;i<Num_calculated;i++) {
        ss.str("");
        ss << "POCC_" << setfill('0') << setw(2) <<(i+1);
        ssxstr_sorted << AFLOWIN_SEPARATION_LINE<< endl;
        ssxstr_sorted << _VASP_POSCAR_MODE_EXPLICIT_START_ << ss.str() << endl;  //CO20200624
        ssxstr_sorted << pcalc.getUniqueSuperCell(i); //vecgroupxstr_sorted.at(i);
        ssxstr_sorted << _VASP_POSCAR_MODE_EXPLICIT_STOP_ << ss.str() << endl; //CO20200624
        ssxstr_sorted << AFLOWIN_SEPARATION_LINE<< endl;
      }
    }

    file_aflowin << ssxstr_sorted.str();
    file_aflowin.close();
  } 
  else {
    aus << "00000  MESSAGE POCC input file already created " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    return FALSE;
  }
  return TRUE;
}

// ***************************************************************************
// void remove_duplicates(vector<T>& vec)
// ***************************************************************************
template <typename T>
void remove_duplicates(vector<T>& vec) {
  sort(vec.begin(), vec.end());
  vec.erase(unique(vec.begin(), vec.end()), vec.end());
}

str_num_data OptimizePoccValue(double dvalue, double tol) {
  str_num_data optimized_dvalue;

  if(tol<1E-6) {
    tol = 0.001;
  }

  int Nmax=100;
  vector<double> vec_double;
  for (int i=1; i<=Nmax; i++) {
    for (int j=i+1; j<=Nmax;j++) {
      double tmp= (1.0*i)/j;
      vec_double.push_back(tmp);
    }
  }
  vec_double.push_back(1.0);
  sort(vec_double.begin(), vec_double.end());
  remove_duplicates(vec_double);

  vector<double> double_range;
  for (uint i=0; i<vec_double.size();i++) {
    double num_diff = abs(vec_double.at(i)-dvalue);
    if(num_diff<=tol) {
      double_range.push_back(vec_double.at(i));
    }
  }

  vector<int> frac_tmp;
  vector<str_num_data> vec_str_num_data;
  str_num_data tmp_str_num_data;
  for (uint i=0; i<double_range.size();i++) {
    double db_tmp =  double_range.at(i);
    tmp_str_num_data = double2str_num_data(db_tmp);
    vec_str_num_data.push_back(tmp_str_num_data);
  }

  sort(vec_str_num_data.begin(), vec_str_num_data.end(), sortdenom); //sorting according to denom
  optimized_dvalue = vec_str_num_data.at(0); 


  return optimized_dvalue;
}

// ***************************************************************************
// void UpdateXstructure(xstructure &b)
// ***************************************************************************
void UpdateXstr_comp_each_type(xstructure &b) {
  //This function is only used to update partial occupation value i.e., comp_each_type
  vector<double> poccaus; // partial occupation local host
  for (uint i=0; i<b.atoms.size();i++) {
    double tmp=b.atoms.at(i).partial_occupation_value;
    poccaus.push_back(tmp);
  }

  b.comp_each_type.clear();
  for(uint i=0;i<b.num_each_type.size();i++) {
    b.comp_each_type.push_back(0.0);
  }

  for(uint i=0;i<b.atoms.size();i++) {
    if(poccaus.at(i)<1.0-1e-5) {
      b.atoms.at(i).partial_occupation_flag=TRUE;
      b.atoms.at(i).partial_occupation_value=poccaus.at(i);
      b.comp_each_type.at(b.atoms.at(i).type)+=b.atoms.at(i).partial_occupation_value;
    } else {
      b.atoms.at(i).partial_occupation_flag=FALSE;
      b.atoms.at(i).partial_occupation_value=1.0;
      b.comp_each_type.at(b.atoms.at(i).type)+=b.atoms.at(i).partial_occupation_value;
    }
  }
}


// ***************************************************************************
// void UpdateXstr(xstructure &xstr_orig, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
void UpdateXstr(xstructure &xstr_orig, ofstream &FileMESSAGE, _aflags &aflags) {
  double epsilon = 1E-2;
  //Update paritial occupation value
  if(xstr_orig.partial_occupation_flag) {
    ostringstream aus;
    aus << "0000 MESSAGE    \"partial_occupation_site_tol\" is " << xstr_orig.partial_occupation_site_tol << " " << Message(__AFLOW_FILE__,aflags); //CO20180409
    if(xstr_orig.partial_occupation_site_tol>1E-6) { //CO20180409
      epsilon = xstr_orig.partial_occupation_site_tol; //CO20180409
    }
    else {
      aus << "0000 MESSAGE    \"partial_occupation_site_tol\" is not set " << Message(__AFLOW_FILE__,aflags); //CO20180409
      aus << "0000 MESSAGE    Default value (1E-2) of \"partial_occupation_site_tol\" is set " << Message(__AFLOW_FILE__,aflags); //CO20180409
    }
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
    double partial_value_tmp;
    for (uint i=0; i<xstr_orig.atoms.size();i++) {
      partial_value_tmp = xstr_orig.atoms.at(i).partial_occupation_value;
      str_num_data tmp;
      tmp = OptimizePoccValue(partial_value_tmp, epsilon);
      xstr_orig.atoms.at(i).partial_occupation_value = tmp.decimal;
      UpdateXstr_comp_each_type(xstr_orig);
    }
  }
}

// ***************************************************************************
// Function pocc::HNFCELL
// ***************************************************************************
void OptimizeXstr(xstructure &a, ofstream &FileMESSAGE, _aflags &aflags) {
  //This function optimize partial values of xstructure
  ostringstream aus;
  aus <<"0000 MESSAGE    Optimizing POSCAR " << Message(__AFLOW_FILE__,aflags);
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
  UpdateXstr(a, FileMESSAGE, aflags);
  aus <<"0000 MESSAGE    Printing optimized POSCAR " << Message(__AFLOW_FILE__,aflags);
  aus << AFLOWIN_SEPARATION_LINE << endl;
  aus << a;
  aus << AFLOWIN_SEPARATION_LINE << endl;
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
}

// ***************************************************************************
// void string_replaceAll(std::string& str, const std::string& from, const std::string& to)
// ***************************************************************************
void string_replaceAll(std::string& str, const std::string& from, const std::string& to) 
{
  size_t start_pos = 0;
  while((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
  }
}

// ***************************************************************************
//  bool is_number(const std::string& s)
// ***************************************************************************
bool is_number(const std::string& s) {
  std::string::const_iterator it = s.begin();
  while (it != s.end() && std::isdigit(*it)) {
    ++it;
  }
  return !s.empty() && it == s.end();
}

// ***************************************************************************
//  double CalculateDistanceXstructure(xstructure& xstr, int i, int j)
// ***************************************************************************
double CalculateDistanceXstructure(xstructure& xstr, int i, int j) {
  double d1=xstr.atoms.at(i).cpos(1)-xstr.atoms.at(j).cpos(1);
  double d2=xstr.atoms.at(i).cpos(2)-xstr.atoms.at(j).cpos(2);
  double d3=xstr.atoms.at(i).cpos(3)-xstr.atoms.at(j).cpos(3);
  double d=sqrt(d1*d1+d2*d2+d3*d3);
  return d;
}

// ***************************************************************************
// vector<vector<int> > CalculateNumberMultiOccupied(xstructure& xstr)
// ***************************************************************************
vector<vector<int> > NormalisedNumberXstructure(xstructure& xstr) {
  //Format Xstructure, if two different atoms occupy the same site, then they will
  //be put in the same vector
  vector<vector<int> > a; 
  bool insert_flag;
  vector<int> temp;
  temp.push_back(0); //Put the first number of xstr
  a.push_back(temp);
  for (unsigned int i=1; i<xstr.atoms.size();i++) {
    int k1=i;
    insert_flag = false;
    for (unsigned int j=0; j<a.size();j++) {
      int k2=a.at(j).at(0);
      if(CheckOneSite(xstr,k1,k2)) {
        insert_flag = true;
        a.at(j).push_back(k1);
        break;
      }
    }
    if(! insert_flag) {
      temp.clear();
      temp.push_back(k1);
      a.push_back(temp);
    }
  }
  return a;
}

// ***************************************************************************
// bool CheckPartialOccupation(xstructure& xstr)
// ***************************************************************************
vector<vector<int> > CalculateLableXstructure(xstructure& xstr) {
  vector<vector<int> > a=NormalisedNumberXstructure(xstr);
  vector<vector<int> > NormalisedLable;
  vector<int> LabelTmp;
  for (unsigned int i=0; i<a.size();i++) {
    LabelTmp.clear();
    for (unsigned int j=0; j<a.at(i).size();j++) {
      int k=a.at(i).at(j);
      LabelTmp.push_back(xstr.atoms.at(k).type);
    }
    if(CheckVacancyOnOnesite(xstr,i)) {
      int NumVacancy=xstr.num_each_type.size();
      LabelTmp.push_back(NumVacancy);
    }
    //LabelTmp.push_back(-1); //Add "-1" to terminate
    NormalisedLable.push_back(LabelTmp);
  }
  return NormalisedLable;
}
// ***************************************************************************
// bool CheckPartialOccupation(xstructure& xstr)
// ***************************************************************************
bool CheckPartialOccupation(xstructure& xstr) {
  //if it is a partially occupied, then return true, else return false
  bool RunFLAG=false;
  if(xstr.partial_occupation_flag) {
    for (unsigned int i=0; i<xstr.atoms.size(); i++) {
      if(xstr.atoms.at(i).partial_occupation_flag) {
        RunFLAG=true;
      }
    }
  }
  return RunFLAG;
}

// ***************************************************************************
// bool CheckOneSite(xstructure& xstr, int i, int j)
// ***************************************************************************
bool CheckOneSite(xstructure& xstr, int i, int j) {
  //if xstr.atoms.at(i) and xstr.atoms.at(j) belong to the same site, return true, else return false
  double epsilon=1E-3;
  bool RunFLAG=false; 
  if(CalculateDistanceXstructure(xstr, i, j)<epsilon) {
    RunFLAG=true;
  }
  return RunFLAG;
}

// ***************************************************************************
// bool CheckMultiOccupied(xstructure& xstr)
// ***************************************************************************
bool CheckMultiOccupied(xstructure& xstr) {
  //if it is double occupied, then return true, else return false
  return(NormalisedNumberXstructure(xstr).size()<xstr.atoms.size()); 
}

// ***************************************************************************
// double CalculatePartialValueOfVacancy(xstructure& xstr, unsigned int k)
// ***************************************************************************
double CalculatePartialValueOfVacancy(xstructure& xstr, unsigned int k) {
  double SumTmp=0.0;
  vector<vector<int> > a=NormalisedNumberXstructure(xstr);
  if(k>=a.size()) {
    throw aurostd::xerror(__AFLOW_FILE__,"CalculatePartialValueOfVacancy():","k>=a.size()",_INDEX_MISMATCH_); //CO20200624
  }
  if(a.at(k).size()==1) {
    int n_xstr_orig=a.at(k).at(0);
    SumTmp=xstr.atoms.at(n_xstr_orig).partial_occupation_value;
  } else {
    SumTmp=0.0;
    for (unsigned int i=0; i<a.at(k).size();i++) {
      int j=a.at(k).at(i);
      SumTmp+=xstr.atoms.at(j).partial_occupation_value;
    }
  }
  return(1.0-SumTmp);
}
// ***************************************************************************
// bool CheckVacancyOnOnesite(xstructure& xstr, int k)
// ***************************************************************************
bool CheckVacancyOnOnesite(xstructure& xstr, unsigned int k) {
  //Note that this k is the number is the normalised xstructure, not orignal xstructure
  bool RunFLAG=false;
  vector<vector<int> > a=NormalisedNumberXstructure(xstr);
  if(k>=a.size()) {
    throw aurostd::xerror(__AFLOW_FILE__,"CheckVacancyOnOnesite():","k>=a.size()",_INDEX_MISMATCH_); //CO20200624
  }

  if(a.at(k).size()==1) {
    int n_xstr_orig=a.at(k).at(0);
    if(xstr.atoms.at(n_xstr_orig).partial_occupation_value<1.0) {
      RunFLAG=true;
    }
  }  else {
    double SumTmp=0.0;
    for (unsigned int i=0; i<a.at(k).size();i++) {
      int j=a.at(k).at(i);
      SumTmp+=xstr.atoms.at(j).partial_occupation_value;
    }
    if(SumTmp<1.0) {
      RunFLAG=true;
    }
  }
  return RunFLAG;
}
// ***************************************************************************
//bool CleanVacancy(xstructure& xstr) 
// ***************************************************************************
xstructure CleanVacancy(xstructure& xstr) 
{
  //This routine assumes that the vacancy is at the last type of the atoms
  xstructure CleanedXstr;
  int Nions=xstr.atoms.size();
  int Ntype=xstr.num_each_type.size();
  int NumLastType=xstr.num_each_type.at(Ntype-1);
  for (int i=(Nions-1); i>=(Nions-NumLastType);i--) {
    xstr.RemoveAtom(i);
  } 
  CleanedXstr=xstructure(xstr);
  return CleanedXstr;
}

// ***************************************************************************
//bool CheckVacancy(xstructure& xstr) 
// ***************************************************************************
bool CheckVacancy(xstructure& xstr) 
{
  //if xstr has a vacancy, then return true, else return false
  bool RunFLAG=false;
  if(CheckPartialOccupation(xstr)) {
    if(!CheckMultiOccupied(xstr)) {
      RunFLAG=true;
    }
    else {
      vector<vector<int> > a=NormalisedNumberXstructure(xstr);
      for (uint i=0;i<a.size();i++) {
        for (uint j=0; j<a.at(i).size();j++) {
          cout << a.at(i).at(j) << " ";
        }
        cout << endl;
      }

      vector<double> SumOccupation;
      double SumTmp;
      for (unsigned int i=0; i<a.size();i++) {
        SumTmp=0.0;
        for (unsigned int j=0; j<a.at(i).size(); j++) {
          int k=a.at(i).at(j);
          SumTmp+=xstr.atoms.at(k).partial_occupation_value;
        }
        SumOccupation.push_back(SumTmp);
      }
      for(unsigned int i=0; i<SumOccupation.size();i++) { 
        if((1-SumOccupation.at(i))>1E-4) {
          RunFLAG=true;
        }
        else if(SumOccupation.at(i)>1) {
          throw aurostd::xerror(__AFLOW_FILE__,"CheckVacancy():","total occupation on site exceeds 1.0",_INDEX_MISMATCH_); //CO20200624
        }
      }
    }
  }
  return RunFLAG;
}

// ***************************************************************************
// vector<vector<int> > CalculatecRange(xstructure& xstr)
// ***************************************************************************
vector<vector<int> > CalculatecRange(xstructure& xstr) {
  int Ntype;
  if(CheckVacancy(xstr)) {
    Ntype=xstr.num_each_type.size()+1;
  }
  else {
    Ntype=xstr.num_each_type.size();
  }
  int nDFull=NormalisedNumberXstructure(xstr).size(); 

  vector<vector<int> > XstrcRange;
  vector<int> cRangeTmp;
  if(CheckPartialOccupation(xstr)) {
    double ctmp=0.0;
    vector<double>  ConUniqueSite;
    for (int i=0; i<(Ntype-1);i++) {
      ctmp+=xstr.comp_each_type.at(i);
      ConUniqueSite.push_back(xstr.comp_each_type.at(i)/nDFull);
    }
    ConUniqueSite.push_back((nDFull-ctmp)/nDFull);

    vector<int> FracTmp; //Midlle Step
    vector<int> ConFrac; //Get the fraction for each unique site
    vector<vector<int> > AllConFrac; //Combine all the fraction into vector
    for (unsigned int i=0; i<ConUniqueSite.size();i++) {
      FracTmp=Decimal2Fraction(ConUniqueSite.at(i));
      ConFrac=GetFraction(FracTmp);
      AllConFrac.push_back(ConFrac);
    }

    vector<int> denominatortmp;
    for (uint i=0; i<AllConFrac.size();i++) {
      denominatortmp.push_back(AllConFrac.at(i).at(1));
    }
    int Xstrlcm=CalculateLcmofVector(denominatortmp);

    int tmp;
    for (uint i=0; i<AllConFrac.size();i++) {
      cRangeTmp.clear();
      tmp=Xstrlcm/AllConFrac.at(i).at(1)*AllConFrac.at(i).at(0);
      cRangeTmp.push_back(tmp);
      cRangeTmp.push_back(tmp);
      cRangeTmp.push_back(Xstrlcm);
      XstrcRange.push_back(cRangeTmp);
    }
    //The second algorithm
    //int sum_num=0;
    //for (uint i=0; i<XstrcRange.size();i++) {
    //sum_num+=XstrcRange.at(i).at(0);
    //}
    //cRangeTmp.clear();
    //cRangeTmp.push_back(XstrcRange.at(0).at(2)-sum_num);
    //cRangeTmp.push_back(XstrcRange.at(0).at(2)-sum_num);
    //cRangeTmp.push_back(XstrcRange.at(0).at(2));
    //XstrcRange.push_back(cRangeTmp);

    int sumtmp=0;
    for (uint i=0; i<XstrcRange.size();i++) {
      sumtmp+=XstrcRange.at(i).at(0); 
    }

    if(sumtmp!=XstrcRange.at(0).at(2)) { //sumtmp2 does not contain vacancy; sumtmp contains vacancy
      int sumtmp2=0;
      int k=XstrcRange.size();
      for (int i=0; i<(k-1);i++) {
        sumtmp2+=XstrcRange.at(i).at(0);
      }
      XstrcRange.at(k-1).at(0)=XstrcRange.at(0).at(2)-sumtmp2;
      XstrcRange.at(k-1).at(1)=XstrcRange.at(0).at(2)-sumtmp2;
    }
  }
  else {
    for (unsigned int i=0; i<xstr.num_each_type.size();i++) {
      cRangeTmp.clear();
      int k=xstr.num_each_type.at(i);
      int Nions=xstr.atoms.size();
      cRangeTmp.push_back(k);
      cRangeTmp.push_back(k);
      cRangeTmp.push_back(Nions);
      XstrcRange.push_back(cRangeTmp);
    }
  }
  return XstrcRange;
}

// ***************************************************************************
// int CalculateLcmofVector(vector<int> Num)
// ***************************************************************************
int CalculateLcmofVector(vector<int> Num) {
  vector<int> vlcm;
  int lcmtmp=0;
  if(Num.size()<2) {
    throw aurostd::xerror(__AFLOW_FILE__,"CalculateLcmofVector():","cannot find LCM: fewer than 2 inputs entered",_INPUT_NUMBER_); //CO20200624
  }
  lcmtmp=lcm(Num.at(0),Num.at(1));
  vlcm.push_back(lcmtmp);
  for (unsigned int i=1; i<(Num.size()-1);i++) {
    lcmtmp=lcm(Num.at(i),Num.at(i+1));
    vlcm.push_back(lcmtmp);
  }
  lcmtmp = max(vlcm);
  return lcmtmp;
}

// ***************************************************************************
// int lcm (int a, int b)
// ***************************************************************************
int lcm (int a, int b) {
  return a/gcd(a,b)*b;
}
// ***************************************************************************
// int gcd(int a, int b)
// ***************************************************************************
int gcd(int a, int b) {
  if(a%b==0) {
    return b;
  }
  else {
    return gcd(b, a%b);
  }
}
// ***************************************************************************
// vector<int> Decimal2Fraction(double a)
// ***************************************************************************
vector<int> Decimal2Fraction(double a) {
  //This function converts decimal number to fraction number (not real),
  //You need to call GetFraction to get the real fraction
  vector<int> IntList;
  vector<int> FracTmp;
  int b;
  double DecNum=a; 
  //Debug
  double FracTmpDouble;
  double NumDiff=1;
  if(abs(a-int(a))<Tol_DecimalPart) {
    IntList.push_back(0);
    IntList.push_back(int(a));
  }
  else {
    while(abs(a)>Tol_DecimalPart && NumDiff>Tol_FracDeciDiff) {
      b=int(1.0/a);
      if(b>MaxDenominator) {
        break;
      }
      IntList.push_back(b);
      a=1.0/a-b;
      FracTmp=GetFraction(IntList);
      FracTmpDouble=double (FracTmp.at(0))/FracTmp.at(1);
      NumDiff=abs(DecNum-FracTmpDouble);
    }
  }
  return IntList;
}

// ***************************************************************************
// vector<int> GetFraction(vector<int> IntList)
// ***************************************************************************
vector<int> GetFraction(vector<int> IntList) {
  //This function read the vector generated by Decimal2Fraction(), and generate real fraction
  vector<int> Fraction;
  int LenghIntList=IntList.size();
  int tmp;
  int fac=1;
  int Denom;
  Denom=IntList.at(LenghIntList-1);
  for (int i=LenghIntList-2; i>=0; i--) {
    tmp=Denom;
    Denom=fac+Denom*IntList.at(i);
    fac=tmp;
  }
  Fraction.push_back(fac);
  Fraction.push_back(Denom);
  return Fraction;
}

// ***************************************************************************
// vector<int> double2fraction(double a)
// ***************************************************************************
vector<int> double2fraction(double a) {
  vector<int> IntList = Decimal2Fraction(a);
  vector<int> Fraction = GetFraction(IntList);
  return Fraction;
}

// ***************************************************************************
// bool sortdenom (str_num_data str_i, str_num_data str_j)
// ***************************************************************************
bool sortdenom (str_num_data str_i, str_num_data str_j) {
  return (str_i.denom<str_j.denom); //Ascending sort according to denom 
}


// ***************************************************************************
// str_num_data double2str_num_data(double a)
// ***************************************************************************
str_num_data double2str_num_data(double a) { 
  str_num_data result;
  vector<int> fraction = double2fraction(a);

  result.decimal = a;
  result.fac = fraction.at(0);
  result.denom = fraction.at(1);
  return result;
}

// ***************************************************************************
// pocc::POSCAR2GULP(istream& input)
// ***************************************************************************
namespace pocc {
  void POSCAR2GULP(istream& input) {
    xstructure xstr; 
    xstr.clear(); //DX20191220 - uppercase to lowercase clear
    xstr=xstructure(input, IOVASP_POSCAR);
    xstr.ReScale(1.0);
    string GulpInput;
    GulpInput.clear();

    vector<string> AtomSpecies;
    aurostd::string2tokens(xstr.SpeciesString(), AtomSpecies, " ");
    if(!xstr.atoms.at(0).name_is_given) {
      cerr << "Warnning!!! Atom Species are not assigned! They are defaulted to H." << endl;
      AtomSpecies.clear();
      for (uint i=0; i<xstr.num_each_type.size();i++) {
        AtomSpecies.push_back("H");
      }
    }

    GulpInput=pocc::POSCAR2GulpInput(xstr, AtomSpecies);
    cout << GulpInput ;
  }
}

// ***************************************************************************
// pocc::POSCAR2GulpInput(xstructure & xstr)
// ***************************************************************************
namespace pocc {
  string POSCAR2GulpInput(xstructure& xstr) {
    vector<string> AtomSpecies;
    aurostd::string2tokens(xstr.SpeciesString(), AtomSpecies, " ");
    if(!xstr.atoms.at(0).name_is_given) {
      cerr << "Warnning!!! Atom Species are not assigned! They are defaulted to H." << endl;
      AtomSpecies.clear();
      for (uint i=0; i<xstr.num_each_type.size();i++) {
        AtomSpecies.push_back("H");
      }
    }

    string GulpInput = pocc::POSCAR2GulpInput(xstr, AtomSpecies);
    return GulpInput;
  }
}

// ***************************************************************************
// pocc::POSCAR2GulpInput(xstructure & xstr)
// ***************************************************************************
namespace pocc {
  string POSCAR2GulpInput(xstructure& xstr, vector<string> AtomSpecies) {
    if(xstr.num_each_type.size()!=AtomSpecies.size()) {
      throw aurostd::xerror(__AFLOW_FILE__,"POSCAR2GulpInput():","xstr.num_each_type.size()!=AtomSpecies.size()",_INDEX_MISMATCH_); //CO20200624
    }
    //Assign names to xstructure 
    int iatom=0;
    for(uint itype=0; itype<xstr.num_each_type.size();itype++) {
      for (int j=0; j<xstr.num_each_type.at(itype);j++) {
        xstr.atoms.at(iatom).name=AtomSpecies.at(itype);
        xstr.atoms.at(iatom).name_is_given=TRUE;
        iatom++;
      }
    }

    stringstream  oss;
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(10);
    oss << "single dist comp conp" << endl;
    oss << "cell" << endl;
    oss << xstr.a << '\t' << xstr.b << '\t'  << xstr.c << '\t' << xstr.alpha << '\t' << xstr.beta << '\t' << xstr.gamma << endl;
    oss << "fractional" << endl;
    oss.setf(ios_base::fixed,ios_base::floatfield);
    for (uint i=0; i < xstr.atoms.size(); i++) {
      oss <<xstr.atoms.at(i).name << '\t' << xstr.atoms.at(i).fpos(1) << '\t' << xstr.atoms.at(i).fpos(2) << '\t' << xstr.atoms.at(i).fpos(3) << endl;
    }
    oss << "#" << endl;
    oss << "species    " << xstr.num_each_type.size() << endl;

    string SpeciesTmp;
    for (uint i=0; i<AtomSpecies.size();i++) {
      SpeciesTmp.clear();
      SpeciesTmp=AtomSpecies.at(i);
      oss << pocc::ReturnAtomSpecies(SpeciesTmp) << endl;
    }
    oss << "#" << endl;
    oss << "uff1 kcal" << endl;
    for (uint i=0; i<AtomSpecies.size();i++) {
      SpeciesTmp.clear();
      SpeciesTmp=AtomSpecies.at(i);
      oss << pocc::ReturnAtomSpeciesPotential(SpeciesTmp) << endl;
    }
    return oss.str();
  }
}

// ***************************************************************************
// double pocc::CalculateEnergyUsingGulp(xstructure & xstr)
// ***************************************************************************
namespace pocc {
  double CalculateEnergyUsingGulp(xstructure & xstr, vector<string> AtomSpecies) {
    double TotalEnergy=0.0;
    string XstrGulpInput;
    string GulpOut;
    string GulpTmp=aurostd::TmpFileCreate("GulpTmp");

    if(AtomSpecies.at(0).compare("unknown")==0) {
      throw aurostd::xerror(__AFLOW_FILE__,"pocc::CalculateEnergyUsingGulp():","atom species are not assigned",_INPUT_MISSING_); //CO20200624
    }

    XstrGulpInput=pocc::POSCAR2GulpInput(xstr, AtomSpecies);
    aurostd::string2file(XstrGulpInput,GulpTmp);

    string RunGulp="cat "+ GulpTmp +" | gulp";
    GulpOut=aurostd::execute2string(RunGulp);
    aurostd::execute("rm " + GulpTmp);

    stringstream ss_gulp_out;
    ss_gulp_out << GulpOut;

    string anchor_word_Total_energy="Total lattice energy";
    string line;
    vector<string> vgulp_output_line, vgulp_output;
    while(getline(ss_gulp_out, line)) {
      if(line.find(anchor_word_Total_energy) !=string::npos) {
        vgulp_output.push_back(line);
      }
    }
    aurostd::string2tokens(vgulp_output.at(0), vgulp_output_line, " ");
    TotalEnergy=aurostd::string2utype<double>(vgulp_output_line.at(4));
    return TotalEnergy;
  }
} // namespace pocc

// ***************************************************************************
// double pocc::CalculateEnergyUsingGulp(xstructure & xstr)
// ***************************************************************************
namespace pocc {
  double CalculateEnergyUsingGulp(xstructure& xstr) {
    double TotalEnergy=0.0;
    string XstrGulpInput;
    string GulpOut;
    string GulpTmp=aurostd::TmpFileCreate("GulpTmp");
    vector<string> AtomSpecies;
    AtomSpecies.clear();
    aurostd::string2tokens(xstr.SpeciesString(), AtomSpecies, " ");

    XstrGulpInput=pocc::POSCAR2GulpInput(xstr, AtomSpecies);
    aurostd::string2file(XstrGulpInput,GulpTmp);

    string RunGulp="cat "+ GulpTmp +" | gulp";
    GulpOut=aurostd::execute2string(RunGulp);
    aurostd::execute("rm " + GulpTmp);

    stringstream ss_gulp_out;
    ss_gulp_out << GulpOut;

    string anchor_word_Total_energy="Total lattice energy";
    string line;
    vector<string> vgulp_output_line, vgulp_output;
    while(getline(ss_gulp_out, line)) {
      if(line.find(anchor_word_Total_energy) !=string::npos) {
        vgulp_output.push_back(line);
      }
    }
    aurostd::string2tokens(vgulp_output.at(0), vgulp_output_line, " ");
    TotalEnergy=aurostd::string2utype<double>(vgulp_output_line.at(4));
    return TotalEnergy;
  }
} // namespace pocc


// ***************************************************************************
// double pocc::CalculateEnergyUsingGulp(xstructure & xstr)
// ***************************************************************************
namespace pocc {
  vector<double> CalculateEnergyUsingGulp(vector<xstructure>& vxstr) {
    vector<double> Etot_vxstr; 
    for (uint i=0; i<vxstr.size();i++) {
      double dbtmp=0.0;
      xstructure xstr_tmp;
      xstr_tmp = vxstr.at(i);
      dbtmp=pocc::CalculateEnergyUsingGulp(xstr_tmp);
      Etot_vxstr.push_back(dbtmp);
    }
    return  Etot_vxstr;
  }
} // namespace pocc

void PrintGulpEnergies(vector<xstructure>& vxstr) {
  vector<double> etot_vxstr = pocc::CalculateEnergyUsingGulp(vxstr);
  for (uint i=0;i<etot_vxstr.size();i++) {
    ostringstream oss;
    oss.setf(std::ios::fixed, std::ios::floatfield);
    oss.precision(12);
    oss << "No. " << i+1 << " " << setw(15) << etot_vxstr.at(i) << endl;
    cout << oss.str();
  }
}

// ***************************************************************************
// bool sortenergy (strno_energy str_i, strno_energy str_j)
// ***************************************************************************
bool sortenergy (strno_energy str_i, strno_energy str_j) {
  return (str_i.energy<str_j.energy);
}

// ***************************************************************************
// vector<int> pocc::SortGroupXstr(vector<xstructure> groupxstr, vector<string> AtomSpecies)
// ***************************************************************************
namespace pocc {
  vector<xstructure> SortGroupXstr(vector<xstructure> groupxstr, vector<string> AtomSpecies) {
    vector<strno_energy> vec_group_strno_energy;
    vector<int> number;

    cerr.precision(8);
    cerr.setf(ios_base::fixed,ios_base::floatfield);

    double TotalEnergyTmp;
    for (uint i=0; i<groupxstr.size();i++) {
      strno_energy strno_energy_tmp;
      TotalEnergyTmp=pocc::CalculateEnergyUsingGulp(groupxstr.at(i),AtomSpecies);
      cerr << "Estimating Total Energy: Str # " <<  setfill('0') << setw(6) << i+1 << " : " << TotalEnergyTmp << " eV" << endl;
      strno_energy_tmp.energy=TotalEnergyTmp;
      strno_energy_tmp.number=i;
      vec_group_strno_energy.push_back(strno_energy_tmp);
    }
    sort (vec_group_strno_energy.begin(), vec_group_strno_energy.end(), sortenergy);

    cerr << "Sorting energy: " << endl;
    for (uint i=0; i<vec_group_strno_energy.size();i++) {
      cerr << "No. " << setfill('0') << setw(6) << i+1 << "   Str # " << setfill('0') << setw(6) << vec_group_strno_energy.at(i).number << " : " << vec_group_strno_energy.at(i).energy << " eV" << endl; 
      number.push_back(vec_group_strno_energy.at(i).number);
    }

    vector<xstructure> SortedGroupXstr;

    for (uint i=0; i<number.size();i++) {
      int k=number.at(i);
      SortedGroupXstr.push_back(groupxstr.at(k));
    }
    return SortedGroupXstr;
  }
} // namespace pocc

// ***************************************************************************
// vector<int> pocc::SortGroupXstrUFFEnergy(vector<xstructure> groupxstr, vector<string> AtomSpecies)
// ***************************************************************************
namespace pocc {
  vector<xstructure> SortGroupXstrUFFEnergy(vector<xstructure> groupxstr) {
    vector<strno_energy> vec_group_strno_energy;
    vector<int> number;

    cout.precision(8);
    cout.setf(ios_base::fixed,ios_base::floatfield);

    double TotalEnergyTmp;
    for (uint i=0; i<groupxstr.size();i++) {
      strno_energy strno_energy_tmp;
      TotalEnergyTmp=pocc::CalculateUFFEnergy(groupxstr.at(i));
      cout << "Estimating Total Energy: Str # " <<  setfill('0') << setw(6) << i+1 << " : " << TotalEnergyTmp << " eV" << endl;
      strno_energy_tmp.energy=TotalEnergyTmp;
      strno_energy_tmp.number=i;
      vec_group_strno_energy.push_back(strno_energy_tmp);
    }
    sort (vec_group_strno_energy.begin(), vec_group_strno_energy.end(), sortenergy);

    cout << "Sorting energy: " << endl;
    for (uint i=0; i<vec_group_strno_energy.size();i++) {
      cout << "No. " << setfill('0') << setw(6) << i+1 << "   Str # " << setfill('0') << setw(6) << vec_group_strno_energy.at(i).number << " : " << vec_group_strno_energy.at(i).energy << " eV" << endl; 
      number.push_back(vec_group_strno_energy.at(i).number);
    }

    vector<xstructure> SortedGroupXstr;

    for (uint i=0; i<number.size();i++) {
      int k=number.at(i);
      SortedGroupXstr.push_back(groupxstr.at(k));
    }
    return SortedGroupXstr;
  }
} // namespace pocc

// ***************************************************************************
// void combine(vector<int> &range, vector<int> &cur, vector<vector<int> > &final_result, int start, int depth)
// ***************************************************************************
// subroutine of combine
void combine(vector<int> &range, vector<int> &cur, vector<vector<int> > &final_result, int start, int depth) {
  if(depth == (int) cur.size()) {
    vector<int> tmp;
    for (vector<int>::iterator i = cur.begin(); i != cur.end(); ++i) {
      tmp.push_back(*i);
    }   
    final_result.push_back(tmp);
    return;
  }   

  for (int i = start; i != (int) range.size(); ++i) {
    cur[depth] = range[i];
    combine(range, cur, final_result, i + 1, depth + 1); 
  }   
}

// ***************************************************************************
// void combine(vector<int> &range, vector<vector<int> > &final_result, int n)
// ***************************************************************************
void combine(vector<int> &range, vector<vector<int> > &final_result, int n) {
  //Combine function, choose all the combinations of n numbers from range
  if(n > (int) range.size()) {
    throw aurostd::xerror(__AFLOW_FILE__,"combine():","n > (int) range.size()",_INPUT_ILLEGAL_); //CO20200624
  }   
  vector<int> result(n);
  combine(range, result, final_result, 0, 0); 
}

// ***************************************************************************
// vector<xmatrix<double> > CalculateHNF(xstructure a, int n)
// ***************************************************************************
vector<xmatrix<double> > CalculateHNF(xstructure str, int n) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  //xstructure str=a; //CO, already copies
  str.BringInCell();

  // GET SYMMETRY
  str.FixLattices();
  xstructure non_pocc_xstr=NormalizeXstructure(str);
  non_pocc_xstr.CalculateSymmetryPointGroup();
  //str.CalculateSymmetryPointGroup();  //CO, cannot calculate symmetry of POCC structure as atoms are on top of each other
  // str.CalculateSymmetryFactorGroup();

  // PHYSICAL REVIEW B 77, 224115 2008
  // Algorithm for generating derivative structures
  // Gus L. W. Hart and Rodney W. Forcade
  xmatrix<double> HNF(3,3);
  vector<xmatrix<double> > vHNF;
  xmatrix<double> A(3,3),B(3,3),Bi(3,3),Bj(3,3),R(3,3),H(3,3);A=trasp(str.lattice);
  vector<xmatrix<double> > vB;    // contains lattices of the potential xstructures
  vector<xmatrix<double> > sHNF;  // contains only the supercells HNF that are uniques
  long long int i=0;
  bool found=FALSE;
  for(int a=1;a<=n;a++) {
    for(int c=1;c<=n/a;c++) {
      for(int f=1;f<=n/a/c;f++) {
        if(a*c*f==n) {
          for(int b=0;b<c;b++) {
            for(int d=0;d<f;d++) {
              for(int e=0;e<f;e++) {
                i++;
                HNF(1,1)=a;HNF(1,2)=0;HNF(1,3)=0;
                HNF(2,1)=b;HNF(2,2)=c;HNF(2,3)=0;
                HNF(3,1)=d;HNF(3,2)=e;HNF(3,3)=f;
                vHNF.push_back(HNF);
              }
            }
          }
        }
      }
    }
  }

  for(i=0;i<(int) vHNF.size();i++) {
    if(LDEBUG) { //CO20180409
      cerr << "HNF MAT " << endl;
      cerr << vHNF.at(i) << endl;
    }
    Bi=A*vHNF.at(i);
    found=FALSE;
    if(LDEBUG) { //CO20180409
      cerr << "SUPERLATTICE " << endl;
      cerr << Bi << endl;
    }
    for(uint istr=0;istr<vB.size()&&found==FALSE;istr++)  {                              // cycle through the unique vBs
      Bj=vB.at(istr);
      if(LDEBUG) { //CO20180409
        cerr << "SUPERLATTICE OLD " << endl;
        cerr << Bj << endl;
      }
      for(uint pgroup=0;pgroup<non_pocc_xstr.pgroup.size()&&found==FALSE;pgroup++) {                  // cycle through the pgroup of str
        R=trasp(non_pocc_xstr.pgroup.at(pgroup).Uc); 
        if(LDEBUG) { //CO20180409
          cerr << "POINT GROUP " << endl;
          cerr << R << endl;
        }
        H=aurostd::inverse(Bj)*aurostd::inverse(R)*Bi;
        if(LDEBUG) { //CO20180409
          cerr << "H" << endl;
          cerr << H << endl;
        }
        if(aurostd::isinteger(H)) {found=TRUE;}
      }
    }
    if(found==FALSE) { // not found, then plug
      if(LDEBUG) { //CO20180409
        cerr << "----------------------" << endl;
        cerr << "FOUND " << endl;
        cerr << vHNF.at(i) << endl;
        cerr << "----------------------" << endl;
      }
      vB.push_back(Bi);
      sHNF.push_back(vHNF.at(i));
    }
  }
  //debug 
  return sHNF;
}

// ***************************************************************************
// xstructure XstrSubstitute(xstructure &a, int n, string b)
// ***************************************************************************
xstructure XstrSubstitute(xstructure &a, int n, string b) {
  xstructure xstrb=a;

  _atom atom_replaced=xstrb.atoms.at(n); //changed the name of the substituted atom
  atom_replaced.name=b;

  xstrb.RemoveAtom(n);    //Remove this atom
  xstrb.AddAtom_POCC(atom_replaced); //Add the substituted atom, note the added atom always appear in the last 
  return xstrb;
}

// ***************************************************************************
// bool sort_function_struct_num_name (atom_number_name i, atom_number_name j)
// ***************************************************************************
bool sort_function_struct_num_name (atom_number_name i, atom_number_name j) {
  return (i.number>j.number);
}

// ***************************************************************************
// vector<atom_number_name> SortMax2Min(vector<atom_number_name>& a)
// ***************************************************************************
vector<atom_number_name> SortMax2Min(vector<atom_number_name>& a) {
  sort(a.begin(), a.end(), sort_function_struct_num_name);
  return a;
}

// ***************************************************************************
// bool myfunction (int i,int j) 
// ***************************************************************************
bool myfunction_int (int i,int j) 
{return (i>j);}

// ***************************************************************************
// vector<int> SortMax2Min(vector<int> a)
// ***************************************************************************
vector<int> SortMax2Min(vector<int>& a) {
  //sort from max to min
  sort(a.begin(), a.end(), myfunction_int);
  return a;
}

// ***************************************************************************
// bool myfunction_double (double i,double j) 
// ***************************************************************************
bool myfunction_double (double i,double j) {
  return (i>j);
}

// ***************************************************************************
// vector<double> SortMax2Min(vector<double>& a)
// ***************************************************************************
vector<double> SortMax2Min(vector<double>& a) {
  //sort from max to min
  sort(a.begin(), a.end(), myfunction_double);
  return a;
}

// ***************************************************************************
// xstructure XstrSubstitute(xstructure &a, vector<int> vec_n, string b)
// ***************************************************************************
xstructure XstrSubstitute(xstructure &a, vector<int> vec_n, string b) {
  xstructure xstr_tmp=a;
  vector<int> vec_n_sorted = SortMax2Min(vec_n);

  //Firstly, we add all the atoms we need
  for (uint i=0; i<vec_n_sorted.size();i++) {
    int k= vec_n_sorted.at(i);
    _atom atom_replaced = xstr_tmp.atoms.at(k);
    atom_replaced.name = b;
    xstr_tmp.AddAtom_POCC(atom_replaced);
  }
  //Now we remove all the atoms that we do not need
  for (uint i=0; i<vec_n_sorted.size();i++) {
    int k= vec_n_sorted.at(i);
    xstr_tmp.RemoveAtom(k);
  }
  int ntype = xstr_tmp.num_each_type.size();
  xstr_tmp.SpeciesSwap(0,ntype-1);   //Force to re-assign
  xstr_tmp.SpeciesPutAlphabetic();
  xstr_tmp.MakeBasis();
  return xstr_tmp;
}

// ***************************************************************************
// xstructure XstrSubstitute(xstructure &a, vector<int> vec_n, vector<string> b)
// ***************************************************************************
xstructure XstrSubstitute(xstructure &a, vector<int> vec_n, vector<string> b) {
  if(vec_n.size()!=b.size()) { 
    throw aurostd::xerror(__AFLOW_FILE__,"XstrSubstitute():","vec_n.size()!=b.size()",_INDEX_MISMATCH_); //CO20200624
  }
  vector<atom_number_name> vec_num_name;
  atom_number_name tmp;
  for (uint i=0; i<vec_n.size();i++) {
    tmp.number = vec_n.at(i);
    tmp.name = b.at(i);
    vec_num_name.push_back(tmp);
  }
  vector<atom_number_name> vec_num_name_new = SortMax2Min(vec_num_name);
  vector<int> vec_num_sorted;
  vector<string> vec_str_sorted;
  for (uint i=0; i<vec_num_name_new.size();i++) {
    vec_num_sorted.push_back(vec_num_name_new.at(i).number);
    vec_str_sorted.push_back(vec_num_name_new.at(i).name);
  }

  xstructure xstr_tmp=a;
  //Firstly, we add all the atoms we need
  for (uint i=0; i<vec_num_sorted.size();i++) {
    int k= vec_num_sorted.at(i);
    _atom atom_replaced = xstr_tmp.atoms.at(k);
    atom_replaced.name = vec_str_sorted.at(i);
    xstr_tmp.AddAtom_POCC(atom_replaced);
  }
  //Now we remove all the atoms that we do not need
  for (uint i=0; i<vec_num_sorted.size();i++) {
    int k= vec_num_sorted.at(i);
    xstr_tmp.RemoveAtom(k);
  }
  xstr_tmp.SpeciesPutAlphabetic();
  int ntype = xstr_tmp.num_each_type.size();
  xstr_tmp.SpeciesSwap(0,ntype-1);   //Force to re-assign
  xstr_tmp.SpeciesPutAlphabetic();
  xstr_tmp.MakeBasis();

  //Clean Vacancies
  int natoms = xstr_tmp.atoms.size();
  if(xstr_tmp.atoms.at(natoms-1).name == name_vacancy) {
    return CleanVacancy(xstr_tmp);
  }
  return xstr_tmp;
}

// ***************************************************************************
// xstructure NormalizeXstructure(xstructure a)
// ***************************************************************************
xstructure NormalizeXstructure(xstructure a) {
  xstructure b=a;
  vector<vector<int> > NormalisedNumber=NormalisedNumberXstructure(a);

  vector<int> removed_number; 
  for(uint i=0; i<NormalisedNumber.size(); i++) {
    if(NormalisedNumber.at(i).size() >= 2) {
      for (uint j=1; j<NormalisedNumber.at(i).size();j++) { //Keep the first one
        int k=NormalisedNumber.at(i).at(j);
        removed_number.push_back(k);
      }
    }
  }
  vector<int> sorted_removed_number = SortMax2Min(removed_number);
  for (uint i=0; i<sorted_removed_number.size();i++) {
    int k=sorted_removed_number.at(i);
    b.RemoveAtom(k);
  }
  return b;
}

// ***************************************************************************
// vector<vector<int> > CalculateXstrNumberSupercell(xstructure a, int n)
// ***************************************************************************
vector<vector<int> > CalculateXstrNumberSupercell(xstructure a, int n) {
  //This function is used to produce the coresponding number in the supercell for the special atom in xstructure a
  vector<vector<int> > number_supercell;
  vector<int> number_tmp;

  for (int i=0; i< (int) a.atoms.size(); i++) {
    number_tmp.clear();
    for (int j=i*n; j<(i+1)*n; j++) {
      number_tmp.push_back(j);
    }
    number_supercell.push_back(number_tmp);
  }
  return number_supercell;
}

// ***************************************************************************
// int hnf_double2int(double a)
// ***************************************************************************
int hnf_double2int(double a) {   
  double adiff=a-int(a);
  if(adiff >=0 && adiff < 0.5 ) {
    return int(a);
  }
  else {
    return int(a)+1;
  }
}

// ***************************************************************************
// bool CheckDoubleOccupied(xstructure xstr_orig)
// ***************************************************************************
bool CheckDoubleOccupied(xstructure xstr_orig) {
  //DoubleOccupied means two atoms occupy one position or one atom and one vacancy occupy one position
  bool FLAG_double_occupied = false;
  vector<vector<int> > num_norm_xstr = NormalisedNumberXstructure(xstr_orig);
  vector<int> vec_one_site_num;
  for (uint i=0; i<num_norm_xstr.size();i++) {
    if(num_norm_xstr.at(i).size()==1 && CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_double_occupied = true;}
    if(num_norm_xstr.at(i).size()==2 && !CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_double_occupied = true;}
  }
  return FLAG_double_occupied;
}

// ***************************************************************************
// bool CheckTripleOccupied(xstructure xstr_orig)
// ***************************************************************************
bool CheckTripleOccupied(xstructure xstr_orig) {
  //TripleOccupied means three atoms occupy one position or two atom and one vacancy occupy one position
  bool FLAG_triple_occupied = false;
  vector<vector<int> > num_norm_xstr = NormalisedNumberXstructure(xstr_orig);
  vector<int> vec_one_site_num;
  for (uint i=0; i<num_norm_xstr.size();i++) {
    if(num_norm_xstr.at(i).size()==2 && CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_triple_occupied = true;}
    if(num_norm_xstr.at(i).size()==3 && !CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_triple_occupied = true;}
  }
  return FLAG_triple_occupied;
}

// ***************************************************************************
// bool CheckFourfoldOccupied(xstructure xstr_orig)
// ***************************************************************************
bool CheckFourfoldOccupied(xstructure xstr_orig) {
  //FourfoldOccupied means three atoms occupy one position or two atom and one vacancy occupy one position
  bool FLAG_fourfold_occupied = false;
  vector<vector<int> > num_norm_xstr = NormalisedNumberXstructure(xstr_orig);
  vector<int> vec_one_site_num;
  for (uint i=0; i<num_norm_xstr.size();i++) {
    if(num_norm_xstr.at(i).size()==3 && CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_fourfold_occupied = true;}
    if(num_norm_xstr.at(i).size()==4 && !CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_fourfold_occupied = true;}
  }
  return FLAG_fourfold_occupied;
}

// ***************************************************************************
// bool CheckFivefoldOccupied(xstructure xstr_orig)
// ***************************************************************************
bool CheckFivefoldOccupied(xstructure xstr_orig) {
  //FourfoldOccupied means three atoms occupy one position or two atom and one vacancy occupy one position
  bool FLAG_fivefold_occupied = false;
  vector<vector<int> > num_norm_xstr = NormalisedNumberXstructure(xstr_orig);
  vector<int> vec_one_site_num;
  for (uint i=0; i<num_norm_xstr.size();i++) {
    if(num_norm_xstr.at(i).size()==4 && CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_fivefold_occupied = true;}
    if(num_norm_xstr.at(i).size()==5 && !CheckVacancyOnOnesite(xstr_orig, i)) {FLAG_fivefold_occupied = true;}
  }
  return FLAG_fivefold_occupied;
}


// ***************************************************************************
// void CombineAll(const vector<vector<vector<int> > > &allVecs, size_t vecIndex, vector<vector<int> > &intSoFar, vector<vector<vector<int> > > &results)
// ***************************************************************************
void CombineAll(const vector<vector<vector<int> > > &allVecs, size_t vecIndex, vector<vector<int> > &intSoFar, vector<vector<vector<int> > > &results) {
  //function to combine all the vectors (choose one from each one)
  //Call method: CombineAll(d, 0, e, f); d initial vector; e temporary vector, f final result
  if(vecIndex >= allVecs.size())
  {   
    vector<vector<int> > str_tmp;
    for(uint i=0;i<intSoFar.size();i++) { 
      str_tmp.push_back(intSoFar.at(i));
    }   
    results.push_back(str_tmp);
    return;
  }   
  for (size_t i=0; i<allVecs.at(vecIndex).size(); i++) {
    vector<vector<int> > tmp=intSoFar;
    tmp.push_back(allVecs.at(vecIndex).at(i));
    CombineAll(allVecs, vecIndex+1, tmp, results);
  }   
}

// ***************************************************************************
// vector<string> CalculateSecondAtomicNameSupercell(xstructure xstr_orig, int n)
// ***************************************************************************
vector<string> CalculateSecondAtomicNameSupercell(xstructure xstr_orig, int n) {
  vector<string> AtomicNameSupercell;
  //    if(CheckDoubleOccupied(xstr_orig)) {  //[CO20200106 - close bracket for indenting]}
  xstructure xstr_norm = NormalizeXstructure(xstr_orig); //Format POSCAR
  vector<vector<int> > num_norm_xstr = NormalisedNumberXstructure(xstr_orig); //Format normalized POSCAR's number using 2D vector
  vector<vector<int> > vec_num_supercell=CalculateXstrNumberSupercell(xstr_norm,n); //Produce supercell's number using 2D vector

  for (uint i=0; i<xstr_norm.atoms.size();i++) {
    if(xstr_norm.atoms.at(i).partial_occupation_flag) {
      //if this is partial occupation
      if(num_norm_xstr.at(i).size()==1) {
        //This means this site only contains vacancy
        for (uint j=0;j<vec_num_supercell.at(i).size();j++) {
          AtomicNameSupercell.push_back(name_vacancy);
        }
      }
      else {
        //This means site site contains another element
        int int_tmp = num_norm_xstr.at(i).at(1); //Read the atomic number in xstr_orig
        string tmp_second_element_name = xstr_orig.atoms.at(int_tmp).name; //Read the name
        for (uint j=0;j<vec_num_supercell.at(i).size();j++) {
          AtomicNameSupercell.push_back(tmp_second_element_name);
        }
      }
    }
    else {
      for (uint j=0;j<vec_num_supercell.at(i).size();j++) {
        string tmp_name=xstr_norm.atoms.at(i).name;
        AtomicNameSupercell.push_back(tmp_name);
      }
    }
  }
  return AtomicNameSupercell;
}


// ***************************************************************************
// vector<vector<int> > GenerateSupercellAtomNumber(xstructure xstr_orig, int n)
// ***************************************************************************
vector<vector<int> > GenerateSupercellAtomNumber(xstructure xstr_orig, int n) {
  //This function generates numbers of all the atoms which will be replaced
  vector<vector<vector<int> > > vec_vec_vec_num_vacancy;
  vector<vector<vector<int> > > vec_two_num_vacancy;

  vector<vector<int> > num_norm_xstr = NormalisedNumberXstructure(xstr_orig); //format like 0 1 2; 3

  xstructure xstr_norm = NormalizeXstructure(xstr_orig);
  vector<vector<int> > atom_num_supercell = CalculateXstrNumberSupercell(xstr_norm, n); //format like 0 1 2 3; 4 5 6 7
  //    vector<vector<int> > atom_num_supercell = CalculateXstrNumberSupercell(xstr_orig, n); //format like 0 1 2 3; 4 5 6 7
  //    We just need to store the atoms' number using normalized xstructure

  //if(CheckDoubleOccupied(xstr_orig)) {  //[CO20200106 - close bracket for indenting]}
  for (uint i=0; i<xstr_norm.atoms.size();i++) {
    if(xstr_norm.atoms.at(i).partial_occupation_flag) {
      int num_vacancy= hnf_double2int((1-xstr_norm.atoms.at(i).partial_occupation_value)*n);
      vector<int> vec_num_part=atom_num_supercell.at(i);
      vector<vector<int> > num_com_vacancy;
      num_com_vacancy.clear();
      combine(vec_num_part, num_com_vacancy, num_vacancy); 
      //generate all the combinations and store them into 3D vectors num_com_vacancy
      vec_vec_vec_num_vacancy.push_back(num_com_vacancy);
      //push the 2D vector into 3D vector according to the atoms' number
    }
  }

  vector<vector<vector<int> > > final_num_result;
  vector<vector<int> >  tmp_vect;
  CombineAll(vec_vec_vec_num_vacancy, 0, tmp_vect, final_num_result);
  //choose only one combination (vector) from each atom's site, and recombine them as one entire combination (2D vector) for one structure
  //All the possible combinations are stored as 3D vector final_num_result
  vector<vector<int> > formated_final_num;

  //Extract the actural atom's number from the 3D vector, and store them into 2D vector formated_final_num
  vector<int> tmp; 
  for (uint i=0; i<final_num_result.size();i++) {
    tmp.clear();
    for (uint j=0; j<final_num_result.at(i).size();j++) {
      for (uint k=0; k<final_num_result.at(i).at(j).size();k++) {
        tmp.push_back(final_num_result.at(i).at(j).at(k));
      }
    }
    formated_final_num.push_back(tmp);
  }
  return formated_final_num;
}


// ***************************************************************************
// vector<xstructure> AssignPartialOccupation(xstructure xstr_supercell, xstructure xstr_orig, int n)
// ***************************************************************************
vector<xstructure> AssignPartialOccupation(xstructure xstr_supercell, xstructure xstr_orig, int n) {
  //This function use two of the most important function GenerateSupercellAtomNumber(xstr_orig, n) & CalculateSecondAtomicNameSupercell(xstr_orig, n) 

  vector<xstructure> xstr_final_list_supercell;

  //Produce atom number which needs to be replaced
  vector<vector<int> > tmp_atom_num_supercell = GenerateSupercellAtomNumber(xstr_orig, n);
  //DEBUG
  //Produce coressponding atomic names which will replace original atom
  vector<string> name_str=CalculateSecondAtomicNameSupercell(xstr_orig, n);

  for (uint i=0; i<tmp_atom_num_supercell.size();i++) {
    xstructure xstr_tmp=xstr_supercell;
    vector<int> atom_number = tmp_atom_num_supercell.at(i);
    vector<int> atom_number_sorted = SortMax2Min(tmp_atom_num_supercell.at(i));

    //Add all the needed atoms 
    for (uint i=0; i<atom_number_sorted.size();i++) {
      int k=atom_number_sorted.at(i);
      _atom atom_replaced = xstr_tmp.atoms.at(k);
      atom_replaced.name = name_str.at(k);
      xstr_tmp.AddAtom_POCC(atom_replaced);
    }

    //Remove all the additional atoms
    for (uint i=0; i<atom_number_sorted.size();i++) {
      int k=atom_number_sorted.at(i);
      xstr_tmp.RemoveAtom(k);
    }

    //Force to re-assign the structure 
    int ntype = xstr_tmp.num_each_type.size();
    xstr_tmp.SpeciesSwap(0,ntype-1);   //Force to re-assign
    xstr_tmp.SpeciesPutAlphabetic();
    xstr_tmp.MakeBasis();

    int LastAtom = xstr_tmp.atoms.size()-1;
    if(xstr_tmp.atoms.at(LastAtom).name == name_vacancy) { //this is safe 
      xstr_final_list_supercell.push_back(CleanVacancy(xstr_tmp));
    }
    else {
      xstr_final_list_supercell.push_back(xstr_tmp);
    }
  }
  return xstr_final_list_supercell;
}

// ***************************************************************************
// int CalculateSupercellSize(xstructure xstr_orig)
// ***************************************************************************
int CalculateSupercellSize(xstructure xstr_orig) {
  int n;
  if(!xstr_orig.partial_occupation_HNF) {
    cerr << "partial_occupation_HNF is not found!" << endl;
    cerr << "Calculating supercell size using occupation value ..." << endl; 
    vector<int> vec_denom;
    for (uint i=0; i<xstr_orig.atoms.size();i++) {
      double partial_value_tmp = xstr_orig.atoms.at(i).partial_occupation_value;
      str_num_data tmp;
      tmp = double2str_num_data(partial_value_tmp);
      vec_denom.push_back(tmp.denom);
    }
    n = CalculateLcmofVector(vec_denom);
  }
  else {
    n = xstr_orig.partial_occupation_HNF;
  }
  return n;
}

// ***************************************************************************
// xstructure AssignAlphabeticNameXstr(xstructure& xstr)
// ***************************************************************************
xstructure AssignAlphabeticNameXstr(xstructure& xstr) {
  vector<string> names;
  names.push_back("As"); names.push_back("Ba"); names.push_back("Cl"); names.push_back("Dy");
  names.push_back("Eu"); names.push_back("Fe"); names.push_back("Ge"); names.push_back("He");

  xstructure xstr_new = xstr;
  int iatom=0;
  for (uint itype=0; itype<xstr.num_each_type.size();itype++) {
    string species = string(names.at(xstr_new.atoms.at(iatom).type));
    xstr_new.species.at(itype)=species;
    for (int j=0;j<xstr.num_each_type.at(itype);j++) {
      xstr_new.atoms.at(iatom).name=species;
      xstr_new.atoms.at(iatom).CleanName();
      xstr_new.atoms.at(iatom).CleanSpin();
      xstr_new.atoms.at(iatom).name_is_given=TRUE;
      iatom++;
    }
  }
  return xstr_new;
}

// ***************************************************************************
// xstructure CleanNameXstr(xstructure& xstr)
// ***************************************************************************
xstructure CleanNameXstr(xstructure& xstr) {
  xstructure xstr_new = xstr;
  for (uint i=0; i<xstr_new.atoms.size();i++) {
    xstr_new.atoms.at(i).name_is_given=false;
  }
  return xstr_new;
}

// ***************************************************************************
// xstructure AssignNameXstr(xstructure& xstr, vector<string> names)
// ***************************************************************************
xstructure AssignNameXstr(xstructure& xstr, vector<string> names) {
  //Assign names to xstr
  xstructure xstr_new = xstr;
  int iatom=0;
  for (uint itype=0; itype<xstr.num_each_type.size();itype++) {
    string species = string(names.at(xstr_new.atoms.at(iatom).type));
    xstr_new.species.at(itype)=species;
    for (int j=0;j<xstr.num_each_type.at(itype);j++) {
      xstr_new.atoms.at(iatom).name=species;
      xstr_new.atoms.at(iatom).CleanName();
      xstr_new.atoms.at(iatom).CleanSpin();
      xstr_new.atoms.at(iatom).name_is_given=TRUE;
      iatom++;
    }
  }
  return xstr_new;
}


// ***************************************************************************
// vector<xstructure> AssignNameXstr(vector<xstructure>& vxstr, vector<string> names)
// ***************************************************************************
vector<xstructure> AssignNameXstr(vector<xstructure>& vxstr, vector<string> names) {
  vector<xstructure> vxstr_new;
  xstructure xstr_tmp;
  for (uint i=0; i<vxstr.size(); i++) {
    xstr_tmp = AssignNameXstr(vxstr.at(i), names);
    vxstr_new.push_back(xstr_tmp);
  }
  return vxstr_new;
}


// ***************************************************************************
// pocc::HNFCELL
// ***************************************************************************
namespace pocc {
  void HNFCELL(istream& input) {
    xstructure xstr_ori(input,IOVASP_AUTO);
    vector<xstructure> vxstr = Partial2Supercell(xstr_ori);
  }
}

// ***************************************************************************
// vector<xstructure> Partial2Supercell(xstructure xstr)
// ***************************************************************************
vector<xstructure> Partial2Supercell(xstructure xstr_ori) {
  xstructure xstr = AssignAlphabeticNameXstr(xstr_ori);
  string str_species_ori;
  if(xstr_ori.atoms.at(0).name_is_given) {str_species_ori = xstr_ori.SpeciesString();}
  else {str_species_ori = "As Ba Cl Dy Eu Fe Ga He In K";}
  vector<string> vxstr_species_ori;
  aurostd::string2tokens(str_species_ori, vxstr_species_ori, " ");

  ofstream FileMESSAGE;
  FileMESSAGE.open("LOG.POCC");
  _aflags aflags;

  int nHNF = InitializeXstr(xstr, vxstr_species_ori, FileMESSAGE, aflags);
  vector<vector<xstructure> > vv_xstr;
  if(CheckFivefoldOccupied(xstr)) { //fivefold or quadruple
    vv_xstr = Partial2Xstr_Fivefold_Occupied(xstr, nHNF, FileMESSAGE, aflags);
  }
  else if(CheckFourfoldOccupied(xstr)) { //fourfold or quadruple
    vv_xstr = Partial2Xstr_Fourfold_Occupied(xstr, nHNF, FileMESSAGE, aflags);
  }
  else if(CheckTripleOccupied(xstr)) { //Triple
    vv_xstr = Partial2Xstr_Triple_Occupied(xstr, nHNF, FileMESSAGE, aflags);
  }
  else if(CheckDoubleOccupied(xstr)) {  //Double
    vv_xstr = Partial2Xstr(xstr, nHNF, FileMESSAGE, aflags);
  }

  vector<xstructure> vxstr;
  for (uint i=0; i<vv_xstr.size();i++) {
    for (uint j=0; j<vv_xstr.at(i).size();j++) {
      xstructure xstr_tmp = vv_xstr.at(i).at(j);
      xstructure xstr_tmp_name = AssignNameXstr(xstr_tmp, vxstr_species_ori);
      xstr_tmp_name.SpeciesPutAlphabetic();
      vxstr.push_back(xstr_tmp_name);
    }
  }
  vector<xstructure> vxstr_sorted = pocc::SortGroupXstrUFFEnergy(vxstr);

  ostringstream aus;
  aus << AFLOWIN_SEPARATION_LINE << endl;
  aus << "0000 MESSAGE    Printing Derivative Supercells ... " << Message(__AFLOW_FILE__,aflags) << endl;
  for (uint i=0; i<vxstr_sorted.size();i++) {
    aus << AFLOWIN_SEPARATION_LINE << endl;
    aus << "[VASP_POSCAR_MODE_EXPLICIT]START.HNFCELL_" << i+1 << endl;
    aus << vxstr_sorted.at(i);
    aus << "[VASP_POSCAR_MODE_EXPLICIT]STOP.HNFCELL_" << i+1 << endl;
    aus << AFLOWIN_SEPARATION_LINE << endl;
    aus.flush();
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
    aus.str("");
  }
  FileMESSAGE.close();
  return vxstr_sorted;
}

// ***************************************************************************
// vector<xstructure> CalculateInitialSupercell(xstructure xstr, int n, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
vector<xstructure> CalculateInitialSupercell(xstructure xstr, int n, ofstream &FileMESSAGE, _aflags &aflags) {
  XHOST.QUIET=FALSE;
  ostringstream aus;
  aus << "0000 MESSAGE    Calculating HNF index ... " << Message(__AFLOW_FILE__,aflags) << endl;
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
  aus.str("");
  vector<xmatrix<double> > sHNF=CalculateHNF(xstr, n);
  vector<xstructure> vxstr_init;
  xstructure b_tmp;
  aus << "0000 MESSAGE    Generating supercells ... " << Message(__AFLOW_FILE__,aflags) << endl;
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
  aus.str("");
  stringstream ss;
  aus << AFLOWIN_SEPARATION_LINE << endl;
  for (int i=0; i<(int) sHNF.size();i++) {
    b_tmp=GetSuperCell(xstr,trasp(sHNF.at(i)));//when generates structure, using trasp
    b_tmp.title+="  [HNF("+aurostd::utype2string(n)+")="+aurostd::utype2string(i+1)+"/"+aurostd::utype2string(sHNF.size())+"=";
    string str_sHNF = " " 
      + aurostd::utype2string(sHNF.at(i)(1,1))+" " 
      + aurostd::utype2string(sHNF.at(i)(1,2))+" "
      + aurostd::utype2string(sHNF.at(i)(1,3))+"; "
      + aurostd::utype2string(sHNF.at(i)(2,1))+" "
      + aurostd::utype2string(sHNF.at(i)(2,2))+" "
      + aurostd::utype2string(sHNF.at(i)(2,3))+"; "
      + aurostd::utype2string(sHNF.at(i)(3,1))+" "
      + aurostd::utype2string(sHNF.at(i)(3,2))+" "
      + aurostd::utype2string(sHNF.at(i)(3,3))+"]";
    b_tmp.title+=str_sHNF;
    b_tmp.partial_occupation_flag=FALSE; //Not partial occupation anymore
    vxstr_init.push_back(b_tmp);
    aus << "[VASP_POSCAR_MODE_EXPLICIT]START.HNF_" << n << "_" << i+1 << "_" << sHNF.size() << endl;
    aus << CleanNameXstr(b_tmp);
    aus << "[VASP_POSCAR_MODE_EXPLICIT]STOP.HNF_" << n << "_" << i+1 << "_" << sHNF.size() << endl;
    aus << AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
    aus.str("");
  }
  aus << "0000 MESSAGE    Total number of initial supercells is " << vxstr_init.size() <<" " << Message(__AFLOW_FILE__,aflags) << endl;
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
  aus.str("");

  return vxstr_init;
}

// ***************************************************************************
// int InitializeXstr(xstructure &xstr, vector<string> vxstr_species_ori, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
int InitializeXstr(xstructure &xstr, vector<string> vxstr_species_ori, ofstream &FileMESSAGE, _aflags &aflags) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  //This part is based on "pflow::HNFTOL" 
  ostringstream oss;
  xstr.BringInCell();
  //xstr.CalculateSymmetryPointGroup(); //CO, cannot calculate symmetry of a structure with overlapping atoms
  xstr.ReScale(1.0);

  uint digits=6, digits1=4, digits2=2*digits+11;

  if(xstr.partial_occupation_HNF) {return xstr.partial_occupation_HNF;} //if HNF exists, no need to optimize pocc values
  double tolerance=DEFAULT_POCC_SITE_TOL; //DEFAULT_PARTIAL_OCCUPATION_TOLERANCE;
  if(LDEBUG) { //CO20180409
    cerr << "default " << DEFAULT_POCC_SITE_TOL << endl; //DEFAULT_PARTIAL_OCCUPATION_TOLERANCE << endl;
    cerr << xstr.partial_occupation_site_tol << endl;
  }
  if(xstr.partial_occupation_site_tol>0) {tolerance=xstr.partial_occupation_site_tol;}
  oss << "[AFLOW] hnf_tolerance=" << tolerance << endl;
  oss << "Printing Input PARTCAR ... " << Message(__AFLOW_FILE__,aflags) << endl;
  oss << AFLOWIN_SEPARATION_LINE << endl;    // --------------------------------
  oss << "[VASP_POSCAR_MODE_EXPLICIT]START" << endl;
  oss << AssignNameXstr(xstr, vxstr_species_ori);// << endl;
  oss << "[VASP_POSCAR_MODE_EXPLICIT]STOP" << endl;
  oss << AFLOWIN_SEPARATION_LINE << endl;    // --------------------------------
  oss << "Optimizing HNF number ... " << endl;
  oss << AFLOWIN_SEPARATION_LINE << endl;    // --------------------------------
  oss.precision(digits);
  double error=1.0,eps;
  int HNFi=0;
  vector<double> error_iatom;  vector<double> effective_pocc_iatom;  vector<uint> M_iatom;
  for(uint iatom=0;iatom<xstr.atoms.size();iatom++) {error_iatom.push_back(1.0);}
  for(uint iatom=0;iatom<xstr.atoms.size();iatom++) {effective_pocc_iatom.push_back(0.0);}
  for(uint iatom=0;iatom<xstr.atoms.size();iatom++) {M_iatom.push_back(0);}

  oss << aurostd::PaddedPRE(aurostd::utype2string(0),digits1) << "  " << "--------" << " ";
  for(uint iatom=0;iatom<xstr.atoms.size();iatom++) {
    if(xstr.atoms.at(iatom).partial_occupation_value<1.0) {
      oss << "| "+aurostd::PaddedPOST("iatom="+aurostd::utype2string(iatom+1)+"/"+aurostd::utype2string(xstr.atoms.size()),digits2) << " " ;  //CO20170629
    }
  }
  oss << " | " << "error" << endl;
  if(LDEBUG) { //CO20180409
    cerr << "ERROR " << error << endl;
    cerr << "TOL " << tolerance << endl;
  }
  while(error>tolerance) {
    HNFi++;
    error=1.0/HNFi;
    oss << aurostd::PaddedPRE(aurostd::utype2string(HNFi),digits1) << "  " << error << " ";//endl;
    for(uint iatom=0;iatom<xstr.atoms.size();iatom++) {
      error_iatom.at(iatom)=1.0;
      M_iatom.at(iatom)=0;
      for(int j=0;j<=HNFi;j++) {
        eps=abs((double) j/HNFi-xstr.atoms.at(iatom).partial_occupation_value);
        if(eps<error_iatom.at(iatom)) {
          error_iatom.at(iatom)=eps;
          M_iatom.at(iatom)=(uint) j;
          effective_pocc_iatom.at(iatom)=(double) M_iatom.at(iatom)/HNFi;
        }
      }
    }

    for(uint iatom=0;iatom<xstr.atoms.size();iatom++) {
      if(xstr.atoms.at(iatom).partial_occupation_value<1.0) {
        stringstream aus;aus.precision(digits);aus.setf(std::ios::fixed,std::ios::floatfield);
        aus << effective_pocc_iatom.at(iatom) << "," << error_iatom.at(iatom) << "," << M_iatom.at(iatom) << "/" << HNFi;
        oss << "| " << aurostd::PaddedPOST(aus.str(),digits2) << " " ;
      }
    }
    error=aurostd::max(error_iatom);
    oss << " | " << error << endl;
  }
  oss << AFLOWIN_SEPARATION_LINE << endl;

  for(uint iatom=0;iatom<xstr.atoms.size();iatom++) {
    if(xstr.atoms.at(iatom).partial_occupation_value<1.0) {xstr.atoms.at(iatom).partial_occupation_value = effective_pocc_iatom.at(iatom) ;}
  }
  UpdateXstr_comp_each_type(xstr); //update partial occupation
  oss << "Printing optimized input PARTCAR ..." << Message(__AFLOW_FILE__,aflags) << endl;
  oss << AFLOWIN_SEPARATION_LINE << endl;
  oss << "[VASP_POSCAR_MODE_EXPLICIT]START" << endl;
  oss << AssignNameXstr(xstr, vxstr_species_ori);// << endl;
  oss << "[VASP_POSCAR_MODE_EXPLICIT]STOP" << endl;
  oss << AFLOWIN_SEPARATION_LINE << endl;

  aurostd::PrintMessageStream(FileMESSAGE, oss,XHOST.QUIET);
  int nHNF=HNFi;
  return nHNF;
}

// ***************************************************************************
//vector<vector<xstructure> > Partial2Xstr_Fivefold_Occupied(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
vector<vector<xstructure> > Partial2Xstr_Fivefold_Occupied(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags) {
  XHOST.QUIET= TRUE;
  //convert fivefold into fourfold firstly
  vector<vector<int> > NormalisedNumber=NormalisedNumberXstructure(xstr);
  xstructure xstr_norm = NormalizeXstructure(xstr);
  xstructure xstr_tri = xstr;
  for (uint i=0;i<xstr_norm.atoms.size();i++) {
    if(NormalisedNumber.at(i).size()==5) {
      int k=NormalisedNumber.at(i).at(4);  //begin from 0, so the 5th atom should be at(4)
      xstr_tri.RemoveAtom(k);
    }
  }
  //Call fourfold routine
  vector<vector<xstructure> > vec_vec_xstr = Partial2Xstr_Fourfold_Occupied(xstr_tri, nHNF, FileMESSAGE, aflags);
  ostringstream aus;

  vector<vector<vector<int> > > vec_vec_vec_num_vacancy;
  vector<vector<int> > atom_num_supercell = CalculateXstrNumberSupercell(xstr_norm, nHNF);
  vector<string>  element_name; //replacing original name
  for (uint i=0; i<xstr_norm.atoms.size();i++) {
    if((NormalisedNumber.at(i).size()==5) || 
        ((NormalisedNumber.at(i).size()==4)&& CheckVacancyOnOnesite(xstr,i)))
    { //tripled occupied
      int k1= NormalisedNumber.at(i).at(0);
      int k2= NormalisedNumber.at(i).at(1);
      int k3= NormalisedNumber.at(i).at(2);
      int k4= NormalisedNumber.at(i).at(3);
      double sum_first_four_pocc_value = xstr.atoms.at(k1).partial_occupation_value
        + xstr.atoms.at(k2).partial_occupation_value 
        + xstr.atoms.at(k3).partial_occupation_value
        + xstr.atoms.at(k4).partial_occupation_value;
      int num_vacancy= hnf_double2int((1-sum_first_four_pocc_value)*nHNF);
      for (int k=0; k<num_vacancy; k++) {
        if(NormalisedNumber.at(i).size()==5) { //If it contains four elements
          int k5=NormalisedNumber.at(i).at(4); //get the the name of the 5th element
          element_name.push_back(xstr.atoms.at(k5).name);
        }
        else { //if it contains three elements and one vacancy
          element_name.push_back(name_vacancy);
        }
      }
      vector<int> vec_num_part=atom_num_supercell.at(i);
      //get rid of the numbers of the first three atom
      //the 4th and 5th atoms will be thought as one entry, so we need to remove the first three
      float first_three_pocc_value = xstr.atoms.at(k1).partial_occupation_value 
        + xstr.atoms.at(k2).partial_occupation_value
        + xstr.atoms.at(k3).partial_occupation_value;
      int num_first_three = hnf_double2int(first_three_pocc_value*nHNF);
      vector<int> vec_num_part_new; 
      for (uint i=num_first_three; i<vec_num_part.size();i++) {
        vec_num_part_new.push_back(vec_num_part.at(i));
      }
      vector<vector<int> > num_com_vacancy; num_com_vacancy.clear();
      combine(vec_num_part_new, num_com_vacancy, num_vacancy);
      vec_vec_vec_num_vacancy.push_back(num_com_vacancy);
    }
  }
  vector<vector<vector<int> > > final_num_result;
  vector<vector<int> >  tmp_vect;
  vector<vector<int> > formated_final_num;

  CombineAll(vec_vec_vec_num_vacancy, 0, tmp_vect, final_num_result);
  vector<int> tmp;
  for (uint i=0; i<final_num_result.size();i++) {
    tmp.clear();
    for (uint j=0; j<final_num_result.at(i).size();j++) {
      for (uint k=0; k<final_num_result.at(i).at(j).size();k++) {
        tmp.push_back(final_num_result.at(i).at(j).at(k));
      }
    }
    formated_final_num.push_back(tmp);
  }


  vector<vector<xstructure> > vec_vec_xstr_fivefold;
  int total=0;
  for (uint k=0; k<vec_vec_xstr.size();k++) {  //2-D vector
    vector<xstructure> vec_xstr=vec_vec_xstr.at(k);
    vector<xstructure> vec_new_xstr;
    for (uint i=0; i<vec_xstr.size();i++) {
      xstructure xstr_tmp = vec_xstr.at(i);
      for (uint j=0; j<formated_final_num.size();j++) {
        total++;
        vector<int> vec_int = formated_final_num.at(j);
        xstructure new_xstr = XstrSubstitute(xstr_tmp, vec_int, element_name);
        vec_new_xstr.push_back(new_xstr);
      }
    }
    vec_vec_xstr_fivefold.push_back(vec_new_xstr);
  }

  vector<vector<xstructure> > vec_vec_xstr_final;
  for (uint i=0; i<vec_vec_xstr_fivefold.size();i++) {
    vector<xstructure> vec_new_xstr = vec_vec_xstr_fivefold.at(i);
    vector<xstructure> vxstr_final = RemoveEquivalentXstr(vec_new_xstr, FileMESSAGE, aflags); 
    vec_vec_xstr_final.push_back(vxstr_final);
  }

  return vec_vec_xstr_final;
}

// ***************************************************************************
// vector<vector<xstructure> > Partial2Xstr_Fourfold_Occupied(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
vector<vector<xstructure> > Partial2Xstr_Fourfold_Occupied(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags) {
  XHOST.QUIET= TRUE;
  //convert fourfold into triple
  vector<vector<int> > NormalisedNumber=NormalisedNumberXstructure(xstr);
  xstructure xstr_norm = NormalizeXstructure(xstr);
  xstructure xstr_tri = xstr;
  for (uint i=0;i<xstr_norm.atoms.size();i++) {
    if(NormalisedNumber.at(i).size()==4) { //if the 4th element is vacancy, then it is fake triple one
      int k=NormalisedNumber.at(i).at(3);
      xstr_tri.RemoveAtom(k);
    }
  }

  vector<vector<xstructure> > vec_vec_xstr = Partial2Xstr_Triple_Occupied(xstr_tri, nHNF, FileMESSAGE, aflags);
  ostringstream aus;

  vector<vector<vector<int> > > vec_vec_vec_num_vacancy;
  vector<vector<int> > atom_num_supercell = CalculateXstrNumberSupercell(xstr_norm, nHNF);
  vector<string>  element_name; //replacing original name
  for (uint i=0; i<xstr_norm.atoms.size();i++) {
    if((NormalisedNumber.at(i).size()==4) || 
        ((NormalisedNumber.at(i).size()==3)&& CheckVacancyOnOnesite(xstr,i)))
    { //tripled occupied
      int k1= NormalisedNumber.at(i).at(0);
      int k2= NormalisedNumber.at(i).at(1);
      int k3= NormalisedNumber.at(i).at(2);
      double sum_first_three_pocc_value = xstr.atoms.at(k1).partial_occupation_value
        + xstr.atoms.at(k2).partial_occupation_value + xstr.atoms.at(k3).partial_occupation_value;
      int num_vacancy= hnf_double2int((1-sum_first_three_pocc_value)*nHNF);
      for (int k=0; k<num_vacancy; k++) {
        if(NormalisedNumber.at(i).size()==4) { //If it contains four elements
          int k4=NormalisedNumber.at(i).at(3); //get the the name of  4th element
          element_name.push_back(xstr.atoms.at(k4).name);
        }
        else { //if it contains three elements and one vacancy
          element_name.push_back(name_vacancy);
        }
      }
      vector<int> vec_num_part=atom_num_supercell.at(i);
      //get rid of the numbers of the first two atom
      float first_two_pocc_value = xstr.atoms.at(k1).partial_occupation_value + xstr.atoms.at(k2).partial_occupation_value;
      int num_first_two = hnf_double2int(first_two_pocc_value*nHNF);
      vector<int> vec_num_part_new; 
      for (uint i=num_first_two; i<vec_num_part.size();i++) {
        vec_num_part_new.push_back(vec_num_part.at(i));
      }
      vector<vector<int> > num_com_vacancy; num_com_vacancy.clear();
      combine(vec_num_part_new, num_com_vacancy, num_vacancy);
      vec_vec_vec_num_vacancy.push_back(num_com_vacancy);
    }
  }
  vector<vector<vector<int> > > final_num_result;
  vector<vector<int> >  tmp_vect;
  vector<vector<int> > formated_final_num;

  CombineAll(vec_vec_vec_num_vacancy, 0, tmp_vect, final_num_result);
  vector<int> tmp;
  for (uint i=0; i<final_num_result.size();i++) {
    tmp.clear();
    for (uint j=0; j<final_num_result.at(i).size();j++) {
      for (uint k=0; k<final_num_result.at(i).at(j).size();k++) {
        tmp.push_back(final_num_result.at(i).at(j).at(k));
      }
    }
    formated_final_num.push_back(tmp);
  }


  vector<vector<xstructure> > vec_vec_xstr_fourfold;
  int total=0;
  for (uint k=0; k<vec_vec_xstr.size();k++) {  //2-D vector
    vector<xstructure> vec_xstr=vec_vec_xstr.at(k);
    vector<xstructure> vec_new_xstr;
    for (uint i=0; i<vec_xstr.size();i++) {
      xstructure xstr_tmp = vec_xstr.at(i);
      for (uint j=0; j<formated_final_num.size();j++) {
        total++;
        vector<int> vec_int = formated_final_num.at(j);
        xstructure new_xstr = XstrSubstitute(xstr_tmp, vec_int, element_name);
        vec_new_xstr.push_back(new_xstr);
      }
    }
    vec_vec_xstr_fourfold.push_back(vec_new_xstr);
  }

  vector<vector<xstructure> > vec_vec_xstr_final;
  for (uint i=0; i<vec_vec_xstr_fourfold.size();i++) {
    vector<xstructure> vec_new_xstr = vec_vec_xstr_fourfold.at(i);
    vector<xstructure> vxstr_final = RemoveEquivalentXstr(vec_new_xstr, FileMESSAGE, aflags); 
    vec_vec_xstr_final.push_back(vxstr_final);
  }
  return vec_vec_xstr_final;
}

// ***************************************************************************
// vector<vector<xstructure> > Partial2Xstr_Triple_Occupied(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
vector<vector<xstructure> > Partial2Xstr_Triple_Occupied(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags) {
  XHOST.QUIET=TRUE;
  vector<vector<int> > NormalisedNumber=NormalisedNumberXstructure(xstr);
  xstructure xstr_norm = NormalizeXstructure(xstr);
  //convert triple into double, in fact, it is not very necessary here
  xstructure xstr_tri = xstr;
  for (uint i=0;i<xstr_norm.atoms.size();i++) {
    if(NormalisedNumber.at(i).size()==3) { //if the 3th element is vacancy, then it is fake triple one
      int k=NormalisedNumber.at(i).at(2);
      xstr_tri.RemoveAtom(k);
    }
  }

  vector<vector<xstructure> > vec_vec_xstr = Partial2Xstr(xstr_tri, nHNF, FileMESSAGE, aflags);
  ostringstream aus;

  vector<vector<vector<int> > > vec_vec_vec_num_vacancy;
  vector<vector<int> > atom_num_supercell = CalculateXstrNumberSupercell(xstr_norm, nHNF);
  vector<string>  element_name; //replacing original name
  for (uint i=0; i<xstr_norm.atoms.size();i++) {
    if((NormalisedNumber.at(i).size()==3) || 
        ((NormalisedNumber.at(i).size()==2)&& CheckVacancyOnOnesite(xstr,i)))
    { //tripled occupied
      int k1= NormalisedNumber.at(i).at(0);
      int k2= NormalisedNumber.at(i).at(1);
      double sum_first_two_pocc_value = xstr.atoms.at(k1).partial_occupation_value
        + xstr.atoms.at(k2).partial_occupation_value;
      int num_vacancy= hnf_double2int((1-sum_first_two_pocc_value)*nHNF);
      for (int k=0; k<num_vacancy; k++) {
        if(NormalisedNumber.at(i).size()==3) { //If it contains three elements
          int k3=NormalisedNumber.at(i).at(2);
          element_name.push_back(xstr.atoms.at(k3).name);
        }
        else { //if it contains two elements and one vacancy
          element_name.push_back(name_vacancy);
        }
      }
      vector<int> vec_num_part=atom_num_supercell.at(i);
      //get rid of the numbers of the first atom
      int num_first_one = hnf_double2int(xstr.atoms.at(k1).partial_occupation_value*nHNF);
      vector<int> vec_num_part_new; 
      for (uint i=num_first_one; i<vec_num_part.size();i++) {
        vec_num_part_new.push_back(vec_num_part.at(i));
      }
      vector<vector<int> > num_com_vacancy; num_com_vacancy.clear();
      combine(vec_num_part_new, num_com_vacancy, num_vacancy);
      vec_vec_vec_num_vacancy.push_back(num_com_vacancy);
    }
  }
  vector<vector<vector<int> > > final_num_result;
  vector<vector<int> >  tmp_vect;
  vector<vector<int> > formated_final_num;

  CombineAll(vec_vec_vec_num_vacancy, 0, tmp_vect, final_num_result);
  vector<int> tmp;
  for (uint i=0; i<final_num_result.size();i++) {
    tmp.clear();
    for (uint j=0; j<final_num_result.at(i).size();j++) {
      for (uint k=0; k<final_num_result.at(i).at(j).size();k++) {
        tmp.push_back(final_num_result.at(i).at(j).at(k));
      }
    }
    formated_final_num.push_back(tmp);
  }


  vector<vector<xstructure> > vec_vec_xstr_triple;
  int total=0;
  for (uint k=0; k<vec_vec_xstr.size();k++) {  //2-D vector
    vector<xstructure> vec_xstr=vec_vec_xstr.at(k);
    vector<xstructure> vec_new_xstr;
    for (uint i=0; i<vec_xstr.size();i++) {
      xstructure xstr_tmp = vec_xstr.at(i);
      for (uint j=0; j<formated_final_num.size();j++) {
        total++;
        vector<int> vec_int = formated_final_num.at(j);
        xstructure new_xstr = XstrSubstitute(xstr_tmp, vec_int, element_name);
        vec_new_xstr.push_back(new_xstr);
      }
    }
    vec_vec_xstr_triple.push_back(vec_new_xstr);
  }

  vector<vector<xstructure> > vec_vec_xstr_final;
  for (uint i=0; i<vec_vec_xstr_triple.size();i++) {
    vector<xstructure> vec_new_xstr = vec_vec_xstr_triple.at(i);
    vector<xstructure> vxstr_final = RemoveEquivalentXstr(vec_new_xstr, FileMESSAGE, aflags); 
    vec_vec_xstr_final.push_back(vxstr_final);
  }
  return vec_vec_xstr_final;
}


// ***************************************************************************
// vector<vector<xstructure> > Partial2Xstr(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
vector<vector<xstructure> > Partial2Xstr(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags) {
  //doubly occupied
  xstr.ReScale(1.0);
  xstr.iomode=IOVASP_POSCAR;

  XHOST.QUIET= FALSE;
  ostringstream aus;
  aus << "0000 MESSAGE    POCC running " << Message(__AFLOW_FILE__,aflags) << endl;
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
  aus.str("");
  xstructure xstr_c = NormalizeXstructure(xstr);  

  int n = nHNF;
  vector<xstructure> vxstr_init =  CalculateInitialSupercell(xstr_c, n, FileMESSAGE, aflags);
  aus << "0000 MESSAGE    Generating derivate supercells ..." << Message(__AFLOW_FILE__,aflags) << endl;
  aus << "0000 MESSAGE    Please be patient ..." << Message(__AFLOW_FILE__,aflags) << endl;
  aus << AFLOWIN_SEPARATION_LINE << endl;
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);

  vector<xstructure> vxstr_mid; vxstr_mid.clear(); 
  vector<vector<xstructure> > vv_xstr; vv_xstr.clear();
  for (uint i=0; i<vxstr_init.size();i++) {
    xstructure xstr_tmp =vxstr_init.at(i);
    vxstr_mid=AssignPartialOccupation(xstr_tmp, xstr, n);
    vv_xstr.push_back(vxstr_mid);
  }

  vector<xstructure> vxstr_mid2; vxstr_mid2.clear(); 
  vector<vector<xstructure> > vv_xstr2; vv_xstr2.clear();
  int total=0;
  for (uint i=0; i<vv_xstr.size();i++) {
    for (uint j=0; j<vv_xstr.at(i).size();j++) {
      total++;
      aus << AFLOWIN_SEPARATION_LINE << endl;
      aus << "[VASP_POSCAR_MODE_EXPLICIT]START.HNF_" << i+1 << "_" << j+1 << "_" << total << endl;
      //string str_title = " HNF_" + aurostd::double2string(i+1) + "_" + aurostd::double2string(j+1) + "_" + aurostd::double2string(total);
      //vv_xstr.at(i).at(j).title+=str_title;
      aus << vv_xstr.at(i).at(j);
      vxstr_mid2.push_back(vv_xstr.at(i).at(j));
      aus << "[VASP_POSCAR_MODE_EXPLICIT]STOP.HNF_" << i+1 << "_" << j+1 << "_" << total << endl;
      //aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
      aus.str("");
    }
    vv_xstr2.push_back(vxstr_mid2);
    vxstr_mid2.clear();
  }

  vector<vector<xstructure> > vv_xstr_final; vv_xstr_final.clear();
  for (uint i=0;i<vv_xstr2.size();i++) {
    vector<xstructure> vxstr_tmp = vv_xstr2.at(i);
    vector<xstructure> vxstr_nq2 = RemoveEquivalentXstr(vxstr_tmp, FileMESSAGE, aflags); //
    vv_xstr_final.push_back(vxstr_nq2);
  } 
  return vv_xstr_final;
}

// ***************************************************************************
// vector<xstructure> CalculatePrimitiveCell(vector<xstructure> &vec_xstr)
// ***************************************************************************
vector<xstructure> CalculatePrimitiveCell(vector<xstructure> &vec_xstr) {
  cerr << "Calculating Pritimive cells" << endl;
  vector<xstructure> vec_xstr_sp;
  vec_xstr_sp.clear();
  xstructure xstr_tmp_sp, a;
  for (uint i=0; i<vec_xstr.size();i++) {
    a = vec_xstr.at(i);
    xstr_tmp_sp = GetStandardPrimitive(a);
    vec_xstr_sp.push_back(xstr_tmp_sp);
    cerr << i+1 << " is finished!" << endl;
  }
  return vec_xstr_sp;
}
// ***************************************************************************
// bool comparison_atom_fpos(_atom atom1, _atom atom2)
// ***************************************************************************
bool comparison_atom_fpos(_atom atom1, _atom atom2) {
  double _EQUAL_DOUBLE = 1E-6;
  if( abs(atom1.fpos[1] - atom2.fpos[1]) > _EQUAL_DOUBLE ) {
    return (atom1.fpos[1] < atom2.fpos[1]);
  } else if( abs(atom1.fpos[2] - atom2.fpos[2]) > _EQUAL_DOUBLE ) {
    return (atom1.fpos[2] < atom2.fpos[2]);
  } else if( abs(atom1.fpos[3] - atom2.fpos[3]) > _EQUAL_DOUBLE ) {
    return (atom1.fpos[3] < atom2.fpos[3]);
  } else {
    return true;
  }
}

// ***************************************************************************
// pflow::Sort_atom_fpos(xstructure &xstr)
// ***************************************************************************
namespace pflow {
  void Sort_atom_fpos(xstructure &xstr) {
    xstr.BringInCell();
    xstr.SetCoordinates(_COORDS_FRACTIONAL_);
    int begin = 0;
    for (uint i=0; i <xstr.num_each_type.size(); i++) {
      int j= xstr.num_each_type.at(i);
      sort(xstr.atoms.begin()+begin, xstr.atoms.begin()+begin+j, comparison_atom_fpos);
      begin += j;
    }
  }
}

// ***************************************************************************
// bool comparison_atom_cpos(_atom atom1, _atom atom2)
// ***************************************************************************
bool comparison_atom_cpos(_atom atom1, _atom atom2) {
  double _EQUAL_DOUBLE = 1E-6;
  if( abs(atom1.cpos[1] - atom2.cpos[1]) > _EQUAL_DOUBLE ) {
    return (atom1.cpos[1] < atom2.cpos[1]);
  } else if( abs(atom1.cpos[2] - atom2.cpos[2]) > _EQUAL_DOUBLE ) {
    return (atom1.cpos[2] < atom2.cpos[2]);
  } else if( abs(atom1.cpos[3] - atom2.cpos[3]) > _EQUAL_DOUBLE ) {
    return (atom1.cpos[3] < atom2.cpos[3]);
  } else {
    return true;
  }
}

// ***************************************************************************
// pflow::Sort_atom_cpos(xstructure &xstr)
// ***************************************************************************
namespace pflow {
  void Sort_atom_cpos(xstructure &xstr) {
    xstr.BringInCell();
    xstr.SetCoordinates(_COORDS_CARTESIAN_);
    int begin = 0;
    for (uint i=0; i <xstr.num_each_type.size(); i++) {
      int j= xstr.num_each_type.at(i);
      sort(xstr.atoms.begin()+begin, xstr.atoms.begin()+begin+j, comparison_atom_cpos);
      begin += j;
    }
  }
}

// ***************************************************************************
// bool LatticeCompare(xstructure xstr1, xstructure xstr2)
// ***************************************************************************
bool LatticeCompare(xstructure xstr1, xstructure xstr2) {
  bool RUN_FLAG = false;
  double tmp1=0, tmp2=0, sum=0, epsilon=1E-4;
  tmp1 = abs(xstr1.a-xstr2.a) + abs(xstr1.b-xstr2.b) + abs(xstr1.c-xstr2.c);
  tmp2 = abs(xstr1.alpha-xstr2.alpha) + abs(xstr1.beta-xstr2.beta) + abs(xstr1.gamma-xstr2.gamma);
  sum = tmp1 +tmp2;  //Make sure sum is absolute value
  if(sum <epsilon) {RUN_FLAG = true;}
  return RUN_FLAG;
}

// ***************************************************************************
// bool CoordCompare(xstructure xstr1, xstructure xstr2)
// ***************************************************************************
bool CoordCompare(xstructure xstr1, xstructure xstr2) {
  bool RUN_FLAG = false;
  double sum=0, epsilon=1E-4;
  if(!LatticeCompare(xstr1, xstr2)) {return RUN_FLAG;}
  else  {
    if(xstr1.atoms.size()==xstr2.atoms.size()) {
      for (uint i=0;i<xstr1.atoms.size();i++) {
        double b1=0, b2=0, b3=0;
        b1= xstr1.atoms.at(i).fpos[1] - xstr2.atoms.at(i).fpos[1];
        b2= xstr1.atoms.at(i).fpos[2] - xstr2.atoms.at(i).fpos[2];
        b3= xstr1.atoms.at(i).fpos[3] - xstr2.atoms.at(i).fpos[3];
        sum += (b1*b1+b2*b2+b3*b3);
      }
      if(sum <epsilon) {RUN_FLAG = true;}
    }
  }
  return RUN_FLAG;
}

// ***************************************************************************
// vector<xstr_atom_number> CalculateXstrNumEachType(xstructure xstr)
// ***************************************************************************
vector<xstr_atom_number> CalculateXstrNumEachType(xstructure xstr) {
  vector<xstr_atom_number> vec_xstr_number;
  int begin=0;
  for (uint i=0; i<xstr.num_each_type.size();i++) {
    int k=xstr.num_each_type.at(i);
    xstr_atom_number atom_number_tmp;
    atom_number_tmp.number = k;
    for (int j=0; j<k; j++) {
      atom_number_tmp.vec_atom_number.push_back(begin);
      begin++;
    }
    vec_xstr_number.push_back(atom_number_tmp);
  }
  return vec_xstr_number;
}

bool sort_atom_number(xstr_atom_number str1, xstr_atom_number str2) {
  return(str1.number<str2.number);
}

// ***************************************************************************
// vector<strno_energy> xstructure2strno_energy(vector<xstructure> vec_xstr)
// ***************************************************************************
vector<xstr_energy> xstructure2xstr_energy(vector<xstructure> vec_xstr) {
  vector<xstr_energy> vxstr_energy;
  xstr_energy xstr_energy_tmp;
  for (uint i=0; i<vec_xstr.size();i++) {
    //xstructure xstr_tmp = vec_xstr.at(i);
    //xstr_energy_tmp.xstr = xstr_tmp;
    xstr_energy_tmp.xstr = vec_xstr.at(i);
    double TotalEnergyTmp;
    if(CALCULATE_ENERGY_UFF) {  
      //TotalEnergyTmp = pocc::CalculateUFFEnergy(xstr_tmp);
      TotalEnergyTmp = pocc::CalculateUFFEnergy(vec_xstr.at(i));
    }
    else {
      //TotalEnergyTmp = pocc::CalculateEnergyUsingGulp(xstr_tmp);
      TotalEnergyTmp = pocc::CalculateEnergyUsingGulp(vec_xstr.at(i));
    }
    xstr_energy_tmp.energy = TotalEnergyTmp;
    vxstr_energy.push_back(xstr_energy_tmp);
  }
  return vxstr_energy;
}

// ***************************************************************************
// vector<xstructure> RemoveEquivalentXstr(vector<xstructure> vec_xstr, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
vector<xstructure> RemoveEquivalentXstr(vector<xstructure> vec_xstr, ofstream &FileMESSAGE, _aflags &aflags) {
  vector<xstr_energy> vec_xstr_energy = xstructure2xstr_energy(vec_xstr);
  const double epsilon_etot = 1E-6;
  vector<xstr_energy> vec_xstr_energy_final;
  vec_xstr_energy_final.push_back(vec_xstr_energy.at(0));
  for (uint i=1; i<vec_xstr_energy.size();i++) {
    uint j;
    for (j=0; j<vec_xstr_energy_final.size();j++) {
      if(abs(vec_xstr_energy.at(i).energy-vec_xstr_energy_final.at(j).energy)<epsilon_etot) { //if they have same energy
        break;
      }
    }
    if(j==vec_xstr_energy_final.size()) {
      vec_xstr_energy_final.push_back(vec_xstr_energy.at(i));
    }
  }

  //cerr << "KY TEST" << endl;
  vector<int> vecDG; vecDG.clear(); //Degenracy
  for (uint i=0;i<vec_xstr_energy_final.size();i++) {
    int DGi =0;
    for (uint j=0; j<vec_xstr_energy.size();j++) {
      if(abs(vec_xstr_energy_final.at(i).energy- vec_xstr_energy.at(j).energy)<epsilon_etot) {DGi++;}
    }
    vecDG.push_back(DGi);
  }

  int number = vec_xstr_energy.size() - vec_xstr_energy_final.size();
  ostringstream aus;
  aus << "0000 MESSAGE    " << number << " equivalent structures are removed!" << Message(__AFLOW_FILE__,aflags) << endl;
  aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
  aus.str("");
  vector<xstructure> vec_xstr_final;
  for (uint i=0;i<vec_xstr_energy_final.size();i++) {
    xstructure xstr_tmp = vec_xstr_energy_final.at(i).xstr;
    string str_title = " DG=" + aurostd::utype2string(vecDG.at(i));
    xstr_tmp.title += str_title;
    vec_xstr_final.push_back(xstr_tmp);  //xstr_energy, structure, which contains xstr, and energy
  }
  return vec_xstr_final;
}

// ***************************************************************************
// pocc::DIFF(xstructure xstr1, xstructure xstr2)
// ***************************************************************************
namespace pocc {
  bool DIFF_UFFENERGY(xstructure xstr1, xstructure xstr2) {
    const double  energy_eps_tol_ =  1E-10;
    double e1 = pocc::CalculateUFFEnergy(xstr1);
    double e2 = pocc::CalculateUFFEnergy(xstr2);
    double e_delta = abs(e1-e2);
    if(e_delta < energy_eps_tol_) {return true;}
    else {return false;}
  }
} // namespace pocc

// ***************************************************************************
// pocc::DIFF(xstructure xstr1, xstructure xstr2)
// ***************************************************************************
namespace pocc {
  bool DIFF(xstructure xstr1, xstructure xstr2) {
    bool RUN_FLAG = false;

    if(pocc::DIFF_LEVEL1(xstr1, xstr2)) { RUN_FLAG = true; return RUN_FLAG;}
    //Highest lelvel
    if(pocc::DIFF_LEVEL4(xstr1, xstr2)) { RUN_FLAG = true; return RUN_FLAG;}
    if(pocc::DIFF_LEVEL4(xstr2, xstr1)) { RUN_FLAG = true; return RUN_FLAG;}
    return RUN_FLAG;
  }
} // namespace pocc

// ***************************************************************************
// pocc::DIFF_LEVEL1(xstructure xstr1, xstructure xstr2)
// ***************************************************************************
namespace pocc {
  bool DIFF_LEVEL1(xstructure xstr1, xstructure xstr2) {
    //pocc::DIFF_Level: 
    //Only compare coordinates, assume names_is_given
    bool RUN_FLAG = false;
    xstr1.SpeciesPutAlphabetic();
    xstr2.SpeciesPutAlphabetic();
    xstr1.ShiftOriginToAtom(0);
    pflow::Sort_atom_cpos(xstr1);
    for (int i=0; i<xstr2.num_each_type.at(0);i++) {
      xstr2.ShiftOriginToAtom(i);
      xstr2.BringInCell();
      pflow::Sort_atom_cpos(xstr2);
      if(CoordCompare(xstr1, xstr2)) {
        RUN_FLAG = true;
        break;
      }
    }
    return RUN_FLAG;
  }
} // namespace pocc

// ***************************************************************************
// pocc::DIFF_LEVEL4(xstructure xstr1, xstructure xstr2)
// ***************************************************************************
namespace pocc {
  bool DIFF_LEVEL4(xstructure xstr1, xstructure xstr2) {
    //pocc::DIFF_Level: 
    //This is the highest level to compare two POSCARs
    bool RUN_FLAG = false;
    xstructure xstr1_sp, xstr2_sp;
    xstr1_sp = GetStandardPrimitive(xstr1);
    xstr2_sp = GetStandardPrimitive(xstr2);
    xstr1_sp.ShiftOriginToAtom(0);
    xstr1_sp.BringInCell();
    pflow::Sort_atom_cpos(xstr1_sp);

    vector<xstructure> vec_rot = pocc::GenerateRotatedXstr(xstr2_sp);
    for (uint i=0; i<vec_rot.size();i++) {
      xstructure xstr2_sp_rot_tmp = vec_rot.at(i);
      for (uint j=0; j<xstr2_sp_rot_tmp.atoms.size();j++) {
        xstr2_sp_rot_tmp.ShiftOriginToAtom(j);
        xstr2_sp.BringInCell();
        pflow::Sort_atom_cpos(xstr2_sp_rot_tmp);
        if(CoordCompare(xstr1_sp, xstr2_sp_rot_tmp)) {
          RUN_FLAG = true;
          break;
        }
      }
    }
    return RUN_FLAG;
  }
} // namespace pocc

// ***************************************************************************
// vector<xstructure> pocc::GenerateRotatedXstr(xstructure xstr)
// ***************************************************************************
namespace pocc {
  vector<xstructure> GenerateRotatedXstr(xstructure xstr) {
    //only including rotations, because origin shift is alreay included in pocc::DIFF
    vector<xstructure> vec_xstr_rotated;
    xstructure xstr_tmp=xstr;

    xstr.CalculateSymmetryPointGroup();

    for(uint i=0; i<xstr.pgroup.size(); i++) {
      for (uint j=0; j<xstr.atoms.size();j++) {
        xstr_tmp.atoms.at(j).fpos = xstr.pgroup.at(i).Uf*xstr.atoms.at(j).fpos;
      }
      vec_xstr_rotated.push_back(xstr_tmp); 
    }
    return vec_xstr_rotated;
  }
} // namespace pocc

// ***************************************************************************
// pocc::DIFF(string options)
// ***************************************************************************
namespace pocc {
  void DIFF(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " options=[" << options << "]" << endl;}
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;}
    if(LDEBUG) {for(uint i=0;i<tokens.size();i++) cerr << __AFLOW_FUNC__ << " tokens.at(i)=" << tokens.at(i) << endl;}

    if(tokens.size()<2) {
      init::ErrorOption(options,__AFLOW_FUNC__,"aflow --diff=POSCAR1,POSCAR2");
    }
    if(!aurostd::FileExist(tokens.at(0))) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"file not found: "+tokens.at(0),_FILE_CORRUPT_);} //CO20200624
    if(!aurostd::FileExist(tokens.at(1))) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"file not found: "+tokens.at(1),_FILE_CORRUPT_);} //CO20200624
    xstructure xstr1, xstr2;
    xstr1 = xstructure(tokens.at(0), IOAFLOW_AUTO);
    xstr2 = xstructure(tokens.at(1), IOAFLOW_AUTO);
    if(tokens.at(0)==tokens.at(1)) {
      cerr << "Do not try to fool me" << endl;
      cerr << "But I will still continue for test!" << endl;
    }
    cerr << "Now running void pocc::DIFF(string options) " << endl;

    if(!pocc::DIFF_UFFENERGY(xstr1, xstr2)) {
      cerr << "Not equivalent!" << endl;
    }
    else if(pocc::DIFF(xstr1, xstr2)) {
      cerr << "Equivalent!" << endl;
    }
    else {
      cerr << "Not equivalent!" << endl;
    }
  }
} // namespace pocc

// **************************************************************************
// xstructure::AddAtom_POCC
// **************************************************************************
// This adds an atom to the structure.
// When this code was written, AddAtom added atoms to the end of the structure, 
// majority of the code intact, this variant of AddAtom was added. 

void xstructure::AddAtom_POCC(const _atom& atom) {
  bool LDEBUG=(FALSE || XHOST.DEBUG); 
  _atom btom;btom=atom;

  // check that this atom is not already present
  xvector<double> a1(3),a2(3),a3(3),aijk(3);
  a1=lattice(1);a2=lattice(2);a3=lattice(3);

  bool FOUND_POSITION=FALSE;
  for(uint iat=0;iat<atoms.size()&&FOUND_POSITION==FALSE;iat++) {
    if(atoms.at(iat).type==atom.type && atoms.at(iat).name==atom.name) {
      for(int i=-1;i<=1&&FOUND_POSITION==FALSE;i++) {
        for(int j=-1;j<=1&&FOUND_POSITION==FALSE;j++) {
          for(int k=-1;k<=1&&FOUND_POSITION==FALSE;k++) {
            aijk[1]=i;aijk[2]=j;aijk[3]=k;
            //	if(aurostd::modulus(atoms.at(iat).cpos-(((double)i)*a1+((double)j)*a2+((double)k)*a3+atom.cpos))<0.1) FOUND_POSITION=TRUE;
            if(aurostd::modulus(atoms.at(iat).fpos-(aijk+atom.fpos))<0.01) {FOUND_POSITION=TRUE;}
          }
        }
      }
    }
  }
  if(FOUND_POSITION==TRUE) {return;} // found no need to add it further

  // now found that it does not exist check type
  //  cerr << "AddAtom new atom" << endl;
  bool FOUND_SPECIES=FALSE;
  uint species_position=0;
  for(uint isp=0;isp<species.size()&&FOUND_SPECIES==FALSE;isp++) {
    if(atom.name==species.at(isp)) {FOUND_SPECIES=TRUE;species_position=isp;}
  }

  if(FOUND_SPECIES==FALSE) {
    if(LDEBUG) {cerr << "AddAtom new_specie=" << atom.name << endl;}
    num_each_type.push_back(1);
    comp_each_type.push_back(atom.partial_occupation_value);
    species.push_back(atom.name); // cerr << "AddAtom=" << atom.name << endl;
    species_pp.push_back(atom.name); // cerr << "AddAtom=" << atom.name << endl;
    species_pp_type.push_back(""); // cerr << "AddAtom=" << atom.name << endl;
    species_pp_version.push_back(""); // cerr << "AddAtom=" << atom.name << endl;
    species_pp_ZVAL.push_back(0.0); // cerr << "AddAtom=" << atom.name << endl;
    species_pp_vLDAU.push_back(deque<double>()); // cerr << "AddAtom=" << atom.name << endl;
    species_volume.push_back(GetAtomVolume(atom.name)); // cerr << "AddAtom=" << atom.name << endl;
    species_mass.push_back(GetAtomMass(atom.name)); // cerr << "AddAtom=" << atom.name << endl;
  } else {
    // cerr << num_each_type.size() << " " <<  btom.type << endl;
    // cerr << comp_each_type.size() << " " <<  btom.type << endl;
    if(LDEBUG) {cerr << "AddAtom increasing species_position " << species_position << endl;}
    num_each_type.at(species_position)++;
    comp_each_type.at(species_position)+=atom.partial_occupation_value;
  }
  if(btom.name_is_given) {
    btom.CleanName();
    btom.CleanSpin();
  }

  // OLD STYLE
  atoms.push_back(btom);  MakeBasis(); return;

  //// NEW STYLE
  //bool found=FALSE;
  //if(0)  for(uint iat=0;iat<atoms.size()&&!found;iat++) {
  //if(iat<atoms.size()-1) {
  //if(atoms.at(iat).type==btom.type && atoms.at(iat+1).type!=btom.type) {
  ////	if(LDEBUG)
  //cerr << "HERE1 iat=" << iat << "  atoms.at(iat).type=" << atoms.at(iat).type << "  btom.type=" << btom.type << endl;//" atoms.begin()=" <<  long(atoms.begin()) << endl;
  //atoms.insert(iat+atoms.begin()+1,btom); // potential problem  with CAST
  //found=TRUE;
  //}
  //}
  //}
  //if(1) {
  //std::deque<_atom>::iterator it=atoms.begin();
  //for(uint iat=0;iat<atoms.size()&&!found;iat++,it++) {
  //if(iat<atoms.size()-1) {
  ////	cerr << "HERE0 iat=" << iat << "  atoms.at(iat).type=" << atoms.at(iat).type << "  btom.type=" << btom.type << endl;
  //if(atoms.at(iat).type==btom.type && atoms.at(iat+1).type!=btom.type) {
  ////	if(LDEBUG)
  ////	  cerr << "HERE1 iat=" << iat << "  atoms.at(iat).type=" << atoms.at(iat).type << "  btom.type=" << btom.type << endl;//" atoms.begin()=" <<  long(atoms.begin()) << endl;
  //atoms.insert(it+1,btom);  // it is iterator, fine for insert.
  //found=TRUE;
  //}
  //}
  //}
  //}
  //// if never found add at the end
  //if(!found) atoms.push_back(btom);
  //// 
  ////  atoms.push_back(btom);
  //// done
  //MakeBasis(); // need to update NUMBER and BASIS
}

// ***************************************************************************
// pocc::POSCAR2ENUM(istream& input)
// ***************************************************************************
namespace pocc {
  void POSCAR2ENUM(istream& input) {

    xstructure xstr_in; //, xstr_out;
    xstr_in.clear(); //DX20191220 - uppercase to lowercase clear
    xstr_in=xstructure(input, IOVASP_POSCAR);
    xstr_in.ReScale(1.00000);

    stringstream oss;
    oss.str("");
    ofstream FileMESSAGE;
    FileMESSAGE.open("LOG.ENUM.INPUT");
    _aflags aflags;

    ostringstream aus;
    aus <<"0000 MESSAGE    Printing input POSCAR " << Message(__AFLOW_FILE__,aflags);
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
    aus << xstr_in;
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);

    pocc::POSCAR2ENUM(xstr_in, oss, FileMESSAGE, aflags);
    aurostd::stringstream2file(oss,"struct_enum.in");
    FileMESSAGE.close();
  }
}

// ***************************************************************************
// pocc::POSCAR2ENUM(xstructure &a, stringstream &oss, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
namespace pocc {
  void POSCAR2ENUM(xstructure &a, stringstream &oss, ofstream &FileMESSAGE, _aflags &aflags) {
    xstructure xstr_in=a;

    string str_species_ori = xstr_in.SpeciesString();
    vector<string> vxstr_species_ori;
    aurostd::string2tokens(str_species_ori, vxstr_species_ori, " ");

    int nHNF = InitializeXstr(xstr_in, vxstr_species_ori, FileMESSAGE, aflags);
    int nMin = nHNF;
    int nMax = nHNF;

    //OptimizeXstr(xstr_in, FileMESSAGE, aflags);

    oss.precision(6);
    oss.setf(ios_base::fixed,ios_base::floatfield);
    oss << xstr_in.title << endl;
    oss << "bulk" << endl;
    oss <<xstr_in.lattice(1,1)<<"  "<<xstr_in.lattice(1,2) <<"  "<<xstr_in.lattice(1,3)<<endl;
    oss <<xstr_in.lattice(2,1)<<"  "<<xstr_in.lattice(2,2) <<"  "<<xstr_in.lattice(2,3)<<endl;
    oss <<xstr_in.lattice(3,1)<<"  "<<xstr_in.lattice(3,2) <<"  "<<xstr_in.lattice(3,3)<<endl;

    vector<vector<int> > NormalisedNumber=NormalisedNumberXstructure(xstr_in);
    int nDFull=NormalisedNumber.size();
    int Ntype;
    if(CheckVacancy(xstr_in)) {Ntype=xstr_in.num_each_type.size()+1;}
    else {Ntype=xstr_in.num_each_type.size();}

    oss << Ntype << " -nary case" << endl;
    oss << nDFull << "  #Number of points in the multilattice" << endl;

    vector<vector<int> > LabelXstr=CalculateLableXstructure(xstr_in);
    for (unsigned int i=0; i<NormalisedNumber.size();i++) {
      int j=NormalisedNumber.at(i).at(0);
      oss << xstr_in.atoms.at(j).cpos[1] << "   " << xstr_in.atoms.at(j).cpos[2] << "   " << xstr_in.atoms.at(j).cpos[3] << "   ";
      for (unsigned int k=0; k<LabelXstr.at(i).size();k++) { 
        oss << LabelXstr.at(i).at(k);
        if(LabelXstr.at(i).size()-k-1) {
          oss << "/";
        }
      }
      oss << endl;
    }

    vector<vector<int> > cRangeFrac=CalculatecRange(xstr_in);
    int NumCell=cRangeFrac.at(0).at(2)/nDFull;
    if(NumCell>nMax) {nMax=NumCell;}
    if(!CheckPartialOccupation(xstr_in)) {nMax=1;}
    oss << nMin << "  " << nMax << "  # Starting and ending cell sizes for search" << endl;
    oss << "1E-6   # Epsilon (finite precision parameter)" << endl;
    oss << "full list of labelings" << endl;
    oss << "# Concentration ranges" << endl;
    for (int i=0; i<Ntype;i++) {
      oss << cRangeFrac.at(i).at(0) << " " << cRangeFrac.at(i).at(1) << " " << cRangeFrac.at(i).at(2) << endl;
    }

    ostringstream aus;
    aus.str("");
    aus <<"0000 MESSAGE    Printing \'struct_enum.in\' file " << Message(__AFLOW_FILE__,aflags);
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
    aus << oss.str();
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
  }
}

// ***************************************************************************
// pocc::POSCAR2ENUM(xstructure &a)
// ***************************************************************************
namespace pocc {
  string POSCAR2ENUM(xstructure &a) {
    stringstream oss;
    oss.str("");
    xstructure xstr_in=a;
    string str_species_ori = xstr_in.SpeciesString();
    vector<string> vxstr_species_ori;
    aurostd::string2tokens(str_species_ori, vxstr_species_ori, " ");

    ofstream FileMESSAGE;
    _aflags aflags;
    int nHNF = InitializeXstr(xstr_in, vxstr_species_ori, FileMESSAGE, aflags);
    int nMin = nHNF;
    int nMax = nHNF;

    //OptimizeXstr(xstr_in, FileMESSAGE, aflags);

    oss.precision(6);
    oss.setf(ios_base::fixed,ios_base::floatfield);
    oss << xstr_in.title << endl;
    oss << "bulk" << endl;
    oss <<xstr_in.lattice(1,1)<<"  "<<xstr_in.lattice(1,2) <<"  "<<xstr_in.lattice(1,3)<<endl;
    oss <<xstr_in.lattice(2,1)<<"  "<<xstr_in.lattice(2,2) <<"  "<<xstr_in.lattice(2,3)<<endl;
    oss <<xstr_in.lattice(3,1)<<"  "<<xstr_in.lattice(3,2) <<"  "<<xstr_in.lattice(3,3)<<endl;

    vector<vector<int> > NormalisedNumber=NormalisedNumberXstructure(xstr_in);
    int nDFull=NormalisedNumber.size();
    int Ntype;
    if(CheckVacancy(xstr_in)) {Ntype=xstr_in.num_each_type.size()+1;}
    else {Ntype=xstr_in.num_each_type.size();}

    oss << Ntype << " -nary case" << endl;
    oss << nDFull << "  #Number of points in the multilattice" << endl;

    vector<vector<int> > LabelXstr=CalculateLableXstructure(xstr_in);
    for (unsigned int i=0; i<NormalisedNumber.size();i++) {
      int j=NormalisedNumber.at(i).at(0);
      oss << xstr_in.atoms.at(j).cpos[1] << "   " << xstr_in.atoms.at(j).cpos[2] << "   " << xstr_in.atoms.at(j).cpos[3] << "   ";
      for (unsigned int k=0; k<LabelXstr.at(i).size();k++) { 
        oss << LabelXstr.at(i).at(k);
        if(LabelXstr.at(i).size()-k-1) {
          oss << "/";
        }
      }
      oss << endl;
    }

    vector<vector<int> > cRangeFrac=CalculatecRange(xstr_in);
    int NumCell=cRangeFrac.at(0).at(2)/nDFull;
    if(NumCell>nMax) {nMax=NumCell;}
    if(!CheckPartialOccupation(xstr_in)) {nMax=1;}
    oss << nMin << "  " << nMax << "  # Starting and ending cell sizes for search" << endl;
    oss << "1E-6   # Epsilon (finite precision parameter)" << endl;
    oss << "full list of labelings" << endl;
    oss << "# Concentration ranges" << endl;
    for (int i=0; i<Ntype;i++) {
      oss << cRangeFrac.at(i).at(0) << " " << cRangeFrac.at(i).at(1) << " " << cRangeFrac.at(i).at(2) << endl;
    }
    return oss.str();
  }
}


// ***************************************************************************
// pocc::MultienumPrintAllXstr(istream& input)
// ***************************************************************************
namespace pocc {
  bool MultienumPrintAllXstr(istream& input) {
    xstructure xstr(input,IOVASP_AUTO);
    xstr.ReScale(1.00000);

    ofstream FileMESSAGE;
    FileMESSAGE.open("LOG.ENUM");
    _aflags aflags;
    ostringstream aus;
    aus << "0000 MESSAGE    Printing input POSCAR " << Message(__AFLOW_FILE__,aflags);
    aus << AFLOWIN_SEPARATION_LINE<< endl;
    aus << xstr;
    aus << AFLOWIN_SEPARATION_LINE<< endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);

    //OptimizeXstr(xstr, FileMESSAGE, aflags);

    stringstream ss;
    stringstream oss;
    vector<xstructure> groupxstr;
#ifdef _AFLOW_GUS_POCC_
    groupxstr=pocc::MultienumGenerateXstr(xstr, FileMESSAGE, aflags);
#endif
    for (uint i=0; i<groupxstr.size();i++) {
      ss.str("");
      ss << setfill('0') << setw(6) <<(i+1);
      oss << AFLOWIN_SEPARATION_LINE<< endl;
      oss << _VASP_POSCAR_MODE_EXPLICIT_START_ << ss.str() << endl;  //CO20200624
      oss << groupxstr.at(i);
      oss << _VASP_POSCAR_MODE_EXPLICIT_STOP_ << ss.str() << endl; //CO20200624
      oss << AFLOWIN_SEPARATION_LINE<< endl;
    }
    aus << "0000 MESSAGE    Printing derivate POSCARs " << Message(__AFLOW_FILE__,aflags);
    aus << oss.str();
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
    FileMESSAGE.close();
    return TRUE;
  }
}

// ***************************************************************************
// vector<xstructure> pocc::MultienumGenerateXstr(xstructure& xstr, ofstream &FileMESSAGE, _aflags &aflags)
// ***************************************************************************
extern "C" {
#ifdef _AFLOW_GUS_POCC_
  void __aflow_call_enum_MOD_aflow_pass_parameter(double *parLV, int *nDFull,double *dFull,int *rdFull,int *cdFull,int *k, int *nMin, int *nMax, char *pLatTyp, const double *eps, bool *full, int *labelFull, int *rlabelFull, int *clabelFull, int *digitFull, int *ndigitFull, int *equivalencies, int *nequivalencies,bool *conc_check,int *cRange, int *rcRange, int *ccRange);
  void __aflow_call_makestr_MOD_makestr(int *strNi, int *strNf);
#endif
}

#ifdef _AFLOW_GUS_POCC_
namespace pocc {
  vector<xstructure> MultienumGenerateXstr(xstructure& xstr, ofstream &FileMESSAGE, _aflags &aflags) {
    const double eps=1E-6;
    vector<string> vectitle;
    string title;
    aurostd::string2tokens(xstr.title, vectitle, " ");
    title=vectitle.at(0);

    string str_species_ori = xstr.SpeciesString();
    vector<string> vxstr_species_ori;
    aurostd::string2tokens(str_species_ori, vxstr_species_ori, " ");

    int nHNF = InitializeXstr(xstr, vxstr_species_ori,FileMESSAGE, aflags);
    int nMin = nHNF;
    int nMax = nHNF;


    ostringstream aus;
    aus << "0000 MESSAGE    Print multienum input file " << Message(__AFLOW_FILE__,aflags);
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aus << pocc::POSCAR2ENUM(xstr);
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    double parLV[3*3];
    parLV[0]=xstr.lattice(1,1); parLV[1]=xstr.lattice(1,2); parLV[2]=xstr.lattice(1,3);
    parLV[3]=xstr.lattice(2,1); parLV[4]=xstr.lattice(2,2); parLV[5]=xstr.lattice(2,3);
    parLV[6]=xstr.lattice(3,1); parLV[7]=xstr.lattice(3,2); parLV[8]=xstr.lattice(3,3);

    vector<vector<int> > NormalisedNumber=NormalisedNumberXstructure(xstr);
    int nDFull=NormalisedNumber.size();

    int Ntype;
    bool VacancyFLAG=CheckVacancy(xstr);
    if(VacancyFLAG) {Ntype=xstr.num_each_type.size()+1;} //Ntype means the number of types of atoms
    else {Ntype=xstr.num_each_type.size();}
    int k=Ntype;

    char pLatTyp[1];
    pLatTyp[0]='B';

    bool full=TRUE;
    bool conc_check=TRUE;

    int nequivalencies=nDFull;
    int equivalencies[nDFull];
    for (int i=0; i<nDFull; i++) {
      equivalencies[i]=i+1;
    }

    int rdFull=3;                  //rdFull means the row number of dFull
    int cdFull=nDFull;             //cdFull means the column number of dFull
    double dFull[rdFull*cdFull];   //dFull is the data array storing the atom sites
    for (unsigned int i=0; i<NormalisedNumber.size();i++) {
      int j=NormalisedNumber.at(i).at(0);
      dFull[3*i+0]=xstr.atoms.at(j).cpos[1];
      dFull[3*i+1]=xstr.atoms.at(j).cpos[2];
      dFull[3*i+2]=xstr.atoms.at(j).cpos[3];
    }

    unsigned int MaxNum=Ntype;
    vector<vector<int> > LabelXstr=CalculateLableXstructure(xstr);
    vector<vector<int> > NewLabelXstr; //This label fills the data array
    vector<int> LabelTmp;
    for (unsigned int i=0; i<LabelXstr.size();i++) {
      LabelTmp.clear();
      for (unsigned int j=0; j<LabelXstr.at(i).size();j++) {
        LabelTmp.push_back(LabelXstr.at(i).at(j));
      }
      if(LabelXstr.at(i).size()<MaxNum) {
        for (unsigned int k=0; k<(MaxNum-LabelXstr.at(i).size());k++) {
          LabelTmp.push_back(-1);
        }
      }
      NewLabelXstr.push_back(LabelTmp);
    }
    //-----------------------------------------------------------------------
    //Convert the multi dimensional array into one dimensional
    vector<int> OneDimensionLabelXstr;
    for (unsigned int i=0; i<NewLabelXstr.size();i++) {
      for (unsigned int j=0; j<NewLabelXstr.at(i).size();j++) {
        OneDimensionLabelXstr.push_back(NewLabelXstr.at(i).at(j));
      }
    }
    //-----------------------------------------------------------------------
    int clabelFull=nDFull;    
    int rlabelFull=MaxNum;         //rlabelFull means number of columns, different from c++
    int LabelSize=clabelFull*rlabelFull;
    int labelFull[LabelSize];
    for (int i=0;i<LabelSize;i++) {
      labelFull[i]=OneDimensionLabelXstr.at(i);
    }

    int ndigitFull=nDFull;
    int digitFull[nDFull];
    for(int i=0; i<nDFull; i++) {
      digitFull[i]=LabelXstr.at(i).size();
    }

    //Note the format of cRange is different from that of dFull, if you take a look at its fortran code
    vector<vector<int> > cRangeFrac=CalculatecRange(xstr);
    int ccRange=3;
    int rcRange=Ntype;
    int cRange[rcRange*ccRange];
    int cRange2[rcRange*ccRange];
    for (int i=0; i<Ntype;i++) {
      cRange2[i*3+0]=cRangeFrac.at(i).at(0);
      cRange2[i*3+1]=cRangeFrac.at(i).at(1);
      cRange2[i*3+2]=cRangeFrac.at(i).at(2);
    }
    for( int i=0; i<3; i++) {
      for (int j=0; j<Ntype; j++) {
        cRange[i*Ntype+j]=cRange2[j*3+i];
      }
    }
    int NumCell=cRangeFrac.at(0).at(2)/nDFull;
    if(NumCell>nMax) {nMax=NumCell;}
    if(!CheckPartialOccupation(xstr)) {nMax=1;}
    //[CO20191112 - OBSOLETE]char *cur_dir_name = getcwd(NULL, 0);
    string cur_dir_name=aurostd::getPWD();
    string tmpdir=aurostd::TmpDirectoryCreate("PartialOccupation");
    //[CO20191112 - OBSOLETE]char new_dir_name[1024];
    //[CO20191112 - OBSOLETE]strcpy(new_dir_name, tmpdir.c_str());
    //[CO20191112 - OBSOLETE]chdir(new_dir_name);
    chdir(tmpdir.c_str());  //CO20191112

    __aflow_call_enum_MOD_aflow_pass_parameter(parLV, &nDFull,dFull,&rdFull,&cdFull,&k, &nMin, &nMax, pLatTyp,&eps, &full, labelFull,&rlabelFull,&clabelFull,digitFull,&ndigitFull,equivalencies,&nequivalencies,&conc_check,cRange,&rcRange,&ccRange);

    //debug
    //--------------------------------------------------------------------------------
    //Call makestr
    stringstream str_out;
    string filename_str_out="struct_enum.out";
    if(!aurostd::FileExist(filename_str_out)) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"file not found: struct_enum.out",_FILE_CORRUPT_); //CO20200624
    }
    string line,lastline;
    file2stringstream(filename_str_out, str_out);
    while(!str_out.eof()) {
      line.clear();
      getline(str_out,line);
      if(line.size()>0) {
        lastline=line;
      }
    }
    int strNi;
    int strNf;
    vector<string> vec_str_out;
    aurostd::string2tokens(lastline, vec_str_out, " ");
    if(is_number(vec_str_out.at(0))) {
      strNi=1;
      strNf=aurostd::string2utype<int>(vec_str_out.at(0));
    }
    else {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"multienum did not generate a valid struct_enum.out file",_FILE_CORRUPT_); //CO20200624
    }

    ////output
    aus << "0000 MESSAGE    Print multienum output file " << Message(__AFLOW_FILE__,aflags);
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    aus << str_out.str();
    aus <<AFLOWIN_SEPARATION_LINE << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    //output

    __aflow_call_makestr_MOD_makestr(&strNi, &strNf);

    stringstream ss;
    xstructure xstr_cleaned;
    vector<xstructure> groupxstr, groupxstr_names;
    for (int i=1; i<=strNf;i++) {
      ss.str("");
      ss << setfill('0') << setw(6) << i;
      string strnumber=ss.str();

      xstructure xstr_partocc="xstr"+strnumber;
      string vaspname="vasp."+strnumber;
      xstr_partocc=xstructure(vaspname, IOVASP_POSCAR);
      string_replaceAll(xstr_partocc.title,"MULTIENUM", title); //Replace the name of POSCAR using the actual title
      //if xstr has vacancy, then we need to clean them
      if(VacancyFLAG) {
        xstr_cleaned=CleanVacancy(xstr_partocc);
        groupxstr.push_back(xstr_cleaned);
      }
      else {
        groupxstr.push_back(xstr_partocc);
      }
    }
    //To make sure safely remove the directory, we divide it into two steps
    //Clean all the files in the temporary directory
    stringstream ss_cmd;
    ss_cmd << "rm -rf " << tmpdir << endl;
    aurostd::execute(ss_cmd);

    chdir(cur_dir_name.c_str());  //CO20191112
    return groupxstr;
  }
}
#endif

// ***************************************************************************
// pocc::MULTIENUM(vector<string> argv,istream& input)
// ***************************************************************************
// pocc::MultienumFilterXstr(vector<xstructure> groupxstr, vector<string> AtomSpecies) 
namespace pocc {
  bool MultienumPrintSortedXstr(istream& input) {
    xstructure xstr; 
    xstr.clear(); //DX20191220 - uppercase to lowercase clear
    xstr=xstructure(input, IOVASP_POSCAR);

    vector<string> AtomSpecies;
    aurostd::string2tokens(xstr.SpeciesString(), AtomSpecies, " ");

    if(!xstr.atoms.at(0).name_is_given) {
      throw aurostd::xerror(__AFLOW_FILE__,"pocc::MultienumPrintSortedXstr():","atom species are not assigned",_INPUT_ILLEGAL_); //CO20200624
    }

    ofstream FileMESSAGE;
    FileMESSAGE.open("LOG.ENUM.SORTED");
    _aflags aflags;
    ostringstream aus;
    OptimizeXstr(xstr, FileMESSAGE, aflags);

    vector<xstructure> groupxstr;
#ifdef _AFLOW_GUS_POCC_
    groupxstr=pocc::MultienumGenerateXstr(xstr, FileMESSAGE, aflags);
#endif

    vector<xstructure> sortedgroupxstr;
    //sortedgroupxstr=pocc::SortGroupXstr(groupxstr, AtomSpecies);

    vector<xstructure> groupxstr_names = AssignNameXstr(groupxstr, AtomSpecies);
    sortedgroupxstr=pocc::SortGroupXstrUFFEnergy(groupxstr_names);

    //vector<xstructure> vxstr_final = sortedgroupxstr;
    vector<xstructure> vxstr_final = RemoveEquivalentXstr(sortedgroupxstr, FileMESSAGE, aflags);
    //Format xstructure
    vector<xstructure> vxstr_final_alphabetic;
    for (uint i=0; i<vxstr_final.size();i++) {
      xstructure xstr_tmp = vxstr_final.at(i);
      vector<string> vxstr_tmp;
      aurostd::string2tokens(xstr_tmp.SpeciesString(), vxstr_tmp, " ");
      xstructure xstr2 = AssignNameXstr(xstr_tmp, vxstr_tmp);
      xstr2.SpeciesPutAlphabetic();
      vxstr_final_alphabetic.push_back(xstr2);
    }

    stringstream ss;
    stringstream oss;
    for (uint i=0; i<vxstr_final_alphabetic.size();i++) {
      ss.str("");
      ss << setfill('0') << setw(6) <<(i+1);
      oss << AFLOWIN_SEPARATION_LINE<< endl;
      oss << _VASP_POSCAR_MODE_EXPLICIT_START_ << ss.str() << endl;  //CO20200624
      oss << vxstr_final_alphabetic.at(i);
      oss << _VASP_POSCAR_MODE_EXPLICIT_STOP_ << ss.str() << endl; //CO20200624
      oss << AFLOWIN_SEPARATION_LINE<< endl;
    }
    aus << "0000 MESSAGE    Printing sorted derivate POSCARs " << Message(__AFLOW_FILE__,aflags);
    aus << oss.str();
    aurostd::PrintMessageStream(FileMESSAGE, aus,XHOST.QUIET);
    /*
    //This test proves that hfncell program can exactly give same number supercells with multienum
    //if they are different, this is due to the energy tollerance in gulp.
    cout << "total number " << vxstr_final.size() << endl;
    */
    FileMESSAGE.close();

    return TRUE;
  }
}
// ***************************************************************************

// ***************************************************************************
// void pocc::POCC_SetOptions(string options, string& directory, double& T, double& DOS_Emin, double& DOS_Emax, double& DOSSCALE)
// ***************************************************************************
namespace pocc {
  void POCC_SetOptions(string options, string& directory, double& T, double& DOS_Emin, double& DOS_Emax, double& DOSSCALE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) {
      cerr << "pocc::POCC_SetOptions: BEGIN" << endl;
    }
    directory="./"; DOS_Emin=DEFAULT_DOS_EMIN/2.0; DOS_Emax=DEFAULT_DOS_EMAX/2.0; DOSSCALE=DEFAULT_DOS_SCALE; // some defaults  
    T = 300.0;  //300K room temperature default
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");

    if(LDEBUG) {
      cerr << "pocc::POCC_SetOptions: options=[" << options << "]" << endl;
    }
    if(LDEBUG) {
      cerr << "pocc::POCC_SetOptions: tokens.size()=" << tokens.size() << endl;
    }
    if(LDEBUG) {
      for(uint i=0;i<tokens.size();i++) cerr << "pocc::POCC_SetOptions: tokens.at(i)=" << tokens.at(i) << endl;
    }

    //Usage: program 0; --pocc_dos 1; directory 2; temperature 3, DOS_Emin 4; DOS_Emax 5; DOSSCALE 6

    if(tokens.size()>=1) {
      directory = tokens.at(0); 
    }
    if(tokens.size()>=2) {
      T = aurostd::string2utype<double>(tokens.at(1)); 
    }
    if(tokens.size()>=3) {
      DOS_Emin = aurostd::string2utype<double>(tokens.at(2)); 
    }
    if(tokens.size()>=4) {
      DOS_Emax = aurostd::string2utype<double>(tokens.at(3));
    }
    if(tokens.size()>=5) {
      DOSSCALE = aurostd::string2utype<double>(tokens.at(4));
    }

    // [OBSOLETE] if(directory=="") directory="./";
    // [OBSOLETE] if(T<0.01) T=0.01;
    // [OBSOLETE] if(DOS_Emin<=-999.9) DOS_Emin=DEFAULT_DOS_EMIN;
    // [OBSOLETE] if(DOS_Emax<=-999.9) DOS_Emax=DEFAULT_DOS_EMAX;
    // [OBSOLETE] if(DOSSCALE<=-999.9) DOSSCALE=DEFAULT_DOS_SCALE;

    if(LDEBUG) {
      cerr << "pocc::POCC_SetOptions: directory=[" << directory << "]" << endl;
    }
    if(LDEBUG) {
      cerr << "pocc::POCC_SetOptions: T=" << T << endl;
    }
    if(LDEBUG) {
      cerr << "pocc::POCC_SetOptions: DOS_Emin=" << DOS_Emin << endl;
    }
    if(LDEBUG) {
      cerr << "pocc::POCC_SetOptions: DOS_Emax=" << DOS_Emax << endl;
    }
    if(LDEBUG) {
      cerr << "pocc::POCC_SetOptions: DOSSCALE=" << DOSSCALE << endl;
    }

    if(LDEBUG) {
      cerr << "pocc::POCC_SetOptions: END" << endl;
    }
  }
} // namespace pocc

// ***************************************************************************
// pocc::POCC_DOS
// ***************************************************************************
namespace pocc {
  void POCC_DOS(ostream &oss,string options) {
    string directory;
    double T, DOS_Emin, DOS_Emax, DOSSCALE; 
    pocc::POCC_SetOptions(options,directory,T,DOS_Emin,DOS_Emax,DOSSCALE);
    oss << "pocc::POCC_DOS is working on " << directory  << endl;
    oss << "Temperature is " << T << "K" << endl;
    pocc::POCC_GENERATE_OUTPUT(directory,T,DOS_Emin,DOS_Emax,DOSSCALE);
  }
} // namespace pocc

// ***************************************************************************
// pocc::POCC_Already_Calculated_Input(const string& str_AflowIn)
// ***************************************************************************
namespace pocc {
  bool POCC_Already_Calculated_Input(const string& str_AflowIn,const string& dir) {
    return aurostd::substring2bool(str_AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]START.POCC") || aurostd::EFileExist(dir+"/"+POCC_FILE_PREFIX+POCC_ALL_SUPERCELLS_FILE);
  }
} // namespace pocc

// ***************************************************************************
// void pocc::POCC_CalculateDistribution(vector<double>& vdelta_toten, vector<double>& vprob, const double& T, const string& NameDist, const vector<int>& vDE)
// ***************************************************************************
namespace pocc {
  void POCC_CalculateDistribution(vector<double>& vdelta_toten, vector<double>& vprob, const double& T, const string& NameDist, const vector<int>& vDE) {
    //Check 
    if(vdelta_toten.size()!=vDE.size()) {
      cerr << "Error!!! Check the number of POSCARs, Degeneracy!" << endl;
    }
    vector<double> vFi; vprob.clear();
    double Fi, DEi;
    for (uint i=0; i<vdelta_toten.size();i++) {
      double delta_toten = vdelta_toten.at(i);
      DEi = vDE.at(i);
      if(NameDist=="B")  {
        Fi = DEi*exp((-1.0)*delta_toten/(KBOLTZEV*T));    // Boltzmann distribution
      }
      if(NameDist=="FD") {
        Fi = DEi*1.0/(exp(delta_toten/(KBOLTZEV*T))+1.0); // Fermi-Dirac distribution
      }
      if(NameDist=="BE") {
        Fi = DEi*1.0/(exp(delta_toten/(KBOLTZEV*T))-1.0); // Bose-Einstein distribution
      }
      vFi.push_back(Fi);
    }
    double sum_vFi = aurostd::sum<double> (vFi);
    for (uint i=0; i<vFi.size();i++) {
      double probi = vFi.at(i)/sum_vFi;
      vprob.push_back(probi);
    }
  }
} // namespace pocc

// ***************************************************************************
// vector<vector<double> > pocc::SpinSplitDOS(const vector<vector<double> >& vva)
// ***************************************************************************
namespace pocc {
  vector<vector<double> > SpinFlipDOS(const vector<vector<double> >& vva) {
    vector<vector<double> > vvb; vvb.clear();
    if(vva.at(0).size()%2==0) {throw aurostd::xerror(__AFLOW_FILE__,"pocc::SpinFlipDOS():","input not spin-polarized",_INPUT_ILLEGAL_);} //CO20200624
    for (uint i=0; i<vva.size();i++) {
      vector<double> vtmp; vtmp.clear();
      vtmp.push_back(vva.at(i).at(0));
      for (uint j=1; j<vva.at(0).size();j=j+2) {
        vtmp.push_back(vva.at(i).at(j+1)*(-1));   //exchange columns 1 and 2 spin_dn becomes spin_up
        vtmp.push_back(vva.at(i).at(j)*(-1));   //spin_up becomes spin_dn
      }
      vvb.push_back(vtmp);
    }
    return vvb;
  }
} // namespace pocc

// ***************************************************************************
// vector<vector<double> > pocc::SpinSplitDOS(const vector<vector<double> >& vva)
// ***************************************************************************
namespace pocc {
  vector<vector<double> > SpinSplitDOS(const vector<vector<double> >& vva) {
    vector<vector<double> > vvb;
    for (uint i=0; i<vva.size();i++) {
      vector<double> vtmp; vtmp.clear();
      vtmp.push_back(vva.at(i).at(0));
      for (uint j=1; j<vva.at(i).size();j++) {
        double tmp = 0.5*vva.at(i).at(j);  
        vtmp.push_back(tmp); vtmp.push_back((-1)*tmp);
      }
      vvb.push_back(vtmp);
    }
    return vvb;
  }
} // namespace pocc


// ***************************************************************************
// vector<vector<double> > pocc::POCC_Formalise(SPIN_FLAG, vweight, TOTALPDOS)
// ***************************************************************************
namespace pocc {
  vector<vector<double> > POCC_Formalise(const bool& SPIN_FLAG, const double& weight, const double& mag, const vector<vector<double> >& TDOS) {
    vector<vector<double> > vva = TDOS;
    vector<vector<double> > vvb;
    //weight
    for (uint i=0; i<vva.size();i++) {
      for (uint j=1; j<vva.at(i).size();j++) {
        vva.at(i).at(j) *= weight;
      }
    }
    //spin
    if(SPIN_FLAG) {
      if(TDOS.at(0).size()==3 || TDOS.at(0).size()==4) { //TDOS or TOTALPDOS, E, s, p, d
        vvb = pocc::SpinSplitDOS(vva);
      } else if(TDOS.at(0).size()==5 && abs(mag) < 1E-6) { //must be TOTALPDOS: E, s, p, d, f, no-spin-polarized
        vvb = pocc::SpinSplitDOS(vva);
      } else if((TDOS.at(0).size()==5) && ((-1)*mag > 1E-3 )) { //if magnetic moment is smaller than -0.001 and size is 5,  must be TDOS, 
        vvb = SpinFlipDOS(vva);
      } else if((TDOS.at(0).size()==7) && ((-1)*mag > 1E-3 )) { //if magnetic moment is smaller than -0.001, size is 7, must be PDOS
        vvb = SpinFlipDOS(vva);
      } else if((TDOS.at(0).size()==9) && ((-1)*mag > 1E-3 )) { //if magnetic moment is smaller than -0.001, size is 7, must be PDOS
        vvb = SpinFlipDOS(vva);
      } else {
        vvb = vva;
      }
    } else {
      vvb = vva;
    }
    return vvb;
  }
} // namespace pocc

// ***************************************************************************
//void ExtracAllPOSCARSFromAflowin(vector<xstructure>& vxstr, const string& str_aflowin)
// ***************************************************************************
void ExtracAllPOSCARSFromAflowin(vector<xstructure>& vxstr, const string& str_aflowin) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  vxstr.clear();
  string POSCAR_START_DELIMITER=_VASP_POSCAR_MODE_EXPLICIT_START_; //CO20200624
  string POSCAR_STOP_DELIMITER=_VASP_POSCAR_MODE_EXPLICIT_STOP_; //CO20200624
  vector<string> vKBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING;
  aurostd::substring2strings(str_aflowin,vKBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING,POSCAR_START_DELIMITER);
  // load up the structures
  for(uint i=0;i<vKBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++) {
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vKBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING[i]=" << endl;cerr << vKBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING[i] << endl;}
    string START="[VASP_POSCAR_MODE_EXPLICIT]START";
    string STOP="[VASP_POSCAR_MODE_EXPLICIT]STOP";
    START=POSCAR_START_DELIMITER+vKBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING[i];
    STOP=POSCAR_STOP_DELIMITER+vKBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING[i];
    stringstream POSCAR;
    POSCAR.str(aurostd::substring2string(str_aflowin,START,STOP,-1));
    //[SD20220520 - OBSOLETE]if(aurostd::substring2bool(str_aflowin,START) && aurostd::substring2bool(str_aflowin,STOP)) {
    //[SD20220520 - OBSOLETE]  aurostd::ExtractToStringstreamEXPLICIT(str_aflowin,POSCAR,START,STOP);
    //[SD20220520 - OBSOLETE]}
    if(!POSCAR.str().empty()) {vxstr.push_back(xstructure(POSCAR,IOVASP_AUTO));}
  }
}

// ***************************************************************************
// void GetDegeneracyFromVPOSCAR(const vector<xstructure>& vxstr, vector<int>& vDE)
// ***************************************************************************
void GetDegeneracyFromVPOSCAR(const vector<xstructure>& vxstr, vector<int>& vDE) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string title="";
  vector<string> vtitle,vtitle2;
  int DEI=1;
  for (uint i=0; i<vxstr.size();i++) {
    title = vxstr.at(i).title;
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " title=\"" << title << "\"" << endl;}
    if(!aurostd::substring2bool(title, "DG=")) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"degeneracy data not found: "+_AFLOWIN_,_FILE_CORRUPT_); //CO20200624
    }
    //CO20180220 - multiply occupied sites will yield titles with more than one DG= value
    vtitle.clear();
    aurostd::string2tokensAdd(title,vtitle," ");
    DEI=1;
    for(uint j=0;j<vtitle.size();j++){
      if(!aurostd::substring2bool(vtitle[j],"DG=")){continue;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " found DG tag: \"" << vtitle[j] << "\"" << endl;}
      vtitle2.clear();
      aurostd::string2tokens(vtitle[j],vtitle2,"=");
      if(vtitle2.size()!=2){
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"cannot read degeneracy value: "+_AFLOWIN_,_FILE_CORRUPT_); //CO20200624
      }
      DEI*=aurostd::string2utype<int>(vtitle2[1]);
    }
    if(LDEBUG){cerr << " DEI[i=" << i << "]=" << DEI << endl;}
    //[OBSOLETE CO20180220]string last_part_title = vtitle.at(vtitle.size()-1);
    //[OBSOLETE CO20180220]vector<string> vtitle2;
    //[OBSOLETE CO20180220]aurostd::string2tokensAdd(last_part_title,vtitle2,"=");
    //[OBSOLETE CO20180220]int DEI;
    //[OBSOLETE CO20180220]DEI=aurostd::string2utype<int>(vtitle2.at(vtitle2.size()-1));
    vDE.push_back(DEI);
  }
}

// ***************************************************************************
// void pocc::POCC_GENERATE_DOSDATA(const string& directory, const double& T, vector<vector<double> >& DOS, vector<double>& POCC_Efermi, double& POCC_vmag, vector<double>& vprob)
// ***************************************************************************
namespace pocc {
  void POCC_GENERATE_DOSDATA(const string& directory, const double& T, vector<vector<double> >& TDOS_ONLY, vector<vector<double> >& PDOS_ONLY, vector<double>& POCC_Efermi, double& POCC_mag, double& Egap_net, vector<double>& Egap, vector<double>& vprob) {
    //Warnning: DOS is absolute value, no shift, and the output DOS, it is format is Energy, s, p, d, f, TDOS, TDOS_sum

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " directory=[" << directory << "]" << endl;}

    string aflowin,MESSAGE="pocc::POCC_DOS ERROR";
    aflowin=string(directory +"/"+_AFLOWIN_);
    if(!aurostd::FileExist(aflowin)) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"file not found: "+aflowin,_FILE_CORRUPT_);} //CO20200624
    string str_AflowIn; aurostd::file2string(aflowin, str_AflowIn);
    bool pocc_already_calculated_input=pocc::POCC_Already_Calculated_Input(str_AflowIn,directory);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " pocc_already_calculated_input=" << pocc_already_calculated_input << endl;}
    if(pocc_already_calculated_input) {
      //Fix Degeneracy, Reading all POSCARs from aflowin
      vector<xstructure> vxstr;
      if(aurostd::EFileExist(directory+"/"+POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE)){ //read from unique structures file
        string unique_structures_file_contents="";
        aurostd::efile2string(directory+"/"+POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE,unique_structures_file_contents);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " unique_structures_file_contents=" << endl;cerr << unique_structures_file_contents << endl;}
        ExtracAllPOSCARSFromAflowin(vxstr,unique_structures_file_contents);
      }else{
        ExtracAllPOSCARSFromAflowin(vxstr,str_AflowIn);
      }
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vxstr.size()=" << vxstr.size() << endl;}
      vector<int> vDE; vDE.clear();
      GetDegeneracyFromVPOSCAR(vxstr, vDE);

      //Fixing Degeneracy can be very easy if we begin from below, but I do not want to repeat calculations, so I fixt it above
      vector<string> vrun, vdoscar_files, voutcar_files;
      string command, aflow_pocc_out;
      command = "grep \"\\[VASP_POSCAR_MODE_EXPLICIT\\]START\\.\" " + aflowin;
      aflow_pocc_out = aurostd::execute2string(command);
      aurostd::StringSubst(aflow_pocc_out,_VASP_POSCAR_MODE_EXPLICIT_START_,""); //replace string  //CO20200624
      aurostd::string2vectorstring(aflow_pocc_out,vrun);

      //Check files
      bool FLAG_ALLFILES_EXIST = true;
      ostringstream aus;
      ofstream FileMESSAGE;
      //double kpt_tol;

      bool found_ext=false;
      string outcar_file,doscar_file;

      for(uint i=0;i<vrun.size();i++) {
        found_ext=false;
        for(uint iext=0;iext<XHOST.vext.size();iext++){
          outcar_file = aurostd::CleanFileName(directory + "/ARUN."+vrun.at(i)+"/OUTCAR.static"+XHOST.vext.at(iext));
          if(aurostd::FileExist(outcar_file)) {
            if(LDEBUG) {
              aus << "00000  MESSAGE POCC OUTCAR file OK: " << outcar_file << " " << endl;
            }
            if(LDEBUG) {
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            voutcar_files.push_back(outcar_file);  // worse
            found_ext=true;
            break;
          }
        }
        if(!found_ext){
          aus << "ERROR" << ": file not found " << aurostd::CleanFileName(directory + "/ARUN."+vrun.at(i)+"/OUTCAR.static") << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          FLAG_ALLFILES_EXIST = false;
        }
      }
      for(uint i=0;i<vrun.size();i++) {
        found_ext=false;
        for(uint iext=0;iext<XHOST.vext.size();iext++){
          doscar_file = aurostd::CleanFileName(directory + "/ARUN."+vrun.at(i)+"/DOSCAR.static"+XHOST.vext.at(iext));
          if(aurostd::FileExist(doscar_file)) {
            if(LDEBUG) {
              aus << "00000  MESSAGE POCC DOSCAR file OK: " << doscar_file << " " << endl;
            }
            if(LDEBUG) {
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            vdoscar_files.push_back(doscar_file);
            found_ext=true;
            break;
          }
        }
        if(!found_ext){
          aus << "ERROR" << ": file not found " << aurostd::CleanFileName(directory + "/ARUN."+vrun.at(i)+"/DOSCAR.static") << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          FLAG_ALLFILES_EXIST = false;
        }
      }
      if(!FLAG_ALLFILES_EXIST) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"input files not found",_FILE_CORRUPT_); //CO20200624
      }

      //ions, total energy, spin flag, magnetic moment
      vector<int> vions; vector<double> vtoten_per_atom; vector<double> vmag; vector<double> vEgap_net;
      vector<vector<double> > vEgap;
      bool POCC_SPIN_FLAG = false;  //default non-spin-polarized
      for (uint i=0; i<vrun.size();i++) {
        if(LDEBUG) {
          cerr << "pocc::POCC_GENERATE_DOSDATA: vrun.at(" << i << ")=" << vrun.at(i) << endl;
        }
        //xstructure xstr=vxstr[i]; //CO20171002 - grab xstr's too
        //CO20171002 - using tolerance from symmetry calc - START
        //if(xstr.CalculateSymmetry()){kpt_tol=xstr.sym_eps;}
        //else {kpt_tol=SYM::defaultTolerance(xstr);}
        //CO20171002 - using tolerance from symmetry calc - STOP

        xOUTCAR outcar_aus;
        if(!outcar_aus.GetPropertiesFile(voutcar_files.at(i))){
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, "xOUTCAR::GetPropertiesFile() failed", _RUNTIME_ERROR_); //CO20200404
        }
        double EFERMI=outcar_aus.Efermi;
        outcar_aus.GetBandGap();
        vions.push_back(outcar_aus.NIONS);
        vtoten_per_atom.push_back(outcar_aus.enthalpy_atom);
        vmag.push_back(outcar_aus.mag_cell);
        string outcar_bands=voutcar_files.at(i);aurostd::StringSubst(outcar_bands,"static","bands");
        if(aurostd::FileExist(outcar_bands)) {
          if(LDEBUG) {
            cerr << "pocc::POCC_GENERATE_DOSDATA: loading outcar=" << outcar_bands << endl;
          }
          outcar_aus.GetPropertiesFile(outcar_bands);
          outcar_aus.GetBandGap(EFERMI);
        } 
        vEgap_net.push_back(outcar_aus.Egap_net);
        //	cerr << "CAMILO ??? outcar_aus.Egap.size()=" << outcar_aus.Egap.size() << endl;
        vEgap.push_back(outcar_aus.Egap);

        if(LDEBUG||1) {
          cerr << "pocc::POCC_GENERATE_DOSDATA: vrun.at(" << i << ")=" << vrun.at(i) << "; DE=" << vDE[i] << "; NIONS=" << vions.at(i) << "; enthalpy_atom=" << vtoten_per_atom.at(i) << "; mag_cell=" << vmag.at(i) << "; Egap_net=" << vEgap_net.at(i) << endl;    //CO - should read out enthalpy/atom from STATIC, not BANDS
        }

        //	if(LDEBUG||1) {
        //	cerr << "pocc::POCC_GENERATE_DOSDATA: vrun.at(" << i << ")=" << vrun.at(i) << "; NIONS=" << outcar_aus.NIONS << "; enthalpy_atom=" << outcar_aus.enthalpy_atom << "; mag_cell=" << outcar_aus.mag_cell << "; Egap_net=" << outcar_aus.Egap_net << endl;   //CO
        //	}
      }

      int min_IONS = min(vions); 
      vector<double> vweight; //different ions in each structure
      for (uint i=0; i<vrun.size();i++) {   // will use this for DOS as it is extensive
        vweight.push_back(double(1.0*min_IONS/vions.at(i))); //1.0 make int to double, otherwise it is zero!
      }

      double min_toten = min(vtoten_per_atom);
      vector<double> vdelta_toten;
      for (uint i=0; i<vtoten_per_atom.size();i++) {
        vdelta_toten.push_back(vtoten_per_atom.at(i) - min_toten);
      }

      // Calculate the probability boltzmann distribution
      vprob.clear(); 
      pocc::POCC_CalculateDistribution(vdelta_toten, vprob, T, "B", vDE); //Boltzmann

      if(0) { 
        cerr << "check vprob: ";
        for (uint i=0; i<vprob.size();i++) {
          cerr << vprob.at(i) << " ";
        }
        cerr << " = " << aurostd::sum(vprob) << endl;
      }

      // vector<double> vmag_abs; 
      //  for (uint i=0; i<vmag.size();i++) {
      //  vmag_abs.push_back(abs(vmag.at(i))); //ferromagnetic alignment
      //  }

      POCC_mag=0.0;
      for (uint i=0;i<vmag.size();i++) {
        POCC_mag+=abs(vmag.at(i))*vprob.at(i); //Calculate average absolute magnetic moment
      }

      Egap_net=0.0;
      for (uint i=0;i<vEgap_net.size();i++) {
        Egap_net+=vEgap_net.at(i)*vprob.at(i); //Calculate average Egap_net
      }

      Egap.clear();
      for(uint j=0;j<vEgap.at(0).size();j++) {
        Egap.push_back(0);
        for (uint i=0;i<vEgap.size();i++) {
          //Calculate average Egap
          Egap.at(j)+=vEgap.at(0).at(j)*vprob.at(i);    //CO - otherwise, it crashes if you have Egap up AND Egap down (magnetized system)
          //Egap.at(j)+=vEgap.at(i).at(j)*vprob.at(i);  //CO
        }
      }

      // Get TDOS & TOTALPDOS DATA
      vector<vector<vector<double> > > POCC_TDOS, POCC_TOTALPDOS;
      for (uint i=0; i<vdoscar_files.size();i++) {
        string doscar_file = vdoscar_files.at(i); stringstream ss_doscar; aurostd::efile2stringstream(doscar_file, ss_doscar); 
        string outcar_file = voutcar_files.at(i); stringstream ss_outcar; aurostd::efile2stringstream(outcar_file, ss_outcar);
        double Efermi; vector<vector<double> > TDOS, TOTALPDOS;
        //CO20180218 - let's not mess around with KY's functions too much
        //for now, assume POCC runs have standard DOSCAR.static (PDOS in it)
        //if no PDOS, exit
        //I will fix later
        if(!(estructure::GET_DOS_DATA(ss_doscar, ss_outcar, Efermi, TDOS, TOTALPDOS) && TOTALPDOS.size()>0)){
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no PDOS found in DOSCAR",_FILE_CORRUPT_); //CO20200624
        }
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " TDOS.size()=" << TDOS.size() << endl;
          for(uint i=0;i<TDOS.size();i++){
            cerr << __AFLOW_FUNC__ << " TDOS[i=" << i << "].size()=" << TDOS[i].size() << endl;
          }
        }
        //format TDOS, if spin and non-spin coexist, then format them into spin
        vector<vector<double> > TDOSf, TOTALPDOSf;
        TDOSf = pocc::POCC_Formalise(POCC_SPIN_FLAG, vweight.at(i), vmag.at(i), TDOS); 
        TOTALPDOSf = pocc::POCC_Formalise(POCC_SPIN_FLAG, vweight.at(i), vmag.at(i), TOTALPDOS);
        POCC_Efermi.push_back(Efermi);
        POCC_TDOS.push_back(aurostd::ShiftFirstColumn(TDOSf, -1*Efermi)); //conver it into absolute value
        POCC_TOTALPDOS.push_back(aurostd::ShiftFirstColumn(TOTALPDOSf, -1*Efermi));
      }

      vector<vector<double> > POCC_TDOS_normalized = aurostd::NormalizeAndSum3DVector(POCC_TDOS, vprob);  //normalize TDOS
      vector<vector<double> > POCC_TOTALPDOS_normalized = aurostd::NormalizeAndSum3DVector(POCC_TOTALPDOS, vprob); //normalize pdos
      TDOS_ONLY = POCC_TDOS_normalized;
      PDOS_ONLY = POCC_TOTALPDOS_normalized;
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " POCC_TDOS_normalized.size()=" << POCC_TDOS_normalized.size() << endl;
        for(uint i=0;i<POCC_TDOS_normalized.size();i++){
          cerr << __AFLOW_FUNC__ << " POCC_TDOS_normalized[i=" << i << "].size()=" << POCC_TDOS_normalized[i].size() << endl;
        }
      }
    }
  }
} // namespace pocc

// ***************************************************************************
// void pocc::POCC_COMBINE_TDOS_PDOS_ONEDOS(const vector<vector<double> >& TDOS, const vector<vector<double> >& PDOS, vector<vector<double> >& DOS)
// ***************************************************************************
namespace pocc {
  void POCC_COMBINE_TDOS_PDOS_ONEDOS(const vector<vector<double> >& TDOS, const vector<vector<double> >& PDOS, vector<vector<double> >& DOS, vector<vector<double> >& DOS_IDOS) {   //CO20190808 - one without IDOS and one with
    if(TDOS.size()!=PDOS.size()) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"TDOS and PDOS have different size",_FILE_CORRUPT_); //CO20200624
    }
    vector<double> vtmp;
    if(TDOS.at(0).size()==3) { //non-spin
      for (uint i=0; i<TDOS.size();i++) {
        vtmp.clear();
        vtmp.push_back(TDOS.at(i).at(0)); //Energy
        for (uint j=1; j<PDOS.at(i).size();j++) {
          vtmp.push_back(PDOS.at(i).at(j)); //s, p, d
        }
        vtmp.push_back(TDOS.at(i).at(1)); //TDOS
        DOS.push_back(vtmp);
        vtmp.push_back(TDOS.at(i).at(2));   //CO20190808 - IDOS
        DOS_IDOS.push_back(vtmp); //CO20190808 - IDOS
      }
    }
    else if(TDOS.at(0).size()==5) { //spin
      for (uint i=0; i<TDOS.size();i++) {
        vtmp.clear();
        vtmp.push_back(TDOS.at(i).at(0)); //Energy
        for (uint j=1; j<PDOS.at(i).size();j++) {
          vtmp.push_back(PDOS.at(i).at(j)); //s, p ,d, f
        }
        vtmp.push_back(TDOS.at(i).at(1)); //TDOS up
        vtmp.push_back(TDOS.at(i).at(2)); //TDOS dn
        DOS.push_back(vtmp);
        vtmp.push_back(TDOS.at(i).at(3));vtmp.push_back(TDOS.at(i).at(4));   //CO20190808 - IDOS up/down
        DOS_IDOS.push_back(vtmp);  //CO20190808 - IDOS up/down
      }
    }
    else {
      cerr << "ERROR!!!!!!!!" << endl;
    }
  }
} // namespace pocc

// ***************************************************************************
// string pocc::POCC_GENERATE_GNUPLOTSCRIPT(vector<vector<double> >& DOS, const string& SystemName, const string& dosdatafile, const double& T, const double& DOS_Emin, double& DOS_Emax, const double& DOSSCALE, const double& DOSMAX)
// ***************************************************************************
namespace pocc {
  string POCC_GENERATE_GNUPLOTSCRIPT(vector<vector<double> >& DOS, const string& SystemName, const string& dosdatafile, const double& T, const double& DOS_Emin, double& DOS_Emax, const double& DOSSCALE, const double& DOSMAX) {
    stringstream  ss_gnu; ss_gnu.str(std::string());
    ss_gnu << "#Generated by AFLOW (Kesong Yang [kesong.yang@gmail.com], 2011, Duke)" << endl;
    ss_gnu << "set term postscript eps enhanced color font \"Times-Roman, 40\" size 18, 10.125" << endl;
    ss_gnu << "set output " << "\"" << SystemName <<"_DOS.eps" << "\"" << endl;
    ss_gnu << "set title \"" << SystemName << "\""<< endl;
    ss_gnu << endl;

    ss_gnu << "#DOS PLOT" << endl;
    ss_gnu << "set xtics 2" << endl;
    ss_gnu << "set ytics" << endl;
    ss_gnu << "set xrange [" << DOS_Emin << ":" << DOS_Emax << "]" << endl;
    if(DOS.at(0).size()==10||(DOS.at(0).size()==12)) {
      if(DOSMAX*DOSSCALE!=0) {
        ss_gnu << "set yrange [" << DOSMAX*DOSSCALE*(-1) <<":"  << DOSMAX*DOSSCALE << "]" << endl;
      } else {
        ss_gnu << "set yrange [0:2]" << endl;
      }
    }
    if(DOS.at(0).size()==5||(DOS.at(0).size()==6)) {
      if(DOSMAX*DOSSCALE!=0) {
        ss_gnu << "set yrange [0:" << DOSMAX*DOSSCALE << "]" << endl;
      }
    }
    ss_gnu << endl;
    ss_gnu << "set label '" << AFLOWLIB_CONSORTIUM_STRING << "' at screen 0.70, 0.02 font \"Times-Roman, 32\"" << endl;
    ss_gnu << "set xlabel 'Energy (eV)' offset graph 0.00" << endl;
    ss_gnu << "set ylabel 'Density of States (States/eV)' offset graph 0.00" << endl;
    ss_gnu << "set label 'E_F' at 0.05, graph 0.95" << endl;
    ss_gnu << "set label 'T = " << aurostd::utype2string(T) << "K' " << "at 0.3, graph 0.95" << endl;
    ss_gnu << "set arrow from 0, 0 to first 0, graph 1 nohead lt 3 lw 1.5" << endl;
    ss_gnu << "set arrow from 0, 0 to first 0, graph 0 nohead lt 3 lw 1.5" << endl;
    ss_gnu << endl;
    ss_gnu << "plot[][] \\" << endl;
    if(DOS.at(0).size()==5) {
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:2 w l lt  2 lw 6 title 's', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:3 w l lt  3 lw 6 title 'p', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:4 w l lt  8 lw 6 title 'd', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:5 w l lt -1 lw 2 title 'Total'" << endl;
    }
    if(DOS.at(0).size()==6) {
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:2 w l lt  2 lw 6 title 's', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:3 w l lt  3 lw 6 title 'p', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:4 w l lt  8 lw 6 title 'd', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:5 w l lt  5 lw 6 title 'f', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:6 w l lt -1 lw 2 title 'Total'" << endl;
    }
    if(DOS.at(0).size()==9) { // Only works for s, p, d and f orbitals; Spin-polarized
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:2 w l lt  2 lw 6 title 's', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:3 w l lt  2 lw 6 notitle, \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:4 w l lt  3 lw 6 title 'p', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:5 w l lt  3 lw 6 notitle, \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:6 w l lt  8 lw 6 title 'd', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:7 w l lt  8 lw 6 notitle, \\"  << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:8 w l lt -1 lw 2 title 'Total', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:9 w l lt -1 lw 2 notitle" << endl;
    }
    if(DOS.at(0).size()==11) { // Only works for s, p, d and f orbitals; Spin-polarized
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:2 w l lt  2 lw 6 title 's', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:3 w l lt  2 lw 6 notitle, \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:4 w l lt  3 lw 6 title 'p', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:5 w l lt  3 lw 6 notitle, \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:6 w l lt  8 lw 6 title 'd', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:7 w l lt  8 lw 6 notitle, \\"  << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:8 w l lt  5 lw 6 title 'f', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:9 w l lt  5 lw 6 notitle, \\"  << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:10 w l lt -1 lw 2 title 'Total', \\" << endl;
      ss_gnu << "\"" << dosdatafile << "\"" << " u 1:11 w l lt -1 lw 2 notitle"<< endl;
    }
    ss_gnu << endl;
    return ss_gnu.str();
  }
} // namespace pocc


// ***************************************************************************
// void pocc::POCC_GENERATE_OUTPUT(const string& directory, const double& T, const double& DOS_Emax, double& DOS_Emin, const double& DOSSCALE)
// ***************************************************************************
namespace pocc {
  void POCC_GENERATE_OUTPUT(const string& directory, const double& T, const double& DOS_Emin, double& DOS_Emax, const double& DOSSCALE) {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    //produce DOS data
    if(!XHOST.is_command("gnuplot")) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"gnuplot is missing",_RUNTIME_INIT_);}  //CO20200624
    if(!XHOST.is_command("convert")) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"convert is missing",_RUNTIME_INIT_);}  //CO20200624
    vector<vector<double> > TDOS_ONLY, PDOS_ONLY, DOS, DOS_IDOS;  vector<double> vEfermi,Egap,vprob; double mag,Egap_net; //CO20190808 - DOS_IDOS
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " starting" << endl;}
    POCC_GENERATE_DOSDATA(directory,T,TDOS_ONLY,PDOS_ONLY,vEfermi,mag,Egap_net,Egap,vprob);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " POCC_GENERATE_DOSDATA done" << endl;}
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " TDOS_ONLY.size()=" << TDOS_ONLY.size() << endl;
      for(uint i=0;i<TDOS_ONLY.size();i++){
        cerr << __AFLOW_FUNC__ << " TDOS_ONLY[i=" << i << "].size()=" << TDOS_ONLY[i].size() << endl;
      }
      cerr << __AFLOW_FUNC__ << " PDOS_ONLY.size()=" << PDOS_ONLY.size() << endl;
      for(uint i=0;i<PDOS_ONLY.size();i++){
        cerr << __AFLOW_FUNC__ << " PDOS_ONLY[i=" << i << "].size()=" << PDOS_ONLY[i].size() << endl;
      }
      cerr << __AFLOW_FUNC__ << " vEfermi.size()=" << vEfermi.size() << endl;
      cerr << __AFLOW_FUNC__ << " Egap.size()=" << Egap.size() << endl;
      cerr << __AFLOW_FUNC__ << " vprob.size()=" << vprob.size() << endl;
    }
    POCC_COMBINE_TDOS_PDOS_ONEDOS(TDOS_ONLY,PDOS_ONLY,DOS,DOS_IDOS);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " POCC_COMBINE_TDOS_PDOS_ONEDOS done" << endl;}
    DOS = aurostd::ShiftFirstColumn(DOS, max(vEfermi)); //shift to Efermi
    DOS_IDOS = aurostd::ShiftFirstColumn(DOS_IDOS, max(vEfermi)); //shift to Efermi
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ShiftFirstColumn done" << endl;}
    double DOSMAX = aurostd::FindMaxIn2DvectorExcept1stColumn(DOS, DOS_Emin, DOS_Emax);

    //write 2D vector into files
    vector<vector<double> >* p_DOS2FILE=&DOS;             //CO20190808 - write out DOS_IDOS if desired
    if(PRINT_AVG_IDOS){p_DOS2FILE=&DOS_IDOS;}             //CO20190808 - write out DOS_IDOS if desired
    const vector<vector<double> >& DOS2FILE=*p_DOS2FILE;  //CO20190808 - write out DOS_IDOS if desired
    string str_dos = aurostd::vector2string(DOS2FILE);    //CO20190808 - write out DOS_IDOS if desired
    string dosdatafile = "DOSDATA_" + aurostd::utype2string(T) + "K";
    aurostd::string2file(str_dos, dosdatafile);

    //generate gnuplot script
    string POCC_NAME = "untitled_" + aurostd::utype2string(T) + "K";
    string gnuplotscript = "GNUPLOT_POCC_" + POCC_NAME + "_DOS.gp";
    string str_gnu = POCC_GENERATE_GNUPLOTSCRIPT(DOS, POCC_NAME, dosdatafile, T, DOS_Emin, DOS_Emax, DOSSCALE, DOSMAX);
    aurostd::string2file(str_gnu, gnuplotscript);

    //run system command to plot
    aurostd::execute(XHOST.command("gnuplot")+" " + gnuplotscript);
    aurostd::execute(XHOST.command("convert")+" -background white ./" + POCC_NAME + "_DOS.eps ./" + POCC_NAME + "_DOS.png");

    if(!XHOST.vflag_control.flag("KEEP::GPL")) {
      aurostd::execute("rm  -f " + gnuplotscript);
      //     aurostd::execute("rm  -f " + dosdatafile); // always keep it
    }
    if(!XHOST.vflag_control.flag("KEEP::EPS")) {
      aurostd::execute("rm  -f " + POCC_NAME + "_DOS.eps");
    } 
  }
} // namespace pocc

// ***************************************************************************
// pocc::POCC_BANDGAP(string options) {
// ***************************************************************************
namespace pocc {
  void POCC_BANDGAP(string options) {
    // modified: Camilo Calderon 04 AUG 2014
    string directory;
    double Temperature, DOS_Emin, DOS_Emax, DOSSCALE; 
    pocc::POCC_SetOptions(options, directory, Temperature, DOS_Emin, DOS_Emax, DOSSCALE);
    vector<vector<double> > DOS, PDOS;  
    vector<double> vEfermi,vprob, Egap; double mag = 0.0, Egap_net=0.0;
    POCC_GENERATE_DOSDATA(directory,Temperature,DOS,PDOS,vEfermi,mag,Egap_net,Egap,vprob);
    // [NON_NECESSARY] DOS = aurostd::ShiftFirstColumn(DOS, max(vEfermi)); //shift to Efermi
    //  cout << "Egap.size()=" << Egap.size() << endl;
    if(Egap.size() == 0) {  // CAMILO, can Egap have size = 0 ?? this is the <xOUTCAR.Egap> with vprob
      cout << "Temperature: " << Temperature 
        << "  Egap_net:  " << Egap_net << endl ;
    }
    if(Egap.size() == 1) {
      cout << "Temperature: " << Temperature 
        << "  Egap up:  " << Egap.at(0) 
        << "  Egap_net:  " << Egap_net << endl ;
    }
    if(Egap.size() == 2) {
      cout << "Temperature: " << Temperature 
        << "  Egap up:  " << Egap.at(0) 
        << "  Egap dn:  " << Egap.at(1) 
        << "  Egap net: " << Egap_net << endl ;
    }
  }
}

// ***************************************************************************
// pocc::POCC_MAG(string options)
// ***************************************************************************
namespace pocc {
  void POCC_MAG(string options) {
    string directory;
    double T, DOS_Emin, DOS_Emax, DOSSCALE; 
    pocc::POCC_SetOptions(options, directory, T, DOS_Emin, DOS_Emax, DOSSCALE);
    vector<vector<double> > DOS, PDOS;  vector<double> vEfermi, vprob,Egap; double mag = 0.0,Egap_net=0.0;
    POCC_GENERATE_DOSDATA(directory,T,DOS,PDOS,vEfermi,mag,Egap_net,Egap,vprob);
    cout << "magnetic moment: "  << mag << endl;
  }
}

// ***************************************************************************
// Bond Class function
// ***************************************************************************

namespace pocc {
  Bond::Bond() {}
  Bond::Bond(const Bond& b) {copy(b);}
  Bond::~Bond() {}
  void Bond::free() {}
  void Bond::copy(const Bond& b) {
    bgn = b.bgn;
    end = b.end;
    length = b.length;
  }

  void Bond::Set(xstructure xstr, _atom atomBGN, _atom atomEND) {
    xstr.ReScale(1.0); //safety
    bgn = atomBGN;
    end = atomEND;
    length = AtomDist(bgn, end);
  }

  const Bond& Bond::operator=(const Bond& b) {
    if(this != &b) {copy(b);}
    return *this;
  }

  bool Bond::operator==(const Bond &other) const {
    bool FLAG = false;
    if((SameAtom(bgn, other.end) && SameAtom(end, other.bgn)) ||
        (SameAtom(bgn, other.bgn) && SameAtom(end, other.end))) {
      FLAG=true;
    }
    return FLAG;
  }

  bool Bond::operator!=(const Bond &other) const {
    return !(*this == other);
  }

  ostream& operator<<(ostream& os ,const Bond& bond) {
    cout << bond.bgn << endl;
    cout << bond.end << endl;
    cout << bond.length << endl;
    return os;
  }
}

// ***************************************************************************
// void pocc::RemoveSameBond(vector<Bond>& Bonds_orig, vector<Bond>& Bonds_new)
// ***************************************************************************
namespace pocc {
  void RemoveSameBond(vector<Bond>& Bonds_orig, vector<Bond>& Bonds_new) {
    Bonds_new.push_back(Bonds_orig.at(0));
    for (uint i=1; i<Bonds_orig.size();i++) {
      uint j;
      for (j=0; j<Bonds_new.size();j++) {
        if(Bonds_orig.at(i)==Bonds_new.at(j)) {
          break;
        }
      }
      if(j==Bonds_new.size()) {
        Bonds_new.push_back(Bonds_orig.at(i));
      }
    }
  }
} // namespace pocc

// ***************************************************************************
// void pocc::SetUFFPara(_atom atomi, _atom atomj, double& R0, double& Kij, double& Xij, double& Dij)
// ***************************************************************************
// Set UFF parameters
namespace pocc {
  void SetUFFPara(_atom atomi, _atom atomj, double& R0, double& Kij, double& Xij, double& Dij) {
    bool LDEBUG=(FALSE || XHOST.DEBUG); //CO20180409
    //double bondorder = 1.0; //single bond
    string atomi_name = atomi.name;
    string atomj_name = atomj.name;
    pocc::UFFPara uffparai; uffparai.GetUFFParameters(atomi_name);
    pocc::UFFPara uffparaj; uffparaj.GetUFFParameters(atomj_name);
    double ri = uffparai.r1, rj = uffparaj.r1; //bond length
    double Zi = uffparai.Z1, Zj = uffparaj.Z1; //effective charge
    double chiI = uffparai.Xi, chiJ = uffparaj.Xi; //GMP electronegativity
    double Xi = uffparai.x1, Xj = uffparaj.x1; //nonbond distance
    double Di = uffparai.D1, Dj = uffparaj.D1; //nonbond energy
    //From Equation 3
    //double rbo = -0.1332*(ri+rj)*log(bondorder); //zero for single bond
    //From Equation 4
    double ren = ri*rj*(pow((sqrt(chiI) - sqrt(chiJ)), 2.0)) / (chiI*ri +chiJ*rj);
    //From Equation 2
    //NOTE: See http://towhee.sourceforge.net/forcefields/uff.html
    //There is a typo in the published paper
    //R0 = ri + rj + rbo - ren;
    R0 = ri + rj - ren;
    Kij = 664.12*Zi*Zj/(R0*R0*R0); //From Equation 6
    Xij = sqrt(Xi*Xj); //From Equation 21b
    Dij = sqrt(Di*Dj); //From Equation 22
    if(LDEBUG) { //CO20180409
      cerr << "uffb.ren " << ren << endl;
      cerr << "uffb.R0 " << R0 << endl;
    }
  }   
} // namespace pocc

// ***************************************************************************
// double pocc::CalculateBondEnergy(xstructure xstr, _atom atomi, _atom atomj)
// ***************************************************************************
namespace pocc {
  double CalculateBondEnergy(xstructure xstr, _atom atomi, _atom atomj) {
    bool LDEBUG=(FALSE || XHOST.DEBUG); //CO20180409
    double r0, kij, xij_tmp, dij_tmp;
    pocc::SetUFFPara(atomi, atomj, r0, kij, xij_tmp, dij_tmp);

    xstr.ReScale(1.0); //Safety
    double rij = AtomDist(atomi, atomj); //Notice here rij is r, while r0 is rij in paper
    if(LDEBUG) { //CO20180409
      cerr << "distij " << rij << endl;
    }
    double delta = rij - r0;
    double delta2 = delta*delta;
    if(LDEBUG) { //CO20180409
      cerr << "uffb.Kij " << kij << endl;
      cerr << "uffb.delta " << delta << endl;
    }
    //From Equation 1a
    double energy = 0.5 * kij * delta2;  //unit kcal/mol
    energy = energy*KCAL_TO_EV;
    return energy;
  }
} // namespace pocc

// ***************************************************************************
// double pocc::CalculateNonBondEnergy(xstructure xstr, _atom atomi, _atom atomj)
// ***************************************************************************
namespace pocc {
  double CalculateNonBondEnergy(xstructure xstr, _atom atomi, _atom atomj) {
    bool LDEBUG=(FALSE || XHOST.DEBUG); //CO20180409
    //van der Waals (Nonbonded interaction)
    double r0_tmp, kij_tmp, Xij, Dij;
    pocc::SetUFFPara(atomi, atomj, r0_tmp, kij_tmp, Xij, Dij);

    xstr.ReScale(1.0); //Safety
    double x = AtomDist(atomi, atomj); //x distance
    double X6 = pow((Xij/x), 6); 
    double X12 = pow((Xij/x), 12);

    if(LDEBUG) { //CO20180409
      pocc::UFFPara uffparai; uffparai.GetUFFParameters(atomi.name);
      pocc::UFFPara uffparaj; uffparaj.GetUFFParameters(atomj.name);
      cerr << "Xi " << uffparai.Xi << endl;
      cerr << "Xj " << uffparaj.Xi << endl;
      cerr << "Xij " << Xij << endl;
      cerr << "distij " << x << endl;
      cerr << "Xij/distij " << Xij/x << endl;
      cerr << "Dij " << Dij << endl;
      cerr << "X12 " << X12 << endl;
      cerr << "X6 " << X6 << endl;
    }

    //From Equation 20
    double energy = Dij*(X12-2*X6);
    if(LDEBUG) { //CO20180409
      cerr << "atomi " << atomi << endl;
      cerr << "atomj " << atomj << endl;
      cerr << "distij " << x << endl;
      cerr << "energy " << energy << endl;
    }
    energy = energy*KCAL_TO_EV;
    return energy;
  }
} // namespace pocc

// ***************************************************************************
// double pocc::CalculateUFFEnergy(xstructure xstr)
// ***************************************************************************
namespace pocc {
  double CalculateUFFEnergy(xstructure xstr) {
    bool LDEBUG=(FALSE || XHOST.DEBUG); //CO20180409
    if(LDEBUG) { //CO20180409
      cerr << xstr << endl;
    }
    vector<Bond>  Bonds;
    vector<Bond>  NonBonds;
    pocc::AnalyzeBonds(xstr, Bonds, NonBonds);
    double totalenergy=0.0;

    //bond interaction
    double bondenergy=0.0;
    if(LDEBUG) { //CO20180409
      cerr << "BONDING LENGTH " << Bonds.size() << endl;
    }
    for (uint i=0; i<Bonds.size();i++) {
      Bond bondi = Bonds.at(i);
      _atom atomi = bondi.bgn;
      _atom atomj = bondi.end;
      if(LDEBUG) { //CO20180409
        cerr << "atom1 " << atomi << endl;
        cerr << "atom2 " << atomj << endl;
      }
      bondenergy += pocc::CalculateBondEnergy(xstr, atomi, atomj);
      if(LDEBUG) { //CO20180409
        cerr << "energy " << pocc::CalculateBondEnergy(xstr, atomi, atomj) << endl;
      }
    }

    if(LDEBUG) { //CO20180409
      cerr << "BONDING ENERGY " << bondenergy << endl;
    }

    //nonbond interaction, van der Waals
    double nonbondenergy=0.0;
    for (uint i=0; i<NonBonds.size();i++) {
      Bond bondi = NonBonds.at(i);
      _atom atomi = bondi.bgn;
      _atom atomj = bondi.end;
      nonbondenergy += pocc::CalculateNonBondEnergy(xstr, atomi, atomj);
    }

    totalenergy = bondenergy + nonbondenergy;
    return (totalenergy);
  }
} // namespace pocc

// ***************************************************************************
// void pocc::ExtractBonds(const xstructure& xstr, deque<deque<_atom> >& neigh_mat_bonded, deque<deque<_atom> >& neigh_mat_nonbonded)
// ***************************************************************************
namespace pocc {
  void ExtractBonds(const xstructure& xstr, deque<deque<_atom> >& neigh_mat_bonded, deque<deque<_atom> >& neigh_mat_nonbonded) {
    bool LDEBUG=(FALSE || XHOST.DEBUG); //CO20180409
    //Extract bonded and nonbonded atoms
    const double radius = 10.00; //12;//8.05323*1.5;//10.00; //CO TEST, this number is very important, and not well explained...
    deque<deque<_atom> > neigh_mat;
    xstructure xstr_tmp = xstr;
    xstr_tmp.GetNeighData(radius, neigh_mat); // radius 12 angstrom //CO20220623 - using new GetNeighData()

    deque<_atom>  atom_tmp_bonded;
    deque<_atom>  atom_tmp_nonbonded;

    if(LDEBUG) { //CO20180409
      cerr << "FULL " << endl;
      for(uint i=0;i<neigh_mat.size();i++){
        cerr << "i=" << i << " " << neigh_mat[i].size() << endl;
      }
    }

    for(uint ia=0; ia<neigh_mat.size();ia++) {
      atom_tmp_bonded.clear(); atom_tmp_nonbonded.clear();
      _atom a = neigh_mat.at(ia).at(0);
      _atom a1 = neigh_mat.at(ia).at(1);
      atom_tmp_bonded.push_back(a); atom_tmp_bonded.push_back(a1); //a, a1 must be bonded
      atom_tmp_nonbonded.push_back(a);

      double dr0 = AtomDist(a,a1);
      for(uint in=2; in<neigh_mat.at(ia).size();in++) {
        _atom an = neigh_mat.at(ia).at(in);
        if(abs(AtomDist(a,an)-dr0)<0.5) {  
          atom_tmp_bonded.push_back(an);
        }
        else {
          atom_tmp_nonbonded.push_back(an);
        }
      }
      neigh_mat_bonded.push_back(atom_tmp_bonded);
      neigh_mat_nonbonded.push_back(atom_tmp_nonbonded);
    }

    if(LDEBUG) { //CO20180409
      cerr << "BOND " << endl;
      for(uint i=0;i<neigh_mat_bonded.size();i++){
        cerr << "i=" << i << " " << neigh_mat_bonded[i].size() << endl;
      }
      cerr << "NONBOND " << endl;
      for(uint i=0;i<neigh_mat_nonbonded.size();i++){
        cerr << "i=" << i << " " << neigh_mat_nonbonded[i].size() << endl;
      }
    }
  }
} // namespace pocc


// ***************************************************************************
// void pocc::AnalyzeBonds(const xstructure& xstr, vector<Bond>& Bonds, vector<Bond>& NonBonds)
// ***************************************************************************
namespace pocc {
  void AnalyzeBonds(const xstructure& xstr, vector<Bond>& Bonds, vector<Bond>& NonBonds) {
    deque<deque<_atom> > neigh_mat_bonded;
    deque<deque<_atom> > neigh_mat_nonbonded;
    pocc::ExtractBonds(xstr, neigh_mat_bonded, neigh_mat_nonbonded);

    //Bonds 
    vector<vector<Bond> > vv_xstrBonds;
    vector<Bond>  v_xstrBonds;
    for (uint i=0; i<neigh_mat_bonded.size();i++) {
      _atom atomI=neigh_mat_bonded.at(i).at(0);
      v_xstrBonds.clear();
      for (uint j=1; j<neigh_mat_bonded.at(i).size();j++) {
        _atom atomJ = neigh_mat_bonded.at(i).at(j);
        Bond bondIJ; bondIJ.Set(xstr, atomI, atomJ);
        v_xstrBonds.push_back(bondIJ);
      }
      vv_xstrBonds.push_back(v_xstrBonds);
    }
    for (uint i=0;i<vv_xstrBonds.size();i++) {
      for (uint j=0; j<vv_xstrBonds.at(i).size();j++) {
        Bonds.push_back(vv_xstrBonds.at(i).at(j));
      }
    }

    //NonBonds
    vector<vector<Bond> > vv_xstrNonBonds;
    vector<Bond>  v_xstrNonBonds;
    for (uint i=0; i<neigh_mat_nonbonded.size();i++) {
      _atom atomI=neigh_mat_nonbonded.at(i).at(0);
      v_xstrNonBonds.clear();
      for (uint j=1; j<neigh_mat_nonbonded.at(i).size();j++) {
        _atom atomJ = neigh_mat_nonbonded.at(i).at(j);
        Bond bondIJ; bondIJ.Set(xstr, atomI, atomJ);
        v_xstrNonBonds.push_back(bondIJ);
      }
      vv_xstrNonBonds.push_back(v_xstrNonBonds);
    }
    for (uint i=0;i<vv_xstrNonBonds.size();i++) {
      for (uint j=0; j<vv_xstrNonBonds.at(i).size();j++) {
        NonBonds.push_back(vv_xstrNonBonds.at(i).at(j));
      }
    }
  }
} // namespace pocc

// ***************************************************************************
// pocc::UFFENERGY(istream& input)
// ***************************************************************************
namespace pocc {
  void UFFENERGY(istream& input) {
    xstructure xstr; 
    xstr.clear(); //DX20191220 - uppercase to lowercase clear
    xstr=xstructure(input, IOVASP_POSCAR);
    cout << "total energy: " <<  pocc::CalculateUFFEnergy(xstr) << endl;
  }
} // namespace pocc

//References: 
//"UFF, a full periodic table force field for molecular mechanics and molecular dynamics simulations"
//A. K. Rappe, C. J. Casewit, K. S. Colwell, W. A. Goddard III, W. M. Skiff, JACS, 114, 10024 (1992)

// ***************************************************************************
// string pocc::ReturnAtomSpecies(string& atom)
// ***************************************************************************
// input for GULP program
namespace pocc {
  string ReturnAtomSpecies(string atom) {
    if(atom.compare("H")==0)    {return  "H  core";}
    else if(atom.compare("He")==0)   {return  "He core";}
    else if(atom.compare("Li")==0)   {return  "Li core";}
    else if(atom.compare("Be")==0)   {return  "Be core";}
    else if(atom.compare("B")==0)    {return  "B  core";}
    else if(atom.compare("C")==0)    {return  "C  core";}
    else if(atom.compare("N")==0)    {return  "N  core";}
    else if(atom.compare("O")==0)    {return  "O  core";}
    else if(atom.compare("F")==0)    {return  "F  core";}
    else if(atom.compare("Ne")==0)   {return  "Ne core";}
    else if(atom.compare("Na")==0)   {return  "Na core";}
    else if(atom.compare("Mg")==0)   {return  "Mg core";}
    else if(atom.compare("Al")==0)   {return  "Al core";}
    else if(atom.compare("Si")==0)   {return  "Si core";}
    else if(atom.compare("P")==0)    {return  "P  core";}
    else if(atom.compare("S")==0)    {return  "S  core";}
    else if(atom.compare("Cl")==0)   {return  "Cl core";}
    else if(atom.compare("Ar")==0)   {return  "Ar core";}
    else if(atom.compare("K")==0)    {return  "K  core";}
    else if(atom.compare("Ca")==0)   {return  "Ca core";}
    else if(atom.compare("Sc")==0)   {return  "Sc core";}
    else if(atom.compare("Ti")==0)   {return  "Ti core";}
    else if(atom.compare("V")==0)    {return  "V  core";}
    else if(atom.compare("Cr")==0)   {return  "Cr core";}
    else if(atom.compare("Mn")==0)   {return  "Mn core";}
    else if(atom.compare("Fe")==0)   {return  "Fe core";}
    else if(atom.compare("Co")==0)   {return  "Co core";}
    else if(atom.compare("Ni")==0)   {return  "Ni core";}
    else if(atom.compare("Cu")==0)   {return  "Cu core";}
    else if(atom.compare("Zn")==0)   {return  "Zn core";}
    else if(atom.compare("Ga")==0)   {return  "Ga core";}
    else if(atom.compare("Ge")==0)   {return  "Ge core";}
    else if(atom.compare("As")==0)   {return  "As core";}
    else if(atom.compare("Se")==0)   {return  "Se core";}
    else if(atom.compare("Br")==0)   {return  "Br core";}
    else if(atom.compare("Kr")==0)   {return  "Kr core";}
    else if(atom.compare("Rb")==0)   {return  "Rb core";}
    else if(atom.compare("Sr")==0)   {return  "Sr core";}
    else if(atom.compare("Y")==0)    {return  "Y  core";}
    else if(atom.compare("Zr")==0)   {return  "Zr core";}
    else if(atom.compare("Nb")==0)   {return  "Nb core";}
    else if(atom.compare("Mo")==0)   {return  "Mo core";}
    else if(atom.compare("Tc")==0)   {return  "Tc core";}
    else if(atom.compare("Ru")==0)   {return  "Ru core";}
    else if(atom.compare("Rh")==0)   {return  "Rh core";}
    else if(atom.compare("Pd")==0)   {return  "Pd core";}
    else if(atom.compare("Ag")==0)   {return  "Ag core";}
    else if(atom.compare("Cd")==0)   {return  "Cd core";}
    else if(atom.compare("In")==0)   {return  "In core";}
    else if(atom.compare("Sn")==0)   {return  "Sn core";}
    else if(atom.compare("Sb")==0)   {return  "Sb core";}
    else if(atom.compare("Te")==0)   {return  "Te core";}
    else if(atom.compare("I")==0)    {return  "I  core";}
    else if(atom.compare("Xe")==0)   {return  "Xe core";}
    else if(atom.compare("Cs")==0)   {return  "Cs core";}
    else if(atom.compare("Ba")==0)   {return  "Ba core";}
    else if(atom.compare("La")==0)   {return  "La core";}
    else if(atom.compare("Ce")==0)   {return  "Ce core";}
    else if(atom.compare("Pr")==0)   {return  "Pr core";}
    else if(atom.compare("Nd")==0)   {return  "Nd core";}
    else if(atom.compare("Pm")==0)   {return  "Pm core";}
    else if(atom.compare("Sm")==0)   {return  "Sm core";}
    else if(atom.compare("Eu")==0)   {return  "Eu core";}
    else if(atom.compare("Gd")==0)   {return  "Gd core";}
    else if(atom.compare("Tb")==0)   {return  "Tb core";}
    else if(atom.compare("Dy")==0)   {return  "Dy core";}
    else if(atom.compare("Ho")==0)   {return  "Ho core";}
    else if(atom.compare("Er")==0)   {return  "Er core";}
    else if(atom.compare("Tm")==0)   {return  "Tm core";}
    else if(atom.compare("Yb")==0)   {return  "Yb core";}
    else if(atom.compare("Lu")==0)   {return  "Lu core";}
    else if(atom.compare("Hf")==0)   {return  "Hf core";}
    else if(atom.compare("Ta")==0)   {return  "Ta core";}
    else if(atom.compare("W")==0)    {return  "W  core";}
    else if(atom.compare("Re")==0)   {return  "Re core";}
    else if(atom.compare("Os")==0)   {return  "Os core";}
    else if(atom.compare("Ir")==0)   {return  "Ir core";}
    else if(atom.compare("Pt")==0)   {return  "Pt core";}
    else if(atom.compare("Au")==0)   {return  "Au core";}
    else if(atom.compare("Hg")==0)   {return  "Hg core";}
    else if(atom.compare("Tl")==0)   {return  "Tl core";}
    else if(atom.compare("Pb")==0)   {return  "Pb core";}
    else if(atom.compare("Bi")==0)   {return  "Bi core";}
    else if(atom.compare("Po")==0)   {return  "Po core";}
    else if(atom.compare("At")==0)   {return  "At core";}
    else if(atom.compare("Rn")==0)   {return  "Rn core";}
    else if(atom.compare("Fr")==0)   {return  "Fr core";}
    else if(atom.compare("Ra")==0)   {return  "Ra core";}
    else if(atom.compare("Ac")==0)   {return  "Ac core";}
    else if(atom.compare("Th")==0)   {return  "Th core";}
    else if(atom.compare("Pa")==0)   {return  "Pa core";}
    else if(atom.compare("U")==0)    {return  "U  core";}
    else if(atom.compare("Np")==0)   {return  "Np core";}
    else if(atom.compare("Pu")==0)   {return  "Pu core";}
    else if(atom.compare("Am")==0)   {return  "Am core";}
    else if(atom.compare("Cm")==0)   {return  "Cm core";}
    else if(atom.compare("Bk")==0)   {return  "Bk core";}
    else if(atom.compare("Cf")==0)   {return  "Cf core";}
    else if(atom.compare("Es")==0)   {return  "Es core";}
    else if(atom.compare("Fm")==0)   {return  "Fm core";}
    else if(atom.compare("Md")==0)   {return  "Md core";}
    else if(atom.compare("No")==0)   {return  "No core";}
    else if(atom.compare("Lr")==0)   {return  "Lr core";}
    else {return "Unknown";}
  }
} // namespace pocc 
// ***************************************************************************
// string pocc::ReturnAtomSpecies(string& atom)
// ***************************************************************************
// Potential for GULP, only consider sp1
namespace pocc {
  string ReturnAtomSpeciesPotential(string atom) {
    if(atom.compare("H")==0)        {return  "H    0.354   180.000 2.886   0.044   12.000 0.7120 0 0.000  0.0  0.0000  4.528";}      
    else if(atom.compare("He")==0)  {return  "He   0.849   90.000  2.362   0.056   15.240 0.0972 4 0.000  0.0  0.0000  9.660";}      
    else if(atom.compare("Li")==0)  {return  "Li   1.336   180.000 2.451   0.025   12.000 1.0255 0 0.000  0.0  0.0000  3.006";}      
    else if(atom.compare("Be")==0)  {return  "Be   1.074   109.470 2.745   0.085   12.000 1.5650 3 0.000  0.0  0.0000  4.877";}      
    else if(atom.compare("B")==0)   {return  "B    0.838   109.470 4.083   0.180   12.052 1.7550 3 0.000  0.0  0.0000  4.068";}      
    else if(atom.compare("C")==0)   {return  "C    0.706   180.000 3.851   0.105   12.730 1.9120 3 2.119  0.0  0.0000  5.343";}      
    else if(atom.compare("N")==0)   {return  "N    0.656   180.000 3.660   0.069   13.407 2.5438 3 0.450  0.0 61.2230  6.899";}      
    else if(atom.compare("O")==0)   {return  "O    0.639   180.000 3.500   0.060   14.085 2.2998 3 0.018  0.0  0.0000  8.741";}      
    else if(atom.compare("F")==0)   {return  "F    0.668   180.000 3.364   0.050   14.762 1.7350 0 0.000  0.0  0.0000 10.874";}      
    else if(atom.compare("Ne")==0)  {return  "Ne   0.920   90.000  3.243   0.042   15.440 0.1944 4 0.000  0.0  0.0000 11.040";}      
    else if(atom.compare("Na")==0)  {return  "Na   1.539   180.000 2.983   0.030   12.000 1.0809 0 0.000  0.0  0.0000  2.843";}      
    else if(atom.compare("Mg")==0)  {return  "Mg   1.421   109.470 3.021   0.111   12.000 1.7866 3 0.000  0.0  0.0000  3.951";}      
    else if(atom.compare("Al")==0)  {return  "Al   1.244   109.470 4.499   0.505   11.278 1.7924 3 0.000  0.0  0.0000  3.041";}      
    else if(atom.compare("Si")==0)  {return  "Si   1.117   109.470 4.295   0.402   12.175 2.3232 3 1.225  0.0  0.0000  4.168";}      
    else if(atom.compare("P")==0)   {return  "P    1.101   93.800  4.147   0.305   13.072 2.8627 3 2.400 22.0 84.4339  5.463";}      
    else if(atom.compare("S")==0)   {return  "S    1.064   92.1000 4.035   0.274   13.969 2.7032 3 0.484  0.0  0.0000  6.928";}      
    else if(atom.compare("Cl")==0)  {return  "Cl   1.044   180.000 3.947   0.227   14.866 2.3484 0 0.000  0.0  0.0000  8.564";}      
    else if(atom.compare("Ar")==0)  {return  "Ar   1.032   90.000  3.868   0.185   15.763 0.2994 4 0.000  0.0  0.0000  9.465";}      
    else if(atom.compare("K")==0)   {return  "K    1.953   180.000 3.812   0.035   12.000 1.1645 0 0.000  0.0  0.0000  2.421";}      
    else if(atom.compare("Ca")==0)  {return  "Ca   1.761   90.000  3.399   0.238   12.000 2.1414 6 0.000  0.0  0.0000  3.231";}      
    else if(atom.compare("Sc")==0)  {return  "Sc   1.513   109.470 3.295   0.019   12.000 2.5924 3 0.000  0.0  0.0000  3.395";}      
    else if(atom.compare("Ti")==0)  {return  "Ti   1.412   109.470 3.175   0.017   12.000 2.6595 3 0.000  0.0  0.0000  3.470";}      
    else if(atom.compare("V")==0)   {return  "V    1.402   109.470 3.144   0.016   12.000 2.6789 3 0.000  0.0  0.0000  3.650";}      
    else if(atom.compare("Cr")==0)  {return  "Cr   1.345   90.000  3.023   0.015   12.000 2.4631 6 0.000  0.0  0.0000  3.415";}      
    else if(atom.compare("Mn")==0)  {return  "Mn   1.382   90.000  2.961   0.013   12.000 2.4301 6 0.000  0.0  0.0000  3.325";}      
    else if(atom.compare("Fe")==0)  {return  "Fe   1.270   109.470 2.912   0.013   12.000 2.4301 3 0.000  0.0  0.0000  3.760";}      
    else if(atom.compare("Co")==0)  {return  "Co   1.241   90.000  2.872   0.014   12.000 2.4301 6 0.000  0.0  0.0000  4.105";}      
    else if(atom.compare("Ni")==0)  {return  "Ni   1.164   90.000  2.834   0.015   12.000 2.4301 4 0.000  0.0  0.0000  4.465";}      
    else if(atom.compare("Cu")==0)  {return  "Cu   1.302   109.470 3.495   0.005   12.000 1.7565 3 0.000  0.0  0.0000  3.729";}      
    else if(atom.compare("Zn")==0)  {return  "Zn   1.193   109.470 2.763   0.124   12.000 1.3084 3 0.000  0.0  0.0000  5.106";}      
    else if(atom.compare("Ga")==0)  {return  "Ga   1.260   109.470 4.383   0.415   11.000 1.8206 3 0.000  0.0  0.0000  2.999";}      
    else if(atom.compare("Ge")==0)  {return  "Ge   1.197   109.470 4.280   0.379   12.000 2.7888 3 0.701  0.0  0.0000  4.051";}      
    else if(atom.compare("As")==0)  {return  "As   1.211   92.100  4.230   0.309   13.000 2.8640 3 1.500 22.0 86.9735  5.188";}      
    else if(atom.compare("Se")==0)  {return  "Se   1.190   90.600  4.205   0.291   14.000 2.7645 3 0.335  0.0  0.0000  6.428";}      
    else if(atom.compare("Br")==0)  {return  "Br   1.192   180.000 4.189   0.251   15.000 2.5186 0 0.000  0.0  0.0000  7.790";}      
    else if(atom.compare("Kr")==0)  {return  "Kr   1.147   90.000  4.141   0.220   16.000 0.4520 4 0.000  0.0  0.0000  8.505";}      
    else if(atom.compare("Rb")==0)  {return  "Rb   2.260   180.000 4.114   0.040   12.000 1.5922 0 0.000  0.0  0.0000  2.331";}      
    else if(atom.compare("Sr")==0)  {return  "Sr   2.052   90.000  3.641   0.235   12.000 2.4486 6 0.000  0.0  0.0000  3.024";}      
    else if(atom.compare("Y")==0)   {return  "Y    1.698   109.470 3.345   0.072   12.000 3.2573 3 0.000  0.0  0.0000  3.830";}      
    else if(atom.compare("Zr")==0)  {return  "Zr   1.564   109.470 3.124   0.069   12.000 3.6675 3 0.000  0.0  0.0000  3.400";}      
    else if(atom.compare("Nb")==0)  {return  "Nb   1.473   109.470 3.165   0.059   12.000 3.6179 3 0.000  0.0  0.0000  3.550";}      
    else if(atom.compare("Mo")==0)  {return  "Mo   1.467   90.000  3.052   0.056   12.000 3.4021 6 0.000  0.0  0.0000  3.465";}      
    else if(atom.compare("Tc")==0)  {return  "Mo   1.322   90.000  2.998   0.048   12.000 3.4021 3 0.000  0.0  0.0000  3.465";}      
    else if(atom.compare("Ru")==0)  {return  "Tc   1.478   90.000  2.963   0.056   12.000 3.4021 6 0.000  0.0  0.0000  3.290";}      
    else if(atom.compare("Rh")==0)  {return  "Ru   1.332   90.000  2.929   0.053   12.000 3.4021 6 0.000  0.0  0.0000  3.575";}      
    else if(atom.compare("Pd")==0)  {return  "Pd   1.338   90.000  2.899   0.048   12.000 3.2077 4 0.000  0.0  0.0000  4.320";}      
    else if(atom.compare("Ag")==0)  {return  "Ag   1.386   180.000 3.148   0.036   12.000 1.9557 1 0.200  0.0  0.0000  4.436";}      
    else if(atom.compare("Cd")==0)  {return  "Cd   1.403   109.470 2.848   0.228   12.000 1.6525 3 0.000  0.0  0.0000  5.034";}      
    else if(atom.compare("In")==0)  {return  "In   1.459   109.470 4.463   0.599   11.000 2.0704 3 0.000  0.0  0.0000  2.997";}      
    else if(atom.compare("Sn")==0)  {return  "Sn   1.398   109.470 4.392   0.567   12.000 2.9608 3 0.199  0.0  0.0000  3.987";}      
    else if(atom.compare("Sb")==0)  {return  "Sb   1.407   91.600  4.420   0.449   13.000 2.7042 3 1.100 22.0 87.7047  4.899";}      
    else if(atom.compare("Te")==0)  {return  "Te   1.386   90.250  4.470   0.398   14.000 2.8821 3 0.300  0.0  0.0000  5.816";}     
    else if(atom.compare("I")==0)   {return  "I    1.382   180.000 4.500   0.339   15.000 2.6537 0 0.000  0.0  0.0000  6.822";}     
    else if(atom.compare("Xe")==0)  {return  "Xe   1.267   90.000  4.404   0.332   12.000 0.5560 4 0.000  0.0  0.0000  7.595";}     
    else if(atom.compare("Cs")==0)  {return  "Cs   2.570   180.000 4.517   0.045   12.000 1.5728 0 0.000  0.0  0.0000  2.183";}     
    else if(atom.compare("Ba")==0)  {return  "Ba   2.277   90.000  3.703   0.364   12.000 2.7266 6 0.000  0.0  0.0000  2.814";}     
    else if(atom.compare("La")==0)  {return  "La   1.943   109.470 3.522   0.017   12.000 3.3049 3 0.000  0.0  0.0000  2.8355";}     
    else if(atom.compare("Ce")==0)  {return  "Ce   1.841   90.000  3.556   0.013   12.000 3.3049 6 0.000  0.0  0.0000  2.774";}     
    else if(atom.compare("Pr")==0)  {return  "Pr   1.823   90.000  3.606   0.010   12.000 3.3049 6 0.000  0.0  0.0000  2.858";}     
    else if(atom.compare("Nd")==0)  {return  "Nd   1.816   90.000  3.575   0.010   12.000 3.3049 6 0.000  0.0  0.0000  2.8685";}     
    else if(atom.compare("Pm")==0)  {return  "Pm   1.801   90.000  3.547   0.009   12.000 3.3049 6 0.000  0.0  0.0000  2.881";}     
    else if(atom.compare("Sm")==0)  {return  "Sm   1.780   90.000  3.520   0.008   12.000 3.3049 6 0.000  0.0  0.0000  2.9115";}     
    else if(atom.compare("Eu")==0)  {return  "Eu   1.771   90.000  3.493   0.008   12.000 3.3049 6 0.000  0.0  0.0000  2.8785";}     
    else if(atom.compare("Gd")==0)  {return  "Gd   1.735   90.000  3.368   0.009   12.000 3.3049 6 0.000  0.0  0.0000  3.1665";}     
    else if(atom.compare("Tb")==0)  {return  "Tb   1.732   90.000  3.451   0.007   12.000 3.3049 6 0.000  0.0  0.0000  3.018";}     
    else if(atom.compare("Dy")==0)  {return  "Dy   1.710   90.000  3.428   0.007   12.000 3.3049 6 0.000  0.0  0.0000  3.0555";}     
    else if(atom.compare("Ho")==0)  {return  "Ho   1.696   90.000  3.409   0.007   12.000 3.4157 6 0.000  0.0  0.0000  3.127";}     
    else if(atom.compare("Er")==0)  {return  "Er   1.673   90.000  3.391   0.007   12.000 3.3049 6 0.000  0.0  0.0000  3.1865";}      
    else if(atom.compare("Tm")==0)  {return  "Tm   1.660   90.000  3.374   0.006   12.000 3.3049 6 0.000  0.0  0.0000  3.2514";}      
    else if(atom.compare("Yb")==0)  {return  "Yb   1.637   90.000  3.355   0.228   12.000 2.6177 6 0.000  0.0  0.0000  3.2889";}      
    else if(atom.compare("Lu")==0)  {return  "Lu   1.671   90.000  3.640   0.041   12.000 3.2709 6 0.000  0.0  0.0000  2.9629";}      
    else if(atom.compare("Hf")==0)  {return  "Hf   1.611   109.470 3.141   0.072   12.000 3.9212 3 0.000  0.0  0.0000  3.7000";}      
    else if(atom.compare("Ta")==0)  {return  "Ta   1.511   109.470 3.170   0.081   12.000 4.0748 3 0.000  0.0  0.0000  5.1000";}      
    else if(atom.compare("W")==0)   {return  "W    1.392   90.000  3.069   0.067   12.000 3.6937 6 0.000  0.0  0.0000  4.6300";}      
    else if(atom.compare("Re")==0)  {return  "Re   1.372   90.000  2.954   0.066   12.000 3.6937 6 0.000  0.0  0.0000  3.9600";}      
    else if(atom.compare("Os")==0)  {return  "Os   1.372   90.000  3.120   0.037   12.000 3.6937 6 0.000  0.0  0.0000  5.1400";}      
    else if(atom.compare("Ir")==0)  {return  "Ir   1.371   90.000  2.840   0.073   12.000 3.7307 6 0.000  0.0  0.0000  5.0000";}      
    else if(atom.compare("Pt")==0)  {return  "Pt   1.364   90.000  2.754   0.080   12.000 3.3817 4 0.000  0.0  0.0000  4.7900";}      
    else if(atom.compare("Au")==0)  {return  "Au   1.262   90.000  3.293   0.039   12.000 2.6255 4 0.000  0.0  0.0000  4.8940";}      
    else if(atom.compare("Hg")==0)  {return  "Hg   1.340   180.000 2.705   0.385   12.000 1.7497 1 0.000  0.0  0.0000  6.2700";}      
    else if(atom.compare("Tl")==0)  {return  "Tl   1.518   120.000 4.347   0.680   11.000 2.0685 3 0.000  0.0  0.0000  3.2000";}      
    else if(atom.compare("Pb")==0)  {return  "Pb   1.459   109.470 4.297   0.663   12.000 2.8461 3 0.100  0.0  0.0000  3.9000";}      
    else if(atom.compare("Bi")==0)  {return  "Bi   1.512   90.000  4.370   0.518   13.000 2.4700 3 1.000 22.0 90.0000  4.6900";}      
    else if(atom.compare("Po")==0)  {return  "Po   1.500   90.000  4.709   0.325   14.000 2.3329 3 0.300  0.0  0.0000  4.2100";}      
    else if(atom.compare("At")==0)  {return  "At   1.545   180.000 4.750   0.284   15.000 2.2357 0 0.000  0.0  0.0000  4.7500";}      
    else if(atom.compare("Rn")==0)  {return  "Rn   1.420   90.000  4.765   0.248   16.000 0.5832 4 0.000  0.0  0.0000  5.3700";}      
    else if(atom.compare("Fr")==0)  {return  "Fr   2.880   180.000 4.900   0.050   12.000 1.8469 0 0.000  0.0  0.0000  2.0000";}      
    else if(atom.compare("Ra")==0)  {return  "Ra   2.512   90.000  3.677   0.404   12.000 2.9161 6 0.000  0.0  0.0000  2.8430";}      
    else if(atom.compare("Ac")==0)  {return  "Ac   1.983   90.000  3.478   0.033   12.000 3.8882 6 0.000  0.0  0.0000  2.8350";}      
    else if(atom.compare("Th")==0)  {return  "Th   1.721   90.000  3.396   0.026   12.000 4.2021 6 0.000  0.0  0.0000  3.1750";}      
    else if(atom.compare("Pa")==0)  {return  "Pa   1.711   90.000  3.424   0.022   12.000 3.8882 6 0.000  0.0  0.0000  2.9850";}      
    else if(atom.compare("U")==0)   {return  "U    1.684   90.000  3.395   0.022   12.000 3.8882 6 0.000  0.0  0.0000  3.3410";}      
    else if(atom.compare("Np")==0)  {return  "Np   1.666   90.000  3.424   0.019   12.000 3.8882 6 0.000  0.0  0.0000  3.5490";}      
    else if(atom.compare("Pu")==0)  {return  "Pu   1.657   90.000  3.424   0.016   12.000 3.8882 6 0.000  0.0  0.0000  3.2430";}      
    else if(atom.compare("Am")==0)  {return  "Am   1.660   90.000  3.381   0.014   12.000 3.8882 6 0.000  0.0  0.0000  2.9895";}      
    else if(atom.compare("Cm")==0)  {return  "Cm   1.801   90.000  3.326   0.013   12.000 3.8882 6 0.000  0.0  0.0000  2.8315";}      
    else if(atom.compare("Bk")==0)  {return  "Bk   1.761   90.000  3.339   0.013   12.000 3.8882 6 0.000  0.0  0.0000  3.1935";}      
    else if(atom.compare("Cf")==0)  {return  "Cf   1.750   90.000  3.313   0.013   12.000 3.8882 6 0.000  0.0  0.0000  3.1970";}      
    else if(atom.compare("Es")==0)  {return  "Es   1.724   90.000  3.299   0.012   12.000 3.8882 6 0.000  0.0  0.0000  3.3330";}      
    else if(atom.compare("Fm")==0)  {return  "Fm   1.712   90.000  3.286   0.012   12.000 3.8882 6 0.000  0.0  0.0000  3.4000";}      
    else if(atom.compare("Md")==0)  {return  "Md   1.689   90.000  3.274   0.011   12.000 3.8882 6 0.000  0.0  0.0000  3.4700";}      
    else if(atom.compare("No")==0)  {return  "No   1.679   90.000  3.248   0.011   12.000 3.8882 6 0.000  0.0  0.0000  3.4750";}      
    else if(atom.compare("Lr")==0)  {return  "Lr   1.698   90.000  3.236   0.011   12.000 3.8882 6 0.000  0.0  0.0000  3.5000";}      
    else {return "Unknown";}                                                                                                    
  }                                                                                                                            
} // namespace pocc                                                                                                            

// ***************************************************************************                                                 
// string pocc::ReturnUFFParameters(string atom)                                                                               
// ***************************************************************************                                                 
// UFF parameters                                                                                                              
// Force field parameters for UFF, the Universal Force Field                                                                   
// Atom   r1 theta0  x1  D1  zeta    Z1  Vi  Uj  Xi  Hard    Radius                                                            
namespace pocc {                                                                                                               
  string ReturnUFFParameters(string atom) {                                                                                    
    if(atom=="H")       {return "H   0.354   180.000 2.886   0.044   12.000  0.712   0.000   0.000   4.528   6.945   0.371   ";}   
    else if(atom=="He") {return "He  0.849   90.000  2.362   0.056   15.240  0.098   0.000   0.000   9.660   14.920  1.300   ";}   
    else if(atom=="Li") {return "Li  1.336   180.000 2.451   0.025   12.000  1.026   0.000   2.000   3.006   2.386   1.557   ";}   
    else if(atom=="Be") {return "Be  1.074   109.470 2.745   0.085   12.000  1.565   0.000   2.000   4.877   4.443   1.240   ";}  
    else if(atom=="B")  {return "B   0.838   109.470 4.083   0.180   12.052  1.755   0.000   2.000   5.110   4.750   0.822   ";}  
    else if(atom=="C")  {return "C   0.706   180.000 3.851   0.105   12.730  1.912   0.000   2.000   5.343   5.063   0.759   ";}  
    else if(atom=="N")  {return "N   0.656   180.000 3.660   0.069   13.407  2.544   0.000   2.000   6.899   5.880   0.715   ";}  
    else if(atom=="O")  {return "O   0.639   180.000 3.500   0.060   14.085  2.300   0.000   2.000   8.741   6.682   0.669   ";}  
    else if(atom=="F")  {return "F   0.668   180.000 3.364   0.050   14.762  1.735   0.000   2.000   10.874  7.474   0.706   ";}  
    else if(atom=="Ne") {return "Ne  0.920   90.000  3.243   0.042   15.440  0.194   0.000   2.000   11.040  10.550  1.768   ";}  
    else if(atom=="Na") {return "Na  1.539   180.000 2.983   0.030   12.000  1.081   0.000   1.250   2.843   2.296   2.085   ";}  
    else if(atom=="Mg") {return "Mg  1.421   109.470 3.021   0.111   12.000  1.787   0.000   1.250   3.951   3.693   1.500   ";}  
    else if(atom=="Al") {return "Al  1.244   109.470 4.499   0.505   11.278  1.792   0.000   1.250   4.060   3.590   1.201   ";}  
    else if(atom=="Si") {return "Si  1.117   109.470 4.295   0.402   12.175  2.323   1.225   1.250   4.168   3.487   1.176   ";}  
    else if(atom=="P")  {return "P   1.101   93.800  4.147   0.305   13.072  2.863   2.400   1.250   5.463   4.000   1.102   ";}  
    else if(atom=="S")  {return "S   1.064   92.1000 4.035   0.274   13.969  2.703   0.000   1.250   6.928   4.486   1.047   ";}  
    else if(atom=="Cl") {return "Cl  1.044   180.000 3.947   0.227   14.866  2.348   0.000   1.250   8.564   4.946   0.994   ";}  
    else if(atom=="Ar") {return "Ar  1.032   90.000  3.868   0.185   15.763  0.300   0.000   1.250   9.465   6.355   2.108   ";}  
    else if(atom=="K")  {return "K   1.953   180.000 3.812   0.035   12.000  1.165   0.000   0.700   2.421   1.920   2.586   ";}  
    else if(atom=="Ca") {return "Ca  1.761   90.000  3.399   0.238   12.000  2.141   0.000   0.700   3.231   2.880   2.000   ";}  
    else if(atom=="Sc") {return "Sc  1.513   109.470 3.295   0.019   12.000  2.592   0.000   0.700   3.395   3.080   1.750   ";}  
    else if(atom=="Ti") {return "Ti  1.412   109.470 3.175   0.017   12.000  2.659   0.000   0.700   3.470   3.380   1.607   ";}  
    else if(atom=="V")  {return "V   1.402   109.470 3.144   0.016   12.000  2.679   0.000   0.700   3.650   3.410   1.470   ";}  
    else if(atom=="Cr") {return "Cr  1.345   90.000  3.023   0.015   12.000  2.463   0.000   0.700   3.415   3.865   1.402   ";}  
    else if(atom=="Mn") {return "Mn  1.382   90.000  2.961   0.013   12.000  2.430   0.000   0.700   3.325   4.105   1.533   ";}  
    else if(atom=="Fe") {return "Fe  1.270   109.470 2.912   0.013   12.000  2.430   0.000   0.700   3.760   4.140   1.393   ";}  
    else if(atom=="Co") {return "Co  1.241   90.000  2.872   0.014   12.000  2.430   0.000   0.700   4.105   4.175   1.406   ";}  
    else if(atom=="Ni") {return "Ni  1.164   90.000  2.834   0.015   12.000  2.430   0.000   0.700   4.465   4.205   1.398   ";}  
    else if(atom=="Cu") {return "Cu  1.302   109.470 3.495   0.005   12.000  1.756   0.000   0.700   4.200   4.220   1.434   ";}  
    else if(atom=="Zn") {return "Zn  1.193   109.470 2.763   0.124   12.000  1.308   0.000   0.700   5.106   4.285   1.400   ";}  
    else if(atom=="Ga") {return "Ga  1.260   109.470 4.383   0.415   11.000  1.821   0.000   0.700   3.641   3.160   1.211   ";}  
    else if(atom=="Ge") {return "Ge  1.197   109.470 4.280   0.379   12.000  2.789   0.701   0.700   4.051   3.438   1.189   ";}  
    else if(atom=="As") {return "As  1.211   92.100  4.230   0.309   13.000  2.864   1.500   0.700   5.188   3.809   1.204   ";}  
    else if(atom=="Se") {return "Se  1.190   90.600  4.205   0.291   14.000  2.764   0.335   0.700   6.428   4.131   1.224   ";}  
    else if(atom=="Br") {return "Br  1.192   180.000 4.189   0.251   15.000  2.519   0.000   0.700   7.790   4.425   1.141   ";}  
    else if(atom=="Kr") {return "Kr  1.147   90.000  4.141   0.220   16.000  0.452   0.000   0.700   8.505   5.715   2.270   ";}  
    else if(atom=="Rb") {return "Rb  2.260   180.000 4.114   0.040   12.000  1.592   0.000   0.200   2.331   1.846   2.770   ";}  
    else if(atom=="Sr") {return "Sr  2.052   90.000  3.641   0.235   12.000  2.449   0.000   0.200   3.024   2.440   2.415   ";}  
    else if(atom=="Y")  {return "Y   1.698   109.470 3.345   0.072   12.000  3.257   0.000   0.200   3.830   2.810   1.998   ";}  
    else if(atom=="Zr") {return "Zr  1.564   109.470 3.124   0.069   12.000  3.667   0.000   0.200   3.400   3.550   1.758   ";}  
    else if(atom=="Nb") {return "Nb  1.473   109.470 3.165   0.059   12.000  3.618   0.000   0.200   3.550   3.380   1.603   ";}  
    else if(atom=="Mo") {return "Mo  1.467   90.000  3.052   0.056   12.000  3.400   0.000   0.200   3.465   3.755   1.530   ";}  
    else if(atom=="Tc") {return "Tc  1.322   90.000  2.998   0.048   12.000  3.400   0.000   0.200   3.290   3.990   1.500   ";}  
    else if(atom=="Ru") {return "Ru  1.478   90.000  2.963   0.056   12.000  3.400   0.000   0.200   3.575   4.015   1.500   ";}  
    else if(atom=="Rh") {return "Rh  1.332   90.000  2.929   0.053   12.000  3.500   0.000   0.200   3.975   4.005   1.509   ";}  
    else if(atom=="Pd") {return "Pd  1.338   90.000  2.899   0.048   12.000  3.210   0.000   0.200   4.320   4.000   1.544   ";}  
    else if(atom=="Ag") {return "Ag  1.386   180.000 3.148   0.036   12.000  1.956   0.000   0.200   4.436   3.134   1.622   ";}  
    else if(atom=="Cd") {return "Cd  1.403   109.470 2.848   0.228   12.000  1.650   0.000   0.200   5.034   3.957   1.600   ";}  
    else if(atom=="In") {return "In  1.459   109.470 4.463   0.599   11.000  2.070   0.000   0.200   3.506   2.896   1.404   ";}  
    else if(atom=="Sn") {return "Sn  1.398   109.470 4.392   0.567   12.000  2.961   0.199   0.200   3.987   3.124   1.354   ";}  
    else if(atom=="Sb") {return "Sb  1.407   91.600  4.420   0.449   13.000  2.704   1.100   0.200   4.899   3.342   1.404   ";}  
    else if(atom=="Te") {return "Te  1.386   90.250  4.470   0.398   14.000  2.882   0.300   0.200   5.816   3.526   1.380   ";}  
    else if(atom=="I")  {return "I   1.382   180.000 4.500   0.339   15.000  2.650   0.000   0.200   6.822   3.762   1.333   ";}  
    else if(atom=="Xe") {return "Xe  1.267   90.000  4.404   0.332   12.000  0.556   0.000   0.200   7.595   4.975   2.459   ";}  
    else if(atom=="Cs") {return "Cs  2.570   180.000 4.517   0.045   12.000  1.573   0.000   0.100   2.183   1.711   2.984   ";}  
    else if(atom=="Ba") {return "Ba  2.277   90.000  3.703   0.364   12.000  2.727   0.000   0.100   2.814   2.396   2.442   ";}  
    else if(atom=="La") {return "La  1.943   109.470 3.522   0.017   12.000  3.300   0.000   0.100   2.836   2.742   2.071   ";}  
    else if(atom=="Ce") {return "Ce  1.841   90.000  3.556   0.013   12.000  3.300   0.000   0.100   2.774   2.692   1.925   ";}  
    else if(atom=="Pr") {return "Pr  1.823   90.000  3.606   0.010   12.000  3.300   0.000   0.100   2.858   2.564   2.007   ";}  
    else if(atom=="Nd") {return "Nd  1.816   90.000  3.575   0.010   12.000  3.300   0.000   0.100   2.869   2.621   2.007   ";}  
    else if(atom=="Pm") {return "Pm  1.801   90.000  3.547   0.009   12.000  3.300   0.000   0.100   2.881   2.673   2.000   ";}  
    else if(atom=="Sm") {return "Sm  1.780   90.000  3.520   0.008   12.000  3.300   0.000   0.100   2.912   2.720   1.978   ";}  
    else if(atom=="Eu") {return "Eu  1.771   90.000  3.493   0.008   12.000  3.300   0.000   0.100   2.879   2.788   2.227   ";}  
    else if(atom=="Gd") {return "Gd  1.735   90.000  3.368   0.009   12.000  3.300   0.000   0.100   3.167   2.975   1.968   ";}  
    else if(atom=="Tb") {return "Tb  1.732   90.000  3.451   0.007   12.000  3.300   0.000   0.100   3.018   2.834   1.954   ";}  
    else if(atom=="Dy") {return "Dy  1.710   90.000  3.428   0.007   12.000  3.300   0.000   0.100   3.056   2.872   1.934   ";}  
    else if(atom=="Ho") {return "Ho  1.696   90.000  3.409   0.007   12.000  3.416   0.000   0.100   3.127   2.891   1.925   ";}  
    else if(atom=="Er") {return "Er  1.673   90.000  3.391   0.007   12.000  3.300   0.000   0.100   3.187   2.915   1.915   ";}  
    else if(atom=="Tm") {return "Tm  1.660   90.000  3.374   0.006   12.000  3.300   0.000   0.100   3.251   2.933   2.000   ";}  
    else if(atom=="Yb") {return "Yb  1.637   90.000  3.355   0.228   12.000  2.618   0.000   0.100   3.289   2.965   2.158   ";}  
    else if(atom=="Lu") {return "Lu  1.671   90.000  3.640   0.041   12.000  3.271   0.000   0.100   2.963   2.463   1.896   ";}  
    else if(atom=="Hf") {return "Hf  1.611   109.470 3.141   0.072   12.000  3.921   0.000   0.100   3.700   3.400   1.759   ";}  
    else if(atom=="Ta") {return "Ta  1.511   109.470 3.170   0.081   12.000  4.075   0.000   0.100   5.100   2.850   1.605   ";}  
    else if(atom=="W")  {return "W   1.392   90.000  3.069   0.067   12.000  3.700   0.000   0.100   4.630   3.310   1.538   ";}  
    else if(atom=="Re") {return "Re  1.372   90.000  2.954   0.066   12.000  3.700   0.000   0.100   3.960   3.920   1.600   ";}  
    else if(atom=="Os") {return "Os  1.372   90.000  3.120   0.037   12.000  3.700   0.000   0.100   5.140   3.630   1.700   ";}  
    else if(atom=="Ir") {return "Ir  1.371   90.000  2.840   0.073   12.000  3.731   0.000   0.100   5.000   4.000   1.866   ";}  
    else if(atom=="Pt") {return "Pt  1.364   90.000  2.754   0.080   12.000  3.382   0.000   0.100   4.790   4.430   1.557   ";}  
    else if(atom=="Au") {return "Au  1.262   90.000  3.293   0.039   12.000  2.625   0.000   0.100   4.894   2.586   1.618   ";}  
    else if(atom=="Hg") {return "Hg  1.340   180.000 2.705   0.385   12.000  1.750   0.000   0.100   6.270   4.160   1.600   ";}  
    else if(atom=="Tl") {return "Tl  1.518   120.000 4.347   0.680   11.000  2.068   0.000   0.100   3.200   2.900   1.530   ";}  
    else if(atom=="Pb") {return "Pb  1.459   109.470 4.297   0.663   12.000  2.846   0.100   0.100   3.900   3.530   1.444   ";}  
    else if(atom=="Bi") {return "Bi  1.512   90.000  4.370   0.518   13.000  2.470   1.000   0.100   4.690   3.740   1.514   ";}  
    else if(atom=="Po") {return "Po  1.500   90.000  4.709   0.325   14.000  2.330   0.300   0.100   4.210   4.210   1.480   ";}  
    else if(atom=="At") {return "At  1.545   180.000 4.750   0.284   15.000  2.240   0.000   0.100   4.750   4.750   1.470   ";}  
    else if(atom=="Rn") {return "Rn  1.420   90.000  4.765   0.248   16.000  0.583   0.000   0.100   5.370   5.370   2.200   ";}  
    else if(atom=="Fr") {return "Fr  2.880   180.000 4.900   0.050   12.000  1.847   0.000   0.000   2.000   2.000   2.300   ";}  
    else if(atom=="Ra") {return "Ra  2.512   90.000  3.677   0.404   12.000  2.920   0.000   0.000   2.843   2.434   2.200   ";}  
    else if(atom=="Ac") {return "Ac  1.983   90.000  3.478   0.033   12.000  3.900   0.000   0.000   2.835   2.835   2.108   ";}  
    else if(atom=="Th") {return "Th  1.721   90.000  3.396   0.026   12.000  4.202   0.000   0.000   3.175   2.905   2.018   ";}  
    else if(atom=="Pa") {return "Pa  1.711   90.000  3.424   0.022   12.000  3.900   0.000   0.000   2.985   2.905   1.800   ";}  
    else if(atom=="U")  {return "U   1.684   90.000  3.395   0.022   12.000  3.900   0.000   0.000   3.341   2.853   1.713   ";}  
    else if(atom=="Np") {return "Np  1.666   90.000  3.424   0.019   12.000  3.900   0.000   0.000   3.549   2.717   1.800   ";}  
    else if(atom=="Pu") {return "Pu  1.657   90.000  3.424   0.016   12.000  3.900   0.000   0.000   3.243   2.819   1.840   ";}  
    else if(atom=="Am") {return "Am  1.660   90.000  3.381   0.014   12.000  3.900   0.000   0.000   2.990   3.004   1.942   ";}  
    else if(atom=="Cm") {return "Cm  1.801   90.000  3.326   0.013   12.000  3.900   0.000   0.000   2.832   3.190   1.900   ";}  
    else if(atom=="Bk") {return "Bk  1.761   90.000  3.339   0.013   12.000  3.900   0.000   0.000   3.194   3.036   1.900   ";}  
    else if(atom=="Cf") {return "Cf  1.750   90.000  3.313   0.013   12.000  3.900   0.000   0.000   3.197   3.101   1.900   ";}  
    else if(atom=="Es") {return "Es  1.724   90.000  3.299   0.012   12.000  3.900   0.000   0.000   3.333   3.089   1.900   ";}  
    else if(atom=="Fm") {return "Fm  1.712   90.000  3.286   0.012   12.000  3.900   0.000   0.000   3.400   3.100   1.900   ";}  
    else if(atom=="Md") {return "Md  1.689   90.000  3.274   0.011   12.000  3.900   0.000   0.000   3.470   3.110   1.900   ";}  
    else if(atom=="No") {return "No  1.679   90.000  3.248   0.011   12.000  3.900   0.000   0.000   3.475   3.175   1.900   ";}  
    else if(atom=="Lw") {return "Lw  1.698   90.000  3.236   0.011   12.000  3.900   0.000   0.000   3.500   3.200   1.900   ";}  
    else {return "Unknown";}                                                                                                      
  }                                                                                                                             
} // namespace pocc                                                                                                             

// ***********************************************************************************************************************************************************
// ***********************************************************************************************************************************************************
// constructor for Class UFFPara                                                                                                
namespace pocc {                                                                                                                
  UFFPara::UFFPara() {
    symbol = "";
    r1 = 0;
    theta0 = 0;
    x1 = 0;
    D1 = 0;
    zeta = 0;
    Z1 = 0;
    Vi = 0;
    Uj = 0;
    Xi = 0;
    hard = 0;
    radius = 0;
  }
} // namespace pocc 

// destructor
namespace pocc {
  void UFFPara::free() {}
  UFFPara::~UFFPara() {free();}
} // namespace pocc 

namespace pocc {
  void UFFPara::GetUFFParameters(string atom) {
    //string strData=ReturnUFFParameters(atom);
    string strData=ReturnUFFParameters(KBIN::VASP_PseudoPotential_CleanName(atom));
    vector<string> UFFData;
    aurostd::string2tokens(strData, UFFData, " ");
    symbol = UFFData.at(0); //Symbol
    r1 = aurostd::string2utype<double>(UFFData.at(1)); 
    theta0 = aurostd::string2utype<double>(UFFData.at(2)); 
    x1 = aurostd::string2utype<double>(UFFData.at(3));  //nonbond distance
    D1 = aurostd::string2utype<double>(UFFData.at(4));  //nonbond energy
    zeta = aurostd::string2utype<double>(UFFData.at(5));  //scale
    Z1 = aurostd::string2utype<double>(UFFData.at(6)); 
    Vi = aurostd::string2utype<double>(UFFData.at(7)); 
    Uj = aurostd::string2utype<double>(UFFData.at(8)); 
    Xi = aurostd::string2utype<double>(UFFData.at(9)); //electronegativity
    hard = aurostd::string2utype<double>(UFFData.at(10)); 
    radius = aurostd::string2utype<double>(UFFData.at(11)); 
  }
} // namespace pocc 

// ***********************************************************************************************************************************************************
// ***********************************************************************************************************************************************************

// ***************************************************************************
// string pocc::ReturnAtomProperties(string& atom)
// ***************************************************************************
// http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)                
// name; symbol; atomic number; atomic mass; atomic radius; pauling electronegativity
namespace pocc {
  string ReturnAtomProperties(string atom) {
    if(atom.compare("H")==0)        {return  "Hydrogen         H     1   1.007940    0.53    2.20     ";}  
    else if(atom.compare("He")==0)  {return  "Helium           He    2   4.002602    0.31    4.16     ";}  
    else if(atom.compare("Li")==0)  {return  "Lithium          Li    3   6.941000    1.67    0.98     ";}  
    else if(atom.compare("Be")==0)  {return  "Beryllium        Be    4   9.012182    1.12    1.57     ";}  
    else if(atom.compare("B")==0)   {return  "Boron            B     5   10.811000   0.87    2.04     ";}  
    else if(atom.compare("C")==0)   {return  "Carbon           C     6   12.010700   0.67    2.55     ";}  
    else if(atom.compare("N")==0)   {return  "Nitrogen         N     7   14.006700   0.56    3.04     ";}  
    else if(atom.compare("O")==0)   {return  "Oxygen           O     8   15.999400   0.48    3.44     ";}  
    else if(atom.compare("F")==0)   {return  "Fluorine         F     9   18.998403   0.42    3.98     ";}  
    else if(atom.compare("Ne")==0)  {return  "Neon             Ne    10  20.179700   0.38    4.79     ";}  
    else if(atom.compare("Na")==0)  {return  "Sodium           Na    11  22.989770   1.90    0.93     ";}  
    else if(atom.compare("Mg")==0)  {return  "Magnesium        Mg    12  24.305000   1.45    1.31     ";}  
    else if(atom.compare("Al")==0)  {return  "Aluminium        Al    13  26.981538   1.18    1.61     ";}  
    else if(atom.compare("Si")==0)  {return  "Silicon          Si    14  28.085500   1.11    1.90     ";}  
    else if(atom.compare("P")==0)   {return  "Phosphorus       P     15  30.973761   0.98    2.19     ";}  
    else if(atom.compare("S")==0)   {return  "Sulphur          S     16  32.065000   0.88    2.58     ";}  
    else if(atom.compare("Cl")==0)  {return  "Chlorine         Cl    17  35.453000   0.79    3.16     ";}  
    else if(atom.compare("Ar")==0)  {return  "Argon            Ar    18  39.948000   0.71    3.24     ";}  
    else if(atom.compare("K")==0)   {return  "Potassium        K     19  39.098300   2.43    0.82     ";}  
    else if(atom.compare("Ca")==0)  {return  "Calcium          Ca    20  40.078000   1.94    1.00     ";}  
    else if(atom.compare("Sc")==0)  {return  "Scandium         Sc    21  44.955910   1.84    1.36     ";}  
    else if(atom.compare("Ti")==0)  {return  "Titanium         Ti    22  47.867000   1.76    1.54     ";}  
    else if(atom.compare("V")==0)   {return  "Vanadium         V     23  50.941500   1.71    1.63     ";}  
    else if(atom.compare("Cr")==0)  {return  "Chromium         Cr    24  51.996100   1.66    1.66     ";}  
    else if(atom.compare("Mn")==0)  {return  "Manganese        Mn    25  54.938049   1.61    1.55     ";}  
    else if(atom.compare("Fe")==0)  {return  "Iron             Fe    26  55.845000   1.56    1.83     ";}  
    else if(atom.compare("Co")==0)  {return  "Cobalt           Co    27  58.933200   1.52    1.88     ";}  
    else if(atom.compare("Ni")==0)  {return  "Nickel           Ni    28  58.693400   1.49    1.91     ";}  
    else if(atom.compare("Cu")==0)  {return  "Copper           Cu    29  63.546000   1.45    1.90     ";}  
    else if(atom.compare("Zn")==0)  {return  "Zinc             Zn    30  65.390000   1.42    1.65     ";}  
    else if(atom.compare("Ga")==0)  {return  "Gallium          Ga    31  69.723000   1.36    1.81     ";}  
    else if(atom.compare("Ge")==0)  {return  "Germanium        Ge    32  72.640000   1.25    2.01     ";}  
    else if(atom.compare("As")==0)  {return  "Arsenic          As    33  74.921600   1.14    2.18     ";}  
    else if(atom.compare("Se")==0)  {return  "Selenium         Se    34  78.960000   1.03    2.55     ";}  
    else if(atom.compare("Br")==0)  {return  "Bromine          Br    35  79.904000   0.94    2.96     ";}  
    else if(atom.compare("Kr")==0)  {return  "Krypton          Kr    36  83.800000   0.88    3.00     ";}  
    else if(atom.compare("Rb")==0)  {return  "Rubidium         Rb    37  85.467800   2.65    0.82     ";}  
    else if(atom.compare("Sr")==0)  {return  "Strontium        Sr    38  87.620000   2.19    0.95     ";}  
    else if(atom.compare("Y")==0)   {return  "Yttrium          Y     39  88.905850   2.12    1.22     ";}  
    else if(atom.compare("Zr")==0)  {return  "Zirconium        Zr    40  91.224000   2.06    1.33     ";}  
    else if(atom.compare("Nb")==0)  {return  "Niobium          Nb    41  92.906380   1.98    1.60     ";}  
    else if(atom.compare("Mo")==0)  {return  "Molybdenum       Mo    42  95.940000   1.90    2.16     ";}  
    else if(atom.compare("Tc")==0)  {return  "Technetium       Mo    43  98.000000   1.83    1.90     ";}  
    else if(atom.compare("Ru")==0)  {return  "Ruthenium        Tc    44  101.070000  1.78    2.20     ";}  
    else if(atom.compare("Rh")==0)  {return  "Rhodium          Ru    45  102.905500  1.73    2.28     ";}  
    else if(atom.compare("Pd")==0)  {return  "Palladium        Pd    46  106.420000  1.69    2.20     ";}  
    else if(atom.compare("Ag")==0)  {return  "Silver           Ag    47  107.868200  1.65    1.93     ";}  
    else if(atom.compare("Cd")==0)  {return  "Cadmium          Cd    48  112.411000  1.61    1.69     ";}  
    else if(atom.compare("In")==0)  {return  "Indium           In    49  114.818000  1.56    1.78     ";}  
    else if(atom.compare("Sn")==0)  {return  "Tin              Sn    50  118.710000  1.45    1.96     ";}  
    else if(atom.compare("Sb")==0)  {return  "Antimony         Sb    51  121.760000  1.33    2.05     ";}  
    else if(atom.compare("Te")==0)  {return  "Tellurium        Te    52  127.600000  1.23    2.10     ";}  
    else if(atom.compare("I")==0)   {return  "Iodine           I     53  126.904470  1.15    2.66     ";}  
    else if(atom.compare("Xe")==0)  {return  "Xenon            Xe    54  131.293000  1.08    2.60     ";}  
    else if(atom.compare("Cs")==0)  {return  "Cesium           Cs    55  132.905450  2.98    0.79     ";}  
    else if(atom.compare("Ba")==0)  {return  "Barium           Ba    56  137.327000  2.53    0.89     ";}  
    else if(atom.compare("La")==0)  {return  "Lanthanium       La    57  138.905500  1.95    1.10     ";}   
    else if(atom.compare("Ce")==0)  {return  "Cerium           Ce    58  140.116000  1.85    1.12     ";}  
    else if(atom.compare("Pr")==0)  {return  "Praseodymium     Pr    59  140.907650  2.47    1.13     ";}  
    else if(atom.compare("Nd")==0)  {return  "Neodymium        Nd    60  144.240000  2.06    1.14     ";}   
    else if(atom.compare("Pm")==0)  {return  "Promethium       Pm    61  145.000000  2.05    1.07     ";}  
    else if(atom.compare("Sm")==0)  {return  "Samarium         Sm    62  150.360000  2.38    1.17     ";}   
    else if(atom.compare("Eu")==0)  {return  "Europium         Eu    63  151.964000  2.31    1.01     ";}    
    else if(atom.compare("Gd")==0)  {return  "Gadolinium       Gd    64  157.250000  2.33    1.20     ";}   
    else if(atom.compare("Tb")==0)  {return  "Terbium          Tb    65  158.925340  2.25    1.10     ";}  
    else if(atom.compare("Dy")==0)  {return  "Dysprosium       Dy    66  162.500000  2.28    1.22     ";}   
    else if(atom.compare("Ho")==0)  {return  "Holmium          Ho    67  164.930320  2.26    1.23     ";}  
    else if(atom.compare("Er")==0)  {return  "Erbium           Er    68  167.259000  2.26    1.24     ";}   
    else if(atom.compare("Tm")==0)  {return  "Thulium          Tm    69  168.934210  2.22    1.25     ";}   
    else if(atom.compare("Yb")==0)  {return  "Ytterbium        Yb    70  173.040000  2.22    1.06     ";}   
    else if(atom.compare("Lu")==0)  {return  "Lutetium         Lu    71  174.967000  2.17    1.27     ";}   
    else if(atom.compare("Hf")==0)  {return  "Hafnium          Hf    72  178.490000  2.08    1.30     ";}   
    else if(atom.compare("Ta")==0)  {return  "Tantalum         Ta    73  180.947900  2.00    1.50     ";}   
    else if(atom.compare("W")==0)   {return  "Tungsten         W     74  183.840000  1.93    2.36     ";}   
    else if(atom.compare("Re")==0)  {return  "Rhenium          Re    75  186.207000  1.88    1.90     ";}   
    else if(atom.compare("Os")==0)  {return  "Osmium           Os    76  190.230000  1.85    2.20     ";}   
    else if(atom.compare("Ir")==0)  {return  "Iridium          Ir    77  192.217000  1.80    2.20     ";}   
    else if(atom.compare("Pt")==0)  {return  "Platinum         Pt    78  195.078000  1.77    2.28     ";}   
    else if(atom.compare("Au")==0)  {return  "Gold             Au    79  196.966550  1.74    2.54     ";}   
    else if(atom.compare("Hg")==0)  {return  "Mercury          Hg    80  200.590000  1.71    2.00     ";}   
    else if(atom.compare("Tl")==0)  {return  "Thallium         Tl    81  204.383300  1.56    1.62     ";}   
    else if(atom.compare("Pb")==0)  {return  "Lead             Pb    82  207.200000  1.54    2.33     ";} 
    else if(atom.compare("Bi")==0)  {return  "Bismuth          Bi    83  208.980380  1.43    2.02     ";}   
    else if(atom.compare("Po")==0)  {return  "Polonium         Po    84  209.000000  1.35    2.00     ";}   
    else if(atom.compare("At")==0)  {return  "Astatine         At    85  210.000000  1.27    2.20     ";}   
    else if(atom.compare("Rn")==0)  {return  "Radon            Rn    86  222.000000  1.20    2.59     ";}   
    else if(atom.compare("Fr")==0)  {return  "Francium         Fr    87  223.000000  2.18    0.70     ";}   
    else if(atom.compare("Ra")==0)  {return  "Radium           Ra    88  226.000000  2.40    0.90     ";}   
    else if(atom.compare("Ac")==0)  {return  "Actinium         Ac    89  227.000000  2.20    1.10     ";}   
    else if(atom.compare("Th")==0)  {return  "Thorium          Th    90  232.038100  1.65    1.30     ";}   
    else if(atom.compare("Pa")==0)  {return  "Protoactinium    Pa    91  231.035880  1.65    1.50     ";}   
    else if(atom.compare("U")==0)   {return  "Uranium          U     92  238.028910  1.42    1.38     ";}   
    //else {return "Hydrogen         H     1   1.007940   0.53    2.20     ";}                                        
    else {return "Unknown";}                                        
  }                                 

  // ***********************************************************************************************************************************************************
  // Class for Atom                 
  Atom::Atom() {                     
    name = "";                      
    symbol = "";                    
    number =0;                      
    mass = 0 ;                      
    radius = 0;                     
    Xi = 0;
  }

  void Atom::free() {}
  Atom::~Atom() {free(); }

  void Atom::GetAtomicProperties(string atom) {
    //string strData=pocc::ReturnAtomProperties(atom);
    string strData=pocc::ReturnAtomProperties(KBIN::VASP_PseudoPotential_CleanName(atom));
    vector<string> atomData;
    aurostd::string2tokens(strData, atomData, " ");
    name = atomData.at(0);
    symbol = atomData.at(1);
    number = aurostd::string2utype<int>(atomData.at(2));
    mass = aurostd::string2utype<double>(atomData.at(3));
    radius = aurostd::string2utype<double>(atomData.at(4));
    Xi = aurostd::string2utype<double>(atomData.at(5)); //pauling electronegativity, Greek chi
  }
} // namespace pocc 

// ***********************************************************************************************************************************************************

#endif     
// ***************************************************************************
// *                                                                         *
// *              AFlow KESONG YANG  Duke University 2010-2011               *
// *                                                                         *
// ***************************************************************************
