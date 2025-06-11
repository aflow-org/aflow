// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// **************************************************************************
// *                                                                        *
// * GNDSTATE STEFANO CURTAROLO MASSACHUSETTS INSTITUTE OF TECHNOLOGY 2003  *
// *                                                                        *
// **************************************************************************
//


  
// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// DEFINITIONS !!!

void FixNameNumbersHQTC(void);
void FixNameNumbersGUS(void);

// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// BOOT !!!

// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// EVENTS DEFINITIONS

double defgamma(double e1,double e2) {
  if (e1<e2) return 1.0;
  else return 0.0;
}

template <class utype> bool iszero(const utype& x) {
#define iszero_threshold 0.00001
  if (abs(x)<iszero_threshold) return TRUE;
  else return FALSE;
}


// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// RANK FUNCTIONAL FUNTIONS
template<class utype> double SpearmanRankOrderVV(const int, const xvector<utype>&, const xvector<utype>&);

template<class utype>
double SpearmanDistanceVV(const int n, const xvector<utype>& x, const xvector<utype>& y) {
  if(n<1) return 0;
  double D=0;
  for (register int i=1;i<=n;i++) {
    for(register int j=1;j<=n;j++) {
      if(x[i]==y[j]) {
        D+=(double) (i-j)*(i-j);
        //      cout <<  x[i] << "  " << i << " " << j << " " << D << endl;
        break;
      }
    }
  }
  return (double) D/n;
}


template<class utype>
double SpearmanRankOrderVV(const int n, const xvector<utype>& x, const xvector<utype>& y) {
  if(n<=1) return 0;
  double D=0;
  for (register int i=1;i<=n;i++) {
    for(register int j=1;j<=n;j++) {
      if(x[i]==y[j]) {
        D+=(double) (i-j)*(i-j);
        //      cout <<  x[i] << "  " << i << " " << j << " " << D << endl;
        break;
      }
    }
  }
  //  cout << D << endl;
  return (double) 1.0 - 6*D/(n*n*n-n);
}


// **************************************************************************
// **************************************************************************
// **************************************************************************
// LIBRARY STUFF


void write_Library_energy_per_atom(const APENNSY_Parameters &params) {
  for(uint i=1;i<=params.Nstructures;i++) {
    for(uint j=1;j<=params.Nalloys;j++)
      printf("Elibrary(%i,%i)=%15.12f; ",i,j,params.Library[i][j].energy_per_atom);
    cout << endl;
  }
}

void write_Library_simulation_time(const APENNSY_Parameters &params) {
  for(uint i=1;i<=params.Nstructures;i++) {
    for(uint j=1;j<=params.Nalloys;j++)
      printf("Library[%i][%i].simulation_time=%15.12f; ",i,j,params.Library[i][j].simulation_time);
    cout << endl;
  }
}



double ProbabilityEformationPositive(const xmatrix<double>& Ef) {
  int i,j,Ntot=0,Npos=0,Nneg=0;
  for(i=Ef.lrows;i<=Ef.urows;i++) {
    for(j=Ef.lcols;j<=Ef.ucols;j++) {
      //  if(i!= A1_HTQC_1 && i!= A1_HTQC_2 && i!= A2_HTQC_1 && i!= A2_HTQC_2 && i!= A3_HTQC_1 && i!= A3_HTQC_2 &&
      //  i!= A4_HTQC_1 && i!= A4_HTQC_2 && i!= A6_HTQC_1 && i!= A6_HTQC_2 && i!= A7_HTQC_1 && i!= A7_HTQC_2)
      {
	Ntot++;
	if(Ef(i,j)<0) Nneg++;
	if(Ef(i,j)>0) Npos++;
      }
    }
  }
  cout << "Ntot=" << Ntot << " Npos=" << Npos << " Nneg=" << Nneg << endl;
  return (double) Npos/Ntot;
}

// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************

string SubStringRemove(string str_original, string subremove) {
  string::size_type sub_size1,sub_size2;
  string subSpre,subSpost,strOUT;
  strOUT=str_original;
  sub_size1=strOUT.find(subremove)+subremove.length();
  sub_size2=(strOUT.substr(sub_size1)).find("\n");
  subSpost=strOUT.substr(sub_size1,sub_size2);
  sub_size1=0;
  sub_size2=(strOUT.substr(sub_size1)).find(subremove);
  subSpre=strOUT.substr(sub_size1,sub_size2);
  strOUT=subSpre+subSpost;
  return strOUT;
}

\
// **************************************************************************
// **************************************************************************
// **************************************************************************
// *                                                                        *
// * GNDSTATE STEFANO CURTAROLO MASSACHUSETTS INSTITUTE OF TECHNOLOGY 2003  *
// *                                                                        *
// **************************************************************************
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
