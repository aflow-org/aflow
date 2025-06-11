// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - Nov 2007

#include "aflow.h"

#define cdebug cerr
#define _EPS_ 0.001
#define _EPS_roundoff_ 1.0e-8
#define _CAGEROUGHNESS_ 1.075

// #define radius  6
// #define coordination_position  8
// #define _cages_position_   9
// #define _cages_irrtype_    7

#define _oss_precision_aflow_defects_   14
#define _oss_short_precision_aflow_defects_ 7

pthread_mutex_t mutex_DEFECTS=PTHREAD_MUTEX_INITIALIZER;
#define _PTHREAD_FLUSH_TIME_ 1

#define string_species "X"

// ****************************************************************************************************
// CLASS acage

// constructor
acage::acage(void) {
  origin_fpos.clear();
  origin_cpos.clear();
  radius=0;
  coordination_position=0;
  cages_position=0;
  cages_irrtype=0;
  atoms.clear();
}
// destructor
acage::~acage() { free(); }

void acage::clear(void) {
  origin_fpos.clear();
  origin_cpos.clear();
  radius=0;
  coordination_position=0;
  cages_position=0;
  cages_irrtype=0;
  atoms.clear();
}

void acage::free() {}

void acage::copy(const acage& b) {
  origin_fpos=b.origin_fpos;
  origin_cpos=b.origin_cpos;
  radius=b.radius;
  coordination_position=b.coordination_position;
  cages_position=b.cages_position;
  cages_irrtype=b.cages_irrtype;
  atoms.clear(); for(uint i=0;i<b.atoms.size();i++) atoms.push_back(b.atoms.at(i));
}

acage::acage(const acage& b) {
  //  free(); // *this=b;
  copy(b);
}


// ***************************************************************************
// prototypes

// bool GetSphereFromFourPoints(xvector<double>& orig,double& radius,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3,const xvector<double>& v4);
// bool GetCircumCircleeFromThreePoints(xvector<double>& orig,double& radius,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
// bool GetCircleFromTwoPoints(xvector<double>& orig,double& radius,const xvector<double>& v1,const xvector<double>& v2);
// bool FPositionInsideCell(const xvector<double>& r);
// bool EmptySphere(const deque<_atom>& grid_atoms,const xvector<double>& origin_cpos,const double& radius);
// bool EmptySphere(const xstructure& str,const xvector<double>& origin_cpos,const double& radius);
// uint CoordinationPoint(const deque<_atom>& atoms,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin,const double& rmax);
// uint CoordinationPoint(const deque<_atom>& atoms,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin);
// uint CoordinationPoint(const xstructure& str,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin,const double& rmax);
// uint CoordinationPoint(const xstructure& str,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin);
// bool AddCageToCages(const xstructure& str,const xvector<double>& origin_cpos,const xvector<double>& origin_fpos,const double& radius,
// 		    const int& cage_points_type,const double& roughness,
// 		    vector<acage>& cages,vector<acage>& cagesX,
// 		    const bool& oss1write,ostream& oss1, const bool& oss2write,ostream& oss2);
// uint GetCages4(const xstructure& str,const double& roughness,vector<acage>& cages,vector<acage>& cages4,
// 	      const bool& oss1write,ostream& oss1, const bool& oss2write,ostream& oss2);
// uint GetCages3(const xstructure& str,const double& roughness,vector<acage>& cages,vector<acage>& cages3,
// 	      const bool& oss1write,ostream& oss1, const bool& oss2write,ostream& oss2);
// uint GetCages2(const xstructure& str,const double& roughness,vector<acage>& cages,vector<acage>& cages2,
// 	      const bool& oss1write,ostream& oss1, const bool& oss2write,ostream& oss2);
// bool GetCages(const xstructure& _str,_aflags& aflags,
// 	      vector<acage>& cagesirreducible,vector<acage>& cagesreducible,vector<acage>& cages4,
// 	      vector<acage>& cages3,vector<acage>& cages2,const double& _roughness,const bool& osswrite,ostream& oss);

// ***************************************************************************
// implementations
bool GetSphereFromFourPoints(xvector<double>& orig,double& radius,
    const xvector<double>& v1,const xvector<double>& v2,
    const xvector<double>& v3,const xvector<double>& v4) {
  double eps=_EPS_;
  xmatrix<double> A(5,5);
  A[1][1]=1.0;A[1][2]=1.0;A[1][3]=1.0;A[1][4]=0;A[1][5]=1.0;
  A[2][1]=v1[1]*v1[1]+v1[2]*v1[2]+v1[3]*v1[3];A[2][2]=v1[1];A[2][3]=v1[2];A[2][4]=v1[3];A[2][5]=1;
  A[3][1]=v2[1]*v2[1]+v2[2]*v2[2]+v2[3]*v2[3];A[3][2]=v2[1];A[3][3]=v2[2];A[3][4]=v2[3];A[3][5]=1;
  A[4][1]=v3[1]*v3[1]+v3[2]*v3[2]+v3[3]*v3[3];A[4][2]=v3[1];A[4][3]=v3[2];A[4][4]=v3[3];A[4][5]=1;
  A[5][1]=v4[1]*v4[1]+v4[2]*v4[2]+v4[3]*v4[3];A[5][2]=v4[1];A[5][3]=v4[2];A[5][4]=v4[3];A[5][5]=1;

  double M11=minordeterminant(A,1,1);
  if(abs(M11)<eps) return FALSE;
  double M12=minordeterminant(A,1,2);
  double M13=minordeterminant(A,1,3);
  double M14=minordeterminant(A,1,4);
  double M15=minordeterminant(A,1,5);

  orig[1]= (M12/M11)/2.0;
  orig[2]=-(M13/M11)/2.0;
  orig[3]= (M14/M11)/2.0;
  radius=sqrt(orig[1]*orig[1]+orig[2]*orig[2]+orig[3]*orig[3]-M15/M11);

  double d1=modulus(v1-orig);if(d1<radius-eps || d1>radius+eps) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Sphere Error [d1]",_RUNTIME_ERROR_);} //CO20200624
  double d2=modulus(v2-orig);if(d2<radius-eps || d2>radius+eps) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Sphere Error [d2]",_RUNTIME_ERROR_);} //CO20200624
  double d3=modulus(v3-orig);if(d3<radius-eps || d3>radius+eps) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Sphere Error [d3]",_RUNTIME_ERROR_);} //CO20200624
  double d4=modulus(v4-orig);if(d4<radius-eps || d4>radius+eps) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Sphere Error [d4]",_RUNTIME_ERROR_);} //CO20200624
  return TRUE;
}

bool GetCircumCircleeFromThreePoints(xvector<double>& orig,double& radius,
    const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3) {
  double eps=_EPS_;
  xmatrix<double> A(4,4);
  A[1][1]=0;                                   A[1][2]=0;     A[1][3]=    0; A[1][4]=    0;
  A[2][1]=v1[1]*v1[1]+v1[2]*v1[2]+v1[3]*v1[3]; A[2][2]=v1[1]; A[2][3]=v1[2]; A[2][4]=v1[3];
  A[3][1]=v2[1]*v2[1]+v2[2]*v2[2]+v2[3]*v2[3]; A[3][2]=v2[1]; A[3][3]=v2[2]; A[3][4]=v2[3];
  A[4][1]=v3[1]*v3[1]+v3[2]*v3[2]+v3[3]*v3[3]; A[4][2]=v3[1]; A[4][3]=v3[2]; A[4][4]=v3[3];

  double M11= minordeterminant(A,1,1);
  if(abs(M11)<eps) return FALSE;
  double M12= minordeterminant(A,1,2);
  double M13= minordeterminant(A,1,3);
  double M14= minordeterminant(A,1,4);

  orig[1]= (M12/M11)/2.0;
  orig[2]=-(M13/M11)/2.0;
  orig[3]= (M14/M11)/2.0;
  radius=modulus(v1-orig);
  if(abs(modulus(v1-orig)-modulus(v2-orig))<eps && abs(modulus(v2-orig)-modulus(v3-orig))<eps) return TRUE;
  stringstream message;
  message << modulus(v1-orig) << " " << modulus(v2-orig) << " " << modulus(v3-orig) << " ";
  return FALSE; //CO20200624 - previously an exit
}

bool GetCircleFromTwoPoints(xvector<double>& orig,double& radius,const xvector<double>& v1,const xvector<double>& v2) {
  double eps=_EPS_;
  if(modulus(v2-v1)<eps) return FALSE;
  orig=(v2+v1)/2.0;
  radius=modulus(v1-v2)/2;
  return TRUE;
}

bool FPositionInsideCell(const xvector<double>& r) {
  if(r[1]>=0.0 && r[1]<1.0 && r[2]>=0.0 && r[2]<1.0 && r[3]>=0.0 && r[3]<1.0) return TRUE;
  else return FALSE;
}

bool EmptySphere(const deque<_atom>& grid_atoms,const xvector<double>& origin_cpos,const double& radius) {
  double eps=_EPS_;
  bool empty=TRUE;
  //  for(uint iat=0;iat<grid_atoms.size()&&empty==TRUE;iat++)
  for(uint iat=0;iat<grid_atoms.size();iat++) {
    //    if(abs(grid_atoms.at(iat).cpos[1]-1.728644884961154)<0.01)
    //     //   0.000000000099805   2.261436255651232
    //     cerr << grid_atoms.at(iat).cpos << "      " << grid_atoms.at(iat).fpos << endl;
    if(modulus(grid_atoms.at(iat).cpos-origin_cpos)<radius-5*eps)
      empty=FALSE;
  }
  return empty;
}

bool EmptySphere(const xstructure& str,const xvector<double>& origin_cpos,const double& radius) {
  if(!str.grid_atoms_calculated) {throw aurostd::xerror(__AFLOW_FILE__,"EmptySphere():","str.grid_atoms must be calculated",_INPUT_MISSING_);} //CO20200624
  return EmptySphere(str.grid_atoms,origin_cpos,radius);
}

uint CoordinationPoint(const deque<_atom>& atoms,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin,const double& rmax) {
  uint out=0;
  ratoms.clear();
  double eps=_EPS_,dist;
  for(uint iat=0;iat<atoms.size();iat++) {
    dist=modulus(atoms.at(iat).cpos-point);
    if(dist>=rmin-eps && dist<=rmax+eps) {
      out++; ratoms.push_back(atoms.at(iat));
    }
  }
  return out;
}

uint CoordinationPoint(const deque<_atom>& atoms,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin) {
  return CoordinationPoint(atoms,ratoms,point,rmin/_CAGEROUGHNESS_,rmin*_CAGEROUGHNESS_);
}

uint CoordinationPoint(const xstructure& str,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin,const double& rmax) {
  if(!str.grid_atoms_calculated) {throw aurostd::xerror(__AFLOW_FILE__,"CoordinationPoint():","str.grid_atoms must be calculated",_INPUT_MISSING_);} //CO20200624
  return CoordinationPoint(str.grid_atoms,ratoms,point,rmin,rmax);
}

uint CoordinationPoint(const xstructure& str,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin) {
  return CoordinationPoint(str.grid_atoms,ratoms,point,rmin/_CAGEROUGHNESS_,rmin*_CAGEROUGHNESS_);
}


bool AddCageToCages(const xstructure& str,const xvector<double>& origin_cpos,const xvector<double>& origin_fpos,const double& radius,
    const int& cage_points_type,const double& roughness,
    vector<acage>& cages,vector<acage>& cagesX,
    const bool& oss1write,ostream& oss1, const bool& oss2write,ostream& oss2,int ithread) {
  acage cage;
  double eps=_EPS_;
  bool found=FALSE;
  string species=string_species;
  for(uint ig=0;ig<cages.size()&&!found;ig++)
    found=(modulus(origin_cpos-cages.at(ig).origin_cpos)<eps);
  if(found==FALSE) {
    cage.origin_fpos=origin_fpos;
    cage.origin_cpos=origin_cpos;
    cage.radius=radius;
    cage.coordination_position=CoordinationPoint(str,cage.atoms,origin_cpos,radius/roughness,radius*roughness);
    cage.cages_position=cage_points_type;
    cage.cages_irrtype=0;  // to be fix by the symmetry reduction
    cages.push_back(cage);
    cagesX.push_back(cage);
    // if(oss1write) oss1 << i << " " << j << " " << k << " " << l << " " << m << endl;
    int ig=cagesX.size()-1;
    // if(oss1write) oss1 << ig << " ";
    if(oss1write) {oss1 << "  "; for(int i=1;i<=3;i++) {if(cagesX.at(ig).origin_fpos[i]>=0.0) oss1 << " "; oss1 << " " << cagesX.at(ig).origin_fpos[i];}}
    if(oss2write) {oss2 << "  "; for(int i=1;i<=3;i++) {if(cagesX.at(ig).origin_fpos[i]>=0.0) oss2 << " "; oss2 << " " << cagesX.at(ig).origin_fpos[i];}}
    // if(oss1write) oss1 << "    i=" << ig;
    oss1.precision(_oss_short_precision_aflow_defects_);  oss2.precision(_oss_short_precision_aflow_defects_);
    if(oss1write) oss1 << "  " << aurostd::PaddedPRE(species,2);
    if(oss2write) oss2 << "  " << aurostd::PaddedPRE(species,2);
    if(oss1write) oss1 << "   " << cagesX.at(ig).radius << " ";
    if(oss2write) oss2 << "   " << cagesX.at(ig).radius << " ";
    oss1.precision(_oss_precision_aflow_defects_);  oss2.precision(_oss_precision_aflow_defects_);
    if(oss1write) oss1 << "  -";// << (int) cagesX.at(ig).cages_irrtype;
    if(oss2write) oss2 << "  -";// << (int) cagesX.at(ig).cages_irrtype;
    if(oss1write) oss1 << "  " << (int) cagesX.at(ig).coordination_position;
    if(oss2write) oss2 << "  " << (int) cagesX.at(ig).coordination_position;
    if(oss1write) oss1 << "  " << (int) cagesX.at(ig).cages_position;
    if(oss2write) oss2 << "  " << (int) cagesX.at(ig).cages_position;
    if(oss1write) {oss1 << "  [";for(uint iat=0;iat<cagesX.at(ig).atoms.size();iat++) oss1 << str.species.at(cagesX.at(ig).atoms.at(iat).type) << "(" << cagesX.at(ig).atoms.at(iat).basis << ")" << (iat<cagesX.at(ig).atoms.size()-1?",":""); oss1 << "]"; }
    if(oss2write) {oss2 << "  [";for(uint iat=0;iat<cagesX.at(ig).atoms.size();iat++) oss2 << str.species.at(cagesX.at(ig).atoms.at(iat).type) << "(" << cagesX.at(ig).atoms.at(iat).basis << ")" << (iat<cagesX.at(ig).atoms.size()-1?",":""); oss2 << "]"; }

    // if(oss1write) oss1 << " " << iat1 << "(" << modulus(origin_cpos-v1) << ") ";
    // if(oss1write) oss1 << " " << iat2 << "(" << modulus(origin_cpos-v2) << ") ";
    // if(oss1write) oss1 << " " << iat3 << "(" << modulus(origin_cpos-v3) << ") ";
    // if(oss1write) oss1 << " " << iat4 << "(" << modulus(origin_cpos-v4) << ") ";
    if(ithread>=0) if(oss1write) oss1 << "  ithread=" << ithread;
    if(ithread>=0) if(oss2write) oss2 << "  ithread=" << ithread;
    if(oss1write) oss1 << endl;
    if(oss2write) oss2 << endl;
    return TRUE;
  } else {
    return FALSE;
  }
}

uint GetCages4(const xstructure& str,const double& roughness,vector<acage>& cages,vector<acage>& cages4,
    const bool& oss1write,ostream& oss1, const bool& oss2write,ostream& oss2) {
  xvector<double> v1(3),v2(3),v3(3),v4(3);
  xvector<double> origin_cpos(3),origin_fpos(3);
  double radius;

  cages4.clear();
  bool found;
  uint gridsize=str.grid_atoms.size();
  //  for(uint iat1=0;iat1<gridsize;iat1++)
  for(uint iat1=0;iat1<str.atoms.size();iat1++) {
    v1=str.grid_atoms.at(iat1).cpos;
    for(uint iat2=iat1+1;iat2<gridsize;iat2++) {
      v2=str.grid_atoms.at(iat2).cpos;
      for(uint iat3=iat2+1;iat3<gridsize;iat3++) {
        v3=str.grid_atoms.at(iat3).cpos;
        for(uint iat4=iat3+1;iat4<gridsize;iat4++) {
          v4=str.grid_atoms.at(iat4).cpos;
          found=GetSphereFromFourPoints(origin_cpos,radius,v1,v2,v3,v4);
          if(found) {
            origin_fpos=C2F(str.lattice,origin_cpos);
            //  origin_fpos=str.c2f*origin_cpos;
            // roundoff(origin_fpos,eps);
            if(FPositionInsideCell(origin_fpos))
              //	      cerr << origin_fpos << "  " << origin_cpos << endl;
              if(EmptySphere(str,origin_cpos,radius))
                AddCageToCages(str,origin_cpos,origin_fpos,radius,4,roughness,cages,cages4,oss1write,oss1,oss2write,oss2,-1);
          }
        }
      }
    }
  }
  sort(cages4.begin(),cages4.end(),_isort_acage_radius());
  return cages4.size();
}

typedef struct {
  _aflags  *paflags;
  int      ITHREAD;
  int      THREADS_MAX;
  xstructure *pstr;
  double    *proughness;
  vector<acage> *pcages;
  vector<acage> *pcages2;
  vector<acage> *pcages3;
  vector<acage> *pcages4;
  bool *poss1write;
  ostream *poss1;
  bool *poss2write;
  ofstream *poss2;
  int      itbusy;
} _threaded_GETCAGES_params;


void *_threaded_GetCages4(void *ptr) {
  _threaded_GETCAGES_params* pparams;
  pparams = (_threaded_GETCAGES_params*) ptr;
  (pparams->itbusy)=TRUE;
  AFLOW_PTHREADS::RUNNING++;
  // CODE BEGIN
  xvector<double> v1(3),v2(3),v3(3),v4(3),origin_cpos(3),origin_fpos(3);
  double eps=_EPS_,radius;
  bool found;  
  // pthread_mutex_lock( &mutex_DEFECTS );
  // cerr << "DEBUG: _threaded_GetCages4 " << (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << endl;
  // pthread_mutex_unlock( &mutex_DEFECTS );
  // for(uint iat1=0+(pparams->ITHREAD);iat1<(*pparams->pstr).grid_atoms.size();iat1+=(pparams->THREADS_MAX))
  for(uint iat1=0;iat1<(*pparams->pstr).atoms.size();iat1++) {
    v1=(*pparams->pstr).grid_atoms[iat1].cpos;
    // for(uint iat2=iat1+1;iat2<(*pparams->pstr).grid_atoms.size();iat2++)
    for(uint iat2=iat1+1+(pparams->ITHREAD);iat2<(*pparams->pstr).grid_atoms.size();iat2+=(pparams->THREADS_MAX)) {
      v2=(*pparams->pstr).grid_atoms[iat2].cpos;
      for(uint iat3=iat2+1;iat3<(*pparams->pstr).grid_atoms.size();iat3++) {
        v3=(*pparams->pstr).grid_atoms[iat3].cpos;
        for(uint iat4=iat3+1;iat4<(*pparams->pstr).grid_atoms.size();iat4++) {
          v4=(*pparams->pstr).grid_atoms[iat4].cpos;
          found=GetSphereFromFourPoints(origin_cpos,radius,v1,v2,v3,v4);
          if(found) {
            origin_fpos=C2F((*pparams->pstr).lattice,origin_cpos);
            roundoff(origin_fpos,eps);
            if(FPositionInsideCell(origin_fpos))
              if(EmptySphere((*pparams->pstr),origin_cpos,radius)) {
                pthread_mutex_lock( &mutex_DEFECTS );
                AddCageToCages((*pparams->pstr),origin_cpos,origin_fpos,radius,4,(*pparams->proughness),
                    (*pparams->pcages),(*pparams->pcages4),
                    (*pparams->poss1write),(*pparams->poss1),
                    (*pparams->poss2write),(*pparams->poss2),
                    (pparams->ITHREAD));	    
                pthread_mutex_unlock( &mutex_DEFECTS );
              }
          }
        }
      }
    }
  }
  // CODE END
  (pparams->itbusy)=FALSE;
  AFLOW_PTHREADS::RUNNING--;
  aurostd::Sleep(_PTHREAD_FLUSH_TIME_);
  return NULL;
}

uint GetCages3(const xstructure& str,const double& roughness,vector<acage>& cages,vector<acage>& cages3,
    const bool& oss1write,ostream& oss1, const bool& oss2write,ostream& oss2) {
  xvector<double> v1(3),v2(3),v3(3),origin_cpos(3),origin_fpos(3);
  //  acage cage;
  double eps=_EPS_;
  double radius;
  cages3.clear();
  bool found;
  // for(uint iat1=0;iat1<str.grid_atoms.size();iat1++)
  for(uint iat1=0;iat1<str.atoms.size();iat1++) {
    v1=str.grid_atoms[iat1].cpos;
    for(uint iat2=iat1+1;iat2<str.grid_atoms.size();iat2++) {
      v2=str.grid_atoms[iat2].cpos;
      for(uint iat3=iat2+1;iat3<str.grid_atoms.size();iat3++) {
        v3=str.grid_atoms[iat3].cpos;
        found=GetCircumCircleeFromThreePoints(origin_cpos,radius,v1,v2,v3);
        if(found) {
          // cerr << iat1 << " " << iat2 << " " << iat3 << endl;
          origin_fpos=C2F(str.lattice,origin_cpos);
          roundoff(origin_fpos,eps);
          if(FPositionInsideCell(origin_fpos))
            if(EmptySphere(str,origin_cpos,radius))
              AddCageToCages(str,origin_cpos,origin_fpos,radius,3,roughness,cages,cages3,oss1write,oss1,oss2write,oss2,-1);
        }
      }
    }
  }
  sort(cages3.begin(),cages3.end(),_isort_acage_radius());
  return cages3.size();
}

void *_threaded_GetCages3(void *ptr) {
  _threaded_GETCAGES_params* pparams;
  pparams = (_threaded_GETCAGES_params*) ptr;
  (pparams->itbusy)=TRUE;
  AFLOW_PTHREADS::RUNNING++;
  // CODE BEGIN
  xvector<double> v1(3),v2(3),v3(3),origin_cpos(3),origin_fpos(3);
  //  acage cage;
  double eps=_EPS_,radius;
  bool found;
  // pthread_mutex_lock( &mutex_DEFECTS );
  // cerr << "DEBUG: _threaded_GetCages3 " << (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << endl;
  // pthread_mutex_unlock( &mutex_DEFECTS );
  // for(uint iat1=0+(pparams->ITHREAD);iat1<(*pparams->pstr).grid_atoms.size();iat1+=(pparams->THREADS_MAX))
  for(uint iat1=0;iat1<(*pparams->pstr).atoms.size();iat1++) {
    v1=(*pparams->pstr).grid_atoms[iat1].cpos;
    // for(uint iat2=iat1+1;iat2<(*pparams->pstr).grid_atoms.size();iat2++)
    for(uint iat2=iat1+1+(pparams->ITHREAD);iat2<(*pparams->pstr).grid_atoms.size();iat2+=(pparams->THREADS_MAX)) {
      v2=(*pparams->pstr).grid_atoms[iat2].cpos;
      for(uint iat3=iat2+1;iat3<(*pparams->pstr).grid_atoms.size();iat3++) {
        v3=(*pparams->pstr).grid_atoms[iat3].cpos;
        found=GetCircumCircleeFromThreePoints(origin_cpos,radius,v1,v2,v3);
        if(found) {
          // cerr << iat1 << " " << iat2 << " " << iat3 << endl;
          origin_fpos=C2F((*pparams->pstr).lattice,origin_cpos);
          roundoff(origin_fpos,eps);
          if(FPositionInsideCell(origin_fpos))
            if(EmptySphere((*pparams->pstr),origin_cpos,radius)) {
              pthread_mutex_lock( &mutex_DEFECTS );
              AddCageToCages((*pparams->pstr),origin_cpos,origin_fpos,radius,3,
                  (*pparams->proughness),(*pparams->pcages),(*pparams->pcages3),
                  (*pparams->poss1write),(*pparams->poss1),
                  (*pparams->poss2write),(*pparams->poss2),
                  (pparams->ITHREAD));
              pthread_mutex_unlock( &mutex_DEFECTS );
            }
        }
      }
    }
  }
  // CODE END
  (pparams->itbusy)=FALSE;
  AFLOW_PTHREADS::RUNNING--;
  aurostd::Sleep(_PTHREAD_FLUSH_TIME_);
  return NULL;
}

uint GetCages2(const xstructure& str,const double& roughness,vector<acage>& cages,vector<acage>& cages2,
    const bool& oss1write,ostream& oss1, const bool& oss2write,ostream& oss2) {
  xvector<double> v1(3),v2(3),origin_cpos(3),origin_fpos(3);
  //  acage cage;
  double eps=_EPS_;
  double radius;
  cages2.clear();
  bool found;
  //  for(uint iat1=0;iat1<str.grid_atoms.size();iat1++)
  for(uint iat1=0;iat1<str.atoms.size();iat1++) {
    v1=str.grid_atoms[iat1].cpos;
    for(uint iat2=iat1+1;iat2<str.grid_atoms.size();iat2++) {
      v2=str.grid_atoms[iat2].cpos;
      found=GetCircleFromTwoPoints(origin_cpos,radius,v1,v2);
      if(found) {
        origin_fpos=C2F(str.lattice,origin_cpos);
        roundoff(origin_fpos,eps);
        if(FPositionInsideCell(origin_fpos))
          if(EmptySphere(str,origin_cpos,radius))
            AddCageToCages(str,origin_cpos,origin_fpos,radius,2,roughness,cages,cages2,oss1write,oss1,oss2write,oss2,-1);
      }
    }
  }
  sort(cages2.begin(),cages2.end(),_isort_acage_radius());
  return cages2.size();
}

void *_threaded_GetCages2(void *ptr) {
  _threaded_GETCAGES_params* pparams;
  pparams = (_threaded_GETCAGES_params*) ptr;
  (pparams->itbusy)=TRUE;
  AFLOW_PTHREADS::RUNNING++;
  // CODE BEGIN
  xvector<double> v1(3),v2(3),origin_cpos(3),origin_fpos(3);
  //  acage cage;
  double eps=_EPS_,radius;
  bool found;
  // pthread_mutex_lock( &mutex_DEFECTS );
  // cerr << "DEBUG: _threaded_GetCages2 " << (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << endl;
  // pthread_mutex_unlock( &mutex_DEFECTS );
  // for(uint iat1=0;iat1<(*pparams->pstr).grid_atoms.size();iat1++)
  for(uint iat1=0+(pparams->ITHREAD);iat1<(*pparams->pstr).grid_atoms.size();iat1+=(pparams->THREADS_MAX)) {
    v1=(*pparams->pstr).grid_atoms[iat1].cpos;
    for(uint iat2=iat1+1;iat2<(*pparams->pstr).grid_atoms.size();iat2++) {
      v2=(*pparams->pstr).grid_atoms[iat2].cpos;
      found=GetCircleFromTwoPoints(origin_cpos,radius,v1,v2);
      if(found) {
        origin_fpos=C2F((*pparams->pstr).lattice,origin_cpos);
        roundoff(origin_fpos,eps);
        if(FPositionInsideCell(origin_fpos))
          if(EmptySphere((*pparams->pstr),origin_cpos,radius)) {
            pthread_mutex_lock( &mutex_DEFECTS );
            AddCageToCages((*pparams->pstr),origin_cpos,origin_fpos,radius,2,(*pparams->proughness),
                (*pparams->pcages),(*pparams->pcages2),
                (*pparams->poss1write),(*pparams->poss1),
                (*pparams->poss2write),(*pparams->poss2),
                (pparams->ITHREAD));
            pthread_mutex_unlock( &mutex_DEFECTS );
          }
      }
    }
  }

  // CODE END
  (pparams->itbusy)=FALSE;
  AFLOW_PTHREADS::RUNNING--;
  aurostd::Sleep(_PTHREAD_FLUSH_TIME_);
  return NULL;
}

//ME20210521 - added ofstream write flag (FFFflag) as function input
bool GetCages(const xstructure& _str,_aflags& aflags,
    vector<acage>& cagesirreducible,vector<acage>& cagesreducible,
    vector<acage>& cages4,vector<acage>& cages3,vector<acage>& cages2,
    double _roughness,bool FFFflag, bool osswrite,ostream& oss) {
  // oss << "CAGES BEGIN" << endl;
  bool Krun=TRUE;
  string species=string_species;

  xstructure str(_str);
  // str=BringInCell(_str);
  str=ReScale(BringInCell(_str),1.0);str.neg_scale=FALSE;
  // str.SetCoordinates(_COORDS_CARTESIAN_);

  // if(1) str.ShiftOriginToAtom(0);
  xmatrix<double> lattice(3,3);                         // lattice
  lattice=(str.lattice);                                // lattice
  xvector<double> a1(3),a2(3),a3(3);                    // a1,a2,a3 are the rows of the lattice matrix
  a1=lattice(1);a2=lattice(2);a3=lattice(3);            // a1,a2,a3 are the rows of the lattice matrix
  xvector<double> v1(3),v2(3),v3(3),v4(3);              // vectors for symmetry search
  xvector<double> spherecenter_cpos(3),spherecenter_fpos(3);
  double radius;
  xvector<double> rrr(3),origin_cpos(3),origin_fpos(3);

  double eps=_EPS_;
  double roughness=_roughness;
  double cagesirreducible_rmax,cagesirreducible_rmin;
  double cagesreducible_rmax,cagesreducible_rmin;
  double cage2_rmax,cage2_rmin;
  double cage3_rmax,cage3_rmin;
  double cage4_rmax,cage4_rmin;

  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(_oss_precision_aflow_defects_);
  ofstream FFF;
  string FileNameCAGE;
  if(FFFflag) FileNameCAGE=DEFAULT_AFLOW_ICAGES_OUT;
  else FileNameCAGE="/dev/null";
  FFF.open(FileNameCAGE.c_str(),std::ios::out);
  FFF.setf(std::ios::fixed,std::ios::floatfield);
  FFF.precision(_oss_precision_aflow_defects_);

  radius=max(modulus(a1),modulus(a2),modulus(a3));
  xvector<int> dims(3);
  dims=LatticeDimensionSphere(lattice,radius);
  dims[1]=1;dims[2]=1;dims[3]=1;
  //  oss << dims << endl;

  acage cage;

  string banner="-------------------------------------------------------------------------------------------";
  //  oss << "------------------------------------------------------------------------------------" << endl;
  bool PFSWRITE=FFFflag;  // ME20210521
  ofstream FileDevNull("/dev/null");
  aflags.QUIET=TRUE;
  bool OSSWRITE=FALSE;// if(!OSSWRITE) PFSWRITE=FALSE;
  oss << banner << endl; // ---------------------------------
  FFF << banner << endl; // ---------------------------------
  oss << aflow::Banner("BANNER_TINY") << endl;
  FFF << aflow::Banner("BANNER_TINY") << endl;
  oss << banner << endl; // ---------------------------------
  FFF << banner << endl; // ---------------------------------
  oss << "CAGES CALCULATION - Stefano Curtarolo   [algorithm: DOI:10.1103/PhysRevB.79.134203]" << endl;
  FFF << "CAGES CALCULATION - Stefano Curtarolo   [algorithm: DOI:10.1103/PhysRevB.79.134203]" << endl;
  if(!AFLOW_PTHREADS::FLAG) oss << "NO PTHREADS" << endl;
  if(!AFLOW_PTHREADS::FLAG) FFF << "NO PTHREADS" << endl;
  if(AFLOW_PTHREADS::FLAG) oss << "THREADS - AFLOW_PTHREADS::MAX_PTHREADS=" << AFLOW_PTHREADS::MAX_PTHREADS << endl;
  if(AFLOW_PTHREADS::FLAG) FFF << "THREADS - AFLOW_PTHREADS::MAX_PTHREADS=" << AFLOW_PTHREADS::MAX_PTHREADS << endl;
  oss << banner << endl; // ---------------------------------
  FFF << banner << endl; // ---------------------------------
  str.sgroup_radius=1.5*RadiusSphereLattice(str.lattice);                                   //CO20171024 - new sym framework
  _kflags kflags; pflow::defaultKFlags4SymCalc(kflags,true);                                //CO20171024 - new sym framework
  pflow::defaultKFlags4SymWrite(kflags,PFSWRITE); kflags.KBIN_SYMMETRY_SGROUP_WRITE=false;  //CO20171024 - new sym framework
  pflow::PerformFullSymmetry(str,FileDevNull,aflags,kflags,OSSWRITE,oss);                   //CO20171024 - new sym framework
  // str.LatticeReduction_avoid=TRUE;                                                       //CO20171024 - new sym framework
  /*SYM::CalculatePointGroup(FileDevNull,str,aflags,PFSWRITE,OSSWRITE,oss);*/ oss << "str.pgroup.size()=" << str.pgroup.size() << endl;                   //CO20171024 - new sym framework
  /*str.BringInCell();*/                                                                                                                                  //CO20171024 - new sym framework
  /*SYM::CalculateSitePointGroup(FileDevNull,str,aflags,PFSWRITE,OSSWRITE,oss);*/ oss << "str.agroup.size()=" << str.agroup.size() << endl;               //CO20171024 - new sym framework
  /*SYM::CalculateFactorGroup(FileDevNull,str,aflags,PFSWRITE,OSSWRITE,oss);*/ oss << "str.fgroup.size()=" << str.fgroup.size() << endl;                  //CO20171024 - new sym framework
  /*SYM::CalculatePointGroupCrystal(FileDevNull,str,aflags,PFSWRITE,OSSWRITE,oss);*/ oss << "str.pgroup_xtal.size()=" << str.pgroup_xtal.size() << endl;  //CO20171024 - new sym framework
  /*str.sgroup_radius=1.5*RadiusSphereLattice(str.lattice);*/                                                                                             //CO20171024 - new sym framework
  /*SYM::CalculateSpaceGroup(FileDevNull,str,aflags,FALSE,OSSWRITE,oss);*/ oss << "str.sgroup.size()=" << str.sgroup.size() << endl;                      //CO20171024 - new sym framework
  oss << banner << endl; // ---------------------------------
  FFF << banner << endl; // ---------------------------------
  oss << "REFERENCE STRUCTURE BELOW" << endl; // STRUCTURE TO USE
  FFF << "REFERENCE STRUCTURE BELOW" << endl; // STRUCTURE TO USE
  oss << banner << endl; // ---------------------------------
  FFF << banner << endl; // ---------------------------------
  oss << str;
  FFF << str;
  oss << banner << endl; // ---------------------------------
  FFF << banner << endl; // ---------------------------------
  if(OSSWRITE) oss << banner << endl; // ---------------------------------
  if(OSSWRITE) oss << "CAGES CALCULATION - Stefano Curtarolo" << endl;
  if(OSSWRITE) FFF << "CAGES CALCULATION - Stefano Curtarolo" << endl;
  if(roughness<0.0) roughness=_CAGEROUGHNESS_;
  oss << "osswrite=" << osswrite << endl;
  oss << "roughness=" << roughness << endl;
  oss << "str.scale=" << str.scale << endl;
  GenerateGridAtoms(str,dims);
  oss << "size of the grid=" << str.grid_atoms.size() << endl;
  //  for(uint i=0;i<4;i++) cout << i << " " << str.grid_atoms.at(i).cpos << endl;

  bool go4=1,go3=1,go2=1;
  string caption="        x(fract)          y(fract)          z(fract)       S   radius(A)   T  C  P   coordination";
  bool oss1write=osswrite;

  AFLOW_PTHREADS::Clean_Threads();                                                       // multithread clean
  _threaded_GETCAGES_params params[MAX_ALLOCATABLE_PTHREADS];                        // multithread

  for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {              // multithread
    // construction of params[i]                                         // multithread
    params[ithread].paflags=&aflags;                                     // multithread
    params[ithread].ITHREAD=ithread;                                     // multithread
    params[ithread].THREADS_MAX=AFLOW_PTHREADS::MAX_PTHREADS;                      // multithread
    params[ithread].pstr=&str;                                           // multithread
    params[ithread].proughness=&roughness;                               // multithread
    params[ithread].pcages=&cagesreducible;                              // multithread
    params[ithread].pcages2=&cages2;                                     // multithread
    params[ithread].pcages3=&cages3;                                     // multithread
    params[ithread].pcages4=&cages4;                                     // multithread
    params[ithread].poss1write=&oss1write;                               // multithread
    params[ithread].poss1=&oss;                                          // multithread
    params[ithread].poss2write=&FFFflag;                                 // multithread
    params[ithread].poss2=&FFF;                                          // multithread
    params[ithread].itbusy=ithread;                                      // multithread
  }                                                                      // multithread

  // SCANNING CAGES 4 POINTS
  if(go4) {
    oss << "SCANNING STABLE CAGES (4 points)" << endl;
    FFF << "SCANNING STABLE CAGES (4 points)" << endl;
    oss << caption << endl;
    FFF << caption << endl;
    if(!AFLOW_PTHREADS::FLAG) {
      cerr << "NO PTHREADS" << endl;
      GetCages4(str,roughness,cagesreducible,cages4,TRUE,oss,FFFflag,FFF);
    }
    if(AFLOW_PTHREADS::FLAG) {
      cages4.clear();
      cerr << "START THREADS " << endl;
      for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {
        AFLOW_PTHREADS::viret[ithread]=pthread_create(&(AFLOW_PTHREADS::vpthread[ithread]),NULL,_threaded_GetCages4,(void*)&params[ithread]);
        //  aurostd::Sleep(10);
      }
      for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++)
        pthread_join((AFLOW_PTHREADS::vpthread[ithread]),NULL);
      sort(cages4.begin(),cages4.end(),_isort_acage_radius());
    }
    cage4_rmax=cages4.at(0).radius;
    cage4_rmin=cages4.at(cages4.size()-1).radius;
    oss.precision(_oss_short_precision_aflow_defects_);
    oss << "UNIQUE=" << cages4.size() << "  -  ";
    oss << "cage4 radius max=" << cage4_rmax << "  -  ";
    oss << "cage4 radius min=" << cage4_rmin << "  -  ";
    oss.precision(_oss_precision_aflow_defects_);
    oss << endl;
    FFF.precision(_oss_short_precision_aflow_defects_);
    FFF << "UNIQUE=" << cages4.size() << "  -  ";
    FFF << "cage4 radius max=" << cage4_rmax << "  -  ";
    FFF << "cage4 radius min=" << cage4_rmin << "  -  ";
    FFF.precision(_oss_precision_aflow_defects_);
    FFF << endl;
  }

  // SCANNING CAGES 3 POINTS
  if(go3) {
    oss << "SCANNING STABLE CAGES (3 points)" << endl;
    FFF << "SCANNING STABLE CAGES (3 points)" << endl;
    oss << caption << endl;
    FFF << caption << endl;
    if(!AFLOW_PTHREADS::FLAG) {
      //  cerr << "NO PTHREADS" << endl;
      GetCages3(str,roughness,cagesreducible,cages3,TRUE,oss,FFFflag,FFF);
    }
    if(AFLOW_PTHREADS::FLAG) {
      cages3.clear();
      //    cerr << "START THREADS3 " << endl;
      for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {
        AFLOW_PTHREADS::viret[ithread]=pthread_create(&(AFLOW_PTHREADS::vpthread[ithread]),NULL,_threaded_GetCages3,(void*)&params[ithread]);
        //  aurostd::Sleep(10);
      }
      for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++)
        pthread_join((AFLOW_PTHREADS::vpthread[ithread]),NULL);
      sort(cages3.begin(),cages3.end(),_isort_acage_radius());
    }
    if(cages3.size()>0) {
      cage3_rmax=cages3.at(0).radius;
      cage3_rmin=cages3.at(cages3.size()-1).radius;
      oss.precision(_oss_short_precision_aflow_defects_);
      oss << "UNIQUE=" << cages3.size() << "  -  ";
      oss << "cage3 radius max=" << cage3_rmax << "  -  ";
      oss << "cage3 radius min=" << cage3_rmin << "  -  ";
      oss.precision(_oss_precision_aflow_defects_);
      oss << endl;
      FFF.precision(_oss_short_precision_aflow_defects_);
      FFF << "UNIQUE=" << cages3.size() << "  -  ";
      FFF << "cage3 radius max=" << cage3_rmax << "  -  ";
      FFF << "cage3 radius min=" << cage3_rmin << "  -  ";
      FFF.precision(_oss_precision_aflow_defects_);
      FFF << endl;
    } else {
      oss << "NONE TO BE ADDED" << endl;
      FFF << "NONE TO BE ADDED" << endl;
    }
  }

  // SCANNING CAGES 2 POINTS
  if(go2) {
    oss << "SCANNING STABLE CAGES (2 points)" << endl;
    FFF << "SCANNING STABLE CAGES (2 points)" << endl;
    oss << caption << endl;
    FFF << caption << endl;
    if(!AFLOW_PTHREADS::FLAG) {
      //  cerr << "NO PTHREADS" << endl;
      GetCages2(str,roughness,cagesreducible,cages2,TRUE,oss,FFFflag,FFF);
    }
    if(AFLOW_PTHREADS::FLAG) {
      cages2.clear();
      //     cerr << "START THREADS2 " << endl;
      for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {
        AFLOW_PTHREADS::viret[ithread]=pthread_create(&(AFLOW_PTHREADS::vpthread[ithread]),NULL,_threaded_GetCages2,(void*)&params[ithread]);
        //  aurostd::Sleep(10);
      }
      for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++)
        pthread_join((AFLOW_PTHREADS::vpthread[ithread]),NULL);
      sort(cages2.begin(),cages2.end(),_isort_acage_radius());
    }
    if(cages2.size()>0) {
      cage2_rmax=cages2.at(0).radius;
      cage2_rmin=cages2.at(cages2.size()-1).radius;
      oss.precision(_oss_short_precision_aflow_defects_);
      oss << "UNIQUE=" << cages2.size() << "  -  ";
      oss << "cage2 radius max=" << cage2_rmax << "  -  ";
      oss << "cage2 radius min=" << cage2_rmin << "  -  ";
      oss.precision(_oss_precision_aflow_defects_);
      oss << endl;
      FFF.precision(_oss_short_precision_aflow_defects_);
      FFF << "UNIQUE=" << cages2.size() << "  -  ";
      FFF << "cage2 radius max=" << cage2_rmax << "  -  ";
      FFF << "cage2 radius min=" << cage2_rmin << "  -  ";
      FFF.precision(_oss_precision_aflow_defects_);
      FFF << endl;
    } else {
      oss << "NONE TO BE ADDED" << endl;
      FFF << "NONE TO BE ADDED" << endl;
    }
  }

  oss << "REDUCIBLE: " << cagesreducible.size() << endl;
  FFF << "REDUCIBLE: " << cagesreducible.size() << endl;
  oss << banner << endl; // ---------------------------------
  FFF << banner << endl; // ---------------------------------
  // ------------------------------
  sort(cagesreducible.begin(),cagesreducible.end(),_isort_acage_radius());

  // ------------------------------
  oss << "PERFORMING REDUCTION: " << endl;
  FFF << "PERFORMING REDUCTION: " << endl;
  oss << caption << endl;
  FFF << caption << endl;
  _atom ratom,iatom,tatom;
  bool cage_found;
  for(uint i=0;i<cagesreducible.size();i++) {
    ratom.cpos=cagesreducible[i].origin_cpos;
    ratom.fpos=C2F(str.lattice,ratom.cpos);
    cage_found=FALSE;
    for(uint ii=0;ii<cagesirreducible.size()&&!cage_found;ii++)
      if(abs(cagesreducible[i].radius-cagesirreducible[ii].radius)<eps) { // must have the same radius, at least
        iatom.cpos=cagesirreducible[ii].origin_cpos;
        iatom.fpos=C2F(str.lattice,iatom.cpos);
        // 	for(uint fg=0;fg<str.fgroup.size()&&!cage_found;fg++)  // FGROUP
        // 	  tatom=SYM::ApplyAtom(iatom,str.fgroup[fg],str,TRUE);  // FGROUP
        // 	  cage_found=(modulus(ratom.fpos-tatom.fpos)<eps);    // FGROUP
        for(uint sg=0;sg<str.sgroup.size()&&!cage_found;sg++) {  // SGROUP
          tatom=SYM::ApplyAtom(iatom,str.sgroup[sg],str);  // SGROUP
          cage_found=(modulus(ratom.cpos-tatom.cpos)<eps); // SGROUP
          if(cage_found) cagesreducible[i].cages_irrtype=cagesirreducible[ii].cages_irrtype;
          // DEBUG
          if(1) if(modulus(ratom.cpos-tatom.cpos)>eps && modulus(ratom.fpos-tatom.fpos)<eps) {
            cerr.setf(std::ios::fixed,std::ios::floatfield);
            cerr.precision(14);  //SC to cut/paste from matlab in format long
            cerr << ratom.fpos[1] << " " << ratom.fpos[2] << " " << ratom.fpos[3] << " " << endl;
            cerr << tatom.fpos[1] << " " << tatom.fpos[2] << " " << tatom.fpos[3] << " " << endl;
            cerr << ratom.cpos[1] << " " << ratom.cpos[2] << " " << ratom.cpos[3] << " " << endl;
            cerr << tatom.cpos[1] << " " << tatom.cpos[2] << " " << tatom.cpos[3] << " " << endl;
            return false; //debug
          }
        }
      }
    if(cage_found==FALSE) {                                     // new irreducible operation, generate and save it
      cagesirreducible.push_back(cagesreducible[i]);
      int ig=cagesirreducible.size()-1;
      cagesirreducible.at(ig).cages_irrtype=ig+1; // set type
      oss << "  "; for(int i=1;i<=3;i++) {if(cagesirreducible.at(ig).origin_fpos[i]>=0.0) oss << " "; oss << " " << cagesirreducible.at(ig).origin_fpos[i];}
      FFF << "  "; for(int i=1;i<=3;i++) {if(cagesirreducible.at(ig).origin_fpos[i]>=0.0) FFF << " "; FFF << " " << cagesirreducible.at(ig).origin_fpos[i];}
      // oss << "    i=" << ig;
      oss.precision(_oss_short_precision_aflow_defects_);
      oss << "  " << aurostd::PaddedPRE(species,2);
      oss << "   " << cagesirreducible.at(ig).radius << " ";
      oss << "  " << (int) cagesirreducible.at(ig).cages_irrtype;
      oss << "  " << (int) cagesirreducible.at(ig).coordination_position;
      oss << "  " << (int) cagesirreducible.at(ig).cages_position;
      {oss << "  [";for(uint iat=0;iat<cagesirreducible.at(ig).atoms.size();iat++) oss << str.species.at(cagesirreducible.at(ig).atoms.at(iat).type) << "(" << cagesirreducible.at(ig).atoms.at(iat).basis << ")" << (iat<cagesirreducible.at(ig).atoms.size()-1?",":""); oss << "]"; }
      oss.precision(_oss_precision_aflow_defects_);
      oss << endl;
      // oss << "*";oss.flush();
      // FFF << "    i=" << ig;
      FFF.precision(_oss_short_precision_aflow_defects_);
      FFF << "  " << aurostd::PaddedPRE(species,2);
      FFF << "   " << cagesirreducible.at(ig).radius << " ";
      FFF << "  " << (int) cagesirreducible.at(ig).cages_irrtype;
      FFF << "  " << (int) cagesirreducible.at(ig).coordination_position;
      FFF << "  " << (int) cagesirreducible.at(ig).cages_position;
      {FFF << "  [";for(uint iat=0;iat<cagesirreducible.at(ig).atoms.size();iat++) FFF << str.species.at(cagesirreducible.at(ig).atoms.at(iat).type) << "(" << cagesirreducible.at(ig).atoms.at(iat).basis << ")" << (iat<cagesirreducible.at(ig).atoms.size()-1?",":""); FFF << "]"; }
      FFF.precision(_oss_precision_aflow_defects_);
      FFF << endl;
      // oss << "*";oss.flush();
    }
  }

  // ------------------------------------------------------------------------------------
  sort(cagesirreducible.begin(),cagesirreducible.end(),_isort_acage_radius());
  // ------------------------------------------------------------------------------------
  oss << "IRREDUCIBLE: " << cagesirreducible.size() << endl;
  FFF << "IRREDUCIBLE: " << cagesirreducible.size() << endl;
  if(1) {
    oss << banner << endl; // ---------------------------------
    FFF << banner << endl; // ---------------------------------
    oss << "IRREDUCIBLE/REDUCIBLE summary: " << cagesirreducible.size() << "/" << cagesreducible.size() << endl;
    FFF << "IRREDUCIBLE/REDUCIBLE summary: " << cagesirreducible.size() << "/" << cagesreducible.size() << endl;
    oss << caption << endl;
    FFF << caption << endl;
    for(uint ig=0;ig<cagesirreducible.size();ig++) {
      oss << "  "; for(int i=1;i<=3;i++) {if(cagesirreducible.at(ig).origin_fpos[i]>=0.0) oss << " "; oss << " " << cagesirreducible.at(ig).origin_fpos[i];}
      FFF << "  "; for(int i=1;i<=3;i++) {if(cagesirreducible.at(ig).origin_fpos[i]>=0.0) FFF << " "; FFF << " " << cagesirreducible.at(ig).origin_fpos[i];}
      // oss << "    i=" << ig;
      oss.precision(_oss_short_precision_aflow_defects_);
      oss << "  " << aurostd::PaddedPRE(species,2);
      oss << "   " << cagesirreducible.at(ig).radius << " ";
      oss << "  " << (int) cagesirreducible.at(ig).cages_irrtype;
      oss << "  " << (int) cagesirreducible.at(ig).coordination_position;
      oss << "  " << (int) cagesirreducible.at(ig).cages_position;
      {oss << "  [";for(uint iat=0;iat<cagesirreducible.at(ig).atoms.size();iat++) oss << str.species.at(cagesirreducible.at(ig).atoms.at(iat).type) << "(" << cagesirreducible.at(ig).atoms.at(iat).basis << ")" << (iat<cagesirreducible.at(ig).atoms.size()-1?",":""); oss << "]"; }
      oss.precision(_oss_precision_aflow_defects_);
      oss << endl;
      FFF.precision(_oss_short_precision_aflow_defects_);
      FFF << "  " << aurostd::PaddedPRE(species,2);
      FFF << "   " << cagesirreducible.at(ig).radius << " ";
      FFF << "  " << (int) cagesirreducible.at(ig).cages_irrtype;
      FFF << "  " << (int) cagesirreducible.at(ig).coordination_position;
      FFF << "  " << (int) cagesirreducible.at(ig).cages_position;
      {FFF << "  [";for(uint iat=0;iat<cagesirreducible.at(ig).atoms.size();iat++) FFF << str.species.at(cagesirreducible.at(ig).atoms.at(iat).type) << "(" << cagesirreducible.at(ig).atoms.at(iat).basis << ")" << (iat<cagesirreducible.at(ig).atoms.size()-1?",":""); FFF << "]"; }
      FFF.precision(_oss_precision_aflow_defects_);
      FFF << endl;
      //   oss << " ("; for(int i=1;i<=3;i++) {if(cagesirreducible.at(ig).origin_cpos[i]>=0.0) oss << " "; oss << " " << cagesirreducible.at(ig).origin_cpos[i];};oss << ")c" << endl;
      //    FFF << " ("; for(int i=1;i<=3;i++) {if(cagesirreducible.at(ig).origin_cpos[i]>=0.0) FFF << " "; FFF << " " << cagesirreducible.at(ig).origin_cpos[i];};FFF << ")c" << endl;

      for(uint jg=0;jg<cagesreducible.size();jg++)
        if(abs(cagesreducible[jg].cages_irrtype-cagesirreducible.at(ig).cages_irrtype)<eps) {
          oss.precision(_oss_short_precision_aflow_defects_);
          oss << "      equiv "; for(int i=1;i<=3;i++) {if(cagesreducible[jg].origin_fpos[i]>=0.0) oss << " "; oss << " " << cagesreducible[jg].origin_fpos[i];}
          oss.precision(_oss_precision_aflow_defects_);
          oss << "   type=" << (int) cagesreducible[jg].cages_irrtype << " ";
          oss << endl;
          FFF.precision(_oss_short_precision_aflow_defects_);
          FFF << "      equiv "; for(int i=1;i<=3;i++) {if(cagesreducible[jg].origin_fpos[i]>=0.0) FFF << " "; FFF << " " << cagesreducible[jg].origin_fpos[i];}
          FFF.precision(_oss_precision_aflow_defects_);
          FFF << "   type=" << (int) cagesreducible[jg].cages_irrtype << " ";
          FFF << endl;
        }
    }
  }

  // for(int i=0;i<grid_origin_cpos.size();i++) delete grid_origin_cpos[i]; grid_origin_cpos.clear();
  // for(int i=0;i<grid_origin_fpos.size();i++) delete grid_origin_fpos[i]; grid_origin_fpos.clear();

  sort(cages2.begin(),cages2.end(),_isort_acage_radius());
  sort(cages3.begin(),cages3.end(),_isort_acage_radius());
  sort(cages4.begin(),cages4.end(),_isort_acage_radius());
  sort(cagesreducible.begin(),cagesreducible.end(),_isort_acage_radius());
  sort(cagesreducible.begin(),cagesreducible.end(),_isort_acage_radius());
  cagesreducible_rmax=cagesreducible.at(0).radius;
  cagesreducible_rmin=cagesreducible.at(cagesreducible.size()-1).radius;
  cagesirreducible_rmax=cagesirreducible.at(0).radius;
  cagesirreducible_rmin=cagesirreducible.at(cagesirreducible.size()-1).radius;

  if(cagesreducible_rmax) {;} // phony
  if(cagesreducible_rmin) {;} // phony
  if(cagesirreducible_rmax) {;} // phony
  if(cagesirreducible_rmin) {;} // phony

  oss << banner << endl; // ---------------------------------
  FFF << banner << endl; // ---------------------------------
  oss << "IRREDUCIBLE=" << cagesirreducible.size() << "   REDUCIBLE=" << cagesreducible.size() << endl;
  FFF << "IRREDUCIBLE=" << cagesirreducible.size() << "   REDUCIBLE=" << cagesreducible.size() << endl;
  oss << banner << endl; // ---------------------------------
  FFF << banner << endl; // ---------------------------------
  oss << aflow::Banner("BANNER_TINY") << endl;
  FFF << aflow::Banner("BANNER_TINY") << endl;
  oss << banner << endl; // ---------------------------------
  FFF << banner << endl; // ---------------------------------

  FFF.close();

  // print something

  stringstream aflowin;

  vector<xstructure> vstr;
  vector<string> vtitle;
  for(uint ig=0;ig<cagesirreducible.size();ig++) {
    xstructure straus(str);straus.write_inequivalent_flag=FALSE;
    stringstream aus;
    aus << " R=" << cagesirreducible.at(ig).radius << " T=" << (int) cagesirreducible.at(ig).cages_irrtype;
    aus << " C=" << (int) cagesirreducible.at(ig).coordination_position << " P=" << (int) cagesirreducible.at(ig).cages_position;
    {aus << " [";for(uint iat=0;iat<cagesirreducible.at(ig).atoms.size();iat++) aus << str.species.at(cagesirreducible.at(ig).atoms.at(iat).type) << "(" << cagesirreducible.at(ig).atoms.at(iat).basis << ")" << (iat<cagesirreducible.at(ig).atoms.size()-1?",":""); aus << "]"; }
    straus.title+=aus.str();

    _atom atom;
    atom=str.atoms.at(0);
    atom.fpos=cagesirreducible.at(ig).origin_fpos;
    atom.cpos=cagesirreducible.at(ig).origin_cpos;
    atom.type=str.species.size();       // to force an extra species
    atom.name="H";atom.cleanname="H";   // atom.name=string_species;atom.cleanname=string_species;
    // ADD atom, put alphabetic
    straus.AddAtom(atom,false); //CO20230319 - add by type
    straus.SpeciesPutAlphabetic();
    // make structure
    vstr.push_back(straus);
    vtitle.push_back("cages_T"+aurostd::utype2string(cagesirreducible.at(ig).cages_irrtype));
  }

  aflowin << AFLOWIN_SEPARATION_LINE << endl;
  aflowin << "[AFLOW] CAGES IN AFLOW.IN FORMAT" << endl;
  aflowin << AFLOWIN_SEPARATION_LINE << endl;
  for(uint ivstr=0;ivstr<vstr.size();ivstr++) {
    aflowin << _VASP_POSCAR_MODE_EXPLICIT_START_P_ << vtitle.at(ivstr) << endl;  //CO20200624
    aflowin << vstr.at(ivstr);
    aflowin << _VASP_POSCAR_MODE_EXPLICIT_STOP_P_ << vtitle.at(ivstr) << endl; //CO20200624
    aflowin << AFLOWIN_SEPARATION_LINE << endl;
  }  
  aflowin << AFLOWIN_SEPARATION_LINE << endl;

  //  cout << aflowin.str();
  if(FFFflag) aurostd::stringstream2file(aflowin,_AFLOWIN_+".part");  // ME20210521
  oss << aflowin.str() << endl; // ---------------------------------
  FFF << aflowin.str() << endl; // ---------------------------------

  return Krun;
}


// ***************************************************************************
// Get Sphere from 4 points
// FROM http://home.att.net/~srschmitt/script_sphere_solver.html
// From analytic geometry, we know that there is a unique sphere that passes through
// four non-coplanar points if, and only if, they are not on the same plane.
// If they are on the same plane, either there are no spheres through the 4 points,
// or an infinite number of them if the 4 points are on a circle. Given 4 points,
//    {x1, y1, z1}, {x2, y2, z2}, {x3, y3, z3}, {x4, y4, z4}
// how does one find the center and radius of a sphere exactly fitting those points?
//  They can be found by solving the following determinant equation:
// | x^2 + y^2 + z^2 	x 	y 	z 	1 |
// | x1^2 + y1^2 + z1^2	x1	y1	z1	1 |
// | x2^2 + y2^2 + z2^2	x2	y2	z2	1 |
// | x3^2 + y3^2 + z3^2	x3	y3	z3	1 |
// | x4^2 + y4^2 + z4^2	x4	y4	z4	1 |
// = 0
// Evaluating the cofactors for the first row of the determinant can give us a solution.
// The determinant equation can be written as an equation of these cofactors:
//
// (x^2 + y^2 + z^2)*M11 - x*M12 + y*M13 - z*M14 + M15 = 0
// This can be converted to the canonical form of the equation of a sphere:
//  x^2 + y^2 + z^2 - (M12/M11)*x + (M13/M11)*y - (M14/M11)*z + M15/M11 = 0
// Completing the squares in x and y and z gives:
//     x0 =  0.5*M12/M11
//  y0 = -0.5*M13/M11
//   z0 =  0.5*M14/M11
//   r0^2 = x0^2 + y0^2 + z0^2 - M15/M11
//
// Note that there is no solution when M11 is equal to zero. In this
// case, the points are not on a sphere; they may all be on a plane or three points may be on a straight line.
//

// ***************************************************************************
// Get Circumcircle from 3 points
// http://mathworld.wolfram.com/Circumcircle.html

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
