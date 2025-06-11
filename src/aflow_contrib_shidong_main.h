// ***************************************************************************
// *                                                                         *
// *               AFlow SHIDONG WANG - Duke University 2010-2011            *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_shidong_main.h
// functions written by
// 2010-2011: shidong.wang@duke.edu

#ifndef _AFLOW_CONTRIB_SHIDONG_MAIN_H_
#define _AFLOW_CONTRIB_SHIDONG_MAIN_H_

// put your prototypes here

#include "aflow_contrib_shidong_cluster_expansion.h"

const int _site_num_low = 1;
const int _site_num_up = 6;
const int _NNNum_low = 1;
const int _NNNum_up = 3;

namespace pflow {
  void AClusterExpansionMethodMain(string options);
  void SQS(string options);
  void Superlattice(string options);
  void Cluster(string options);
}
void CheckAllInputFileExistence(string structure_type);

//void ACEPhaseBoundary(montecarlo phase1, montecarlo phase2);
//void ACEPhaseBoundary(montecarlo phase1, montecarlo phase2, double T, double mu);
//bool CheckDiscontinuity(vector< vector<double > > xy, double x,        double y, double y_dev);
//bool CheckPhaseTransition(vector<_state_status> status_list, int fit_pt_num);
//void PrintOutStateStatus(ostream & os, vector<_state_status> state_list);
//
//double polynomial_interpolation_extrapolation(vector< vector<double> > xy, double x);
//
//void EnsemblePotential(montecarlo phase1, double _T_low, double _T_high,        double _mu_low, double _mu_high);
namespace pflow {
  bool EffectiveMass(vector<string> &argv,string directory,ostream& oss);
}

namespace pflow {
  // animation gif by using jmol and imageMagick
  void JMOLAnimation(istream& input, vector<string> argv);
}

#endif

// ***************************************************************************
// *                                                                         *
// *               AFlow SHIDONG WANG - Duke University 2010-2011            *
// *                                                                         *
// ***************************************************************************
