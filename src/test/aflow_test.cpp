// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xrandom.h"
#include "AUROSTD/aurostd_xvector.h"

#include "structure/aflow_xatom.h"

using aurostd::xmatrix;
using aurostd::xvector;
using std::endl;
using std::ifstream;
using std::istream;
using std::istringstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::string;
using std::stringstream;
using std::vector;

#define NTERM 5
#define SPREAD 0.1
#define NPT 100

// ***************************************************************************

// bool AFLOW_PTHREADS::FLAG;
// int  AFLOW_PTHREADS::MAX_PTHREADS;
// int  AFLOW_PTHREADS::RUNNING;
// pthread_t thread[MAX_ALLOCATABLE_PTHREADS];
// int iret[MAX_ALLOCATABLE_PTHREADS];
// bool thread_busy[MAX_ALLOCATABLE_PTHREADS];

void PERFORM_TESTJ(ostream& oss) {
  // load ICSD
  xmatrix<double> A(5, 5);
  xmatrix<double> v(5, 5);
  xvector<double> d(5);
  A(1, 1) = 1;
  A(1, 2) = 2;
  A(1, 3) = 3;
  A(1, 4) = 4;
  A(1, 5) = 5;
  A(2, 1) = 2;
  A(2, 2) = 2;
  A(2, 3) = 6;
  A(2, 4) = 7;
  A(2, 5) = 8;
  A(3, 1) = 3;
  A(3, 2) = 6;
  A(3, 3) = 3;
  A(3, 4) = 2;
  A(3, 5) = 1;
  A(4, 1) = 4;
  A(4, 2) = 7;
  A(4, 3) = 2;
  A(4, 4) = 4;
  A(4, 5) = 1;
  A(5, 1) = 5;
  A(5, 2) = 8;
  A(5, 3) = 1;
  A(5, 4) = 1;
  A(5, 5) = 5;
  oss << A << endl;
  //   int jacobi(const xmatrix<utype> &ain,xvector<utype> &d,xmatrix<utype> &v) {
  // Computes all eigenvalues and eigenvectors of a real symmetric xmatrix a[1..n][1..n].
  // On output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
  // v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of
  // a. The function returns the number of Jacobi rotations that were required.
  aurostd::jacobi(A, d, v);

  oss << aurostd::trasp(d) << endl;

  oss << v << endl;
}

void PERFORM_PRX(ostream& oss) {
  const vector<string> vsystem{"AgPd", "AgPt", "AuPd", "CdPd", "CdPt", "CdRh", "CoPd", "CoPt", "CrIr", "CrOs", "CrPd", "CrPt", "CrRh", "CuPd", "CuPt", "CuRh", "FeIr", "FePd", "FePt", "FeRh",
                               "HfIr", "HfOs", "HfPd", "HfPt", "HfRh", "HfRu", "HgPd", "HgPt", "HgRh", "IrMn", "IrMo", "IrNb", "IrNi", "IrOs", "IrRe", "IrRh", "IrRu", "IrSc", "IrTa", "IrTc",
                               "IrTi", "IrV",  "IrW",  "IrY",  "IrZn", "IrZr", "MnOs", "MnPd", "MnPt", "MnRh", "MnRu", "MoOs", "MoPd", "MoPt", "MoRh", "MoRu", "NbOs", "NbPd", "NbPt", "NbRh",
                               "NbRu", "NiPd", "NiPt", "OsRe", "OsRh", "OsRu", "OsSc", "OsTa", "OsTc", "OsTi", "OsV",  "OsW",  "OsY",  "OsZr", "PdPt", "PdRe", "PdSc", "PdTa", "PdTc", "PdTi",
                               "PdV",  "PdW",  "PdY",  "PdZn", "PdZr", "PtRe", "PtRh", "PtRu", "PtSc", "PtTa", "PtTc", "PtTi", "PtV",  "PtW",  "PtY",  "PtZn", "PtZr", "ReRh", "ReRu", "RhRu",
                               "RhSc", "RhTa", "RhTc", "RhTi", "RhV",  "RhW",  "RhY",  "RhZn", "RhZr", "RuSc", "RuTa", "RuTc", "RuTi", "RuV",  "RuW",  "RuY",  "RuZn", "RuZr"};

  uint figN = 5;
  for (size_t i = 0; i < vsystem.size(); i += 2) {
    oss << "\\begin{figure*} ";// << endl;
    oss << "\\includegraphics[width=0.93\\linewidth,height=116mm]{shulls/fig_" << vsystem[i] << ".eps} ";// << endl;
    oss << "\\includegraphics[width=0.93\\linewidth,height=116mm]{shulls/fig_" << vsystem.at(i + 1) << ".eps} ";// << endl;
    oss << "\\caption{\\small Convex hulls for the systems " << vsystem[i] << " and " << vsystem.at(i + 1) << ".} ";// << endl;
    oss << "\\label{fig" << figN << "} ";// << endl;
    oss << "\\end{figure*}" << endl;
    figN++;
  }
  oss << endl;
  oss << R"(\def\allfigures{{\cref{)";
  for (uint i = 5; i <= figN; i++) {
    oss << "fig" << i << (i < figN ? "," : "");
  }
  oss << "}}}" << endl;
}

bool isPGM(string element) {
  if (element == "Os" || element == "Ru" || element == "Ir" || element == "Rh" || element == "Pd" || element == "Pt") {
    return true;
  }
  return false;
}

void PERFORM_TEST_ALLOYS(ostream& oss) {
  // load ICSD
  vector<string> vprotos;
  // aurostd::string2vectorstring(aflowlib::PrototypesIcsdReadme(),vprotos);oss << "vprotos.size()=" << vprotos.size() << endl;
  vector<string> vprotos2;
  for (size_t i = 0; i < vprotos.size(); i++) {
    if (aurostd::substring2bool(vprotos[i], "(2)")) {
      vprotos2.push_back(vprotos[i]);
    }
  }
  oss << "vprotos2.size()=" << vprotos2.size() << endl;
  vector<string> vprotos3;
  for (size_t i = 0; i < vprotos.size(); i++) {
    if (aurostd::substring2bool(vprotos[i], "(3)")) {
      vprotos3.push_back(vprotos[i]);
    }
  }
  oss << "vprotos3.size()=" << vprotos3.size() << endl;
  vector<string> vprotos4;
  for (size_t i = 0; i < vprotos.size(); i++) {
    if (aurostd::substring2bool(vprotos[i], "(4)")) {
      vprotos4.push_back(vprotos[i]);
    }
  }
  oss << "vprotos4.size()=" << vprotos4.size() << endl;
  vector<string> vprotos5;
  for (size_t i = 0; i < vprotos.size(); i++) {
    if (aurostd::substring2bool(vprotos[i], "(5)")) {
      vprotos5.push_back(vprotos[i]);
    }
  }
  oss << "vprotos5.size()=" << vprotos5.size() << endl;
  vector<string> vprotos6;
  for (size_t i = 0; i < vprotos.size(); i++) {
    if (aurostd::substring2bool(vprotos[i], "(6)")) {
      vprotos6.push_back(vprotos[i]);
    }
  }
  oss << "vprotos6.size()=" << vprotos6.size() << endl;
  vector<string> velement;
  // aurostd::string2tokens(string("Ag,Au,Cd,Co,Cr,Cu,Fe,Hf,Hg,Ir,La,Mn,Mo,Nb,Ni,Os,Pd,Pt,Re,Rh,Ru,Sc,Ta,Tc,Ti,V,W,Y,Zn,Zr"),velement,","); // Cd is bad
  // aurostd::string2tokens(string("Ag,Au,Co,Cr,Cu,Fe,Hf,Hg,Ir,La,Mn,Mo,Nb,Ni,Os,Pd,Pt,Re,Rh,Ru,Sc,Ta,Tc,Ti,V,W,Y,Zn,Zr"),velement,","); // Hg is bad
  // aurostd::string2tokens(string("Ag,Al,Au,Co,Cr,Cu,Fe,Hf,Ir,La,Mg,Mn,Mo,Nb,Ni,Os,Pd,Pt,Re,Rh,Ru,Sc,Ta,Tc,Ti,V,W,Y,Zn,Zr"),velement,","); //-Cd-Hg +Al+Mg
  // aurostd::string2tokens(string("Ag,Al,Au,Cd,Co,Cr,Cu,Fe,Hf,Hg,Ir,La,Mg,Mn,Mo,Nb,Ni,Os,Pd,Pt,Re,Rh,Ru,Sc,Ta,Tc,Ti,V,W,Y,Zn,Zr"),velement,","); // +Al+Mg
  //   aurostd::string2tokens(string("As,Bi,Ba,Be,Ag,Al,Au,Cd,Co,Cr,Cu,Fe,Ga,Ge,Hf,Hg,Ir,La,Li,Mg,Mn,Mo,Nb,Ni,Os,Pd,Pt,Re,Rh,Ru,Sc,Ta,Tc,Ti,V,W,Y,Zn,Zr"),velement,","); // +Al+Mg
  //  aurostd::string2tokens(string("Ag,Au,Cd,Co,Cr,Cu,Fe,Hf,Hg,Ir,La,Mn,Mo,Nb,Ni,Os,Pd,Pt,Re,Rh,Ru,Sc,Ta,Tc,Ti,V,W,Y,Zn,Zr"),velement,",");
  aurostd::string2tokens(string("Ag,Al,As,Au,B,Ba,Be,Bi,Br,Ca,Cd,Cl,Co,Cr,Cu,Fe,Ga,Ge,Hf,Hg,In,Ir,K,La,Li,Mg,Mn,Mo,Na,Nb,Ni,Os,P,Pb,Pd,Pt,Re,Rh,Ru,Sb,Sc,Se,Si,Sn,Sr,Ta,Tc,Te,Ti,Tl,V,W,Y,Zn,Zr"), velement, ",");
  aurostd::sort(velement);

  //  for(size_t i=0;i<velement.size();i++)
  //     for(size_t j=i+1;j<velement.size();j++)
  //       for(size_t k=j+1;k<velement.size();k++)
  // 	cout << "/home/auro/work/AFLOW3/aflow --terdata " << velement.at(i) << " " << velement.at(j) << " " << velement.at(k) << endl;

  vector<string> tokens;
  vector<string> vspecies;
  bool found = false;

  vector<string> vicsd2;
  for (size_t i = 0; i < velement.size(); i++) {
    for (size_t j = i + 1; j < velement.size(); j++) {  // if(isPGM(velement[i]) || isPGM(velement[j]))
      {
        for (size_t iproto = 0; iproto < vprotos2.size(); iproto++) {
          if (aurostd::substring2bool(vprotos2[iproto], velement[i]) && aurostd::substring2bool(vprotos2[iproto], velement[j])) {
            aurostd::string2tokens(vprotos2[iproto], tokens);
            XATOM_SplitAlloySpecies(string(tokens.at(0)), vspecies);
            if (vspecies.at(0) == velement[i] && vspecies.at(1) == velement[j]) {
              found = false;
              vector<string> tokensicsd;
              for (size_t n = 0; n < vicsd2.size() && !found; n++) {
                if (aurostd::substring2bool(vicsd2[n], velement[i]) && aurostd::substring2bool(vicsd2[n], velement[j])) {
                  aurostd::string2tokens(vicsd2[n], tokensicsd);
                  // cerr << tokens.at(0) << ":" << tokensicsd.at(0) << " - " << tokens.at(2) << ":" << tokensicsd.at(2) << endl;
                  found = found || (tokens.at(0) == tokensicsd.at(0) && tokens.at(2) == tokensicsd.at(2)); // alloy composition and pearson
                }
              }
            }
            if (!found) {
              //   oss << vprotos2[iproto] << endl;
              //      oss << "aflow --missing --aflow_proto=" << velement[i] << velement[j] << "/" << tokens.at(tokens.size()-2) << ".AB" << endl;
              vicsd2.push_back(vprotos2[iproto]);
            }
          }
        }
      }
    }
  }
  oss << "vicsd2.size()=" << vicsd2.size() << endl;

  vector<string> vicsd3;
  for (size_t i = 0; i < velement.size(); i++) {
    for (size_t j = i + 1; j < velement.size(); j++) {
      for (size_t k = j + 1; k < velement.size(); k++) {
        // if(isPGM(velement[i]) || isPGM(velement[j]) || isPGM(velement[k]))
        {
          for (size_t iproto = 0; iproto < vprotos3.size(); iproto++) {
            if (aurostd::substring2bool(vprotos3[iproto], velement[i]) && aurostd::substring2bool(vprotos3[iproto], velement[j]) && aurostd::substring2bool(vprotos3[iproto], velement[k])) {
              aurostd::string2tokens(vprotos3[iproto], tokens);
              XATOM_SplitAlloySpecies(string(tokens.at(0)), vspecies);
              if (vspecies.at(0) == velement[i] && vspecies.at(1) == velement[j] && vspecies.at(2) == velement[k]) {
                found = false;
                vector<string> tokensicsd;
                for (size_t n = 0; n < vicsd3.size() && !found; n++) {
                  if (aurostd::substring2bool(vicsd3[n], velement[i]) && aurostd::substring2bool(vicsd3[n], velement[j]) && aurostd::substring2bool(vicsd3[n], velement[k])) {
                    aurostd::string2tokens(vicsd3[n], tokensicsd);
                    found = found || (tokens.at(0) == tokensicsd.at(0) && tokens.at(2) == tokensicsd.at(2)); // alloy composition and pearson
                  }
                }
              }
              if (!found) {
                //	oss << vprotos3[iproto] << endl;
                oss << "aflow --missing --aflow_proto=" << velement[i] << velement[j] << velement[k] << "/" << tokens.at(tokens.size() - 2) << ".ABC" << endl;
                vicsd3.push_back(vprotos3[iproto]);
              }
            }
          }
        }
      }
    }
  }
  oss << "vicsd3.size()=" << vicsd3.size() << endl;
}

void PERFORM_TEST3(ostream& oss) {
  vector<string> vprotos;
  // FIX aurostd::string2vectorstring(aflowlib::PrototypesIcsdReadme(),vprotos);
  vector<string> vprotos3;
  for (size_t i = 0; i < vprotos.size(); i++) {
    if (aurostd::substring2bool(vprotos[i], "(3)")) {
      vprotos3.push_back(vprotos[i]);
    }
  }

  oss << vprotos3.size() << endl;

  vector<string> velement;
  aurostd::string2tokens(string("Ag,Al,As,Au,B,Ba,Be,Bi,Br,Ca,Cd,Cl,Co,Cr,Cu,Fe,Ga,Ge,Hf,Hg,In,Ir,K,La,Li,Mg,Mn,Mo,Na,Nb,Ni,Os,P,Pb,Pd,Pt,Re,Rh,Ru,Sb,Sc,Se,Si,Sn,Sr,Ta,Tc,Te,Ti,Tl,V,W,Y,Zn,Zr"), velement, ",");
  aurostd::sort(velement);

  vector<string> vicsd;
  for (size_t i = 0; i < velement.size(); i++) {
    for (size_t j = i + 1; j < velement.size(); j++) {
      for (size_t k = j + 1; k < velement.size(); k++) {
        bool found105 = false;
        // the 105 list from Jesus
        if (velement[i] == "Al" && velement[j] == "Au" && velement[k] == "Hf") {
          found105 = true;
        }
        if (velement[i] == "Al" && velement[j] == "Ge" && velement[k] == "Li") {
          found105 = true;
        }
        if (velement[i] == "Al" && velement[j] == "Li" && velement[k] == "Si") {
          found105 = true;
        }
        if (velement[i] == "As" && velement[j] == "Co" && velement[k] == "Hf") {
          found105 = true;
        }
        if (velement[i] == "As" && velement[j] == "Co" && velement[k] == "Ti") {
          found105 = true;
        }
        if (velement[i] == "As" && velement[j] == "Co" && velement[k] == "Zr") {
          found105 = true;
        }
        if (velement[i] == "As" && velement[j] == "Fe" && velement[k] == "Nb") {
          found105 = true;
        }
        if (velement[i] == "As" && velement[j] == "Fe" && velement[k] == "Ta") {
          found105 = true;
        }
        if (velement[i] == "As" && velement[j] == "Ir" && velement[k] == "Ti") {
          found105 = true;
        }
        if (velement[i] == "As" && velement[j] == "Ir" && velement[k] == "Zr") {
          found105 = true;
        }
        if (velement[i] == "As" && velement[j] == "Nb" && velement[k] == "Ru") {
          found105 = true;
        }
        if (velement[i] == "As" && velement[j] == "Ni" && velement[k] == "Sc") {
          found105 = true;
        }
        if (velement[i] == "As" && velement[j] == "Rh" && velement[k] == "Ti") {
          found105 = true;
        }
        if (velement[i] == "As" && velement[j] == "Rh" && velement[k] == "Zr") {
          found105 = true;
        }
        if (velement[i] == "As" && velement[j] == "Ru" && velement[k] == "Ta") {
          found105 = true;
        }
        if (velement[i] == "Ba" && velement[j] == "Bi" && velement[k] == "K") {
          found105 = true;
        }
        if (velement[i] == "Bi" && velement[j] == "Co" && velement[k] == "Hf") {
          found105 = true;
        }
        if (velement[i] == "Bi" && velement[j] == "Co" && velement[k] == "Ti") {
          found105 = true;
        }
        if (velement[i] == "Bi" && velement[j] == "Co" && velement[k] == "Zr") {
          found105 = true;
        }
        if (velement[i] == "Bi" && velement[j] == "Hf" && velement[k] == "Rh") {
          found105 = true;
        }
        if (velement[i] == "Bi" && velement[j] == "Ir" && velement[k] == "Zr") {
          found105 = true;
        }
        if (velement[i] == "Bi" && velement[j] == "Ni" && velement[k] == "Sc") {
          found105 = true;
        }
        if (velement[i] == "Bi" && velement[j] == "Ni" && velement[k] == "Y") {
          found105 = true;
        }
        if (velement[i] == "Bi" && velement[j] == "Pd" && velement[k] == "Sc") {
          found105 = true;
        }
        if (velement[i] == "Bi" && velement[j] == "Rh" && velement[k] == "Ti") {
          found105 = true;
        }
        if (velement[i] == "Bi" && velement[j] == "Rh" && velement[k] == "Zr") {
          found105 = true;
        }
        if (velement[i] == "B" && velement[j] == "Li" && velement[k] == "Si") {
          found105 = true;
        }
        if (velement[i] == "Cd" && velement[j] == "Na" && velement[k] == "P") {
          found105 = true;
        }
        if (velement[i] == "Cl" && velement[j] == "La" && velement[k] == "Se") {
          found105 = true;
        }
        if (velement[i] == "Co" && velement[j] == "Ge" && velement[k] == "Nb") {
          found105 = true;
        }
        if (velement[i] == "Co" && velement[j] == "Ge" && velement[k] == "Ta") {
          found105 = true;
        }
        if (velement[i] == "Co" && velement[j] == "Ge" && velement[k] == "V") {
          found105 = true;
        }
        if (velement[i] == "Co" && velement[j] == "Hf" && velement[k] == "Sb") {
          found105 = true;
        }
        if (velement[i] == "Co" && velement[j] == "Nb" && velement[k] == "Si") {
          found105 = true;
        }
        if (velement[i] == "Co" && velement[j] == "Nb" && velement[k] == "Sn") {
          found105 = true;
        }
        if (velement[i] == "Co" && velement[j] == "Sb" && velement[k] == "Ti") {
          found105 = true;
        }
        if (velement[i] == "Co" && velement[j] == "Sb" && velement[k] == "Zr") {
          found105 = true;
        }
        if (velement[i] == "Co" && velement[j] == "Si" && velement[k] == "Ta") {
          found105 = true;
        }
        if (velement[i] == "Co" && velement[j] == "Sn" && velement[k] == "Ta") {
          found105 = true;
        }
        if (velement[i] == "Co" && velement[j] == "Sn" && velement[k] == "V") {
          found105 = true;
        }
        if (velement[i] == "Fe" && velement[j] == "Ge" && velement[k] == "W") {
          found105 = true;
        }
        if (velement[i] == "Fe" && velement[j] == "Nb" && velement[k] == "Sb") {
          found105 = true;
        }
        if (velement[i] == "Fe" && velement[j] == "Sb" && velement[k] == "Ta") {
          found105 = true;
        }
        if (velement[i] == "Fe" && velement[j] == "Sb" && velement[k] == "V") {
          found105 = true;
        }
        if (velement[i] == "Fe" && velement[j] == "Te" && velement[k] == "Ti") {
          found105 = true;
        }
        if (velement[i] == "Ga" && velement[j] == "Nb" && velement[k] == "Ni") {
          found105 = true;
        }
        if (velement[i] == "Ga" && velement[j] == "Pt" && velement[k] == "Ta") {
          found105 = true;
        }
        if (velement[i] == "Ge" && velement[j] == "Hf" && velement[k] == "Ni") {
          found105 = true;
        }
        if (velement[i] == "Ge" && velement[j] == "Ir" && velement[k] == "Nb") {
          found105 = true;
        }
        if (velement[i] == "Ge" && velement[j] == "Ir" && velement[k] == "Ta") {
          found105 = true;
        }
        if (velement[i] == "Ge" && velement[j] == "Ir" && velement[k] == "V") {
          found105 = true;
        }
        if (velement[i] == "Ge" && velement[j] == "Ni" && velement[k] == "Ti") {
          found105 = true;
        }
        if (velement[i] == "Ge" && velement[j] == "Ni" && velement[k] == "Zr") {
          found105 = true;
        }
        if (velement[i] == "Ge" && velement[j] == "Pd" && velement[k] == "Zr") {
          found105 = true;
        }
        if (velement[i] == "Ge" && velement[j] == "Pt" && velement[k] == "Ti") {
          found105 = true;
        }
        if (velement[i] == "Ge" && velement[j] == "Pt" && velement[k] == "Zr") {
          found105 = true;
        }
        if (velement[i] == "Hf" && velement[j] == "Ir" && velement[k] == "Sb") {
          found105 = true;
        }
        if (velement[i] == "Hf" && velement[j] == "Ni" && velement[k] == "Sn") {
          found105 = true;
        }
        if (velement[i] == "Hf" && velement[j] == "Pd" && velement[k] == "Sn") {
          found105 = true;
        }
        if (velement[i] == "Ir" && velement[j] == "Nb" && velement[k] == "Sn") {
          found105 = true;
        }
        if (velement[i] == "Ir" && velement[j] == "Sn" && velement[k] == "Ta") {
          found105 = true;
        }
        if (velement[i] == "La" && velement[j] == "Pt" && velement[k] == "Sb") {
          found105 = true;
        }
        if (velement[i] == "La" && velement[j] == "Rh" && velement[k] == "Te") {
          found105 = true;
        }
        if (velement[i] == "Li" && velement[j] == "Sb" && velement[k] == "Zn") {
          found105 = true;
        }
        if (velement[i] == "Na" && velement[j] == "P" && velement[k] == "Sr") {
          found105 = true;
        }
        if (velement[i] == "Na" && velement[j] == "Sb" && velement[k] == "Sr") {
          found105 = true;
        }
        if (velement[i] == "Nb" && velement[j] == "Os" && velement[k] == "Sb") {
          found105 = true;
        }
        if (velement[i] == "Nb" && velement[j] == "Rh" && velement[k] == "Sn") {
          found105 = true;
        }
        if (velement[i] == "Nb" && velement[j] == "Ru" && velement[k] == "Sb") {
          found105 = true;
        }
        if (velement[i] == "Ni" && velement[j] == "Pb" && velement[k] == "Zr") {
          found105 = true;
        }
        if (velement[i] == "Ni" && velement[j] == "Sn" && velement[k] == "Ti") {
          found105 = true;
        }
        if (velement[i] == "Ni" && velement[j] == "Sn" && velement[k] == "Zr") {
          found105 = true;
        }
        if (velement[i] == "Os" && velement[j] == "Sb" && velement[k] == "Ta") {
          found105 = true;
        }
        if (velement[i] == "Pb" && velement[j] == "Pd" && velement[k] == "Zr") {
          found105 = true;
        }
        if (velement[i] == "Rh" && velement[j] == "Sn" && velement[k] == "Ta") {
          found105 = true;
        }
        if (velement[i] == "Ru" && velement[j] == "Sb" && velement[k] == "Ta") {
          found105 = true;
        }
        if (velement[i] == "Ru" && velement[j] == "Te" && velement[k] == "Zr") {
          found105 = true;
        }
        if (found105) {
          for (size_t l = 0; l < vprotos3.size(); l++) {
            if (!aurostd::substring2bool(vprotos3[l], "As1Be1Li1_ICSD_100004") && !aurostd::substring2bool(vprotos3[l], "Be1Li1P1_ICSD_42037 ICSD_42037")) {
              if (aurostd::substring2bool(vprotos3[l], velement[i])) {
                if (aurostd::substring2bool(vprotos3[l], velement[j])) {
                  if (aurostd::substring2bool(vprotos3[l], velement[k])) {
                    vector<string> tokens;
                    vector<string> vspecies;
                    aurostd::string2tokens(vprotos3[l], tokens);
                    // oss << tokens.size() << endl;
                    XATOM_SplitAlloySpecies(tokens.at(0), vspecies);
                    // oss << vspecies.size() << endl;
                    if (vspecies.at(0) == velement[i]) {
                      if (vspecies.at(1) == velement[j]) {
                        if (vspecies.at(2) == velement[k]) {
                          bool found = false;
                          vector<string> tokensicsd;
                          for (size_t n = 0; n < vicsd.size() && !found; n++) {
                            if (aurostd::substring2bool(vicsd[n], velement[i])) {
                              if (aurostd::substring2bool(vicsd[n], velement[j])) {
                                if (aurostd::substring2bool(vicsd[n], velement[k])) {
                                  aurostd::string2tokens(vicsd[n], tokensicsd);
                                  if (tokens.at(0) == tokensicsd.at(0) && tokens.at(1) == tokensicsd.at(1) && tokens.at(3) == tokensicsd.at(3)) {
                                    found = true; // alloy composition and pearson
                                  }
                                }
                              }
                            }
                          }
                          if (!found) {
                            oss << velement[i] << velement[j] << velement[k] << " " << vprotos3[l] << endl;
                            vicsd.push_back(vprotos3[l]);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  oss << vicsd.size() << endl;
  vector<string> tokens;
  const vector<string> tokens2;
  for (size_t i = 0; i < vicsd.size(); i++) {
    aurostd::string2tokens(vicsd[i], tokens);
    oss << "/home/auro/work/AFLOW3/aflow --noldau --aflow_proto=" << tokens.at(tokens.size() - 3) << endl;
  }
  for (size_t i = 0; i < vicsd.size(); i++) {
    aurostd::string2tokens(vicsd[i], tokens);
    aurostd::string2tokens(string(tokens.at(tokens.size() - 3)), tokens, "_");
    aurostd::StringSubst(tokens.at(0), "V1", "V_sv");
    aurostd::StringSubst(tokens.at(0), "Y1", "Y_sv");
    aurostd::StringSubst(tokens.at(0), "0", "");
    aurostd::StringSubst(tokens.at(0), "1", "");
    aurostd::StringSubst(tokens.at(0), "2", "");
    aurostd::StringSubst(tokens.at(0), "3", "");
    aurostd::StringSubst(tokens.at(0), "4", "");
    aurostd::StringSubst(tokens.at(0), "5", "");
    aurostd::StringSubst(tokens.at(0), "6", "");
    aurostd::StringSubst(tokens.at(0), "7", "");
    aurostd::StringSubst(tokens.at(0), "8", "");
    aurostd::StringSubst(tokens.at(0), "9", "");

    aurostd::StringSubst(tokens.at(0), "Ba", "Ba_sv");
    aurostd::StringSubst(tokens.at(0), "Be", "Be_sv");
    aurostd::StringSubst(tokens.at(0), "Bi", "Bi_d");
    aurostd::StringSubst(tokens.at(0), "Ca", "Ca_sv");
    aurostd::StringSubst(tokens.at(0), "Cr", "Cr_pv");
    aurostd::StringSubst(tokens.at(0), "Cu", "Cu_pv");
    aurostd::StringSubst(tokens.at(0), "Fe", "Fe_pv");
    aurostd::StringSubst(tokens.at(0), "Ga", "Ga_h");
    aurostd::StringSubst(tokens.at(0), "Ge", "Ge_h");
    aurostd::StringSubst(tokens.at(0), "Hf", "Hf_pv");
    aurostd::StringSubst(tokens.at(0), "In", "In_d");
    aurostd::StringSubst(tokens.at(0), "Li", "Li_sv");
    aurostd::StringSubst(tokens.at(0), "Mg", "Mg_pv");
    aurostd::StringSubst(tokens.at(0), "Mn", "Mn_pv");
    aurostd::StringSubst(tokens.at(0), "Mo", "Mo_pv");
    aurostd::StringSubst(tokens.at(0), "Na", "Na_sv");
    aurostd::StringSubst(tokens.at(0), "Nb", "Nb_sv");
    aurostd::StringSubst(tokens.at(0), "Ni", "Ni_pv");
    aurostd::StringSubst(tokens.at(0), "Os", "Os_pv");
    aurostd::StringSubst(tokens.at(0), "Pb", "Pb_d");
    aurostd::StringSubst(tokens.at(0), "Pd", "Pd_pv");
    aurostd::StringSubst(tokens.at(0), "Re", "Re_pv");
    aurostd::StringSubst(tokens.at(0), "Rh", "Rh_pv");
    aurostd::StringSubst(tokens.at(0), "Ru", "Ru_pv");
    aurostd::StringSubst(tokens.at(0), "Sc", "Sc_sv");
    aurostd::StringSubst(tokens.at(0), "Sr", "Sr_sv");
    aurostd::StringSubst(tokens.at(0), "Ta", "Ta_pv");
    aurostd::StringSubst(tokens.at(0), "Tc", "Tc_pv");
    aurostd::StringSubst(tokens.at(0), "Ti", "Ti_sv");
    aurostd::StringSubst(tokens.at(0), "Tl", "Tl_d");
    aurostd::StringSubst(tokens.at(0), "Zr", "Zr_sv");

    oss << " ./LIB3/LIB/" << tokens.at(0) << "/" << tokens.at(1) << "_" << tokens.at(2) << ".ABC" << endl;
    //  oss << "mkdir  ./LIB3/" << " && " << "mkdir  ./LIB3/LIB/" << " && " << "mkdir  ./LIB3/LIB/" << tokens.at(0) << " && " "mkdir  ./LIB3/LIB/" << tokens.at(0) << " && " << "mkdir  ./LIB3/LIB/" << tokens.at(0) << "/" << tokens.at(1) << "_" << tokens.at(2) << ".ABC" << endl;
  }
  // oss << vprotos3.size() << endl;
}

void pinku_funcs(float x, aurostd::xvector<float>& afunc) {
  int i;
  afunc[1] = 1.0;
  afunc[2] = x;
  for (i = 3; i <= afunc.urows; i++) {
    afunc[i] = sin(i * x);
  }
}

int pinku_main() {
  int i;
  int j;
  float chisq;

  xvector<int> ia(1, NTERM);
  xvector<float> a(1, NTERM);
  xvector<float> x(1, NPT);
  xvector<float> y(1, NPT);
  xvector<float> sig(1, NPT);
  xmatrix<float> covar(1, NTERM, 1, NTERM);

  for (i = 1; i <= NPT; i++) {
    x[i] = 0.1 * i;
    pinku_funcs(x[i], a);
    y[i] = 0.0;
    for (j = 1; j <= NTERM; j++) {
      y[i] += j * a[j];
    }
    y[i] += SPREAD * aurostd::ran0();
    sig[i] = SPREAD;
  }
  for (i = 1; i <= NTERM; i++) {
    ia[i] = 1;
  }
  aurostd::lfit(x, y, sig, a, ia, covar, chisq, pinku_funcs);
  printf("\n%11s %21s\n", "parameter", "uncertainty");
  for (i = 1; i <= NTERM; i++) {
    printf("  a[%1d] = %8.6f %12.6f\n", i, a[i], sqrt(covar[i][i]));
  }
  printf("chi-squared = %12f\n", chisq);
  printf("full covariance matrix\n");
  for (i = 1; i <= NTERM; i++) {
    for (j = 1; j <= NTERM; j++) {
      printf("%12f", covar[i][j]);
    }
    printf("\n");
  }
  printf("\npress RETURN to continue...\n");
  (void) getchar();
  /* Now check results of restricting fit parameters */
  for (i = 2; i <= NTERM; i += 2) {
    ia[i] = 0;
  }
  lfit(x, y, sig, a, ia, covar, chisq, pinku_funcs);
  printf("\n%11s %21s\n", "parameter", "uncertainty");
  for (i = 1; i <= NTERM; i++) {
    printf("  a[%1d] = %8.6f %12.6f\n", i, a[i], sqrt(covar[i][i]));
  }
  printf("chi-squared = %12f\n", chisq);
  printf("full covariance matrix\n");
  for (i = 1; i <= NTERM; i++) {
    for (j = 1; j <= NTERM; j++) {
      printf("%12f", covar[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  return 0;
}

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2024              *
// *                                                                        *
// **************************************************************************

// aurostd::StringSubst(tokens.at(0),"","B_h");
// aurostd::StringSubst(tokens.at(0),"","K_sv");
// aurostd::StringSubst(tokens.at(0),"","V_sv");
// aurostd::StringSubst(tokens.at(0),"","W_pv");
// aurostd::StringSubst(tokens.at(0),"","Y_sv");
