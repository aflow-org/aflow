// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Richard H. Taylor
// UPDATED by David Hicks
// d.hicks@duke.edu

#ifndef _AFLOW_SYMMETRY_SPACEGROUP_ITC_LIBRARY_CPP_
#define _AFLOW_SYMMETRY_SPACEGROUP_ITC_LIBRARY_CPP_
#include <string>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <map>
#include "aflow_symmetry_spacegroup.h"
#include "aflow.h"

using namespace std;

//DX20190215 [OBSOLETE] #ifdef AFLOW_SYMMETRY_MULTITHREADS_ENABLE
//DX20190215 [OBSOLETE] extern thread_local double _SYM_TOL_;
//DX20190215 [OBSOLETE] #else
//DX20190215 [OBSOLETE] extern double _SYM_TOL_;
//DX20190215 [OBSOLETE] #endif

//DX20190905 [OBSOLETE - use explicit division by aurostd::modulus()] // **********************************************************************************************************************
//DX20190905 [OBSOLETE - use explicit division by aurostd::modulus()] // normalize (Overloaded)
//DX20190905 [OBSOLETE - use explicit division by aurostd::modulus()] // **********************************************************************************************************************
//DX20190905 [OBSOLETE - use explicit division by aurostd::modulus()] namespace SYM {
//DX20190905 [OBSOLETE - use explicit division by aurostd::modulus()]   void normalize(xvector<double>& v) {
//DX20190905 [OBSOLETE - use explicit division by aurostd::modulus()]     double mod = aurostd::modulus(v);
//DX20190905 [OBSOLETE - use explicit division by aurostd::modulus()]     for (int i = 1; i <= v.urows; i++) {
//DX20190905 [OBSOLETE - use explicit division by aurostd::modulus()]       v(i) = v(i) / mod;
//DX20190905 [OBSOLETE - use explicit division by aurostd::modulus()]     }
//DX20190905 [OBSOLETE - use explicit division by aurostd::modulus()]   }
//DX20190905 [OBSOLETE - use explicit division by aurostd::modulus()] } //namespace SYM
//// **********************************************************************************************************************
//// normalize (Overloaded)
//// **********************************************************************************************************************
//namespace SYM {
//void normalize(vector<double>& v) {
//double mod = aurostd::modulus(v);
//for (uint i = 0; i < v.size(); i++) {
//v[i] = v[i] / mod;
//}
//}
//} //namespace SYM

// **********************************************************************************************************************
// xvec2xmat (Overloaded)
// **********************************************************************************************************************
namespace SYM {
  xmatrix<double> xvec2xmat(xvector<double> a, xvector<double> b, xvector<double> c) {  //converts three xvectors into xmatrix
    xmatrix<double> out;
    out(1, 1) = a(1);
    out(1, 2) = a(2);
    out(1, 3) = a(3);
    out(2, 1) = b(1);
    out(2, 2) = b(2);
    out(2, 3) = b(3);
    out(3, 1) = c(1);
    out(3, 2) = c(2);
    out(3, 3) = c(3);
    return out;
  }
} //namespace SYM

// **********************************************************************************************************************
// xvect2xmat (Overloaded)
// **********************************************************************************************************************
namespace SYM {
  xmatrix<double> xvec2xmat(vector<xvector<double> > V) {
    uint vectorLength = (uint)V[0].urows;
    xmatrix<double> M(V.size(), V[0].urows);
    for (uint i = 0; i < V.size(); i++) {
      for (uint j = 0; j < vectorLength; j++) {
        M(i + 1, j + 1) = V[i][j + 1];
      }
    }
    return M;
  }
} //namespace SYM

// **********************************************************************************************************************
// xvect2xmat (Overloaded)
// **********************************************************************************************************************
namespace SYM {
  xmatrix<double> xvec2xmat(vector<xvector<double> > V, vector<double> R) {
    uint vectorLength = (uint)(V[0].urows + 1);
    xmatrix<double> M(V.size(), vectorLength);
    for (uint i = 1; i <= V.size(); i++) {
      for (uint j = 1; j <= vectorLength - 1; j++) {
        M(i, j) = V[i - 1][j];
      }
      M(i, vectorLength) = R[i - 1];
    }
    return M;
  }
} //namespace SYM

// **********************************************************************************************************************
// extract_row
// **********************************************************************************************************************
namespace SYM {
  xvector<double> extract_row(xmatrix<double> a, int row) {
    xvector<double> out;
    out(1) = a(row, 1);
    out(2) = a(row, 2);
    out(3) = a(row, 3);
    return out;
  }
} //namespace SYM

// **********************************************************************************************************************
// get_point_in_line
// **********************************************************************************************************************
namespace SYM {
  xvector<double> get_point_in_line(string line, double param) {
    vector<string> linestringvec = splitstring_2(line);
    for (uint i = 0; i < linestringvec.size(); i++) {
      cleanupstring(linestringvec[i]);
    }
    xvector<double> point_out;
    for (uint i = 0; i < linestringvec.size(); i++) {
      double agg = 0;
      vector<sdouble> sd = simplify(linestringvec[i]);
      for (uint j = 0; j < sd.size(); j++) {
        if(sd[j].chr != '\0') {
          agg += param * sd[j].dbl;
        } else {
          agg += sd[j].dbl;
        }
      }
      point_out[i + 1] = agg;
    }
    return point_out;
  }
} //namespace SYM

// **********************************************************************************************************************
// whereischar
// **********************************************************************************************************************
namespace SYM {
  int whereischar(string str, char c) {
    int out = -1;
    for (uint i = 0; i < str.size(); i++) {
      if(str[i] == c)
        out = i;
    }
    return out;
  }
} //namespace SYM

// **********************************************************************************************************************
// get_triplet
// **********************************************************************************************************************
namespace SYM {
  xvector<double> get_triplet(string str) {
    //extracts value inside parentehsis, ignoring everything else
    //and converts these string values to doubles

    ostringstream tstring;

    xvector<double> out;
    uint start = 0;
    for (uint i = 0; i < str.size(); i++) {
      if(str[i] == '(') {
        start = i;
      }
    }
    uint i = start + 1;
    while (str[i + 1] != ')') {
      //tstring << ITCstring[i+1];
      tstring << str[i];
      i++;
    }
    tstring << str[i];
    vector<string> vs = splitstring(tstring.str());

    for (uint i = 0; i < vs.size(); i++)
      //DX20200313 [OBSOLETE] out(i + 1) = frac2dbl(vs[i]);
      out(i + 1) = aurostd::frac2dbl(vs[i]); //DX20200313

    return out;
  }
} //namespace SYM

// **********************************************************************************************************************
// Inversion Class
// **********************************************************************************************************************
//void Inversion::get_inversion(string ITCstring) {
//xvector<double> outpoint;
//uint start = 0;
//for (uint i = 0; i < ITCstring.size(); i++) {
//if(ITCstring[i + 1] == 'I') {
//start = i + 1;
//}
//}
//ITCstring.erase(0, start + 1);
//vector<string> tmpvec = splitstring(ITCstring);
//
//for (uint i = 0; i < tmpvec.size(); i++) {
//outpoint[i + 1] = aurostd::frac2dbl(tmpvec[i]);
//}
//inversion_point = outpoint;
//}
//xvector<double> operator*(Inversion& I, const xvector<double>& P) {
//return 2 * I.inversion_point - P;
//}
//
//xvector<double> operator*(const xvector<double>& P, Inversion& I) {
//return 2 * I.inversion_point - P;
//}
//END INVERSION CLASS
// **********************************************************************************************************************

// **********************************************************************************************************************
// Translation Class
// **********************************************************************************************************************
//void Translation::get_translation(string ITCstring) {
//uint start = 0;
//for (uint i = 0; i < ITCstring.size(); i++) {
//if(ITCstring[i] == '(') {
//start = i;
//}
//}
//ITCstring.erase(0, start + 1);
//uint i = 0;
//ostringstream oss;
//while (ITCstring[i] != ')') {
//oss << ITCstring[i];
//i++;
//}
//vector<string> vecstring = splitstring(oss.str());
//for (uint i = 0; i < vecstring.size(); i++) {
//translation_vec[i + 1] = aurostd::frac2dbl(vecstring[i]);
//}
////cerr << translation_vec << endl;
//}
//
//xvector<double> operator*(Translation& T, const xvector<double>& P) {
//return P + T.translation_vec;
//}
//xvector<double> operator*(const xvector<double>& P, Translation& T) {
//return P + T.translation_vec;
//}
//xvector<double> operator+(Translation& T, const xvector<double>& P) {
//return P + T.translation_vec;
//}
//xvector<double> operator+(const xvector<double>& P, Translation& T) {
//return P + T.translation_vec;
//}
//END TRANSLATION CLASS
// **********************************************************************************************************************


// **********************************************************************************************************************
// Screw Class
// **********************************************************************************************************************
namespace SYM{
  void Screw::get_line_string(string str) {
    linestring = str;
  }
  void Screw::get_trans_comp(xvector<double> trans) {
    T = trans;
  }
  void Screw::get_angle(double theta) {
    angle = theta;
  }
  void Screw::get_A() {
    xmatrix<double> temp(4, 4);
    vector<string> linestringvector = splitstring(linestring);
    xvector<double> shift;  //vector of doubles that occur when you have e.g., x+1/
    ostringstream oss;

    xvector<double> point_1 = get_point_in_line(linestring, 1);
    xvector<double> point_2 = get_point_in_line(linestring, 0);

    xvector<double> dir_vec = point_2 - point_1;
    //cerr << "dir vec: " << dir_vec << endl;
    double n = aurostd::modulus(dir_vec);
    //DX20191028 [OBSOLETE] if(aurostd::abs(n) < 1e-8) { cerr << "SYM::Screw::get_A(): PROBLEM WITH DIR VEC ( Screw::get_A() ) " << endl; }
    //DX20191028 - if modulus is zero, then do no rotate - START
    if(aurostd::abs(n) < _ZERO_TOL_) {
      A = aurostd::eye<double>(4,4);
      return;
    }
    //DX20191028 - if modulus is zero, then do no rotate - END
    dir_vec = 1 / n * dir_vec;
    one_point = point_2;  //a point on the line to pass around
    //double u = dir_vec[1];
    //double v = dir_vec[2];
    //double w = dir_vec[3];
    //double a = point_2[1];
    //double b = point_2[2];
    //double c = point_2[3];
    //cerr << "(1,4) " << (a*(v*v+w*w)-u*(b*v+c*w))*(1-cos(angle))+(b*w-c*v)*sin(angle) << endl;
    //cerr << "shift: " << shift << endl;
    //cerr << "point 1: " << point_1 << endl;
    //Rotation of angle about an arbitrary line: from https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas
    temp(1, 1) = dir_vec[1] * dir_vec[1] + (dir_vec[2] * dir_vec[2] + dir_vec[3] * dir_vec[3]) * cos(angle);
    temp(1, 2) = dir_vec[1] * dir_vec[2] * (1 - cos(angle)) - dir_vec[3] * sin(angle);
    temp(1, 3) = dir_vec[1] * dir_vec[3] * (1 - cos(angle)) + dir_vec[2] * sin(angle);
    temp(1, 4) = (point_1[1] * (dir_vec[2] * dir_vec[2] + dir_vec[3] * dir_vec[3]) - dir_vec[1] * (point_1[2] * dir_vec[2] + point_1[3] * dir_vec[3])) * (1 - cos(angle)) + (point_1[2] * dir_vec[3] - point_1[3] * dir_vec[2]) * sin(angle);
    temp(2, 1) = dir_vec[1] * dir_vec[2] * (1 - cos(angle)) + dir_vec[3] * sin(angle);
    temp(2, 2) = dir_vec[2] * dir_vec[2] + (dir_vec[1] * dir_vec[1] + dir_vec[3] * dir_vec[3]) * cos(angle);
    temp(2, 3) = dir_vec[2] * dir_vec[3] * (1 - cos(angle)) - dir_vec[1] * sin(angle);
    temp(2, 4) = (point_1[2] * (dir_vec[1] * dir_vec[1] + dir_vec[3] * dir_vec[3]) - dir_vec[2] * (point_1[1] * dir_vec[1] + point_1[3] * dir_vec[3])) * (1 - cos(angle)) + (point_1[3] * dir_vec[1] - point_1[1] * dir_vec[3]) * sin(angle);
    temp(3, 1) = dir_vec[1] * dir_vec[3] * (1 - cos(angle)) - dir_vec[2] * sin(angle);
    temp(3, 2) = dir_vec[2] * dir_vec[3] * (1 - cos(angle)) + dir_vec[1] * sin(angle);
    temp(3, 3) = dir_vec[3] * dir_vec[3] + (dir_vec[1] * dir_vec[1] + dir_vec[2] * dir_vec[2]) * cos(angle);
    temp(3, 4) = (point_1[3] * (dir_vec[1] * dir_vec[1] + dir_vec[2] * dir_vec[2]) - dir_vec[3] * (point_1[1] * dir_vec[1] + point_1[2] * dir_vec[2])) * (1 - cos(angle)) + (point_1[1] * dir_vec[2] - point_1[2] * dir_vec[1]) * sin(angle);
    temp(4, 1) = 0;
    temp(4, 2) = 0;
    temp(4, 3) = 0;
    temp(4, 4) = 1;
    A = temp;
  }

  void Screw::get_screw_direct(xvector<double> axis_dir, xvector<double> point, double order) {
    xmatrix<double> temp(4, 4);
    xvector<double> shift;  //vector of doubles that occur when you have e.g., x+1/
    ostringstream oss;
    xvector<double> dir_vec = axis_dir;
    double n = aurostd::modulus(dir_vec);
    angle = 2 * Pi_r / order;
    //DX20191028 [OBSOLETE] if(aurostd::abs(n) < 1e-8) { cerr << "SYM::get_screw_direct: PROBLEM WITH DIR VEC ( Screw::get_screw_direct() )" << endl; }
    //DX20191028 - if modulus is zero, then do no rotate - START
    if(aurostd::abs(n) < _ZERO_TOL_) {
      A = aurostd::eye<double>(4,4);
      return;
    }
    //DX20191028 - if modulus is zero, then do no rotate - END
    dir_vec = 1 / n * dir_vec;
    direction_vector = dir_vec;
    one_point = point;  //a point on the line to pass around
    xvector<double> point_1 = point;
    temp(1, 1) = dir_vec[1] * dir_vec[1] + (dir_vec[2] * dir_vec[2] + dir_vec[3] * dir_vec[3]) * cos(angle);
    temp(1, 2) = dir_vec[1] * dir_vec[2] * (1 - cos(angle)) - dir_vec[3] * sin(angle);
    temp(1, 3) = dir_vec[1] * dir_vec[3] * (1 - cos(angle)) + dir_vec[2] * sin(angle);
    temp(1, 4) = (point_1[1] * (dir_vec[2] * dir_vec[2] + dir_vec[3] * dir_vec[3]) - dir_vec[1] * (point_1[2] * dir_vec[2] + point_1[3] * dir_vec[3])) * (1 - cos(angle)) + (point_1[2] * dir_vec[3] - point_1[3] * dir_vec[2]) * sin(angle);
    temp(2, 1) = dir_vec[1] * dir_vec[2] * (1 - cos(angle)) + dir_vec[3] * sin(angle);
    temp(2, 2) = dir_vec[2] * dir_vec[2] + (dir_vec[1] * dir_vec[1] + dir_vec[3] * dir_vec[3]) * cos(angle);
    temp(2, 3) = dir_vec[2] * dir_vec[3] * (1 - cos(angle)) - dir_vec[1] * sin(angle);
    temp(2, 4) = (point_1[2] * (dir_vec[1] * dir_vec[1] + dir_vec[3] * dir_vec[3]) - dir_vec[2] * (point_1[1] * dir_vec[1] + point_1[3] * dir_vec[3])) * (1 - cos(angle)) + (point_1[3] * dir_vec[1] - point_1[1] * dir_vec[3]) * sin(angle);
    temp(3, 1) = dir_vec[1] * dir_vec[3] * (1 - cos(angle)) - dir_vec[2] * sin(angle);
    temp(3, 2) = dir_vec[2] * dir_vec[3] * (1 - cos(angle)) + dir_vec[1] * sin(angle);
    temp(3, 3) = dir_vec[3] * dir_vec[3] + (dir_vec[1] * dir_vec[1] + dir_vec[2] * dir_vec[2]) * cos(angle);
    temp(3, 4) = (point_1[3] * (dir_vec[1] * dir_vec[1] + dir_vec[2] * dir_vec[2]) - dir_vec[3] * (point_1[1] * dir_vec[1] + point_1[2] * dir_vec[2])) * (1 - cos(angle)) + (point_1[1] * dir_vec[2] - point_1[2] * dir_vec[1]) * sin(angle);
    temp(4, 1) = 0;
    temp(4, 2) = 0;
    temp(4, 3) = 0;
    temp(4, 4) = 1;
    A = temp;
  }

  xvector<double> Screw::return_direction() {
    xvector<double> xvd;
    xvd = direction_vector;
    return xvd;
  }

  xvector<double> Screw::return_point() {
    xvector<double> xvd;
    xvd = one_point;
    return xvd;
  }

  double Screw::return_order() {
    double out = 0;
    out = 2 * Pi_r / angle;
    return out;
  }
}

// END SCREW CLASS
// **********************************************************************************************************************

// **********************************************************************************************************************
// print (Overloaded)
// **********************************************************************************************************************
void print(vector<int> vec) {
  for (unsigned int i = 0; i < vec.size(); i++) {
    cerr << vec.at(i) << " ";
  }
  cerr << endl;
}

// **********************************************************************************************************************
// print (Overloaded)
// **********************************************************************************************************************
void print(const vector<xvector<double> >& points) {
  for (uint i = 0; i < points.size(); i++) {
    cerr << points[i] << endl;
  }
}

// **********************************************************************************************************************
// print (Overloaded)
// **********************************************************************************************************************
void print(const vector<vector<double> >& points) {
  for (uint i = 0; i < points.size(); i++) {
    for (uint j = 0; j < points[i].size(); j++) {
      cerr << points[i][j] << " ";
    }
    cerr << endl;
  }
}

// **********************************************************************************************************************
// print (Overloaded)
// **********************************************************************************************************************
//Transform matrix to go from direct to cart and vice versa.
void print(const vector<xvector<double> >& points, xmatrix<double> T) {
  for (uint i = 0; i < points.size(); i++) {
    cerr << points[i] * T << endl;
  }
}

// **********************************************************************************************************************
// print (Overloaded)
// **********************************************************************************************************************
void print(const vector<xmatrix<double> >& mats) {
  for (uint i = 0; i < mats.size(); i++) {
    cerr << mats[i] << endl;
    cerr << endl;
  }
}

// **********************************************************************************************************************
// print (Overloaded)
// **********************************************************************************************************************
void print(const vector<_atom>& atoms) {
  for (uint i = 0; i < atoms.size(); i++) {
    //ORIG cerr << atoms[i] << endl;
    cerr << atoms[i].fpos << " " << atoms[i].name << endl;
  }
}

// **********************************************************************************************************************
// print (Overloaded)
// **********************************************************************************************************************
void print(const deque<_atom>& atoms) {
  for (uint i = 0; i < atoms.size(); i++) {
    //ORIG cerr << atoms[i] << endl;
    cerr << atoms[i].fpos << " " << atoms[i].name << endl;
  }
}

// **********************************************************************************************************************
// Screw Class Friends
// **********************************************************************************************************************
namespace SYM {
  ostream& operator<<(ostream& output, const Screw& S) {
    //output << "(" <<  S.A << ", " << S.T <<")";
    output << "(" << S.direction_vector << " ) with trans: (" << S.T << " )"
      << " ord: " << 360.0 / (S.angle * (180 / Pi_r)) << " at point: " << S.one_point << endl;
    return output;  // for multiple << operators.
  }

  xvector<double> operator*(Screw& S, const xvector<double>& P) {
    xvector<double> temp(4);
    xvector<double> out;
    temp[1] = P[1];
    temp[2] = P[2];
    temp[3] = P[3];
    temp[4] = 1;
    //S.get_A();//maybe put this in get_angle or somewhere else.
    out[1] = (S.A * temp)[1] + S.T[1];
    out[2] = (S.A * temp)[2] + S.T[2];
    out[3] = (S.A * temp)[3] + S.T[3];
    return out;
  }
  xvector<double> operator*(const xvector<double>& P, Screw& S) {
    xvector<double> temp(4);
    xvector<double> out;
    temp[1] = P[1];
    temp[2] = P[2];
    temp[3] = P[3];
    temp[4] = 1;
    //S.get_A();//maybe put this in get_angle or somewhere else.
    out[1] = (S.A * temp)[1] + S.T[1];
    out[2] = (S.A * temp)[2] + S.T[2];
    out[3] = (S.A * temp)[3] + S.T[3];
    return out;
  }
} //namespace SYM

// **********************************************************************************************************************
// Glide Class
// **********************************************************************************************************************
namespace SYM {
  ostream& operator<<(ostream& output, const Glide& G) {
    if(G.DIRECT == true) {
      output << "(" << G.a << " ) on (" << G.plane_point << " )" << endl;
    } else {
      output << "reflecting matrix: " << G.A << endl;
      output << "translation vector: " << G.T << endl;
    }
    return output;  // for multiple << operators.
  }

  xvector<double> operator*(Glide& G, const xvector<double>& P) {
    stringstream message;
    if(G.DIRECT == true) {
      //cerr << "DIRECT = true" << endl;
      double tol = 0.0001;
      double orient = 1;
      xvector<double> x_0 = G.plane_point;
      double d1 = aurostd::distance(x_0, P); //DX20210422 - use aurostd
      double d2 = aurostd::distance(x_0 + tol * G.a, P); //DX20210422 - use aurostd
      if(d2 == d1) {
        message << "THERE IS A PROBLEM in operator* Glide(DIRECT=true): " << endl;
        message << "d1,d2: " << d1 << " " << d2 << endl;
        message << "x_0: " << x_0 << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_GENERIC_ERROR_);
      }
      if(d2 < d1) {  //Normal not oriented properly, switch sign. (G.a -> -G.a) et.c
        orient = -1;
        //cerr << "Not OP" << endl;
      }
      xvector<double> norm;
      norm = orient * G.a;
      norm = norm / aurostd::scalar_product(norm, norm);
      G.get_A(norm);
      //double m = aurostd::modulus(norm);
      double d = orient * G.d;
      //return (1/(m*m))*(m*G.A*P+G.d*orient*G.a+aurostd::abs(aurostd::scalar_product(orient*G.a,P)-G.d)*orient*G.a)+G.T;
      //return (1/(m*m))*(m*G.A*P+d*norm+aurostd::abs(aurostd::scalar_product(norm,P)-d)*norm)+G.T;
      return G.A * P + d * norm + aurostd::abs(aurostd::scalar_product(norm, P) - d) * norm + G.T;
    }
    if(G.HEX == true) {
      xvector<double> out(3);
      xvector<double> tmp;
      xmatrix<double> identity(3, 3);
      identity(1, 1) = 1;
      identity(2, 2) = 1;
      identity(3, 3) = 1;
      out = G.A * P + G.T + (identity - G.A) * G.a;
      return out;
    } else {
      double tol = 0.0001;
      double orient = 1;
      xvector<double> x_0 = get_random_point_in_plane(G.planestring);
      double d1 = aurostd::distance(x_0, P); //DX20210422 - use aurostd
      double d2 = aurostd::distance(x_0 + tol * G.a, P); //DX20210422 - use aurostd
      if(d2 == d1) { cerr << "THERE IS A PROBLEM: " << endl; }
      if(d2 < d1) {  //Normal not oriented properly, switch sign. (G.a -> -G.a) et.c
        orient = -1;
        //cerr << "Not OP" << endl;
      }
      xvector<double> norm;
      norm = orient * G.a;
      norm = norm / aurostd::scalar_product(norm, norm);
      G.get_A(norm);
      //double m = aurostd::modulus(norm);
      double d = orient * G.d;
      //return (1/(m*m))*(m*G.A*P+G.d*orient*G.a+aurostd::abs(aurostd::scalar_product(orient*G.a,P)-G.d)*orient*G.a)+G.T;
      //return (1/(m*m))*(m*G.A*P+d*norm+aurostd::abs(aurostd::scalar_product(norm,P)-d)*norm)+G.T;
      return G.A * P + d * norm + aurostd::abs(aurostd::scalar_product(norm, P) - d) * norm + G.T;
    }
  }
  xvector<double> operator*(const xvector<double>& P, Glide& G) {
    double tol = 0.0001;
    double orient = 1;
    xvector<double> x_0 = get_random_point_in_plane(G.planestring);
    double d1 = aurostd::distance(x_0, P); //DX20210422 - used aurostd
    double d2 = aurostd::distance(x_0 + tol * G.a, P); //DX20210422 - use aurostd
    if(d2 == d1) { cerr << "THERE IS A PROBLEM: " << endl; }
    if(d2 < d1) {  //Normal not oriented properly, switch sign. (G.a -> -G.a) et.c
      orient = -1;
      //cerr << "Not OP" << endl;
    }
    xvector<double> norm;
    norm = orient * G.a;
    norm = norm / aurostd::scalar_product(norm, norm);
    G.get_A(norm);
    //double m = aurostd::modulus(norm);
    double d = orient * G.d;
    //return (1/(m*m))*(m*G.A*P+G.d*orient*G.a+aurostd::abs(aurostd::scalar_product(orient*G.a,P)-G.d)*orient*G.a)+G.T;
    //return (1/(m*m))*(m*G.A*P+d*norm+aurostd::abs(aurostd::scalar_product(norm,P)-d)*norm)+G.T;
    return G.A * P + d * norm + aurostd::abs(aurostd::scalar_product(norm, P) - d) * norm + G.T;
  }
} //namespace SYM

// END GLIDE CLASS
// **********************************************************************************************************************

// **********************************************************************************************************************
// xb
// **********************************************************************************************************************
void xb() {
  cerr << "************************************" << endl;
}

// **********************************************************************************************************************
// get_angle
// **********************************************************************************************************************
namespace SYM {
  double get_angle(xvector<double> a, xvector<double> b, string c) {
    double ang = acos(aurostd::scalar_product(a, b) / (aurostd::modulus(a) * aurostd::modulus(b)));
    if(c == "deg" || c == "Deg" || c == "d" || c == "D") {
      ang = (180 / Pi_r) * ang;
    }

    return ang;
  }
} //namespace SYM

// **********************************************************************************************************************
// get_random_point_in_plane
// **********************************************************************************************************************
namespace SYM {
  xvector<double> get_random_point_in_plane(string plane) {
    stringstream message;
    vector<string> parametric_plane = splitstring(plane);
    vector<char> params;
    char c;
    bool found_params=false; //DX20210422
    for (uint i = 0; i < parametric_plane.size(); i++) {
      c = whichchar(parametric_plane[i]);
      if(c != '\0') {
        found_params=false; //DX20210422
        for(uint j=0;j<params.size()&&!found_params;j++){ //DX20210422
          if(c == params[j]){ found_params = true; }
        }
        if(!found_params){
          params.push_back(c);
        }
      }
    }
    if(params.size() > 2) {
      message << "Not a valid representation of a plane.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }
    vector<double> x, y;
    x.push_back((double)rand() / RAND_MAX);
    x.push_back((double)rand() / RAND_MAX);
    x.push_back((double)rand() / RAND_MAX);
    y.push_back((double)rand() / RAND_MAX);
    y.push_back((double)rand() / RAND_MAX);
    y.push_back((double)rand() / RAND_MAX);
    xvector<double> point(3);
    char chr;
    double dbl;
    uint j = 0;
    for (uint i = 0; i < parametric_plane.size(); i++) {
      chr = simplify(parametric_plane[i])[0].chr;
      dbl = simplify(parametric_plane[i])[0].dbl;
      if(chr == params[0]) {
        point(i + 1) = (x[j] * dbl);
      } else {
        if(chr == params[1]) {
          point(i + 1) = (y[j] * dbl);
        } else {
          point(i + 1) = dbl;
        }
      }
    }
    return point;
  }
} //namespace SYM

// **********************************************************************************************************************
// a2m3x3
// **********************************************************************************************************************
namespace SYM {
  xmatrix<double> a2m3x3(double* array) {  //convert a 16 length array into 3x3 xmatrix
    xmatrix<double> S(3, 3);
    S(1, 1) = array[0];
    S(1, 2) = array[1];
    S(1, 3) = array[2];
    S(2, 1) = array[4];
    S(2, 2) = array[5];
    S(2, 3) = array[6];
    S(3, 1) = array[8];
    S(3, 2) = array[9];
    S(3, 3) = array[10];
    return S;
  }
} //namespace SYM

// **********************************************************************************************************************
// Glide::get_glide_direct (Overloaded)
// **********************************************************************************************************************
namespace SYM {
  void Glide::get_glide_direct(xvector<double> n, xvector<double> p) {
    DIRECT = true;
    a = 1 / aurostd::modulus(n) * n;
    xvector<double> translation;
    plane_point = p;
    T = translation;
  }
} //namespace SYM

// **********************************************************************************************************************
// Glide::get_glide_direct (Overloaded)
// **********************************************************************************************************************
namespace SYM {
  void Glide::get_glide_direct(xvector<double> n, xvector<double> p, xvector<double> trans) {
    //n = n;
    n = n+n;  //CO20190329 - keep n busy, it's not needed //compiler doesn't like n=n;
    DIRECT = true;
    plane_point = p;
    T = trans;
  }
} //namespace SYM

// **********************************************************************************************************************
// Glide::return_point
// **********************************************************************************************************************
namespace SYM {
  xvector<double> Glide::return_point() {
    return plane_point;
  }
} //namespace SYM

// **********************************************************************************************************************
// Glide::return_direction
// **********************************************************************************************************************
namespace SYM {
  xvector<double> Glide::return_direction() {
    return a;
  }
} //namespace SYM

// **********************************************************************************************************************
// Glide::get_A
// **********************************************************************************************************************
namespace SYM {
  void Glide::get_A(xvector<double> n) {  //n is plane normal
    //cerr << "plane normal: " << n << endl;
    //xmatrix<double> A(3);
    A(1, 1) = n(2) * n(2) + n(3) * n(3);
    A(1, 2) = -n(1) * n(2);
    A(1, 3) = -n(1) * n(3);
    A(2, 1) = -n(1) * n(2);
    A(2, 2) = n(1) * n(1) + n(3) * n(3);
    A(2, 3) = -n(3) * n(2);
    A(3, 1) = -n(3) * n(1);
    A(3, 2) = -n(3) * n(2);
    A(3, 3) = n(1) * n(1) + n(2) * n(2);
    //cerr << "from get_glide: " <<endl << A << endl;

    //return A;
  }
} //namespace SYM

// **********************************************************************************************************************
// reduce_atom_deques
// **********************************************************************************************************************
namespace SYM {
  void reduce_atom_deques(deque<_atom>& expanded, xmatrix<double>& lattice, double& min_dist, double& sym_tol) { //DX20190215 -added sym_tol
    deque<_atom> tmpvec;
    for (uint i = 0; i < expanded.size(); i++) {
      for (int j = 1; j < 4; j++) {
        //DX20190905 [OBSOLETE - no more mod_one] expanded[i].fpos[j] = SYM::mod_one(expanded[i].fpos[j]);
        BringInCellInPlace(expanded[i].fpos[j]); //DX20190905
      }
    }
    //double min_mod = find_min_lattice_vector(xvec2xmat(expanded[0],expanded[1],expanded[2]));
    xmatrix<double> f2c = trasp(lattice);
    //DX20190619 [OBSOLETE] xmatrix<double> c2f = inverse(trasp(lattice));
    bool skew = SYM::isLatticeSkewed(lattice, min_dist, sym_tol); //DX20190215

    for (uint i = 0; i < expanded.size(); i++) {
      if(!SYM::MapAtom(tmpvec, expanded[i], TRUE, lattice, f2c, skew, sym_tol)) {  //CAN I USE JUST 1 HERE //DX20190215 //DX20190619 - lattice and f2c as input
        tmpvec.push_back(expanded[i]);
      }
    }
    expanded.clear();
    for (uint i = 0; i < tmpvec.size(); i++) {
      expanded.push_back(tmpvec[i]);
    }
  }
} //namespace SYM
//DX20190215 [OBSOLETE] namespace SYM {
//DX20190215 [OBSOLETE] #ifdef AFLOW_SYMMETRY_MULTITHREADS_ENABLE
//DX20190215 [OBSOLETE]thread_local vector<xvector<double> > glideplanes;
//DX20190215 [OBSOLETE]thread_local vector<xvector<double> > glideplanes_hex;
//DX20190215 [OBSOLETE]thread_local vector<xvector<double> > glidetrans;
//DX20190215 [OBSOLETE]thread_local vector<xvector<double> > glidetrans_hex;
//DX20190215 [OBSOLETE]thread_local vector<string> glidesymbols;
//DX20190215 [OBSOLETE]thread_local vector<string> glidesymbols_hex;
//DX20190215 [OBSOLETE]#else
//DX20190215 [OBSOLETE]vector<xvector<double> > glideplanes;
//DX20190215 [OBSOLETE]vector<xvector<double> > glideplanes_hex;
//DX20190215 [OBSOLETE]vector<xvector<double> > glidetrans;
//DX20190215 [OBSOLETE]vector<xvector<double> > glidetrans_hex;
//DX20190215 [OBSOLETE]vector<string> glidesymbols;
//DX20190215 [OBSOLETE]vector<string> glidesymbols_hex;
//DX20190215 [OBSOLETE]#endif
//DX20190215 [OBSOLETE]} //namespace SYM
// **********************************************************************************************************************
// initglides
// **********************************************************************************************************************
bool SymmetryInformationITC::initglides() {
  glideplanes.clear();
  glideplanes_hex.clear();
  glidetrans.clear();
  glidetrans_hex.clear();
  glidesymbols.clear();
  glidesymbols_hex.clear();
  xvector<double> a;
  a(1) = 1;
  xvector<double> b;
  b(2) = 1;
  xvector<double> c;
  c(3) = 1;

  xvector<double> tmp;
  //AXIAL GLIDES
  tmp(2) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * a);
  glidesymbols.push_back("a");

  tmp(3) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * a);
  glidesymbols.push_back("a");

  tmp(3) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * b);
  glidesymbols.push_back("b");

  tmp(1) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * b);
  glidesymbols.push_back("b");

  tmp(1) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * c);
  glidesymbols.push_back("c");

  tmp(2) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * c);
  glidesymbols.push_back("c");

  tmp(1) = 1;
  tmp(2) = -1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * c);
  glidesymbols.push_back("c");

  tmp(1) = 1;
  tmp(2) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * c);
  glidesymbols.push_back("c");

  //HEXAGONAL GLIDE PLANES
  tmp(1) = 1;
  tmp(2) = 0;  //HEX
  glideplanes_hex.push_back(tmp);
  tmp.clear();
  glidetrans_hex.push_back(.5 * c);
  glidesymbols_hex.push_back("c");

  tmp(1) = 0;
  tmp(2) = 1;  //HEX
  glideplanes_hex.push_back(tmp);
  tmp.clear();
  glidetrans_hex.push_back(.5 * c);
  glidesymbols_hex.push_back("c");

  tmp(1) = -1;
  tmp(2) = -1;  //HEX
  glideplanes_hex.push_back(tmp);
  tmp.clear();
  glidetrans_hex.push_back(.5 * c);
  glidesymbols_hex.push_back("c");

  tmp(1) = 1;
  tmp(2) = -1;  //HEX
  glideplanes_hex.push_back(tmp);
  tmp.clear();
  glidetrans_hex.push_back(.5 * c);
  glidesymbols_hex.push_back("c");

  tmp(1) = 1;
  tmp(2) = 2;  //HEX
  glideplanes_hex.push_back(tmp);
  tmp.clear();
  glidetrans_hex.push_back(.5 * c);
  glidesymbols_hex.push_back("c");

  tmp(1) = -2;
  tmp(2) = -1;  //HEX
  glideplanes_hex.push_back(tmp);
  tmp.clear();
  glidetrans_hex.push_back(.5 * c);
  glidesymbols_hex.push_back("c");

  //DIAGONAL GLIDES
  tmp(3) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * (a + b));
  glidesymbols.push_back("n");

  tmp(1) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * (b + c));
  glidesymbols.push_back("n");

  tmp(2) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * (a + c));
  glidesymbols.push_back("n");

  tmp(1) = 1;
  tmp(2) = -1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * (a + b + c));
  glidesymbols.push_back("n");

  tmp(2) = 1;
  tmp(3) = -1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * (a + b + c));
  glidesymbols.push_back("n");

  tmp(1) = -1;
  tmp(3) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * (a + b + c));
  glidesymbols.push_back("n");

  tmp(1) = 1;
  tmp(2) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * (-a + b + c));
  glidesymbols.push_back("n");

  tmp(2) = 1;
  tmp(3) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * (a - b + c));
  glidesymbols.push_back("n");

  tmp(1) = 1;
  tmp(3) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.5 * (a + b - c));
  glidesymbols.push_back("n");

  //DIAMOND GLIDES
  tmp(3) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (a + b));
  glidesymbols.push_back("d");

  tmp(3) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (a - b));
  glidesymbols.push_back("d");

  tmp(1) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (b + c));
  glidesymbols.push_back("d");

  tmp(1) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (b - c));
  glidesymbols.push_back("d");

  tmp(2) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (a + c));
  glidesymbols.push_back("d");

  tmp(2) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (-a + c));
  glidesymbols.push_back("d");

  tmp(1) = 1;
  tmp(2) = -1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (a + b + c));
  glidesymbols.push_back("d");

  tmp(1) = 1;
  tmp(2) = -1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (a + b - c));
  glidesymbols.push_back("d");

  tmp(2) = 1;
  tmp(3) = -1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (a + b + c));
  glidesymbols.push_back("d");

  tmp(2) = 1;
  tmp(3) = -1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (-a + b + c));
  glidesymbols.push_back("d");

  tmp(1) = -1;
  tmp(3) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (a + b + c));
  glidesymbols.push_back("d");

  tmp(1) = -1;
  tmp(3) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (a - b + c));
  glidesymbols.push_back("d");

  tmp(1) = 1;
  tmp(2) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (-a + b + c));
  glidesymbols.push_back("d");

  tmp(1) = 1;
  tmp(2) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (-a + b - c));
  glidesymbols.push_back("d");

  tmp(2) = 1;
  tmp(3) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (a - b + c));
  glidesymbols.push_back("d");

  tmp(2) = 1;
  tmp(3) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (-a - b + c));
  glidesymbols.push_back("d");

  tmp(1) = 1;
  tmp(3) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (a + b - c));
  glidesymbols.push_back("d");

  tmp(1) = 1;
  tmp(3) = 1;
  glideplanes.push_back(tmp);
  tmp.clear();
  glidetrans.push_back(.25 * (a - b - c));
  glidesymbols.push_back("d");

  return true;
}
//DX20190215 [OBSOLETE]} //namespace SYM

//sym_mats gives the symmetry matrices. The directions are specified by the map
//DX20190215 [OBSOLETE]namespace SYM {
//DX20190215 [OBSOLETE]  //counting FROM ZERO
//DX20190215 [OBSOLETE]#ifdef AFLOW_SYMMETRY_MULTITHREADS_ENABLE
//DX20190215 [OBSOLETE]  thread_local vector<int> index_cubic;
//DX20190215 [OBSOLETE]  thread_local vector<int> index_hex;  //To be used with the *_hex vectors
//DX20190215 [OBSOLETE]  thread_local vector<int> index_rhom;
//DX20190215 [OBSOLETE]  thread_local vector<int> index_tetr;
//DX20190215 [OBSOLETE]  thread_local vector<int> index_ortho;
//DX20190215 [OBSOLETE]  thread_local vector<int> index_mono_b;
//DX20190215 [OBSOLETE]  thread_local vector<int> index_mono_c;
//DX20190215 [OBSOLETE]  thread_local vector<int> index_tric;
//DX20190215 [OBSOLETE]
//DX20190215 [OBSOLETE]  thread_local vector<xmatrix<double> > sym_mats;
//DX20190215 [OBSOLETE]  thread_local vector<string> symbol;
//DX20190215 [OBSOLETE]  thread_local vector<string> dirparam;
//DX20190215 [OBSOLETE]  thread_local hash sym_mats_direction;
//DX20190215 [OBSOLETE]
//DX20190215 [OBSOLETE]  thread_local vector<xmatrix<double> > sym_mats_hex;
//DX20190215 [OBSOLETE]  thread_local vector<string> symbol_hex;
//DX20190215 [OBSOLETE]  thread_local vector<string> dirparam_hex;
//DX20190215 [OBSOLETE]  thread_local hash sym_mats_direction_hex;
//DX20190215 [OBSOLETE]#else
//DX20190215 [OBSOLETE]  vector<int> index_cubic;
//DX20190215 [OBSOLETE]  vector<int> index_hex;  //To be used with the *_hex vectors
//DX20190215 [OBSOLETE]  vector<int> index_rhom;
//DX20190215 [OBSOLETE]  vector<int> index_tetr;
//DX20190215 [OBSOLETE]  vector<int> index_ortho;
//DX20190215 [OBSOLETE]  vector<int> index_mono_b;
//DX20190215 [OBSOLETE]  vector<int> index_mono_c;
//DX20190215 [OBSOLETE]  vector<int> index_tric;
//DX20190215 [OBSOLETE]
//DX20190215 [OBSOLETE]  vector<xmatrix<double> > sym_mats;
//DX20190215 [OBSOLETE]  vector<string> symbol;
//DX20190215 [OBSOLETE]  vector<string> dirparam;
//DX20190215 [OBSOLETE]  hash sym_mats_direction;
//DX20190215 [OBSOLETE]
//DX20190215 [OBSOLETE]  vector<xmatrix<double> > sym_mats_hex;
//DX20190215 [OBSOLETE]  vector<string> symbol_hex;
//DX20190215 [OBSOLETE]  vector<string> dirparam_hex;
//DX20190215 [OBSOLETE]  hash sym_mats_direction_hex;
//DX20190215 [OBSOLETE]#endif
//DX20190215 [OBSOLETE]} //namespace SYM

// **********************************************************************************************************************
// initsymmats
// **********************************************************************************************************************
//namespace SYM {
bool SymmetryInformationITC::initsymmats() {
  index_cubic.clear();
  index_hex.clear();
  index_rhom.clear();
  index_tetr.clear();
  index_ortho.clear();
  index_mono_b.clear();
  index_mono_c.clear();
  index_tric.clear();
  sym_mats.clear();
  symbol.clear();
  dirparam.clear();
  sym_mats_direction.clear();
  sym_mats_hex.clear();
  symbol_hex.clear();
  dirparam_hex.clear();
  sym_mats_direction_hex.clear();

  for (int i = 0; i < 48; i++) {
    index_cubic.push_back(i);
  }
  for (int i = 0; i < 24; i++) {
    index_hex.push_back(i);
  }

  index_rhom.push_back(0);
  index_rhom.push_back(1);
  index_rhom.push_back(2);
  index_rhom.push_back(3);
  index_rhom.push_back(4);
  index_rhom.push_back(5);
  index_rhom.push_back(10);
  index_rhom.push_back(17);
  index_rhom.push_back(24);
  index_rhom.push_back(31);
  index_rhom.push_back(38);
  index_rhom.push_back(45);

  index_tetr.push_back(0);
  index_tetr.push_back(3);
  index_tetr.push_back(6);
  index_tetr.push_back(9);
  index_tetr.push_back(10);
  index_tetr.push_back(11);
  index_tetr.push_back(12);
  index_tetr.push_back(13);
  index_tetr.push_back(16);
  index_tetr.push_back(17);
  index_tetr.push_back(18);
  index_tetr.push_back(19);
  index_tetr.push_back(20);
  index_tetr.push_back(27);
  index_tetr.push_back(34);
  index_tetr.push_back(41);

  index_ortho.push_back(0);
  index_ortho.push_back(3);
  index_ortho.push_back(6);
  index_ortho.push_back(13);
  index_ortho.push_back(20);
  index_ortho.push_back(27);
  index_ortho.push_back(34);
  index_ortho.push_back(41);

  //MONOCLINIC UNIQUE AXIS C:
  index_mono_c.push_back(0);
  index_mono_c.push_back(3);
  index_mono_c.push_back(6);
  index_mono_c.push_back(13);

  //MONOCLINIC UNIQUE AXIS B:
  index_mono_b.push_back(0);
  index_mono_b.push_back(3);
  index_mono_b.push_back(20);
  index_mono_b.push_back(27);

  index_tric.push_back(0);
  index_tric.push_back(3);

  symbol.push_back("1");
  symbol.push_back("3+");
  symbol.push_back("3-");
  symbol.push_back("-1");
  symbol.push_back("-3+");
  symbol.push_back("-3-");
  symbol.push_back("2");
  symbol.push_back("3+");
  symbol.push_back("3-");
  symbol.push_back("2");
  symbol.push_back("2");  //10
  symbol.push_back("4+");
  symbol.push_back("4-");
  symbol.push_back("m");
  symbol.push_back("-3+");
  symbol.push_back("-3-");
  symbol.push_back("m");
  symbol.push_back("m");
  symbol.push_back("-4+");
  symbol.push_back("-4-");
  symbol.push_back("2");  //20
  symbol.push_back("3+");
  symbol.push_back("3-");
  symbol.push_back("2");   //23
  symbol.push_back("2");   //24
  symbol.push_back("4+");  //25
  symbol.push_back("4-");  //26
  symbol.push_back("m");   //27
  symbol.push_back("-3+");
  symbol.push_back("-3-");
  symbol.push_back("m");    //30
  symbol.push_back("m");    //31
  symbol.push_back("-4+");  //32
  symbol.push_back("-4-");  //33
  symbol.push_back("2");    //34
  symbol.push_back("3+");
  symbol.push_back("3-");
  symbol.push_back("2");  //37
  symbol.push_back("2");  //38
  symbol.push_back("4+");
  symbol.push_back("4-");
  symbol.push_back("m");
  symbol.push_back("-3+");
  symbol.push_back("-3-");
  symbol.push_back("m");
  symbol.push_back("m");
  symbol.push_back("-4+");
  symbol.push_back("-4-");

  dirparam.push_back("     ");
  dirparam.push_back("x,x,x");
  dirparam.push_back("x,x,x");
  dirparam.push_back("0,0,0");
  dirparam.push_back("x,x,x");
  dirparam.push_back("x,x,x");
  dirparam.push_back("0,0,z");
  dirparam.push_back("x,-x,-x");
  dirparam.push_back("x,-x,-x");
  dirparam.push_back("x,x,0");
  dirparam.push_back("x,-x,0");
  dirparam.push_back("0,0,z");
  dirparam.push_back("0,0,z");
  dirparam.push_back("x,y,0");
  dirparam.push_back("x,-x,-x");
  dirparam.push_back("x,-x,-x");
  dirparam.push_back("x,-x,z");
  dirparam.push_back("x,x,z");
  dirparam.push_back("0,0,z");
  dirparam.push_back("0,0,z");
  dirparam.push_back("0,y,0");
  dirparam.push_back("-x,x,-x");
  dirparam.push_back("-x,x,-x");
  dirparam.push_back("x,0,x");
  dirparam.push_back("-x,0,x");
  dirparam.push_back("0,y,0");
  dirparam.push_back("0,y,0");
  dirparam.push_back("x,0,z");
  dirparam.push_back("-x,x,-x");
  dirparam.push_back("-x,x,-x");
  dirparam.push_back("-x,y,x");
  dirparam.push_back("x,y,x");
  dirparam.push_back("0,y,0");
  dirparam.push_back("0,y,0");
  dirparam.push_back("x,0,0");
  dirparam.push_back("-x,-x,x");
  dirparam.push_back("-x,-x,x");
  dirparam.push_back("0,y,y");
  dirparam.push_back("0,y,-y");
  dirparam.push_back("x,0,0");
  dirparam.push_back("x,0,0");
  dirparam.push_back("0,y,z");
  dirparam.push_back("-x,-x,x");
  dirparam.push_back("-x,-x,x");
  dirparam.push_back("x,y,-y");
  dirparam.push_back("x,y,y");
  dirparam.push_back("x,0,0");
  dirparam.push_back("x,0,0");

  //HEX
  symbol_hex.push_back("1");
  symbol_hex.push_back("2");
  symbol_hex.push_back("2");
  symbol_hex.push_back("2");
  symbol_hex.push_back("-1");
  symbol_hex.push_back("m");
  symbol_hex.push_back("m");
  symbol_hex.push_back("m");
  symbol_hex.push_back("3+");
  symbol_hex.push_back("6+");
  symbol_hex.push_back("2");
  symbol_hex.push_back("2");
  symbol_hex.push_back("-3+");
  symbol_hex.push_back("-6+");
  symbol_hex.push_back("m");
  symbol_hex.push_back("m");
  symbol_hex.push_back("3-");
  symbol_hex.push_back("6-");
  symbol_hex.push_back("2");
  symbol_hex.push_back("2");
  symbol_hex.push_back("-3-");
  symbol_hex.push_back("-6-");
  symbol_hex.push_back("m");
  symbol_hex.push_back("m");

  dirparam_hex.push_back("     ");
  dirparam_hex.push_back("0,0,z");
  dirparam_hex.push_back("x,x,0");
  dirparam_hex.push_back("x,-x,0");
  dirparam_hex.push_back("0,0,0");
  dirparam_hex.push_back("x,y,0");
  dirparam_hex.push_back("x,-x,z");
  dirparam_hex.push_back("x,x,z");
  dirparam_hex.push_back("0,0,z");
  dirparam_hex.push_back("0,0,z");
  dirparam_hex.push_back("x,0,0");
  dirparam_hex.push_back("x,2x,0");
  dirparam_hex.push_back("0,0,z");
  dirparam_hex.push_back("0,0,z");
  dirparam_hex.push_back("x,2x,z");
  dirparam_hex.push_back("x,0,z");
  dirparam_hex.push_back("0,0,z");
  dirparam_hex.push_back("0,0,z");
  dirparam_hex.push_back("0,y,0");
  dirparam_hex.push_back("2x,x,0");
  dirparam_hex.push_back("0,0,z");
  dirparam_hex.push_back("0,0,z");
  dirparam_hex.push_back("2x,x,z");
  dirparam_hex.push_back("0,y,z");
  xvector<double> d000;
  xvector<double> d100;
  d100(1) = 1;
  xvector<double> d010;
  d010(2) = 1;
  xvector<double> d001;
  d001(3) = 1;
  xvector<double> d111;
  d111(1) = 1;
  d111(2) = 1;
  d111(3) = 1;
  xvector<double> d1nn;
  d1nn(1) = 1;
  d1nn(2) = -1;
  d1nn(3) = -1;
  xvector<double> d110;
  d110(1) = 1;
  d110(2) = 1;
  d110(3) = 0;
  xvector<double> d1n0;
  d1n0(1) = 1;
  d1n0(2) = -1;
  d1n0(3) = 0;
  xvector<double> dn1n;
  dn1n(1) = -1;
  dn1n(2) = 1;
  dn1n(3) = -1;
  xvector<double> dn01;
  dn01(1) = -1;
  dn01(2) = 0;
  dn01(3) = 1;
  xvector<double> dnn1;
  dnn1(1) = -1;
  dnn1(2) = -1;
  dnn1(3) = 1;
  xvector<double> d101;
  d101(1) = 1;
  d101(2) = 0;
  d101(3) = 1;
  xvector<double> d011;
  d011(1) = 0;
  d011(2) = 1;
  d011(3) = 1;
  xvector<double> d01n;
  d01n(1) = 0;
  d01n(2) = 1;
  d01n(3) = -1;
  xvector<double> d120;
  d120(1) = 1;
  d120(2) = 2;
  d120(3) = 0;
  xvector<double> d210;
  d210(1) = 2;
  d210(2) = 1;
  d210(3) = 0;
  //CUBIC TETRAGONAL ORTHORHOMBIC MONOCLINIC TRICLINIC OR RHOMBOHEDRAL
  sym_mats_direction[1] = d000;
  sym_mats_direction[2] = d111;
  sym_mats_direction[3] = d111;
  sym_mats_direction[4] = d000;
  sym_mats_direction[5] = d111;
  sym_mats_direction[6] = d111;
  sym_mats_direction[7] = d001;
  sym_mats_direction[8] = d1nn;
  sym_mats_direction[9] = d1nn;
  sym_mats_direction[10] = d110;
  sym_mats_direction[11] = d1n0;
  sym_mats_direction[12] = d001;
  sym_mats_direction[13] = d001;
  sym_mats_direction[14] = d001;
  sym_mats_direction[15] = d1nn;
  sym_mats_direction[16] = d1nn;
  sym_mats_direction[17] = d110;
  sym_mats_direction[18] = d1n0;
  sym_mats_direction[19] = d001;
  sym_mats_direction[20] = d001;
  sym_mats_direction[21] = d010;
  sym_mats_direction[22] = dn1n;
  sym_mats_direction[23] = dn1n;
  sym_mats_direction[24] = d101;
  sym_mats_direction[25] = dn01;
  sym_mats_direction[26] = d010;
  sym_mats_direction[27] = d010;
  sym_mats_direction[28] = d010;
  sym_mats_direction[29] = dn1n;
  sym_mats_direction[30] = dn1n;
  sym_mats_direction[31] = d101;
  sym_mats_direction[32] = dn01;
  sym_mats_direction[33] = d010;
  sym_mats_direction[34] = d010;
  sym_mats_direction[35] = d100;
  sym_mats_direction[36] = dnn1;
  sym_mats_direction[37] = dnn1;
  sym_mats_direction[38] = d011;
  sym_mats_direction[39] = d01n;
  sym_mats_direction[40] = d100;
  sym_mats_direction[41] = d100;
  sym_mats_direction[42] = d100;
  sym_mats_direction[43] = dnn1;
  sym_mats_direction[44] = dnn1;
  sym_mats_direction[45] = d011;
  sym_mats_direction[46] = d01n;
  sym_mats_direction[47] = d100;
  sym_mats_direction[48] = d100;
  //HEXAGONAL
  sym_mats_direction_hex[1] = d000;
  sym_mats_direction_hex[2] = d001;
  sym_mats_direction_hex[3] = d110;
  sym_mats_direction_hex[4] = d1n0;
  sym_mats_direction_hex[5] = d000;
  sym_mats_direction_hex[6] = d001;
  sym_mats_direction_hex[7] = d110;
  sym_mats_direction_hex[8] = d1n0;
  sym_mats_direction_hex[9] = d001;
  sym_mats_direction_hex[10] = d001;
  sym_mats_direction_hex[11] = d100;
  sym_mats_direction_hex[12] = d120;
  sym_mats_direction_hex[13] = d001;
  sym_mats_direction_hex[14] = d001;
  sym_mats_direction_hex[15] = d100;
  sym_mats_direction_hex[16] = d120;
  sym_mats_direction_hex[17] = d001;
  sym_mats_direction_hex[18] = d001;
  sym_mats_direction_hex[19] = d010;
  sym_mats_direction_hex[20] = d210;
  sym_mats_direction_hex[21] = d001;
  sym_mats_direction_hex[22] = d001;
  sym_mats_direction_hex[23] = d010;
  sym_mats_direction_hex[24] = d210;
  double ss1[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss2[16] = {0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1};
  double ss3[16] = {0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1};
  double ss4[16] = {-1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss5[16] = {0, 0, -1, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1};
  double ss6[16] = {0, -1, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 1};
  double ss7[16] = {-1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss8[16] = {0, 0, -1, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1};
  double ss9[16] = {0, -1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 1};
  double ss10[16] = {0, 1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss11[16] = {0, -1, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss12[16] = {0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss13[16] = {0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss14[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss15[16] = {0, 0, 1, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1};
  double ss16[16] = {0, 1, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 1};
  double ss17[16] = {0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss18[16] = {0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss19[16] = {0, 1, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss20[16] = {0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss21[16] = {-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss22[16] = {0, 0, 1, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1};
  double ss23[16] = {0, -1, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 1};
  double ss24[16] = {0, 0, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
  double ss25[16] = {0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1};
  double ss26[16] = {0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1};
  double ss27[16] = {0, 0, -1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
  double ss28[16] = {1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss29[16] = {0, 0, -1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1};
  double ss30[16] = {0, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 1};
  double ss31[16] = {0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1};
  double ss32[16] = {0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
  double ss33[16] = {0, 0, -1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
  double ss34[16] = {0, 0, 1, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1};
  double ss35[16] = {1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss36[16] = {0, 0, -1, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1};
  double ss37[16] = {0, 1, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 1};
  double ss38[16] = {-1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1};
  double ss39[16] = {-1, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 1};
  double ss40[16] = {1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 1};
  double ss41[16] = {1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1};
  double ss42[16] = {-1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss43[16] = {0, 0, 1, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1};
  double ss44[16] = {0, -1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1};
  double ss45[16] = {1, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 1};
  double ss46[16] = {1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1};
  double ss47[16] = {-1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1};
  double ss48[16] = {-1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 1};
  double ss49[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss50[16] = {-1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss51[16] = {0, 1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss52[16] = {0, -1, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss53[16] = {-1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss54[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss55[16] = {0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss56[16] = {0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss57[16] = {0, -1, 0, 0, 1, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss58[16] = {1, -1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss59[16] = {1, -1, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss60[16] = {-1, 1, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss61[16] = {0, 1, 0, 0, -1, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss62[16] = {-1, 1, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss63[16] = {-1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss64[16] = {1, -1, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss65[16] = {-1, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss66[16] = {0, 1, 0, 0, -1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss67[16] = {-1, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss68[16] = {1, 0, 0, 0, 1, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss69[16] = {1, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss70[16] = {0, -1, 0, 0, 1, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1};
  double ss71[16] = {1, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  double ss72[16] = {-1, 0, 0, 0, -1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  sym_mats.push_back(SYM::a2m3x3(ss1));
  sym_mats.push_back(SYM::a2m3x3(ss2));
  sym_mats.push_back(SYM::a2m3x3(ss3));
  sym_mats.push_back(SYM::a2m3x3(ss4));
  sym_mats.push_back(SYM::a2m3x3(ss5));
  sym_mats.push_back(SYM::a2m3x3(ss6));
  sym_mats.push_back(SYM::a2m3x3(ss7));
  sym_mats.push_back(SYM::a2m3x3(ss8));
  sym_mats.push_back(SYM::a2m3x3(ss9));
  sym_mats.push_back(SYM::a2m3x3(ss10));
  sym_mats.push_back(SYM::a2m3x3(ss11));
  sym_mats.push_back(SYM::a2m3x3(ss12));
  sym_mats.push_back(SYM::a2m3x3(ss13));
  sym_mats.push_back(SYM::a2m3x3(ss14));
  sym_mats.push_back(SYM::a2m3x3(ss15));
  sym_mats.push_back(SYM::a2m3x3(ss16));
  sym_mats.push_back(SYM::a2m3x3(ss17));
  sym_mats.push_back(SYM::a2m3x3(ss18));
  sym_mats.push_back(SYM::a2m3x3(ss19));
  sym_mats.push_back(SYM::a2m3x3(ss20));
  sym_mats.push_back(SYM::a2m3x3(ss21));
  sym_mats.push_back(SYM::a2m3x3(ss22));
  sym_mats.push_back(SYM::a2m3x3(ss23));
  sym_mats.push_back(SYM::a2m3x3(ss24));
  sym_mats.push_back(SYM::a2m3x3(ss25));
  sym_mats.push_back(SYM::a2m3x3(ss26));
  sym_mats.push_back(SYM::a2m3x3(ss27));
  sym_mats.push_back(SYM::a2m3x3(ss28));
  sym_mats.push_back(SYM::a2m3x3(ss29));
  sym_mats.push_back(SYM::a2m3x3(ss30));
  sym_mats.push_back(SYM::a2m3x3(ss31));
  sym_mats.push_back(SYM::a2m3x3(ss32));
  sym_mats.push_back(SYM::a2m3x3(ss33));
  sym_mats.push_back(SYM::a2m3x3(ss34));
  sym_mats.push_back(SYM::a2m3x3(ss35));
  sym_mats.push_back(SYM::a2m3x3(ss36));
  sym_mats.push_back(SYM::a2m3x3(ss37));
  sym_mats.push_back(SYM::a2m3x3(ss38));
  sym_mats.push_back(SYM::a2m3x3(ss39));
  sym_mats.push_back(SYM::a2m3x3(ss40));
  sym_mats.push_back(SYM::a2m3x3(ss41));
  sym_mats.push_back(SYM::a2m3x3(ss42));
  sym_mats.push_back(SYM::a2m3x3(ss43));
  sym_mats.push_back(SYM::a2m3x3(ss44));
  sym_mats.push_back(SYM::a2m3x3(ss45));
  sym_mats.push_back(SYM::a2m3x3(ss46));
  sym_mats.push_back(SYM::a2m3x3(ss47));
  sym_mats.push_back(SYM::a2m3x3(ss48));
  sym_mats_hex.push_back(SYM::a2m3x3(ss49));
  sym_mats_hex.push_back(SYM::a2m3x3(ss50));
  sym_mats_hex.push_back(SYM::a2m3x3(ss51));
  sym_mats_hex.push_back(SYM::a2m3x3(ss52));
  sym_mats_hex.push_back(SYM::a2m3x3(ss53));
  sym_mats_hex.push_back(SYM::a2m3x3(ss54));
  sym_mats_hex.push_back(SYM::a2m3x3(ss55));
  sym_mats_hex.push_back(SYM::a2m3x3(ss56));
  sym_mats_hex.push_back(SYM::a2m3x3(ss57));
  sym_mats_hex.push_back(SYM::a2m3x3(ss58));
  sym_mats_hex.push_back(SYM::a2m3x3(ss59));
  sym_mats_hex.push_back(SYM::a2m3x3(ss60));
  sym_mats_hex.push_back(SYM::a2m3x3(ss61));
  sym_mats_hex.push_back(SYM::a2m3x3(ss62));
  sym_mats_hex.push_back(SYM::a2m3x3(ss63));
  sym_mats_hex.push_back(SYM::a2m3x3(ss64));
  sym_mats_hex.push_back(SYM::a2m3x3(ss65));
  sym_mats_hex.push_back(SYM::a2m3x3(ss66));
  sym_mats_hex.push_back(SYM::a2m3x3(ss67));
  sym_mats_hex.push_back(SYM::a2m3x3(ss68));
  sym_mats_hex.push_back(SYM::a2m3x3(ss69));
  sym_mats_hex.push_back(SYM::a2m3x3(ss70));
  sym_mats_hex.push_back(SYM::a2m3x3(ss71));
  sym_mats_hex.push_back(SYM::a2m3x3(ss72));
  return true;
}
//} //namespace SYM
//DX20190215 [OBSOLETE]#ifdef AFLOW_SYMMETRY_MULTITHREADS_ENABLE
//DX20190215 [OBSOLETE]thread_local vector<string> sym_ops;
//DX20190215 [OBSOLETE]#else
//DX20190215 [OBSOLETE]vector<string> sym_ops;
//DX20190215 [OBSOLETE]#endif
// **********************************************************************************************************************
// initsymops
// **********************************************************************************************************************
//DX20190215 [OBSOLETE] namespace SYM {
bool SymmetryInformationITC::initsymops() {
  //DX20190215 [OBSOLETE]extern vector<string> sgs;
  vector<string> temp;

  string so1 = "H 1(0 0 0)";
  string so2 = "H 1(0 0 0) \n I 0 0 0";
  string so3 = "H 1(0 0 0) \n 2(0 0 0) 0 y 0 \n ? \n H 1(0 0 0) \n 2(0 0 0) 0 0 z";
  string so4 = "H 1(0 0 0) \n 2(0 1/2 0) 0 y 0 \n ? \n H 1(0 0 0) \n 2(0 0 1/2) 0 0 z";
  string so5 = "H 1(0 0 0) \n 2(0 0 0) 0 y 0 \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 1/2 0) 1/4 y 0 \n ? \n H 1(0 0 0) \n 2(0 0 0) 0 0 z \n Z(0 1/2 1/2) \n t(0 1/2 1/2) \n 2(0 0 1/2) 0 1/4 z";
  string so6 = "H 1(0 0 0) \n g(0 0 0) x 0 z \n ? \n H 1(0 0 0) \n g(0 0 0) x y 0";  //mirror planes
  string so7 = "H 1(0 0 0) \n c x 0 z \n ? \n H 1(0 0 0) \n a x y 0";
  string so8 = "H 1(0 0 0) \n g(0 0 0) x 0 z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n a x 1/4 z \n ? \n H 1(0 0 0) \n g(0 0 0) x y 0 \n Z(1/2 1/2 0) \n t(0 1/2 1/2) \n b x y 1/4";
  string so9 = "H 1(0 0 0) \n c x 0 z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n n(1/2 0 1/2) x 1/4 z \n ? \n H 1(0 0 0) \n a x y 0 \n Z(0 1/2 1/2) \n t(0 1/2 1/2) \n n(1/2 1/2 0) x y 1/4";
  string so10 = "H 1(0 0 0) \n 2(0 0 0) 0 y 0 \n I 0 0 0 \n g(0 0 0) x 0 z \n ? \n H 1(0 0 0) \n 2(0 0 0) 0 0 z \n I 0 0 0 \n g(0 0 0) x y 0";
  string so11 = "H 1(0 0 0) \n 2(0 1/2 0) 0 y 0 \n I 0 0 0 \n g(0 0 0) x 1/4 z \n ? \n H 1(0 0 0) \n 2(0 0 1/2) 0 0 z \n I 0 0 0 \n g(0 0 0) x y 1/4";
  string so12 = "H 1(0 0 0) \n 2(0 0 0) 0 y 0 \n I 0 0 0 \n g(0 0 0) x 0 z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 1/2 0) 1/4 y 0 \n I 1/4 1/4 0 \n a x 1/4 z \n ? \n H 1(0 0 0) \n 2(0 0 0) 0 0 z \n I 0 0 0 \n g(0 0 0) x y 0 \n Z(0 1/2 1/2) \n t(0 1/2 1/2) \n 2(0 0 1/2) 0 1/4 z \n I 0 1/4 1/4 \n b x y 1/4";
  string so13 = "H 1(0 0 0) \n 2(0 0 0) 0 y 1/4 \n I 0 0 0 \n c x 0 z \n ? \n H 1(0 0 0) \n 2(0 0 0) 1/4 0 z \n I 0 0 0 \n a x y 0";
  string so14 = "H 1(0 0 0) \n 2(0 1/2 0) 0 y 1/4 \n I 0 0 0 \n c x 1/4 z \n ? \n H 1(0 0 0) \n 2(0 0 1/2) 1/4 0 z \n I 0 0 0 \n a x y 1/4";
  string so15 = "H 1(0 0 0) \n 2(0 0 0) 0 y 1/4 \n I 0 0 0 \n c x 0 z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 1/2 0) 1/4 y 1/4 \n I 1/4 1/4 0 \n n(1/2 0 1/2) x 1/4 z \n ? \n H 1(0 0 0) \n 2(0 0 0) 1/4 0 z \n I 0 0 0 \n a x y 0 \n Z(0 1/2 1/2) \n t(0 1/2 1/2) \n 2(0 0 1/2) 1/4 1/4 z \n I 0 1/4 1/4 \n n(1/2 1/2 0) x y 1/4";
  string so16 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 0";
  string so17 = "H 1(0 0 0) \n 2(0 0 1/2) 0 0 z \n 2(0 0 0) 0 y 1/4 \n 2(0 0 0) x 0 0";
  string so18 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 1/2 0) 1/4 y 0 \n 2(1/2 0 0) x 1/4 0";
  string so19 = "H 1(0 0 0) \n 2(0 0 1/2) 1/4 0 z \n 2(0 1/2 0) 0 y 1/4 \n 2(1/2 0 0) x 1/4 0";
  string so20 = "H 1(0 0 0) \n 2(0 0 1/2) 0 0 z \n 2(0 0 0) 0 y 1/4 \n 2(0 0 0) x 0 0 \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 1/2) 1/4 1/4 z \n 2(0 1/2 0) 1/4 y 1/4 \n 2(1/2 0 0) x 1/4 0";
  string so21 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 0 \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 0) 1/4 1/4 z \n 2(0 1/2 0) 1/4 y 0 \n 2(1/2 0 0) x 1/4 0";
  string so22 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 0 \n Z(0 1/2 1/2) \n t(0 1/2 1/2) \n 2(0 0 1/2) 0 1/4 z \n 2(0 1/2 0) 0 y 1/4 \n 2(0 0 0) x 1/4 1/4 \n Z(1/2 0 1/2) \n t(1/2 0 1/2) \n 2(0 0 1/2) 1/4 0 z \n 2(0 0 0) 1/4 y 1/4 \n 2(1/2 0 0) x 0 1/4 \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 0) 1/4 1/4 z \n 2(0 1/2 0) 1/4 y 0 \n 2(1/2 0 0) x 1/4 0";
  string so23 = "H 1(0 0 0) \n 2(0 0 z) \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 0 \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 1/2) 1/4 1/4 z \n 2(0 1/2 0) 1/4 y 1/4 \n 2(1/2 0 0) x 1/4 1/4 ";
  string so24 = "H 1(0 0 0) \n 2(0 0 1/2) 1/4 0 z \n 2(0 1/2 0) 0 y 1/4 \n 2(1/2 0 0) x 1/4 0 \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 0) 0 1/4 z \n 2(0 0 0) 1/4 y 0 \n 2(0 0 0) x 0 1/4";
  string so25 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n g(0 0 0) x 0 z \n g(0 0 0) 0 y z";
  string so26 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n c x 0 z \n g(0 0 0) 0 y z";
  string so27 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n c x 0 z \n c 0 y z";
  string so28 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n a x 0 z \n g(0 0 0) 1/4 y z";
  string so29 = "H 1(0 0 0) \n 2(0 0 1/2) 0 0 z \n a x 0 z \n c 1/4 y z";
  string so30 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n c x 1/4 z \n n(0 1/2 1/2) 0 y z";
  string so31 = "H 1(0 0 0) \n 2(0 0 1/2) 1/4 0 z \n n(1/2 0 1/2) x 0 z \n g(0 0 0) 0 y z";
  string so32 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n a x 1/4 z \n b 1/4 y z";
  string so33 = "H 1(0 0 0) \n 2(0 0 1/2) 0 0 z \n a x 1/4 z \n n(0 1/2 1/2) 1/4 y z";
  string so34 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n n(1/2 0 1/2) x 1/4 z \n n(0 1/2 1/2) 1/4 y z";
  string so35 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n g(0 0 0) x 0 z \n g(0 0 0) 0 y z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 0) 1/4 1/4 z \n a x 1/4 z \n b 1/4 y z";
  string so36 = "H 1(0 0 0) \n 2(0 0 1/2) 0 0 z \n c x 0 z \n g(0 0 0) 0 y z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 1/2) 1/4 1/4 z \n n(1/2 0 1/2) x 1/4 z \n b 1/4 y z";
  string so37 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n c x 0 z \n c 0 y z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 0) 1/4 1/4 z \n n(1/2 0 1/2) x 1/4 z \n n(0 1/2 1/2) 1/4 y z ";
  string so38 = "H 1(0 0 0) \n 2 (0 0 0) 0 0 z \n g(0 0 0) x 0 z \n g(0 0 0) 0 y z \n Z(0 1/2 1/2) \n t(0 1/2 1/2) \n 2 (0 0 1/2) 0 1/4 z \n c x 1/4 z \n n(0 1/2 1/2) 0 y z";
  string so39 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n g(0 0 0) x 1/4 z \n b 0 y z \n Z(0 1/2 1/2) \n t(0 1/2 1/2) \n 2(0 0 1/2) 0 1/4 z \n c x 0 z \n c 0 y z";
  string so40 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n a x 0 z \n g(0 0 0) 1/4 y z \n Z(0 1/2 1/2) \n t(0 1/2 1/2) \n 2(0 0 1/2) 0 1/4 z \n n(1/2 0 1/2) x 1/4 z \n n(0 1/2 1/2) 1/4 y z";
  string so41 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n a x 1/4 z \n b 1/4 y z \n Z(0 1/2 1/2) \n t(0 1/2 12) \n 2(0 0 1/2) 0 1/4 z \n n(1/2 0 1/2) x 0 z \n c 1/4 y z";
  string so42 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n g(0 0 0) x 0 z \n g(0 0 0) 0 y z \n Z(0 1/2 1/2) \n t(0 1/2 1/2) \n 2(0 0 1/2) 0 1/4 z \n c x 1/4 z \n n(0 1/2 1/2) 0 y z \n Z(1/2 0 1/2) \n t(1/2 0 1/2) \n 2(0 0 1/2) 1/4 0 z \n n(1/2 0 1/2) x 0 z \n c 1/4 y z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 0) 1/4 1/4 z \n a x 1/4 z \n b 1/4 y z";
  string so43 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n d(1/4 0 1/4) x 1/8 z \n d(0 1/4 1/4) 1/8 y z \n Z(0 1/2 1/2) \n t(0 1/2 1/2) \n 2(0 0 1/2) 0 1/4 z \n d(1/4 0 3/4) x 3/8 z \n d(0 3/4 3/4) 1/8 y z \n Z(1/2 0 1/2) \n t(1/2 0 1/2) \n 2(0 0 1/2) 1/4 0 z \n d(3/4 0 3/4) x 1/8 z \n d(0 1/4 3/4) 3/8 y z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 0) 1/4 1/4 z \n d(3/4 0 1/4) x 3/8 z \n d(0 3/4 1/4) 3/8 y z";
  string so44 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n g(0 0 0) x 0 z \n g(0 0 0) 0 y z \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 1/2) 1/4 1/4 z \n n(1/2 0 1/2) x 1/4 z \n n(0 1/2 1/2) 1/4 y z";
  string so45 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n a x 1/4 z \n b 1/4 y z \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 1/2) 1/4 1/4 z \n c x 0 z \n c 0 y z";
  string so46 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n a x 0 z \n g(0 0 0) 1/4 y z \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 1/2) 1/4 1/4 z \n c x 1/4 z \n n(0 1/2 1/2) 0 y z";
  string so47 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 0 \n I 0 0 0 \n g(0 0 0) x y 0 \n g(0 0 0) x 0 z \n g(0 0 0) 0 y z ";
  string so48 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 0 \n I 1/4 1/4 1/4 \n n(1/2 1/2 0) x y 1/4 \n n(1/2 0 1/2) x 1/4 z \n n(0 1/2 1/2) 1/4 y z \n ? \n H 1(0 0 0) \n 2(0 0 0) 1/4 1/4 z \n 2(0 0 0) 1/4 y 1/4 \n 2(0 0 0) x 1/4 1/4 \n I 0 0 0 \n n(1/2 1/2 0) x y 0 \n n(1/2 0 1/2) x 0 z \n n(0 1/2 1/2) 0 y z";
  string so49 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 0 0) 0 y 1/4 \n 2(0 0 0) x 0 1/4 \n I 0 0 0 \n g(0 0 0) x y 0 \n c x 0 z \n c 0 y z";
  string so50 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 0 \n I 1/4 1/4 0 \n n(1/2 1/2 0) x y 0 \n a x 1/4 z \n b 1/4 y z \n ? \n H 1(0 0 0) \n 2(0 0 0) 1/4 1/4 z \n 2(0 0 0) 1/4 y 0 \n 2(0 0 0) x 1/4 0 \n I 0 0 0 \n n(1/2 1/2 0) x y 0 \n a x 0 z \n b 0 y z";
  string so51 = "H 1(0 0 0) \n 2(0 0 0) 1/4 0 z \n 2(0 0 0) 0 y 0 \n 2(1/2 0 0) x 0 0 \n I 0 0 0 \n a x y 0 \n g(0 0 0) x 0 z \n g(0 0 0) 1/4 y z";
  string so52 = "H 1(0 0 0) \n 2(0 0 0) 1/4 0 z \n 2(0 1/2 0) 1/4 y 1/4 \n 2(0 0 0) x 1/4 1/4 \n I 0 0 0 \n a x y 0 \n n(1/2 0 1/2) x 1/4 z \n n(0 1/2 1/2) 0 y z ";
  string so53 = "H 1(0 0 0) \n 2(0 0 1/2) 1/4 0 z \n 2(0 0 0) 1/4 y 1/4 \n 2(0 0 0) x 0 0 \n I 0 0 0 \n a x y 1/4 \n n(1/2 0 1/2) x 0 z \n g(0 0 0) 0 y z";
  string so54 = "H 1(0 0 0) \n 2(0 0 0) 1/4 0 z \n 2(0 0 0) 0 y 1/4 \n 2(1/2 0 0) x 0 1/4 \n I 0 0 0 \n a x y 0 \n c x 0 z \n c 1/4 y z";
  string so55 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 1/2 0) 1/4 y 0 \n 2(1/2 0 0) x 1/4 0 \n I 0 0 0 \n g(0 0 0) x y 0 \n a x 1/4 z \n b 1/4 y z ";
  string so56 = "H 1(0 0 0) \n 2(0 0 0) 1/4 1/4 z \n 2(0 1/2 0) 0 y 1/4 \n 2(1/2 0 0) x 0 1/4 \n I 0 0 0 \n n(1/2 1/2 0) x y 0 \n c x 1/4 z \n c 1/4 y z";
  string so57 = "H 1(0 0 0) \n 2(0 0 1/2) 0 0 z \n 2(0 1/2 0) 0 y 1/4 \n 2(0 0 0) x 1/4 0 \n I 0 0 0 \n g(0 0 0) x y 1/4 \n c x 1/4 z \n b 0 y z";
  string so58 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 1/2 0) 1/4 y 1/4 \n 2(1/2 0 0) x 1/4 1/4 \n I 0 0 0 \n g(0 0 0) x y 0 \n n(1/2 0 1/2) x 1/4 z \n n(0 1/2 1/2) 1/4 y z";
  string so59 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 1/2 0) 1/4 y 0 \n 2(1/2 0 0) x 1/4 0 \n I 1/4 1/4 0 \n n(1/2 1/2 0) x y 0 \n g(0 0 0) x 0 z \n g(0 0 0) 0 y z \n ? \n H 1(0 0 0) \n 2(0 0 0) 1/4 1/4 z \n 2(0 1/2 0) 0 y 0 \n 2(1/2 0 0) x 0 0 \n I 0 0 0 \n n(1/2 1/2 0) x y 0 \n g(0 0 0) x 1/4 z \n g(0 0 0) 1/4 y z";
  string so60 = "H 1(0 0 0) \n 2(0 0 1/2) 1/4 1/4 z \n 2(0 0 0) 0 y 1/4 \n 2(1/2 0 0) x 1/4 0 \n I 0 0 0 \n n(1/2 1/2 0) x y 1/4 \n c x 0 z \n b 1/4 y z";
  string so61 = "H 1(0 0 0) \n 2(0 0 1/2) 1/4 0 z \n 2(0 1/2 0) 0 y 1/4 \n 2(1/2 0 0) x 1/4 0 \n I 0 0 0 \n a x y 1/4 \n c x 1/4 z \n b 1/4 y z";
  string so62 = "H 1(0 0 0) \n 2(0 0 1/2) 1/4 0 z \n 2(0 1/2 0) 0 y 0 \n 2(1/2 0 0) x 1/4 1/4 \n I 0 0 0 \n a x y 1/4 \n g(0 0 0) x 1/4 z \n n(0 1/2 1/2) 1/4 y z";
  string so63 = "H 1(0 0 0) \n 2(0 0 1/2) 0 0 z \n 2(0 0 0) 0 y 1/4 \n 2(0 0 0) x 0 0 \n I 0 0 0 \n g(0 0 0) x y 1/4 \n c x 0 z \n g(0 0 0) 0 y z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 1/2) 1/4 1/4 z \n 2(0 1/2 0) 1/4 y 1/4 \n 2(1/2 0 0) x 1/4 0 \n I 1/4 1/4 0 \n n(1/2 1/2 0) x y 1/4 \n n(1/2 0 1/2) x 1/4 z \n b 1/4 y z";
  string so64 = "H 1(0 0 0) \n 2(0 0 1/2) 0 1/4 z \n 2(0 1/2 0) 0 y 1/4 \n 2(0 0 0) x 0 0 \n I 0 0 0 \n b x y 1/4 \n c x 1/4 z \n g(0 0 0) 0 y z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 1/2) 1/4 0 z \n 2(0 0 0) 1/4 y 1/4 \n 2(1/2 0 0) x 1/4 0 \n I 1/4 1/4 0 \n a x y 1/4 \n n(1/2 0 1/2) x 0 z \n b 1/4 y z";
  string so65 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 0 \n I 0 0 0 \n g(0 0 0) x y 0 \n g(0 0 0) x 0 z \n g(0 0 0) 0 y z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 0) 1/4 1/4 z \n 2(0 1/2 0) 1/4 y 0 \n 2(1/2 0 0) x 1/4 0 \n I 1/4 1/4 0 \n n(1/2 1/2 0) x y 0 \n a x 1/4 z \n b 1/4 y z";
  string so66 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 0 0) 0 y 1/4 \n 2(0 0 0) x 0 1/4 \n I 0 0 0 \n g(0 0 0) x y 0 \n c x 0 z\n c 0 y z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 0) 1/4 1/4 z \n 2(0 1/2 0) 1/4 y 1/4 \n 2(1/2 0 0) x 1/4 1/4 \n I 1/4 1/4 0 \n n(1/2 1/2 0) x y 0 \n n(1/2 0 1/2) x 1/4 z \n n(0 1/2 1/2) 1/4 y z";
  string so67 = "H 1(0 0 0) \n 2(0 0 0) 0 1/4 z \n 2(0 1/2 0) 0 y 0 \n 2(0 0 0) x 0 0 \n I 0 0 0 \n b x y 0 \n g(0 0 0) x 1/4 z \n g(0 0 0) 0 y z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 0) 1/4 0 z \n 2(0 0 0) 1/4 y 0 \n 2(1/2 0 0) x 1/4 0 \n I 1/4 1/4 0 \n a x y 0 \n a x 0 z \n b 1/4 y z";
  string so68 = "H 1(0 0 0) \n 2(0 0 0) 1/4 1/4 z \n 2(0 0 0) 0 y 0 \n 2(1/2 0 0) x 1/4 0 \n I 0 1/4 1/4 \n a x y 1/4 \n c x 1/4 z \n c 1/4 y z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 0) 0 0 z \n 2(0 1/2 0) 1/4 y 0 \n 2(0 0 0) x 0 0 \n I 1/4 0 1/4 \n b x y 1/4 \n n(1/2 0 1/2) x 0 z \n n(0 1/2 1/2) 0 y z \n ? \n H 1(0 0 0) \n 2(0 0 0) 1/4 0 z \n 2(0 0 0) 0 y 1/4 \n 2(1/2 0 0) x 0 1/4 \n I 0 0 0 \n a x y 0 \n c x 0 z \n c 1/4 y z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 0) 0 1/4 z \n 2(0 1/2 0) 1/4 y 1/4 \n 2(0 0 0) x 1/4 1/4 \n I 1/4 1/4 0 \n b x y 0 \n n(1/2 0 1/2) x 1/4 z \n n(0 1/2 1/2) 0 y z";
  string so69 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 0 \n I 0 0 0 \n g(0 0 0) x y 0 \n g(0 0 0) x 0 z \n g(0 0 0) 0 y z \n Z(0 1/2 1/2) \n t(0 1/2 1/2) \n 2(0 0 1/2) 0 1/4 z \n 2(0 1/2 0) 0 y 1/4 \n 2(0 0 0) x 1/4 1/4 \n I 0 1/4 1/4 \n b x y 1/4 \n c x 1/4 z \n n(0 1/2 1/2) 0 y z \n Z(1/2 0 1/2) \n t(1/2 0 1/2) \n 2(0 0 1/2) 1/4 0 z \n 2(0 0 0) 1/4 y 1/4 \n 2(1/2 0 0) x 0 1/4 \n I 1/4 0 1/4 \n a x y 1/4 \n n(1/2 0 1/2) x 0 z \n c 1/4 y z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 0) 1/4 1/4 z \n 2(0 1/2 0) 1/4 y 0 \n 2(1/2 0 0) x 1/4 0 \n I 1/4 1/4 0 \n n(1/2 1/2 0) x y 0 \n a x 1/4 z \n b 1/4 y z";
  string so70 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 0 \n I 1/8 1/8 1/8 \n d(1/4 1/4 0) x y 1/8 \n d(1/4 0 1/4) x 1/8 z \n d(0 1/4 1/4) 1/8 y z \n Z(0 1/2 1/2) \n t(0 1/2 1/2) \n 2(0 0 1/2) 0 1/4 z \n 2(0 1/2 0) 0 y 1/4 \n 2(0 0 0) x 1/4 1/4 \n I 1/8 3/8 3/8 \n d(1/4 3/4 0) x y 3/8 \n d(1/4 0 3/4) x 3/8 z \n d(0 3/4 3/4) 1/8 y z \n Z(1/2 0 1/2) \n t(1/2 0 1/2) \n 2(0 0 1/2) 1/4 0 z \n 2(0 0 0) 1/4 y 1/4 \n 2(1/2 0 0) x 0 1/4 \n I 3/8 1/8 3/8 \n d(3/4 1/4 0) x y 3/8 \n d(3/4 0 3/4) x 1/8 z \n d(0 1/4 3/4) 3/8 y z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 0) 1/4 1/4 z \n 2(0 1/2 0) 1/4 y 0 \n 2(1/2 0 0) x 1/4 0 \n I 3/8 3/8 1/8 \n d(3/4 3/4 0) x y 1/8 \n d(3/4 0 1/4) x 3/8 z \n d(0 3/4 1/4) 3/8 y z \n ? \n H 1(0 0 0) \n 2(0 0 0) 3/8 3/8 z \n 2(0 0 0) 3/8 y 3/8 \n 2(0 0 0) x 3/8 3/8 \n I 0 0 0 \n d(1/4 1/4 0) x y 0 \n d(1/4 0 1/4) x 0 z \n d(0 1/4 1/4) 0 y z \n Z(0 1/2 1/2) \n t(0 1/2 1/2) \n 2(0 0 1/2) 3/8 1/8 z \n 2(0 1/2 0) 3/8 y 1/8 \n 2(0 0 0) x 1/8 1/8 \n I 0 1/4 1/4 \n d(1/4 3/4 0) x y 1/4 \n d(1/4 0 3/4) x 1/4 z \n d(0 3/4 3/4) 0 y z \n Z(1/2 0 1/2) \n t(1/2 0 1/2) \n 2(0 0 1/2) 1/8 3/8 z \n 2(0 0 0) 1/8 y 1/8 \n 2(1/2 0 0) x 3/8 1/8 \n I 1/4 0 1/4 \n d(3/4 1/4 0) x y 1/4 \n d(3/4 0 3/4) x 0 z \n d(0 1/4 3/4) 1/4 y z \n Z(1/2 1/2 0) \n t(1/2 1/2 0) \n 2(0 0 0) 1/8 1/8 z \n 2(0 1/2 0) 1/8 y 1/8 \n 2(1/2 0 0) x 1/8 3/8 \n I 1/4 1/4 0 \n d(3/4 3/4 0) x y 0 \n d(3/4 0 1/4) x 1/4 z \n d(0 3/4 1/4) 1/4 y z";
  string so71 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 0 \n I 0 0 0 \n g(0 0 0) x y 0 \n g(0 0 0) x 0 z \n g(0 0 0) 0 y z \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 1/2) 1/4 1/4 z \n 2(0 1/2 0) 1/4 y 1/4 \n 2(1/2 0 0) x 1/4 1/4 \n I 1/4 1/4 1/4 \n n(1/2 1/2 0) x y 1/4 \n n(1/2 0 1/2) x 1/4 z \n n(0 1/2 1/2) 1/4 y z";
  string so72 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 2(0 1/2 0) 1/4 y 0 \n 2(1/2 0 0) x 1/4 0 \n I 0 0 0 \n g(0 0 0) x y 0 \n a x 1/4 z \n b 1/4 y z \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 1/2) 1/4 1/4 z \n 2(0 0 0) 0 y 1/4 \n 2(0 0 0) x 0 1/4 \n I 1/4 1/4 1/4 \n n(1/2 1/2 0) x y 1/4 \n c x 0 z \n c 0 y z";
  string so73 = "H 1(0 0 0) \n 2(0 0 1/2) 1/4 0 z \n 2(0 1/2 0) 0 y 1/4 \n 2(1/2 0 0) x 1/4 0 \n I 0 0 0 \n a x y 1/4 \n c x 1/4 z \n b 1/4 y z \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 0) 0 1/4 z \n 2(0 0 0) 1/4 y 0 \n 2(0 0 0) x 0 1/4 \n I 1/4 1/4 1/4 \n b x y 0 \n a x 0 z \n c 0 y z";
  string so74 = "H 1(0 0 0) \n 2(0 0 0) 0 1/4 z \n 2(0 1/2 0) 0 y 0 \n 2(0 0 0) x 0 0 \n I 0 0 0 \n b x y 0 \n g(0 0 0) x 1/4 z \n g(0 0 0) 0 y z \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 1/2) 1/4 0 z \n 2(0 0 0) 1/4 y 1/4 \n 2(1/2 0 0) x 1/4 1/4 \n I 1/4 1/4 1/4 \n a x y 1/4 \n n(1/2 0 1/2) x 0 z \n n(0 1/2 1/2) 1/4 y z";
  //TETRAGONAL
  string so75 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 0) 0 0 z \n -4(0 0 0) 0 0 z";
  string so76 = "H 1(0 0 0) \n 2(0 0 1/2) 0 0 z \n 4(0 0 1/4) 0 0 z \n -4(0 0 3/4) 0 0 z";
  string so77 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 1/2) 0 0 z \n -4(0 0 1/2) 0 0 z";
  string so78 = "H 1(0 0 0) \n 2(0 0 1/2) 0 0 z \n 4(0 0 3/4) 0 0 z \n -4(0 0 1/4) 0 0 z";
  string so79 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 0) 0 0 z \n -4(0 0 0) 0 0 z \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 1/2) 1/4 1/4 z \n 4(0 0 1/2) 0 1/2 z \n -4(0 0 1/2) 1/2 0 z";
  string so80 = "H 1(0 0 0) \n 2(0 0 1/2) 1/4 1/4 z \n 4(0 0 1/4) -1/4 1/4 z \n -4(0 0 3/4) 1/4 -1/4 z \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 0) 0 0 z \n 4(0 0 3/4) 1/4 1/4 z \n -4(0 0 1/4) 1/4 1/4 z";
  string so81 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4 0 0 z; 0 0 0 \n -4 0 0 z; 0 0 0";
  string so82 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4 0 0 z; 0 0 0 \n -4 0 0 z; 0 0 0 \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 1/2) 1/4 1/4 z \n 4 1/2 0 z; 1/2 0 1/4 \n -4 0 1/2 z; 0 1/2 1/4";
  string so83 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 0) 0 0 z \n -4(0 0 0) 0 0 z \n I 0 0 0 \n g(0 0 0) x y 0 \n 4 0 0 z; 0 0 0 \n -4 0 0 z; 0 0 0";
  string so84 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 1/2) 0 0 z \n -4(0 0 1/2) 0 0 z \n I 0 0 0 \n g(0 0 0) x y 0 \n 4 0 0 z; 0 0 1/4 \n -4 0 0 z; 0 0 1/4";
  string so85 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 0) 0 1/2 z \n -4(0 0 0) 1/2 0 z \n I 1/4 1/4 0 \n n(1/2 1/2 0) x y 0 \n 4 0 0 z; 0 0 0 \n -4 0 0 z; 0 0 0 \n ? \n H 1(0 0 0) \n 2(0 0 0) 1/4 1/4 z \n 4(0 0 0) 1/4 1/4 z \n -4(0 0 0) 1/4 1/4 z \n I 0 0 0 \n n(1/2 1/2 0) x y 0 \n 4 1/4 -1/4 z; 1/4 -1/4 0 \n -4 -1/4 1/4 z; -1/4 1/4 0";
  string so86 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 1/2) 0 1/2 z \n -4(0 0 1/2) 1/2 0 z \n I 1/4 1/4 1/4 \n n(1/2 1/2 0) x y 1/4 \n 4 0 0 z; 0 0 0 \n -4 0 0 z; 0 0 0 \n ? \n H 1(0 0 0) \n 2(0 0 0) 1/4 1/4 z \n 4(0 0 1/2) -1/4 1/4 z \n -4(0 0 1/2) 1/4 -1/4 z \n I 0 0 0 \n n(1/2 1/2 0) x y 0 \n 4 1/4 1/4 z; 1/4 1/4 1/4 \n -4 1/4 1/4 z; 1/4 1/4 1/4";
  string so87 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 0) 0 0 z \n -4(0 0 0) 0 0 z \n I 0 0 0 \n g(0 0 0) x y 0 \n 4 0 0 z; 0 0 0 \n -4 0 0 z; 0 0 0 \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 1/2) 1/4 1/4 z \n 4(0 0 1/2) 0 1/2 z \n -4(0 0 1/2) 1/2 0 z \n I 1/4 1/4 1/4 \n n(1/2 1/2 0) x y 1/4 \n 4 1/2 0 z; 1/2 0 1/4 \n -4 0 1/2 z; 0 1/2 1/4";
  string so88 = "H 1(0 0 0) \n 2(0 0 1/2) 1/4 1/4 z \n 4(0 0 1/4) -1/4 1/4 z \n -4(0 0 3/4) 1/4 -1/4 z \n I 0 1/4 1/8 \n a x y 3/8 \n 4 0 0 z;0 0 0 \n -4 0 1/2 z; 0 1/2 1/4 \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 0) 0 0 z \n 4(0 0 3/4) 1/4 1/4 z \n -4(0 0 1/4) 1/4 1/4 z \n I 1/4 0 3/8 \n b x y 1/8 \n 4 1/2 0 z; 1/2 0 1/4 \n -4 0 0 z; 0 0 0";
  string so89 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 0) 0 0 z \n -4(0 0 0) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 0 \n 2(0 0 0) x x 0 \n 2(0 0 0) x -x 0";
  string so90 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 0) 0 1/2 z \n -4(0 0 0) 1/2 0 z \n 2(0 1/2 0) 1/4 y 0 \n 2(1/2 0 0) x 1/4 0 \n 2(0 0 0) x x 0 \n 2(0 0 0) x -x 0";
  string so91 = "H 1(0 0 0) \n 2(0 0 1/2) 0 0 z \n 4(0 0 1/4) 0 0 z \n -4(0 0 3/4) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 1/4 \n 2(0 0 0) x x 3/8 \n 2(0 0 0) x -x 1/8";
  string so92 = "H 1(0 0 0) \n 2(0 0 1/2) 0 0 z \n 4(0 0 1/4) 0 1/2 z \n -4(0 0 3/4) 1/2 0 z \n 2(0 1/2 0) 1/4 y 1/8 \n 2(1/2 0 0) x 1/4 3/8 \n 2(0 0 0) x x 0 \n 2(0 0 0) x -x 1/4";
  string so93 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 1/2) 0 0 z \n -4(0 0 1/2) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 0 \n 2(0 0 0) x x 1/4 \n 2(0 0 0) x -x 1/4";
  string so94 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 1/2) 0 1/2 z \n -4(0 0 1/2) 1/2 0 z \n 2(0 1/2 0) 1/4 y 1/4 \n 2(1/2 0 0) x 1/4 1/4 \n 2(0 0 0) x x 0 \n 2(0 0 0) x -x 0";
  string so95 = "H 1(0 0 0) \n 2(0 0 1/2) 0 0 z \n 4(0 0 3/4) 0 0 z \n -4(0 0 1/4) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 1/4 \n 2(0 0 0) x x 1/8 \n 2(0 0 0) x -x 3/8";
  string so96 = "H 1(0 0 0) \n 2(0 0 1/2) 0 0 z \n 4(0 0 3/4) 0 1/2 z \n -4(0 0 1/4) 1/2 0 z \n 2(0 1/2 0) 1/4 y 3/8 \n 2(1/2 0 0) x 1/4 1/8 \n 2(0 0 0) x x 0 \n 2(0 0 0) x -x 1/4";
  string so97 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 0) 0 0 z \n -4(0 0 0) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 0 \n 2(0 0 0) x x 0 \n 2(0 0 0) x -x 0 \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 1/2) 1/4 1/4 z \n 4(0 0 1/2) 0 1/2 z \n -4(0 0 1/2) 1/2 0 z \n 2(0 1/2 0) 1/4 y 1/4 \n 2(1/2 0 0) x 1/4 1/4 \n 2(1/2 1/2 0) x x 1/4 \n 2(0 0 0) x -x+1/2 1/4";
  string so98 = "H 1(0 0 0) \n 2(0 0 1/2) 1/4 1/4 z \n 4(0 0 1/4) -1/4 1/4 z \n -4(0 0 3/4) 1/4 -1/4 z \n 2(0 0 0) 1/4 y 3/8 \n 2(0 0 0) x 1/4 1/8 \n 2(1/2 1/2 0) x x 1/4 \n 2(0 0 0) x -x 0 \n Z(1/2 1/2 1/2) \n t(1/2 1/2 1/2) \n 2(0 0 0) 0 0 z \n 4(0 0 3/4) 1/4 1/4 z \n -4(0 0 1/4) 1/4 1/4 z \n 2(0 1/2 0) 0 y 1/8 \n 2(1/2 0 0) x 0 3/8 \n 2(0 0 0) x x 0 \n 2(0 0 0) x -x+1/2 1/4";
  string so99 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 0) 0 0 z \n -4(0 0 0) 0 0 z \n g(0 0 0) x 0 z \n g(0 0 0) 0 y z \n g(0 0 0) x -x z \n g(0 0 0) x x z";
  string so100 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 0) 0 0 z \n -4(0 0 0) 0 0 z \n a x 1/4 z \n b 1/4 y z \n g(0 0 0) x+1/2 -x z \n g(1/2 1/2 0) x x z";

  string so123 = "H 1(0 0 0) \n 2(0 0 0) 0 0 z \n 4(0 0 0) 0 0 z \n -4(0 0 0) 0 0 z \n 2(0 0 0) 0 y 0 \n 2(0 0 0) x 0 0 \n 2(0 0 0) x x 0 \n 2(0 0 0) x -x 0 \n I 0 0 0 \n g(0 0 0) x y 0 \n 4 0 0 z; 0 0 0 \n -4 0 0 z; 0 0 0 \n g(0 0 0) x 0 z \n g(0 0 0) 0 y z \n g(0 0 0) x -x z \n g(0 0 0) x x z";
  string so164 = "H 1(0 0 0) \n H 3(0 0 0) 0 0 z \n H -3(0 0 0) 0 0 z \n H 2(0 0 0) x x 0 \n H 2(0 0 0) x 0 0 \n H 2(0 0 0) 0 y 0 \n I 0 0 0 \n H 3 0 0 z; 0 0 0 \n H -3 0 0 z; 0 0 0 \n H g(0 0 0) x -x z \n H g(0 0 0) x 2x z \n H g(0 0 0) 2x x z ";

  //pushback
  sym_ops.push_back(so1);
  sym_ops.push_back(so2);
  sym_ops.push_back(so3);
  sym_ops.push_back(so4);
  sym_ops.push_back(so5);
  sym_ops.push_back(so6);
  sym_ops.push_back(so7);
  sym_ops.push_back(so8);
  sym_ops.push_back(so9);
  sym_ops.push_back(so10);
  sym_ops.push_back(so11);
  sym_ops.push_back(so12);
  sym_ops.push_back(so13);
  sym_ops.push_back(so14);
  sym_ops.push_back(so15);
  sym_ops.push_back(so16);
  sym_ops.push_back(so17);
  sym_ops.push_back(so18);
  sym_ops.push_back(so19);
  sym_ops.push_back(so20);
  sym_ops.push_back(so21);
  sym_ops.push_back(so22);
  sym_ops.push_back(so23);
  sym_ops.push_back(so24);
  sym_ops.push_back(so25);
  sym_ops.push_back(so26);
  sym_ops.push_back(so27);
  sym_ops.push_back(so28);
  sym_ops.push_back(so29);
  sym_ops.push_back(so30);
  sym_ops.push_back(so31);
  sym_ops.push_back(so32);
  sym_ops.push_back(so33);
  sym_ops.push_back(so34);
  sym_ops.push_back(so35);
  sym_ops.push_back(so36);
  sym_ops.push_back(so37);
  sym_ops.push_back(so38);
  sym_ops.push_back(so39);
  sym_ops.push_back(so40);
  sym_ops.push_back(so41);
  sym_ops.push_back(so42);
  sym_ops.push_back(so43);
  sym_ops.push_back(so44);
  sym_ops.push_back(so45);
  sym_ops.push_back(so46);
  sym_ops.push_back(so47);
  sym_ops.push_back(so48);
  sym_ops.push_back(so49);
  sym_ops.push_back(so50);
  sym_ops.push_back(so51);
  sym_ops.push_back(so52);
  sym_ops.push_back(so53);
  sym_ops.push_back(so54);
  sym_ops.push_back(so55);
  sym_ops.push_back(so56);
  sym_ops.push_back(so57);
  sym_ops.push_back(so58);
  sym_ops.push_back(so59);
  sym_ops.push_back(so60);
  sym_ops.push_back(so61);
  sym_ops.push_back(so62);
  sym_ops.push_back(so63);
  sym_ops.push_back(so64);
  sym_ops.push_back(so65);
  sym_ops.push_back(so66);
  sym_ops.push_back(so67);
  sym_ops.push_back(so68);
  sym_ops.push_back(so69);
  sym_ops.push_back(so70);
  sym_ops.push_back(so71);
  sym_ops.push_back(so72);
  sym_ops.push_back(so73);
  sym_ops.push_back(so74);
  sym_ops.push_back(so75);
  sym_ops.push_back(so76);
  sym_ops.push_back(so77);
  sym_ops.push_back(so78);
  sym_ops.push_back(so79);
  sym_ops.push_back(so80);
  sym_ops.push_back(so81);
  sym_ops.push_back(so82);
  sym_ops.push_back(so83);
  sym_ops.push_back(so84);
  sym_ops.push_back(so85);
  sym_ops.push_back(so86);
  sym_ops.push_back(so87);
  sym_ops.push_back(so88);
  sym_ops.push_back(so89);
  sym_ops.push_back(so90);
  sym_ops.push_back(so91);
  sym_ops.push_back(so92);
  sym_ops.push_back(so93);
  sym_ops.push_back(so94);
  sym_ops.push_back(so95);
  sym_ops.push_back(so96);
  sym_ops.push_back(so97);
  sym_ops.push_back(so98);
  sym_ops.push_back(so99);
  sym_ops.push_back(so100);
  sym_ops.push_back(so123);
  sym_ops.push_back(so164);

  //etc

  return true;
}
//} //namespace SYM

// **********************************************************************************************************************
// point_group_library_search
// **********************************************************************************************************************
namespace SYM {
  string point_group_library_search(vector<string> symmetryoperations, string crystalsystem, int centeringops) {
    //symmetry ops should not be multiplied by the centering ops

    if(1 == 0) {  //holder to eliminate compiler warnings
      cerr << centeringops;
    }

    int order = symmetryoperations.size();
    string out = "";
    ostringstream pgroupss;
    int bins[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    for (uint i = 0; i < symmetryoperations.size(); i++) {
      if(symmetryoperations[i] == "-1") { bins[0]++; }
      if(symmetryoperations[i] == "m") { bins[1]++; }
      if(symmetryoperations[i] == "2") { bins[2]++; }
      if(symmetryoperations[i] == "3+" || symmetryoperations[i] == "3-") { bins[3]++; }
      if(symmetryoperations[i] == "4+" || symmetryoperations[i] == "4-") { bins[4]++; }
      if(symmetryoperations[i] == "6+" || symmetryoperations[i] == "6-") { bins[5]++; }
      if(symmetryoperations[i] == "-2") { bins[6]++; }
      if(symmetryoperations[i] == "-3+" || symmetryoperations[i] == "-3-") { bins[7]++; }
      if(symmetryoperations[i] == "-4+" || symmetryoperations[i] == "-4-") { bins[8]++; }
      if(symmetryoperations[i] == "-6+" || symmetryoperations[i] == "-6-") { bins[9]++; }
    }
    for (int i = 0; i < 10; i++) {
      pgroupss << bins[i] << " ";
    }
    string symlist = pgroupss.str();
    //cerr << "symlist: \"" << symlist <<"\"" << endl;
    map<string, map<int, map<string, string> > > pointgroupmap;
    //typedef  map<string,map<int,map<string,string > > >::iterator pit;
    //-1 m 2 3+ 4+ 6+ -2 3+ -4+ -6+
    pointgroupmap["TRICLINIC"][1]["0 0 0 0 0 0 0 0 0 0 "] = "1";
    pointgroupmap["TRICLINIC"][2]["1 0 0 0 0 0 0 0 0 0 "] = "-1";
    pointgroupmap["MONOCLINIC"][2]["0 0 1 0 0 0 0 0 0 0 "] = "2";
    pointgroupmap["MONOCLINIC"][2]["0 1 0 0 0 0 0 0 0 0 "] = "m";
    pointgroupmap["MONOCLINIC"][4]["1 1 1 0 0 0 0 0 0 0 "] = "2/m";
    pointgroupmap["ORTHORHOMBIC"][4]["0 0 3 0 0 0 0 0 0 0 "] = "222";
    pointgroupmap["ORTHORHOMBIC"][4]["0 2 1 0 0 0 0 0 0 0 "] = "mm2";
    pointgroupmap["ORTHORHOMBIC"][8]["1 3 3 0 0 0 0 0 0 0 "] = "mmm";
    pointgroupmap["TETRAGONAL"][4]["0 0 1 0 2 0 0 0 0 0 "] = "4";
    pointgroupmap["TETRAGONAL"][4]["0 0 1 0 0 0 0 0 2 0 "] = "-4";
    pointgroupmap["TETRAGONAL"][8]["1 1 1 0 2 0 0 0 2 0 "] = "4/m";
    pointgroupmap["TETRAGONAL"][8]["0 0 5 0 2 0 0 0 0 0 "] = "422";
    pointgroupmap["TETRAGONAL"][8]["0 4 1 0 2 0 0 0 0 0 "] = "4mm";
    pointgroupmap["TETRAGONAL"][8]["0 2 3 0 0 0 0 0 2 0 "] = "-42m";
    pointgroupmap["TETRAGONAL"][16]["1 5 5 0 2 0 0 0 2 0 "] = "4/mmm";
    pointgroupmap["TRIGONAL"][3]["0 0 0 2 0 0 0 0 0 0 "] = "3";
    pointgroupmap["TRIGONAL"][6]["1 0 0 2 0 0 0 2 0 0 "] = "-3";
    pointgroupmap["TRIGONAL"][6]["0 0 3 2 0 0 0 0 0 0 "] = "32";
    pointgroupmap["TRIGONAL"][6]["0 3 0 2 0 0 0 0 0 0 "] = "3m";
    pointgroupmap["TRIGONAL"][12]["1 3 3 2 0 0 0 2 0 0 "] = "-3m";
    pointgroupmap["HEXAGONAL"][6]["0 1 0 2 0 2 0 0 0 0 "] = "6";   //POSSIBILITY 1: Since mirror=2fold operation in 2-D
    pointgroupmap["HEXAGONAL"][6]["0 0 1 2 0 2 0 0 0 0 "] = "6";   //POSSIBILITY 2: Since mirror=2fold operation in 2-D
    pointgroupmap["HEXAGONAL"][6]["0 1 0 2 0 0 0 0 0 2 "] = "-6";  //POSSIBILITY 1: Since mirror=2fold operation in 2-D
    pointgroupmap["HEXAGONAL"][6]["0 0 1 2 0 0 0 0 0 2 "] = "-6";  //POSSIBILITY 2: Since mirror=2fold operation in 2-D
    pointgroupmap["HEXAGONAL"][12]["1 1 1 2 0 2 0 2 0 2 "] = "6/m";
    pointgroupmap["HEXAGONAL"][12]["0 0 7 2 0 2 0 0 0 0 "] = "622";
    pointgroupmap["HEXAGONAL"][12]["0 6 1 2 0 2 0 0 0 0 "] = "6mm";
    pointgroupmap["HEXAGONAL"][12]["0 4 3 2 0 0 0 0 0 2 "] = "-6m2";
    pointgroupmap["HEXAGONAL"][24]["1 7 7 2 0 2 0 2 0 2 "] = "6/mmm";
    pointgroupmap["CUBIC"][12]["0 0 3 8 0 0 0 0 0 0 "] = "23";
    pointgroupmap["CUBIC"][24]["1 3 3 8 0 0 0 8 0 0 "] = "m-3";
    pointgroupmap["CUBIC"][24]["0 0 9 8 6 0 0 0 0 0 "] = "432";
    pointgroupmap["CUBIC"][24]["0 6 3 8 0 0 0 0 6 0 "] = "-43m";
    pointgroupmap["CUBIC"][48]["1 9 9 8 6 0 0 8 6 0 "] = "m-3m";
    //-1 m 2 3+ 4+ 6+ -2 3+ -4+ -6+

    //for(pit iterator = pointgroupmap.begin(); iterator != pointgroupmap.end();iterator++){
    //  cerr << iterator->first  << endl;
    //}
    //cerr << "Crystal System: " << crystalsystem << endl;
    //cerr << crystalsystem << endl;
    //cerr << order << endl;
    //cerr << symlist << endl;
    //cerr << "Point Group: " << pointgroupmap[crystalsystem][order][symlist] << endl;
    out = pointgroupmap[crystalsystem][order][symlist];
    return out;
  }
} //namespace SYM

// **********************************************************************************************************************
// crystal_system_library_search
// **********************************************************************************************************************
namespace SYM {
  string crystal_system_library_search(vector<string> symmetryoperations) {
    bool LDEBUG = (FALSE || XHOST.DEBUG);

    //symmetry ops should not be multiplied by the centering ops
    //int order = symmetryoperations.size();
    string out = "";
    ostringstream pgroupss;
    int bins[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    for (uint i = 0; i < symmetryoperations.size(); i++) {
      if(symmetryoperations[i] == "-1") { bins[0]++; }
      if(symmetryoperations[i] == "m") { bins[1]++; }
      if(symmetryoperations[i] == "2") { bins[2]++; }
      if(symmetryoperations[i] == "3+" || symmetryoperations[i] == "3-") { bins[3]++; }
      if(symmetryoperations[i] == "4+" || symmetryoperations[i] == "4-") { bins[4]++; }
      if(symmetryoperations[i] == "6+" || symmetryoperations[i] == "6-") { bins[5]++; }
      if(symmetryoperations[i] == "-2") { bins[6]++; }
      if(symmetryoperations[i] == "-3+" || symmetryoperations[i] == "-3-") { bins[7]++; }
      if(symmetryoperations[i] == "-4+" || symmetryoperations[i] == "-4-") { bins[8]++; }
      if(symmetryoperations[i] == "-6+" || symmetryoperations[i] == "-6-") { bins[9]++; }
    }
    for (int i = 0; i < 10; i++) {
      pgroupss << bins[i] << " ";
    }
    string symlist = pgroupss.str();
    //cerr << "-1 m 2 3+ 4+ 6+ -2 3+ -4+ -6+" << endl;
    if(LDEBUG) { cerr << "SYM::crystal_system_library_search: Symmetry list: \"" << symlist << "\"" << endl; }
    map<string, string> pointgroupmap;
    //typedef  map<string,map<int,map<string,string > > >::iterator pit;
    //-1 m 2 3+ 4+ 6+ -2 3+ -4+ -6+
    pointgroupmap["0 0 0 0 0 0 0 0 0 0 "] = "1 TRICLINIC";
    pointgroupmap["1 0 0 0 0 0 0 0 0 0 "] = "-1 TRICLINIC";
    pointgroupmap["0 0 1 0 0 0 0 0 0 0 "] = "2 MONOCLINIC";
    pointgroupmap["0 1 0 0 0 0 0 0 0 0 "] = "m MONOCLINIC";
    pointgroupmap["1 1 1 0 0 0 0 0 0 0 "] = "2/m MONOCLINIC";
    pointgroupmap["0 0 3 0 0 0 0 0 0 0 "] = "222 ORTHORHOMBIC";
    pointgroupmap["0 2 1 0 0 0 0 0 0 0 "] = "mm2 ORTHORHOMBIC";
    pointgroupmap["1 3 3 0 0 0 0 0 0 0 "] = "mmm ORTHORHOMBIC";
    pointgroupmap["0 0 1 0 2 0 0 0 0 0 "] = "4 TETRAGONAL";
    pointgroupmap["0 0 1 0 0 0 0 0 2 0 "] = "-4  TETRAGONAL";
    pointgroupmap["1 1 1 0 2 0 0 0 2 0 "] = "4/m  TETRAGONAL";
    pointgroupmap["0 0 5 0 2 0 0 0 0 0 "] = "422  TETRAGONAL";
    pointgroupmap["0 4 1 0 2 0 0 0 0 0 "] = "4mm  TETRAGONAL";
    pointgroupmap["0 2 3 0 0 0 0 0 2 0 "] = "-42m  TETRAGONAL";
    pointgroupmap["1 5 5 0 2 0 0 0 2 0 "] = "4/mmm  TETRAGONAL";
    pointgroupmap["0 0 0 2 0 0 0 0 0 0 "] = "3 TRIGONAL";
    pointgroupmap["1 0 0 2 0 0 0 2 0 0 "] = "-3 TRIGONAL";
    pointgroupmap["0 0 3 2 0 0 0 0 0 0 "] = "32 TRIGONAL";
    pointgroupmap["0 3 0 2 0 0 0 0 0 0 "] = "3m TRIGONAL";
    pointgroupmap["1 3 3 2 0 0 0 2 0 0 "] = "-3m TRIGONAL";
    pointgroupmap["0 1 0 2 0 2 0 0 0 0 "] = "6 HEXAGONAL";   //POSSIBILITY 1: Since mirror=2fold operation in 2-D
    pointgroupmap["0 0 1 2 0 2 0 0 0 0 "] = "6 HEXAGONAL";   //POSSIBILITY 2: Since mirror=2fold operation in 2-D
    pointgroupmap["0 1 0 2 0 0 0 0 0 2 "] = "-6 HEXAGONAL";  //POSSIBILITY 1: Since mirror=2fold operation in 2-D
    pointgroupmap["0 0 1 2 0 0 0 0 0 2 "] = "-6 HEXAGONAL";  //POSSIBILITY 2: Since mirror=2fold operation in 2-D
    pointgroupmap["1 1 1 2 0 2 0 2 0 2 "] = "6/m HEXAGONAL";
    pointgroupmap["0 0 7 2 0 2 0 0 0 0 "] = "622 HEXAGONAL";
    pointgroupmap["0 6 1 2 0 2 0 0 0 0 "] = "6mm HEXAGONAL";
    pointgroupmap["0 4 3 2 0 0 0 0 0 2 "] = "-6m2 HEXAGONAL";
    pointgroupmap["1 7 7 2 0 2 0 2 0 2 "] = "6/mmm HEXAGONAL";
    pointgroupmap["0 0 3 8 0 0 0 0 0 0 "] = "23 CUBIC";
    pointgroupmap["1 3 3 8 0 0 0 8 0 0 "] = "m-3  CUBIC";
    pointgroupmap["0 0 9 8 6 0 0 0 0 0 "] = "432  CUBIC";
    pointgroupmap["0 6 3 8 0 0 0 0 6 0 "] = "-43m  CUBIC";
    pointgroupmap["1 9 9 8 6 0 0 8 6 0 "] = "m-3m  CUBIC";
    //-1 m 2 3+ 4+ 6+ -2 3+ -4+ -6+

    //for(pit iterator = pointgroupmap.begin(); iterator != pointgroupmap.end();iterator++){
    //  cerr << iterator->first  << endl;
    //}
    //cerr << "Crystal System: " << crystalsystem << endl;
    //cerr << crystalsystem << endl;
    //cerr << order << endl;
    //cerr << symlist << endl;
    //cerr << "Point Group: " << pointgroupmap[crystalsystem][order][symlist] << endl;

    map<string, string>::iterator it = pointgroupmap.find(symlist);
    //cerr << "it: " << it->first << endl;
    if(it != pointgroupmap.end()) {
      out = pointgroupmap[symlist];
      //cerr <<"OUT: " << out << endl;
      return out;
    } else {
      if(LDEBUG) { cerr << "SYM::crystal_system_library_search: WARNING: Point groups found do not match any possible point group sets." << endl; }
      out = "REDO";
      return out;
    }
  }
} //namespace SYM

// **********************************************************************************************************************
// get_random_double
// **********************************************************************************************************************
namespace SYM {
  double get_random_double(double min, double max) {
    //char sz[64];
    double lf = (max - min) * ((double)rand() / (double)RAND_MAX) + min;
    //sprintf(sz,"%..2lf\n",lf); //TRUNCATE TO 4 DECIMAL PLACES.
    //double lf2 = atof(sz);
    //return lf2;
    return lf;
  }
} //namespace SYM

// **********************************************************************************************************************
// spacegroup_to_lattice
// **********************************************************************************************************************
namespace SYM {
  xmatrix<double> spacegroup_to_lattice(int spacegroup) {
    xmatrix<double> L;
    char lattice = '\0';
    //char centering;

    if(spacegroup <= 2) {
      lattice = 'a';
    } else if(spacegroup <= 15 && spacegroup > 2) {
      lattice = 'm';
    } else if(spacegroup <= 74 && spacegroup > 15) {
      lattice = 'o';
    } else if(spacegroup <= 142 && spacegroup > 74) {
      lattice = 't';
    } else if(spacegroup <= 194 && spacegroup > 142) {
      lattice = 'h';
    } else if(spacegroup <= 230 && spacegroup > 194) {
      lattice = 'c';
    }

    if(lattice == 'c') {
      L(1, 1) = 1.0;
      L(1, 2) = 0;
      L(1, 3) = 0;
      L(2, 1) = 0;
      L(2, 2) = 1.0;
      L(2, 3) = 0;
      L(3, 1) = 0;
      L(3, 2) = 0;
      L(3, 3) = 1.0;
    }
    if(lattice == 'h') {
      L(1, 1) = 1.0;
      L(1, 2) = 0;
      L(1, 3) = 0;
      L(2, 1) = -1.0 / 2.0;
      L(2, 2) = sqrt(3) / 2;
      L(2, 3) = 0;
      L(3, 1) = 0;
      L(3, 2) = 0;
      L(3, 3) = 1.76;
    }
    if(lattice == 't') {
      L(1, 1) = 2.2;
      L(1, 2) = 0;
      L(1, 3) = 0;
      L(2, 1) = 0;
      L(2, 2) = 2.2;
      L(2, 3) = 0;
      L(3, 1) = 0;
      L(3, 2) = 0;
      L(3, 3) = 1.66;
    }
    if(lattice == 'o') {
      L(1, 1) = 1.2;
      L(1, 2) = 0;
      L(1, 3) = 0;
      L(2, 1) = 0;
      L(2, 2) = 1.6;
      L(2, 3) = 0;
      L(3, 1) = 0;
      L(3, 2) = 0;
      L(3, 3) = 2.8;
    }
    if(lattice == 'm') {  //monoclinic unique axis b (to be consistent with wyckoff positions in initsgs)
      L(1, 1) = 1.56;
      L(1, 2) = 0;
      L(1, 3) = 0;
      L(2, 1) = 0;
      L(2, 2) = 1.0;
      L(2, 3) = 0.0;
      L(3, 1) = -0.344;
      L(3, 2) = 0;
      L(3, 3) = .97;
    }
    if(lattice == 'a') {
      L(1, 1) = 1.56;
      L(1, 2) = 0;
      L(1, 3) = 0;
      L(2, 1) = -0.344;
      L(2, 2) = .97;
      L(2, 3) = 0;
      L(3, 1) = .2;
      L(3, 2) = .2;
      L(3, 3) = 1.0;
    }
    //int Csgs[] = {5,8,9,12,13,14,15,20,21,35,36,37,38,39,40,41,63,64,65,66,67,68,146,148,155,160,161,166,167};
    //vector<int> C (Csgs, Csgs + sizeof(Csgs) / sizeof(int) ); //iterator constructor
    //
    //int Fsgs[] = {22,42,43,69,70,119,120,121,122,196,202,203,209,210,216,219,225,228};
    //vector<int> F (Fsgs, Fsgs + sizeof(Fsgs) / sizeof(int) ); //iterator constructor
    return L;
  }
} //namespace SYM

// **********************************************************************************************************************
// spacegroups_CS_LS
// **********************************************************************************************************************
namespace SYM {
  vector<int> spacegroups_CS_LS(string crystalsystem, char latticesystem) {
    vector<int> spacegroups;
    if(latticesystem == 'a') {
      spacegroups.push_back(1);
      spacegroups.push_back(2);
    }
    if(latticesystem == 'm') {
      spacegroups.push_back(3);
      spacegroups.push_back(15);
    }
    if(latticesystem == 'o') {
      spacegroups.push_back(16);
      spacegroups.push_back(74);
    }
    if(latticesystem == 't') {
      spacegroups.push_back(75);
      spacegroups.push_back(142);
    }
    if(latticesystem == 'r') {
      spacegroups.push_back(143);
      spacegroups.push_back(167);
    }
    if(latticesystem == 'h' && crystalsystem == "Trigonal") {
      spacegroups.push_back(143);
      spacegroups.push_back(167);
    }
    if(latticesystem == 'h' && crystalsystem == "Hexagonal") {
      spacegroups.push_back(168);
      spacegroups.push_back(194);
    }
    if(latticesystem == 'c') {
      spacegroups.push_back(195);
      spacegroups.push_back(230);
    }
    return spacegroups;
  }
} //namespace SYM

namespace SYM {
  ostream& operator<<(ostream& output, const symop& a) {
    double tol = 1e-8;
    ostringstream diross;
    string direction_string;
    vector<string> shift_string;
    vector<string> screwglide_string;
    xvector<double> zero;
    for (uint i = 1; i < 4; i++) {
      if(aurostd::abs(a.shift(i)) < tol) {
        shift_string.push_back("");
      } else {
        shift_string.push_back(aurostd::dbl2frac(a.shift(i))); //DX20190724 - now namespace aurostd
      }
    }
    for (uint i = 1; i < 4; i++) {
      if(aurostd::abs(a.screwglide(i)) < tol) {
        screwglide_string.push_back("0");
      } else {
        screwglide_string.push_back(aurostd::dbl2frac(a.screwglide(i))); //DX20190724 - now namespace aurostd
      }
    }
    if(!havechar(a.symbol, 'm') && !havechar(a.symbol, 'n') && !havechar(a.symbol, 'd') && !havechar(a.symbol, 'a') && !havechar(a.symbol, 'c') && !havechar(a.symbol, 'b')) {
      if(a.direction == dir(1, 0, 0)) {
        diross << "x" << shift_string[0] << ",";
        diross << "0" << shift_string[1] << ",";
        diross << "0" << shift_string[2];
      }
      if(a.direction == dir(0, 1, 0)) {
        diross << "0" << shift_string[0] << ",";
        diross << "y" << shift_string[1] << ",";
        diross << "0" << shift_string[2];
      }
      if(a.direction == dir(0, 0, 1)) {
        diross << "0" << shift_string[0] << ",";
        diross << "0" << shift_string[1] << ",";
        diross << "z" << shift_string[2];
      }
      if(a.direction == dir(1, 1, 0)) {
        diross << "x" << shift_string[0] << ",";
        diross << "x" << shift_string[1] << ",";
        diross << "0" << shift_string[2];
      }
      if(a.direction == dir(1, 1, 1)) {
        diross << "x" << shift_string[0] << ",";
        diross << "x" << shift_string[1] << ",";
        diross << "x" << shift_string[2];
      }
      if(a.direction == dir(1, -1, 0)) {
        diross << "x" << shift_string[0] << ",";
        diross << "-x" << shift_string[1] << ",";
        diross << "0" << shift_string[2];
      }
      if(a.direction == dir(-1, 0, 1)) {
        diross << "-x" << shift_string[0] << ",";
        diross << "0" << shift_string[1] << ",";
        diross << "x" << shift_string[2];
      }

    } else {
      if(a.direction == dir(1, 0, 0)) {
        diross << "0" << shift_string[0] << ",";
        diross << "y" << shift_string[1] << ",";
        diross << "z" << shift_string[2];
      }
      if(a.direction == dir(0, 1, 0)) {
        diross << "x" << shift_string[0] << ",";
        diross << "0" << shift_string[1] << ",";
        diross << "z" << shift_string[2];
      }
      if(a.direction == dir(0, 0, 1)) {
        diross << "x" << shift_string[0] << ",";
        diross << "y" << shift_string[1] << ",";
        diross << "0" << shift_string[2];
      }
      if(a.direction == dir(1, 1, 0)) {
        diross << "x" << shift_string[0] << ",";
        diross << "-x" << shift_string[1] << ",";
        diross << "z" << shift_string[2];
      }
      if(a.direction == dir(1, -1, 0)) {
        diross << "x" << shift_string[0] << ",";
        diross << "x" << shift_string[1] << ",";
        diross << "z" << shift_string[2];
      }
      if(a.direction == dir(-1, 0, 1)) {
        diross << "x" << shift_string[0] << ",";
        diross << "y" << shift_string[1] << ",";
        diross << "x" << shift_string[2];
      }
    }

    if(!havechar(a.symbol, '1') && a.symbol[0] != '-') {
      if(a.screwglide != zero) {
        output << a.symbol << "(" << screwglide_string[0] << "," << screwglide_string[1] << "," << screwglide_string[2] << ") " << diross.str();
      } else {
        output << a.symbol << " " << diross.str();
      }
    } else {
      if(a.symbol == "-1") {
        output << a.symbol << " " << screwglide_string[0] << "," << screwglide_string[1] << "," << screwglide_string[2];
      } else if(a.symbol == "1") {
        output << a.symbol << " " << aurostd::dbl2frac(a.shift(1)) << "," << aurostd::dbl2frac(a.shift(2)) << "," << aurostd::dbl2frac(a.shift(3)); //DX20190724 - now namespace aurostd
      } else {
        output << a.symbol << " " << diross.str() << ";" << screwglide_string[0] << "," << screwglide_string[1] << "," << screwglide_string[2];
      }
    }

    return output;
  }
} //namespace SYM

// **********************************************************************************************************************
// generatorindex
// **********************************************************************************************************************
//returns vector of integers that indicate where the generator operations are in complete list of symmetry operations
//given the space group number (sg)
namespace SYM {
  vector<int> generatorindex(int sg) {
    vector<int> indx;
    if(sg == 1) {
      indx.push_back(1);
    }
    if(sg >= 2 && sg <= 9) {
      indx.push_back(2);
    }
    if(sg >= 10 && sg <= 46) {
      indx.push_back(2);
      indx.push_back(3);
    }
    if(sg >= 47 && sg <= 74) {
      indx.push_back(2);
      indx.push_back(3);
      indx.push_back(5);
    }
    if(sg >= 75 && sg <= 82) {
      indx.push_back(2);
      indx.push_back(3);
    }
    if(sg >= 83 && sg <= 122) {
      indx.push_back(2);
      indx.push_back(3);
      indx.push_back(5);
    }
    if(sg >= 123 && sg <= 142) {
      indx.push_back(2);
      indx.push_back(3);
      indx.push_back(5);
      indx.push_back(9);
    }
    if(sg >= 143 && sg <= 146) {
      indx.push_back(2);
    }
    if(sg >= 147 && sg <= 161) {
      indx.push_back(2);
      indx.push_back(4);
    }
    if(sg >= 162 && sg <= 167) {
      indx.push_back(2);
      indx.push_back(4);
      indx.push_back(7);
    }
    if(sg >= 168 && sg <= 174) {
      indx.push_back(2);
      indx.push_back(4);
    }
    if(sg >= 175 && sg <= 190) {
      indx.push_back(2);
      indx.push_back(4);
      indx.push_back(7);
    }
    if(sg >= 191 && sg <= 194) {
      indx.push_back(2);
      indx.push_back(4);
      indx.push_back(7);
      indx.push_back(13);
    }
    if(sg >= 195 && sg <= 199) {
      indx.push_back(2);
      indx.push_back(3);
      indx.push_back(5);
    }
    if(sg >= 200 && sg <= 220) {
      indx.push_back(2);
      indx.push_back(3);
      indx.push_back(5);
      indx.push_back(13);
    }
    if(sg >= 221 && sg <= 230) {
      indx.push_back(2);
      indx.push_back(3);
      indx.push_back(5);
      indx.push_back(13);
      indx.push_back(25);
    }
    return indx;
  }
} //namespace SYM

//DX20190215 - moved higher up -START
// **********************************************************************************************************************
// symop::clear
// **********************************************************************************************************************
namespace SYM {
  void symop::clear() {
    symbol = "";
    direction.clear();
    screwglide.clear();
    shift.clear();
  }
} //namespace SYM
//DX20190215 - moved higher up -END

//DX20190215 [OBSOLETE] namespace SYM {
//DX20190215 [OBSOLETE] #ifdef AFLOW_SYMMETRY_MULTITHREADS_ENABLE
//DX20190215 [OBSOLETE]   thread_local vector<vector<symop> > generators;
//DX20190215 [OBSOLETE]   thread_local vector<int> sgindex;
//DX20190215 [OBSOLETE] #else
//DX20190215 [OBSOLETE]   vector<vector<symop> > generators;
//DX20190215 [OBSOLETE]   vector<int> sgindex;
//DX20190215 [OBSOLETE] #endif
//DX20190215 [OBSOLETE] }
// **********************************************************************************************************************
// initgenerators
// **********************************************************************************************************************
bool SymmetryInformationITC::initgenerators(string axis_cell) {
  string setting_1 = aurostd::utype2string<int>(SG_SETTING_1); //DX20190131
  string setting_2 = aurostd::utype2string<int>(SG_SETTING_2); //DX20190131
  string setting_anrl = aurostd::utype2string<int>(SG_SETTING_ANRL); //DX20190131

  xvector<double> d000;
  xvector<double> d100;
  d100(1) = 1;
  xvector<double> d010;
  d010(2) = 1;
  xvector<double> d001;
  d001(3) = 1;
  xvector<double> d111;
  d111(1) = 1;
  d111(2) = 1;
  d111(3) = 1;
  xvector<double> d1nn;
  d1nn(1) = 1;
  d1nn(2) = -1;
  d1nn(3) = -1;
  xvector<double> d110;
  d110(1) = 1;
  d110(2) = 1;
  d110(3) = 0;
  xvector<double> d1n0;
  d1n0(1) = 1;
  d1n0(2) = -1;
  d1n0(3) = 0;
  xvector<double> dn1n;
  dn1n(1) = 1;
  dn1n(2) = 1;
  dn1n(3) = 1;
  xvector<double> dn01;
  dn01(1) = -1;
  dn01(2) = 0;
  dn01(3) = 1;
  xvector<double> dnn1;
  dnn1(1) = -1;
  dnn1(2) = -1;
  dnn1(3) = 1;
  xvector<double> d101;
  d101(1) = 1;
  d101(2) = 0;
  d101(3) = 1;
  xvector<double> d011;
  d011(1) = 0;
  d011(2) = 1;
  d011(3) = 1;
  xvector<double> d01n;
  d01n(1) = 0;
  d01n(2) = 1;
  d01n(3) = -1;
  xvector<double> dn10;
  dn10(1) = -1;
  dn10(2) = 1;
  dn10(3) = 0;
  xvector<double> d120;
  d120(1) = 1;
  d120(2) = 2;
  d120(3) = 0;
  xvector<double> d210;
  d210(1) = 2;
  d210(2) = 1;
  d210(3) = 0;

  bool includecentering = false;
  SYM::symop tmp;
  vector<SYM::symop> tmpv;

  generators.clear();
  sgindex.clear();
  //////////////////
  //SPACE GROUP 1
  //////////////////
  tmp.symbol = "1"; //DX20180614 - added so that sg=1 is not a special case in Wyckoff determination
  tmp.direction = 0.0 * d000; //DX20180614 - added so that sg=1 is not a special case in Wyckoff determination
  tmp.screwglide = 0.0 * d000; //DX20180614 - added so that sg=1 is not a special case in Wyckoff determination
  tmp.shift = 0.0 * d000; //DX20180614 - added so that sg=1 is not a special case in Wyckoff determination
  tmpv.push_back(tmp); //DX20180614 - added so that sg=1 is not a special case in Wyckoff determination
  tmp.clear(); //DX20180614 - added so that sg=1 is not a special case in Wyckoff determination
  generators.push_back(tmpv);
  tmp.clear();
  tmpv.clear();
  sgindex.push_back(1);
  //////////////////
  //SPACE GROUP 2
  //////////////////
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(2);
  //////////////////
  //SPACE GROUP 3
  //////////////////
  //**************
  //DX20180816 - added cell choice options - START
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl){
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(3);
  }
  else if(axis_cell == "m1" || axis_cell==setting_2){
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = 0.0*d000;
    tmp.shift = 0.0*d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(3);
  }
  else { 
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(3);
  }
  //DX20180816 - added cell choice options - END
  //////////////////
  //SPACE GROUP 4
  //////////////////
  //**************
  //DX20180816 - added cell choice options - START
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl){
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = .5 * d010;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(4);
  }
  else if(axis_cell == "m1" || axis_cell==setting_2){
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = .5*d001;
    tmp.shift = 0.0*d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(4);
  }
  else {
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = .5 * d010;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(4);
  }
  //DX20180816 - added cell choice options - END
  //////////////////
  //SPACE GROUP 5
  //////////////////
  if(axis_cell == "b1" || axis_cell == "b2" || axis_cell == "b3") {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d110;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(5);
  }
  else if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180817 - added option 2
    //**************
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d011;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(5);
  }
  //**************
  //CELL CHOICE 2
  //**************
  else if(axis_cell == "m2") {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d101;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(5);
  }
  //**************
  //CELL CHOICE 3
  //**************
  else if(axis_cell == "m3") {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d111;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(5);
  } else {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d110;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(5);
  }
  //////////////////
  //SPACE GROUP 6
  //////////////////
  //**************
  //DX20180816 - added cell choice options - START
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl){
    tmp.symbol = "m";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(6);
  }
  else if(axis_cell == "m1" || axis_cell==setting_2){
    tmp.symbol = "m";
    tmp.direction = d001;
    tmp.screwglide = 0.0*d000;
    tmp.shift = 0.0*d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(6);
  }
  else {
    tmp.symbol = "m";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(6);
  }
  //DX20180816 - added cell choice options - END
  //////////////////
  //SPACE GROUP 7
  //////////////////
  //**************
  //CELL CHOICE 1
  if(axis_cell == "b1" || axis_cell == "b2" || axis_cell == "b3" || axis_cell==setting_1 || axis_cell==setting_anrl) { //DX20180816 - added option 1
    tmp.symbol = "c";
    tmp.direction = d010;
    tmp.screwglide = 0.5 * d001;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(7);
  }
  //**************
  //**************
  //CELL CHOICE 1
  else if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180816 - added option 2
    tmp.symbol = "a";
    tmp.direction = d001;
    tmp.screwglide = 0.5 * d100;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(7);
  }
  //**************
  //CELL CHOICE 2
  else if(axis_cell == "m2") {
    tmp.symbol = "n";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(7);
  }
  //**************
  //CELL CHOICE 3
  else if(axis_cell == "m3") {
    tmp.symbol = "b";
    tmp.direction = d001;
    tmp.screwglide = 0.5 * d010;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(7);
  } else {
    tmp.symbol = "c";
    tmp.direction = d010;
    tmp.screwglide = 0.5 * d001;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(7);
  }
  //////////////////
  //SPACE GROUP 8
  //////////////////
  //CELL CHOICE 1
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl) { //DX20180816 - added option 1
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d110;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "m";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(8);
  } else if(axis_cell == "b2") {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d011;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "m";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(8);
  } else if(axis_cell == "b3") {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d111;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "m";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(8);
  }

  //**************
  //CELL CHOICE 1
  else if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180816 - added option 2
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d011;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "m";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(8);
  }
  //////////////////
  //**************
  //CELL CHOICE 2
  else if(axis_cell == "m2") {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d101;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "m";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(8);
  }

  //////////////////
  //**************
  //CELL CHOICE 3
  else if(axis_cell == "m3") {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d111;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "m";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(8);
  } else {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d110;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "m";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(8);
  }

  //////////////////
  //SPACE GROUP 9
  //////////////////
  //CELL CHOICE 1
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl) { //DX20180816 - added option 1
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d110;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "c";
    tmp.direction = d010;
    tmp.screwglide = 0.5 * d001;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(9);
  } else if(axis_cell == "b2") {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d011;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "c";
    tmp.direction = d010;
    tmp.screwglide = 0.5 * d001;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(9);
  } else if(axis_cell == "b3") {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d111;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "c";
    tmp.direction = d010;
    tmp.screwglide = 0.5 * d001;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(9);
  }

  //////////////////
  else if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180816 - added option 2
    //**************
    //CELL CHOICE 1
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d011;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "a";
    tmp.direction = d001;
    tmp.screwglide = 0.5 * d100;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(9);
  } else if(axis_cell == "m2") {
    //////////////////
    //**************
    //CELL CHOICE 2
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d101;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "n";
    tmp.direction = d001;
    tmp.screwglide = 0.5 * d110;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(9);
  } else if(axis_cell == "m3") {
    //////////////////
    //**************
    //CELL CHOICE 3
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d111;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "b";
    tmp.direction = d001;
    tmp.screwglide = 0.5 * d010;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(9);
  } else {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d110;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "c";
    tmp.direction = d010;
    tmp.screwglide = 0.5 * d001;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(9);
  }
  //////////////////
  //SPACE GROUP 10
  //////////////////
  if(axis_cell == "b1" || axis_cell == "b2" || axis_cell == "b3" || axis_cell==setting_1 || axis_cell==setting_anrl) { //DX20180816 - added option 1
    //**************
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(10);
  }
  if(axis_cell == "m1" || axis_cell == "m2" || axis_cell == "m3" || axis_cell==setting_2) { //DX20180816 - added option 2
    //**************
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(10);
  } else {
    //**************
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();  //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(10);
  }

  //////////////////
  //SPACE GROUP 11
  //////////////////
  if(axis_cell == "b1" || axis_cell == "b2" || axis_cell == "b3" || axis_cell==setting_1 || axis_cell==setting_anrl) { //DX20180816 - added option 1
    //**************
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = .5 * d010;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(11);
  }

  if(axis_cell == "m1" || axis_cell == "m2" || axis_cell == "m3" || axis_cell==setting_2) { //DX20180816 - added option 2
    //**************
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = .5 * d001;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(11);
  } else {
    //**************
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = .5 * d010;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(11);
  }

  //////////////////
  //SPACE GROUP 12
  //////////////////
  //**************
  //CELL CHOICE 1
  if(axis_cell == "b1" || axis_cell == "b2" || axis_cell == "b3" || axis_cell==setting_1 || axis_cell==setting_anrl) { //DX20180816 - added option 1
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d110;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(12);
  }
  //////////////////

  //**************
  //CELL CHOICE 1
  if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180816 - added option 2
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d011;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(12);
  }
  //////////////////
  //CELL CHOICE 2
  if(axis_cell == "m2") {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d101;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(12);
  }
  //////////////////
  //CELL CHOICE 3
  if(axis_cell == "m3") {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d111;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(12);
  } else {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d110;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(12);
  }
  //////////////////
  //SPACE GROUP 13
  //////////////////
  //CELL CHOICE 1
  if(axis_cell == "b1" || axis_cell == "b2" || axis_cell == "b3" || axis_cell==setting_1 || axis_cell==setting_anrl) { //DX20180816 - added option 1
    //**************
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .25 * d001;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(13);
  }
  //CELL CHOICE 1
  //**************
  if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180816 - added option 2
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .25 * d100;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(13);
  }
  if(axis_cell == "m2") {
    //////////////////
    //CELL CHOICE 2
    //**************
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .25 * d010;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(13);
  }
  //////////////////
  //CELL CHOICE 3
  //**************
  if(axis_cell == "m3") {
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .25 * d110;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(13);
  } else {
    //**************
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .25 * d001;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(13);
  }
  //////////////////
  //SPACE GROUP 14
  //////////////////
  //CELL CHOICE 1
  //**************
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl) { //DX20180816 - added option 1
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = .5 * d010;
    tmp.shift = .25 * d001;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(14);
  }
  //**************
  else if(axis_cell == "b2") {
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = .5 * d010;
    tmp.shift = .25 * d101;  //??? right direction???
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(14);
  }

  //      //**************
  else if(axis_cell == "b3") {
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = .5 * d010;
    tmp.shift = .25 * d100;  //??? right direction???
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(14);
  }

  //CELL CHOICE 1
  //**************
  else if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180816 - added option 2
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = .5 * d001;
    tmp.shift = .25 * d100;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(14);
  }
  //////////////////
  //CELL CHOICE 2
  //**************
  else if(axis_cell == "m2") {
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = .5 * d001;
    tmp.shift = .25 * d110;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(14);
  }
  //      //////////////////
  //      //CELL CHOICE 3
  //      //**************
  else if(axis_cell == "m3") {
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = .5 * d001;
    tmp.shift = .25 * d010;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(14);
  } else {
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = .5 * d010;
    tmp.shift = .25 * d001;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(14);
  }
  //////////////////
  //SPACE GROUP 15
  //////////////////
  //CELL CHOICE 1
  //**************
  if(axis_cell == "b1" || axis_cell == "b2" || axis_cell == "b3" || axis_cell==setting_1 || axis_cell==setting_anrl) { //DX20180816 - added option 1
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d110;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .25 * d001;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(15);
  }

  //CELL CHOICE 1
  //**************
  if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180816 - added option 2
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d011;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .25 * d100;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(15);
  }
  if(axis_cell == "m2") {
    //////////////////
    //CELL CHOICE 2
    //**************
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d101;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .25 * d110;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(15);
  }
  if(axis_cell == "m3") {
    //////////////////
    //CELL CHOICE 3
    //**************
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d111;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "2";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .25 * d010;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(15);
  } else {
    if(includecentering == true) {
      tmp.symbol = "1";
      tmp.direction = 0.0 * d000;
      tmp.screwglide = 0.0 * d000;
      tmp.shift = .5 * d110;
      tmpv.push_back(tmp);
      tmp.clear();
    }
    //**************
    tmp.symbol = "2";
    tmp.direction = d010;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .25 * d001;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "-1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = 0.0 * d000;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    generators.push_back(tmpv);
    tmpv.clear();
    sgindex.push_back(15);
  }
  /////////////////////////////////
  //ORTHORHOMBIC
  //ORTHORHOMBIC
  //ORTHORHOMBIC
  /////////////////////////////////
  //SPACE GROUP 16
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(16);
  /////////////////////////////////
  //SPACE GROUP 17
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(17);
  /////////////////////////////////
  //SPACE GROUP 18
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(18);
  /////////////////////////////////
  //SPACE GROUP 19
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(19);
  /////////////////////////////////
  //SPACE GROUP 20
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d110;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(20);
  /////////////////////////////////
  //SPACE GROUP 21
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d110;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(21);
  /////////////////////////////////
  //SPACE GROUP 22
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(22);
  /////////////////////////////////
  //SPACE GROUP 23
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(23);
  /////////////////////////////////
  //SPACE GROUP 24
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(24);
  /////////////////////////////////
  //SPACE GROUP 25
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(25);
  /////////////////////////////////
  //SPACE GROUP 26
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(26);
  /////////////////////////////////
  //SPACE GROUP 27
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(27);
  /////////////////////////////////
  //SPACE GROUP 28
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "a";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d100;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(28);
  /////////////////////////////////
  //SPACE GROUP 29
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "a";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d100;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(29);
  /////////////////////////////////
  //SPACE GROUP 30
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(30);
  /////////////////////////////////
  //SPACE GROUP 31
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "n";
  tmp.direction = d010;
  tmp.screwglide = .5 * d101;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(31);
  /////////////////////////////////
  //SPACE GROUP 32
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "a";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d100;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(32);
  /////////////////////////////////
  //SPACE GROUP 33
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "a";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d100;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(33);
  /////////////////////////////////
  //SPACE GROUP 34
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "n";
  tmp.direction = d010;
  tmp.screwglide = .5 * d101;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(34);
  /////////////////////////////////
  //SPACE GROUP 35
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d110;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(35);
  /////////////////////////////////
  //SPACE GROUP 36
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d110;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(36);
  /////////////////////////////////
  //SPACE GROUP 37
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d110;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(37);
  /////////////////////////////////
  //SPACE GROUP 38
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(38);
  /////////////////////////////////
  //SPACE GROUP 39
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(39);
  /////////////////////////////////
  //SPACE GROUP 40
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "a";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d100;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(40);
  /////////////////////////////////
  //SPACE GROUP 41
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "a";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d100;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(41);
  /////////////////////////////////
  //SPACE GROUP 42
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(42);
  /////////////////////////////////
  //SPACE GROUP 43
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "d";
  tmp.direction = d010;
  tmp.screwglide = .25 * d101;
  tmp.shift = (1.0 / 8.0) * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(43);
  /////////////////////////////////
  //SPACE GROUP 44
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(44);
  /////////////////////////////////
  //SPACE GROUP 45
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "a";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d100;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(45);
  /////////////////////////////////
  //SPACE GROUP 46
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "a";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d100;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(46);
  /////////////////////////////////
  //SPACE GROUP 47
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(47);
  /////////////////////////////////
  //SPACE GROUP 48 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.25 * d111;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(48);

  //**************
  //    ORIGIN CHOICE 2  tmp.symbol = "2";
  //    ORIGIN CHOICE 2  tmp.direction = d001;
  //    ORIGIN CHOICE 2  tmp.screwglide = 0.0*d000;
  //    ORIGIN CHOICE 2  tmp.shift = .25*d110;
  //    ORIGIN CHOICE 2  tmpv.push_back(tmp);
  //    ORIGIN CHOICE 2  tmp.clear();
  //    ORIGIN CHOICE 2  //**************
  //    ORIGIN CHOICE 2  tmp.symbol = "2";
  //    ORIGIN CHOICE 2  tmp.direction = d010;
  //    ORIGIN CHOICE 2  tmp.screwglide = 0.0*d000;
  //    ORIGIN CHOICE 2  tmp.shift = .25*d101;
  //    ORIGIN CHOICE 2  tmpv.push_back(tmp);
  //    ORIGIN CHOICE 2  tmp.clear();
  //    ORIGIN CHOICE 2  //**************
  //    ORIGIN CHOICE 2  tmp.symbol = "-1";
  //    ORIGIN CHOICE 2  tmp.direction = 0.0*d000;
  //    ORIGIN CHOICE 2  tmp.screwglide = 0.0*d000;
  //    ORIGIN CHOICE 2  tmp.shift = 0.0*d000;
  //    ORIGIN CHOICE 2  tmpv.push_back(tmp);
  //    ORIGIN CHOICE 2  tmp.clear();
  //    ORIGIN CHOICE 2  //**************
  //    ORIGIN CHOICE 2  generators.push_back(tmpv);
  //    ORIGIN CHOICE 2  tmpv.clear();
  //    ORIGIN CHOICE 2  sgindex.push_back(48);
  /////////////////////////////////
  //SPACE GROUP 49
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(49);
  /////////////////////////////////
  //SPACE GROUP 50 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(50);
  /////////////////////////////////
  //SPACE GROUP 51
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(51);
  /////////////////////////////////
  //SPACE GROUP 52
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d101;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(52);
  /////////////////////////////////
  //SPACE GROUP 53
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d101;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(53);
  /////////////////////////////////
  //SPACE GROUP 54
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(54);
  /////////////////////////////////
  //SPACE GROUP 55
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(55);
  /////////////////////////////////
  //SPACE GROUP 56
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(56);
  /////////////////////////////////
  //SPACE GROUP 57
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(57);
  /////////////////////////////////
  //SPACE GROUP 58
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d101;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(58);
  /////////////////////////////////
  //SPACE GROUP 59 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = 0.25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(59);
  /////////////////////////////////
  //SPACE GROUP 60
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(60);
  /////////////////////////////////
  //SPACE GROUP 61
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(61);
  /////////////////////////////////
  //SPACE GROUP 62
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(62);
  /////////////////////////////////
  //SPACE GROUP 63
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d110;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(63);
  /////////////////////////////////
  //SPACE GROUP 64
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d110;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(64);
  /////////////////////////////////
  //SPACE GROUP 65
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d110;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(65);
  /////////////////////////////////
  //SPACE GROUP 66
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d110;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(66);
  /////////////////////////////////
  //SPACE GROUP 67
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d110;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(67);
  /////////////////////////////////
  //SPACE GROUP 68 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d110;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.25 * d011;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(68);
  /////////////////////////////////
  //SPACE GROUP 69
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(69);
  /////////////////////////////////
  //SPACE GROUP 70 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = (1.0 / 8.0) * d111;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(70);
  /////////////////////////////////
  //SPACE GROUP 71
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(71);
  /////////////////////////////////
  //SPACE GROUP 72
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(72);
  /////////////////////////////////
  //SPACE GROUP 73
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(73);
  /////////////////////////////////
  //SPACE GROUP 74
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(74);
  ////////////TETRAGONAL/////////////////
  ////////////TETRAGONAL/////////////////
  ////////////TETRAGONAL/////////////////
  ////////////TETRAGONAL/////////////////
  /////////////////////////////////
  //SPACE GROUP 75
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(75);
  /////////////////////////////////
  //SPACE GROUP 76
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .25 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(76);
  /////////////////////////////////
  //SPACE GROUP 77
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(77);
  /////////////////////////////////
  //SPACE GROUP 78
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .75 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(78);
  /////////////////////////////////
  //SPACE GROUP 79
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(79);
  /////////////////////////////////
  //SPACE GROUP 80
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .25 * d001;
  tmp.shift = SYM::dir(-1.0 / 4.0, 1.0 / 4.0, 0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(80);
  /////////////////////////////////
  //SPACE GROUP 81
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(81);
  /////////////////////////////////
  //SPACE GROUP 82
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(82);
  /////////////////////////////////
  //SPACE GROUP 83
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(83);
  /////////////////////////////////
  //SPACE GROUP 84
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(84);
  /////////////////////////////////
  //SPACE GROUP 85 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .5 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(85);
  /////////////////////////////////
  //SPACE GROUP 86 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.5 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.25 * d111;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(86);
  /////////////////////////////////
  //SPACE GROUP 87
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(87);
  /////////////////////////////////
  //SPACE GROUP 88 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .25 * d001;
  tmp.shift = SYM::dir(-1.0 / 4.0, 1.0 / 4.0, 0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = SYM::dir(0.0, 1.0 / 4.0, 1.0 / 8.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(88);
  /////////////////////////////////
  //SPACE GROUP 88
  /////////////////////////////////
  //**************
  //  if(includecentering == true){
  //    tmp.symbol = "1";
  //    tmp.direction = 0.0*d000;
  //    tmp.screwglide = 0.0*d000;
  //    tmp.shift = .5*d111;
  //    tmpv.push_back(tmp);
  //    tmp.clear();
  //  }
  //  //**************
  //  tmp.symbol = "2";
  //  tmp.direction = d001;
  //  tmp.screwglide = .5*d001;
  //  tmp.shift = .25*d010;
  //  tmpv.push_back(tmp);
  //  tmp.clear();
  //  //**************
  //  tmp.symbol = "4+";
  //  tmp.direction = d001;
  //  tmp.screwglide = .25*d001;
  //  tmp.shift = SYM::dir(-1.0/2.0,1.0/4.0,0);
  //  tmpv.push_back(tmp);
  //  tmp.clear();
  //  //**************
  //  tmp.symbol = "-1";
  //  tmp.direction = 0.0*d000;
  //  tmp.screwglide = 0.0*d000;
  //  tmp.shift = 0.0*d000;
  //  tmpv.push_back(tmp);
  //  tmp.clear();
  //  //**************
  //  generators.push_back(tmpv);
  //  tmpv.clear();
  //  sgindex.push_back(88);
  /////////////////////////////////
  //SPACE GROUP 89
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(89);
  /////////////////////////////////
  //SPACE GROUP 90
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .5 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(90);
  /////////////////////////////////
  //SPACE GROUP 91
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .25 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(91);
  /////////////////////////////////
  //SPACE GROUP 92
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .25 * d001;
  tmp.shift = .5 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = SYM::dir(1.0 / 4.0, 0, 1.0 / 8.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(92);
  /////////////////////////////////
  //SPACE GROUP 93
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(93);
  /////////////////////////////////
  //SPACE GROUP 94
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .5 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d101;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(94);
  /////////////////////////////////
  //SPACE GROUP 95
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .75 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(95);
  /////////////////////////////////
  //SPACE GROUP 96
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .5 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .75 * d001;
  tmp.shift = .5 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = SYM::dir(1.0 / 4.0, 0, 3.0 / 8.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(96);
  /////////////////////////////////
  //SPACE GROUP 97
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(97);
  /////////////////////////////////
  //SPACE GROUP 98
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .25 * d001;
  tmp.shift = .25 * dn10;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = SYM::dir(1.0 / 4.0, 0, 3.0 / 8.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(98);
  /////////////////////////////////
  //SPACE GROUP 99
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(99);
  /////////////////////////////////
  //SPACE GROUP 100
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "a";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d100;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(100);
  /////////////////////////////////
  //SPACE GROUP 101
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(101);
  /////////////////////////////////
  //SPACE GROUP 102
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .5 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "n";
  tmp.direction = d010;
  tmp.screwglide = .5 * d101;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(102);
  /////////////////////////////////
  //SPACE GROUP 103
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(103);
  /////////////////////////////////
  //SPACE GROUP 104
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "n";
  tmp.direction = d010;
  tmp.screwglide = .5 * d101;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(104);
  /////////////////////////////////
  //SPACE GROUP 105
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(105);
  /////////////////////////////////
  //SPACE GROUP 106
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "a";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d100;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(106);
  /////////////////////////////////
  //SPACE GROUP 107
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(107);
  /////////////////////////////////
  //SPACE GROUP 108
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(108);
  /////////////////////////////////
  //SPACE GROUP 109
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .25 * d001;
  tmp.shift = .25 * dn10;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(109);
  /////////////////////////////////
  //SPACE GROUP 110
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .25 * d001;
  tmp.shift = .25 * dn10;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(110);
  /////////////////////////////////
  //SPACE GROUP 111
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(111);
  /////////////////////////////////
  //SPACE GROUP 112
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(112);
  /////////////////////////////////
  //SPACE GROUP 113
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(113);
  /////////////////////////////////
  //SPACE GROUP 114
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d101;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(114);
  /////////////////////////////////
  //SPACE GROUP 115
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(115);
  /////////////////////////////////
  //SPACE GROUP 116
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(116);
  /////////////////////////////////
  //SPACE GROUP 117
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "a";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d100;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(117);
  /////////////////////////////////
  //SPACE GROUP 118
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "n";
  tmp.direction = d010;
  tmp.screwglide = .5 * d101;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(118);
  /////////////////////////////////
  //SPACE GROUP 119
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(119);
  /////////////////////////////////
  //SPACE GROUP 120
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d010;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(120);
  /////////////////////////////////
  //SPACE GROUP 121
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(121);
  /////////////////////////////////
  //SPACE GROUP 122
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = SYM::dir(1.0 / 4.0, 0, 3.0 / 8.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(122);
  /////////////////////////////////
  //SPACE GROUP 123
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(123);
  /////////////////////////////////
  //SPACE GROUP 124
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(124);
  /////////////////////////////////
  //SPACE GROUP 125 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(125);
  /////////////////////////////////
  //SPACE GROUP 126 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.25 * d111;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(126);
  /////////////////////////////////
  //SPACE GROUP 127
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(127);
  /////////////////////////////////
  //SPACE GROUP 128
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d101;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(128);
  /////////////////////////////////
  //SPACE GROUP 129 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .5 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = 0.25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(129);
  /////////////////////////////////
  //SPACE GROUP 130 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .5 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d101;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(130);
  /////////////////////////////////
  //SPACE GROUP 131
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(131);
  /////////////////////////////////
  //SPACE GROUP 132
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(132);
  /////////////////////////////////
  //SPACE GROUP 133 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .5 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.25 * d111;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(133);
  /////////////////////////////////
  //SPACE GROUP 134 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.5 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.25 * d111;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(134);
  /////////////////////////////////
  //SPACE GROUP 135
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(135);
  /////////////////////////////////
  //SPACE GROUP 136
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .5 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d101;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(136);
  /////////////////////////////////
  //SPACE GROUP 137
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.5 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = 0.25 * d101;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.25 * d111;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(137);
  /////////////////////////////////
  //SPACE GROUP 138 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.5 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.25 * d111;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(138);
  /////////////////////////////////
  //SPACE GROUP 139
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(139);
  /////////////////////////////////
  //SPACE GROUP 140
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(140);
  /////////////////////////////////
  //SPACE GROUP 141 (ORIGIN CHOICE 1)
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .25 * d001;
  tmp.shift = SYM::dir(-1.0 / 4.0, 1.0 / 4.0, 0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = SYM::dir(1.0 / 4.0, 0.0, 3.0 / 8.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = SYM::dir(0.0, 1.0 / 4.0, 1.0 / 8.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(141);
  /////////////////////////////////
  //SPACE GROUP 142
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "4+";
  tmp.direction = d001;
  tmp.screwglide = .25 * d001;
  tmp.shift = SYM::dir(-1.0 / 4.0, 1.0 / 2.0, 0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(142);
  ////////TRIGONAL/////////////
  ////////TRIGONAL/////////////
  ////////TRIGONAL/////////////
  ////////TRIGONAL/////////////
  /////////////////////////////////
  //SPACE GROUP 143
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(143);
  /////////////////////////////////
  //SPACE GROUP 144
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = (1.0 / 3.0) * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(144);
  /////////////////////////////////
  //SPACE GROUP 145
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = (2.0 / 3.0) * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(145);
  /////////////////////////////////
  //SPACE GROUP 146
  /////////////////////////////////
  //CELL CHOICE 1
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(146);
  //**************
  //////////////////
  //CELL CHOICE 2
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = SYM::dir(2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(146);
  /////////////////////////////////
  //SPACE GROUP 147
  /////////////////////////////////
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(147);
  /////////////////////////////////
  //SPACE GROUP 148
  /////////////////////////////////
  //CELL CHOICE 1
  //**************
  //////////////////
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(148);
  //////////////////
  //CELL CHOICE 2
  //////////////////
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = SYM::dir(2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(148);
  /////////////////////////////////
  //SPACE GROUP 149
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d1n0;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(149);
  /////////////////////////////////
  //SPACE GROUP 150
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(150);
  /////////////////////////////////
  //SPACE GROUP 151
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = (1.0 / 3.0) * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d1n0;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = (1.0 / 3.0) * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(151);
  /////////////////////////////////
  //SPACE GROUP 152
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = (1.0 / 3.0) * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(152);
  /////////////////////////////////
  //SPACE GROUP 153
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = (2.0 / 3.0) * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d1n0;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = (1.0 / 6.0) * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(153);
  /////////////////////////////////
  //SPACE GROUP 154
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = (2.0 / 3.0) * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(154);
  /////////////////////////////////
  //SPACE GROUP 155
  /////////////////////////////////
  //////////////////
  //CELL CHOICE 1: RHOMBOHEDRAL
  //////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = dn01;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(155);
  //////////////////
  //CELL CHOICE 2:HEXAGONAL
  //////////////////
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = SYM::dir(2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(155);
  /////////////////////////////////
  //SPACE GROUP 156
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(156);
  /////////////////////////////////
  //SPACE GROUP 157
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d1n0;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(157);
  /////////////////////////////////
  //SPACE GROUP 158
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d110;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(158);
  /////////////////////////////////
  //SPACE GROUP 159
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d1n0;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(159);
  /////////////////////////////////
  //SPACE GROUP 160
  /////////////////////////////////
  //////////////////
  //CELL CHOICE 1: RHOMBOHEDRAL
  //////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = dn01;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(160);
  //////////////////
  //CELL CHOICE 2: HEXAGONAL
  //////////////////
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = SYM::dir(2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(160);
  /////////////////////////////////
  //SPACE GROUP 161
  /////////////////////////////////
  //////////////////
  //CELL CHOICE 1: RHOMBOHEDRAL
  //////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "n";
  tmp.direction = dn01;
  tmp.screwglide = .5 * d111;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(161);
  //////////////////
  //CELL CHOICE 2: HEXAGONAL
  //////////////////
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = SYM::dir(2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(161);
  /////////////////////////////////
  //SPACE GROUP 162
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d1n0;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(162);
  /////////////////////////////////
  //SPACE GROUP 163
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d1n0;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(163);
  /////////////////////////////////
  //SPACE GROUP 164
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(164);
  /////////////////////////////////
  //SPACE GROUP 165
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(165);
  /////////////////////////////////
  //SPACE GROUP 166
  /////////////////////////////////
  //CELL CHOICE 1:RHOMBOHEDRAL
  ////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = dn01;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(166);
  ////////////
  //CELL CHOICE 2:HEXAGONAL
  ////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = SYM::dir(2 / 3, 1 / 3, 1 / 3);
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(166);
  /////////////////////////////////
  //SPACE GROUP 167
  /////////////////////////////////
  ////////////
  //CELL CHOICE 1:RHOMBOHEDRAL
  ////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = dn01;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = SYM::dir(.5, .25, 0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(167);
  ////////////
  //CELL CHOICE 2:HEXAGONAL
  ////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = SYM::dir(2 / 3, 1 / 3, 1 / 3);
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(167);
  /////////////HEXAGONAL/////////////
  /////////////HEXAGONAL/////////////
  /////////////HEXAGONAL/////////////
  /////////////HEXAGONAL/////////////
  /////////////////////////////////
  //SPACE GROUP 168
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(168);
  /////////////////////////////////
  //SPACE GROUP 169
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = (1.0 / 3.0) * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(169);
  /////////////////////////////////
  //SPACE GROUP 170
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = (2.0 / 3.0) * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(170);
  /////////////////////////////////
  //SPACE GROUP 171
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = (2.0 / 3.0) * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(171);
  /////////////////////////////////
  //SPACE GROUP 172
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = (1.0 / 3.0) * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(172);
  /////////////////////////////////
  //SPACE GROUP 173
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(173);
  /////////////////////////////////
  //SPACE GROUP 174
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(174);
  /////////////////////////////////
  //SPACE GROUP 175
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(175);
  /////////////////////////////////
  //SPACE GROUP 176
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(176);
  /////////////////////////////////
  //SPACE GROUP 177
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(177);
  /////////////////////////////////
  //SPACE GROUP 178
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = (1.0 / 3.0) * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = (1 / 6) * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(178);
  /////////////////////////////////
  //SPACE GROUP 179
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = (2.0 / 3.0) * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = (1.0 / 3.0) * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(179);
  /////////////////////////////////
  //SPACE GROUP 180
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = (2.0 / 3.0) * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = (1.0 / 3.0) * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(180);
  /////////////////////////////////
  //SPACE GROUP 181
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = (1.0 / 3.0) * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = (1.0 / 6.0) * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(181);
  /////////////////////////////////
  //SPACE GROUP 182
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(182);
  /////////////////////////////////
  //SPACE GROUP 183
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(183);
  /////////////////////////////////
  //SPACE GROUP 184
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d110;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(184);
  /////////////////////////////////
  //SPACE GROUP 185
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d110;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(185);
  /////////////////////////////////
  //SPACE GROUP 186
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(186);
  /////////////////////////////////
  //SPACE GROUP 187
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(187);
  /////////////////////////////////
  //SPACE GROUP 188
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "c";
  tmp.direction = d110;
  tmp.screwglide = 0.5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(188);
  /////////////////////////////////
  //SPACE GROUP 189
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(189);
  /////////////////////////////////
  //SPACE GROUP 190
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(190);
  /////////////////////////////////
  //SPACE GROUP 191
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(191);
  /////////////////////////////////
  //SPACE GROUP 192
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(192);
  /////////////////////////////////
  //SPACE GROUP 193
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(193);
  /////////////////////////////////
  //SPACE GROUP 194
  /////////////////////////////////
  //**************
  tmp.symbol = "3+";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(194);
  ////////////CUBIC////////////////
  ////////////CUBIC////////////////
  ////////////CUBIC////////////////
  ////////////CUBIC////////////////
  /////////////////////////////////
  //SPACE GROUP 195
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(195);
  /////////////////////////////////
  //SPACE GROUP 196
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(196);
  /////////////////////////////////
  //SPACE GROUP 197
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = d001;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(197);
  /////////////////////////////////
  //SPACE GROUP 198
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(198);
  /////////////////////////////////
  //SPACE GROUP 199
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(199);
  /////////////////////////////////
  //SPACE GROUP 200
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(200);
  /////////////////////////////////
  //SPACE GROUP 201
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d101;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(201);
  /////////////////////////////////
  //SPACE GROUP 202
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(202);
  /////////////////////////////////
  //SPACE GROUP 203
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = (3.0 / 8.0) * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = (3.0 / 8.0) * d101;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(203);
  /////////////////////////////////
  //SPACE GROUP 204
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(204);
  /////////////////////////////////
  //SPACE GROUP 205
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(205);
  /////////////////////////////////
  //SPACE GROUP 206
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(206);
  /////////////////////////////////
  //SPACE GROUP 207
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(207);
  /////////////////////////////////
  //SPACE GROUP 208
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = .5 * d110;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(208);
  /////////////////////////////////
  //SPACE GROUP 209
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(209);
  /////////////////////////////////
  //SPACE GROUP 210
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = .5 * d110;
  tmp.shift = SYM::dir(0, -.25, 3.0 / 8.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(210);
  /////////////////////////////////
  //SPACE GROUP 211
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(211);
  /////////////////////////////////
  //SPACE GROUP 212
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = .5 * d110;
  tmp.shift = SYM::dir(0, .25, 3.0 / 8.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(212);
  /////////////////////////////////
  //SPACE GROUP 213
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = .5 * d110;
  tmp.shift = SYM::dir(0, -.25, 1.0 / 8.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(213);
  /////////////////////////////////
  //SPACE GROUP 214
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = .5 * d110;
  tmp.shift = SYM::dir(0, -.25, 1.0 / 8.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(214);
  /////////////////////////////////
  //SPACE GROUP 215
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d1n0;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(215);
  /////////////////////////////////
  //SPACE GROUP 216
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d1n0;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(216);
  /////////////////////////////////
  //SPACE GROUP 217
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "m";
  tmp.direction = d1n0;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(217);
  /////////////////////////////////
  //SPACE GROUP 218
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "n";
  tmp.direction = d1n0;
  tmp.screwglide = .5 * d111;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(218);
  /////////////////////////////////
  //SPACE GROUP 219
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "n";
  tmp.direction = d1n0;
  tmp.screwglide = .5 * d111;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(219);
  /////////////////////////////////
  //SPACE GROUP 220
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = .25 * d100;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "d";
  tmp.direction = d1n0;
  tmp.screwglide = .25 * d111;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(220);
  /////////////////////////////////
  //SPACE GROUP 221
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(221);
  /////////////////////////////////
  //SPACE GROUP 222
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d101;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(222);
  /////////////////////////////////
  //SPACE GROUP 223
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = .5 * d110;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(223);
  /////////////////////////////////
  //SPACE GROUP 224
  /////////////////////////////////
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d110;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = .25 * d101;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = .5 * d110;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(224);
  /////////////////////////////////
  //SPACE GROUP 225
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(225);
  /////////////////////////////////
  //SPACE GROUP 226
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = .5 * d110;
  tmp.shift = .25 * d001;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(226);
  /////////////////////////////////
  //SPACE GROUP 227
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = SYM::dir(3.0 / 8.0, 1.0 / 8.0, 0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = SYM::dir(1.0 / 8.0, 0, 3.0 / 8.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = .5 * d110;
  tmp.shift = .25 * SYM::dir(0, -1.0, 1.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(227);
  /////////////////////////////////
  //SPACE GROUP 228
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d011;
    tmpv.push_back(tmp);
    tmp.clear();
    //**************
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d101;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = SYM::dir(1.0 / 8.0, 3.0 / 8.0, 0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = SYM::dir(3.0 / 8.0, 0, 1.0 / 8.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = .5 * d110;
  tmp.shift = -.25 * d010;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(228);
  /////////////////////////////////
  //SPACE GROUP 229
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(229);
  /////////////////////////////////
  //SPACE GROUP 230
  /////////////////////////////////
  //**************
  if(includecentering == true) {
    tmp.symbol = "1";
    tmp.direction = 0.0 * d000;
    tmp.screwglide = 0.0 * d000;
    tmp.shift = .5 * d111;
    tmpv.push_back(tmp);
    tmp.clear();
  }
  //**************
  tmp.symbol = "2";
  tmp.direction = d001;
  tmp.screwglide = .5 * d001;
  tmp.shift = SYM::dir(1.0 / 4.0, 0, 0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d010;
  tmp.screwglide = .5 * d010;
  tmp.shift = SYM::dir(0, 0, 1.0 / 4.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "3+";
  tmp.direction = d111;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "2";
  tmp.direction = d110;
  tmp.screwglide = .5 * d110;
  tmp.shift = SYM::dir(0, -.25, 1.0 / 8.0);
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  tmp.symbol = "-1";
  tmp.direction = 0.0 * d000;
  tmp.screwglide = 0.0 * d000;
  tmp.shift = 0.0 * d000;
  tmpv.push_back(tmp);
  tmp.clear();
  //**************
  generators.push_back(tmpv);
  tmpv.clear();
  sgindex.push_back(230);

  return true;
}

//DX20190215 [OBSOLETE] #ifdef AFLOW_SYMMETRY_MULTITHREADS_ENABLE
//DX20190215 [OBSOLETE]   thread_local vector<string> gl_sgs;
//DX20190215 [OBSOLETE] #else
//DX20190215 [OBSOLETE]   vector<string> gl_sgs;
//DX20190215 [OBSOLETE] #endif
bool SymmetryInformationITC::initsgs(string axis_cell) {

  string setting_1 = aurostd::utype2string<int>(SG_SETTING_1);
  string setting_2 = aurostd::utype2string<int>(SG_SETTING_2);
  string setting_anrl = aurostd::utype2string<int>(SG_SETTING_ANRL);
  //DX20190215 [OBSOLETE]extern vector<string> sgs;
  vector<string> sgs;

  string sg1 = "Space Group 1\n (0,0,0) \n 1 a 1 (x,y,z) \n ";

  string sg2 = "Space Group 2\n  (0,0,0) \n 2 i 1 (x,y,z) (-x,-y,-z) \n 1 h -1 (.5,.5,.5) \n 1 g -1 (0,.5,.5) \n 1 f -1 (.5,0.,.5) \n 1 e -1 (.5,.5,0) \n 1 d -1 (.5,0,0) \n 1 c -1 (0,.5,0) \n 1 b -1 (0,0,.5) \n 1 a -1 (0,0,0) \n ";

  //DX20180818 - added unique axis choices - START
  string sg3 = "";
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl){
    sg3 = "Space Group 3\n (0,0,0) \n 2 e 1 (x,y,z) (-x,y,-z) \n 1 d 2 (.5,y,.5) \n 1 c 2 (.5,y,0) \n 1 b 2 (0,y,.5) \n 1 a 2 (0,y,0) \n ";
  } else if(axis_cell == "m1" || axis_cell==setting_2){
    sg3 = "Space Group 3\n (0,0,0) \n 2 e 1 (x,y,z) (-x,-y,z) \n 1 d 2 (.5,.5,z) \n 1 c 2 (0,.5,z) \n 1 b 2 (.5,0,z) \n 1 a 2 (0,0,z) \n ";
  } else {
    sg3 = "Space Group 3\n (0,0,0) \n 2 e 1 (x,y,z) (-x,y,-z) \n 1 d 2 (.5,y,.5) \n 1 c 2 (.5,y,0) \n 1 b 2 (0,y,.5) \n 1 a 2 (0,y,0) \n ";
  }
  //DX20180818 - added unique axis choices - END

  //DX20180818 - added unique axis choices - START
  string sg4 = "";
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl){
    sg4 = "Space Group 4\n  (0,0,0) \n 2 a 1 (x,y,z) (-x,y+.5,-z) \n ";
  } else if(axis_cell == "m1" || axis_cell==setting_2){
    sg4 = "Space Group 4\n  (0,0,0) \n 2 a 1 (x,y,z) (-x,-y,z+.5) \n ";
  } else {
    sg4 = "Space Group 4\n  (0,0,0) \n 2 a 1 (x,y,z) (-x,y+.5,-z) \n ";
  }
  //DX20180818 - added unique axis choices - END

  string sg5 = "";
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl) {
    sg5 = "Space Group 5b1\n (0,0,0) (.5,.5,0) \n 4 c 1 (x,y,z) (-x,y,-z) \n 2 b 2 (0,y,.5) \n 2 a 2 (0,y,0) \n ";
  } else if(axis_cell == "b2") {
    sg5 = "Space Group 5b2\n (0,0,0) (0,.5,.5) \n 4 c 1 (x,y,z) (-x,y,-z) \n 2 b 2 (.5,y,.5) \n 2 a 2 (0,y,0) \n ";
  } else if(axis_cell == "b3") {
    sg5 = "Space Group 5b3\n (0,0,0) (.5,.5,.5) \n 4 c 1 (x,y,z) (-x,y,-z) \n 2 b 2 (.5,y,0) \n 2 a 2 (0,y,0) \n ";
  } else if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180816 - added axis-c as option 2
    sg5 = "Space Group 5c1\n (0,0,0) (0,.5,.5) \n 4 c 1 (x,y,z) (-x,-y,z) \n 2 b 2 (.5,0,z) \n 2 a 2 (0,0,z) \n ";
  } else if(axis_cell == "m2") {
    sg5 = "Space Group 5c2\n (0,0,0) (.5,0,.5) \n 4 c 1 (x,y,z) (-x,-y,z) \n 2 b 2 (.5,.5,z) \n 2 a 2 (0,0,z) \n ";
  } else if(axis_cell == "m3") {
    sg5 = "Space Group 5c3\n (0,0,0) (.5,.5,.5) \n 4 c 1 (x,y,z) (-x,-y,z) \n 2 b 2 (0,.5,z) \n 2 a 2 (0,0,z) \n ";
  } else {
    sg5 = "Space Group 5\n (0,0,0) (.5,.5,0) \n 4 c 1 (x,y,z) (-x,y,-z) \n 2 b 2 (0,y,.5) \n 2 a 2 (0,y,0) \n ";
  }

  //DX20180818 - added unique axis choices - START
  string sg6 = "";
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl){
    sg6 = "Space Group 6\n (0,0,0) \n 2 c 1 (x,y,z) (x,-y,z) \n 1 b m (x,.5,z) \n 1 a m (x,0,z) \n ";
  } else if(axis_cell == "m1" || axis_cell==setting_2){
    sg6 = "Space Group 6\n (0,0,0) \n 2 c 1 (x,y,z) (x,y,-z) \n 1 b m (x,y,.5) \n 1 a m (x,y,0) \n ";
  } else {
    sg6 = "Space Group 6\n (0,0,0) \n 2 c 1 (x,y,z) (x,-y,z) \n 1 b m (x,.5,z) \n 1 a m (x,0,z) \n ";
  }
  //DX20180818 - added unique axis choices - END

  string sg7 = "";
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl) {
    sg7 = "Space Group 7b1\n (0,0,0) \n 2 a 1 (x,y,z) (x,-y,z+.5) \n ";
  } else if(axis_cell == "b2") {
    sg7 = "Space Group 7b2\n (0,0,0) \n 2 a 1 (x,y,z) (x+.5,-y,z+.5) \n ";
  } else if(axis_cell == "b3") {
    sg7 = "Space Group 7b3\n (0,0,0) \n 2 a 1 (x,y,z) (x+.5,-y,z) \n ";
  } else if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180816 - added axis-c as option 2
    sg7 = "Space Group 7c1\n (0,0,0) \n 2 a 1 (x,y,z) (x+.5,y,-z) \n ";
  } else if(axis_cell == "m2") {
    sg7 = "Space Group 7c2\n (0,0,0) \n 2 a 1 (x,y,z) (x+.5,y+.5,-z) \n ";
  } else if(axis_cell == "m3") {
    sg7 = "Space Group 7c3\n (0,0,0) \n 2 a 1 (x,y,z) (x,y+.5,-z) \n ";
  } else {
    sg7 = "Space Group 7\n (0,0,0) \n 2 a 1 (x,y,z) (x,-y,z+.5) \n ";
  }

  string sg8 = "";
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl) {
    sg8 = "Space Group 8\n (0,0,0) (.5,.5,0) \n 4 b 1 (x,y,z) (x,-y,z) \n 2 a m (x,0,z) \n ";
  } else if(axis_cell == "b2") {
    sg8 = "Space Group 8\n (0,0,0) (0,.5,.5) \n 4 b 1 (x,y,z) (x,-y,z) \n 2 a m (x,0,z) \n ";
  } else if(axis_cell == "b3") {
    sg8 = "Space Group 8\n (0,0,0) (.5,.5,.5) \n 4 b 1 (x,y,z) (x,-y,z) \n 2 a m (x,0,z) \n ";
  } else if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180816 - added axis-c as option 2
    sg8 = "Space Group 8\n (0,0,0) (0,.5,.5) \n 4 b 1 (x,y,z) (x,y,-z) \n 2 a m (x,y,0) \n ";
  } else if(axis_cell == "m2") {
    sg8 = "Space Group 8\n (0,0,0) (.5,0,.5) \n 4 b 1 (x,y,z) (x,y,-z) \n 2 a m (x,y,0) \n ";
  } else if(axis_cell == "m3") {
    sg8 = "Space Group 8\n (0,0,0) (.5,.5,.5) \n 4 b 1 (x,y,z) (x,y,-z) \n 2 a m (x,y,0) \n ";
  } else {
    sg8 = "Space Group 8\n (0,0,0) (.5,.5,0) \n 4 b 1 (x,y,z) (x,-y,z) \n 2 a m (x,0,z) \n ";
  }

  string sg9 = "";
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl) {
    sg9 = "Space Group 9b1\n (0,0,0) (.5,.5,0) \n 4 a 1 (x,y,z) (x,-y,z+.5) \n ";
  } else if(axis_cell == "b2") {
    sg9 = "Space Group 9b2\n (0,0,0) (0,.5,.5) \n 4 a 1 (x,y,z) (x+.5,-y,z+.5) \n ";
  } else if(axis_cell == "b3") {
    sg9 = "Space Group 9b3\n (0,0,0) (.5,.5,.5) \n 4 a 1 (x,y,z) (x+.5,-y,z) \n ";
  } else if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180816 - added axis-c as option 2
    sg9 = "Space Group 9c1\n (0,0,0) (0,.5,.5) \n 4 a 1 (x,y,z) (x+.5,y,-z) \n ";
  } else if(axis_cell == "m2") {
    sg9 = "Space Group 9c2\n (0,0,0) (.5,0,.5) \n 4 a 1 (x,y,z) (x+.5,y+.5,-z) \n ";
  } else if(axis_cell == "m3") {
    sg9 = "Space Group 9c3\n (0,0,0) (.5,.5,.5) \n 4 a 1 (x,y,z) (x,y+.5,-z) \n ";
  } else {
    sg9 = "Space Group 9\n (0,0,0) (.5,.5,0) \n 4 a 1 (x,y,z) (x,-y,z+.5) \n ";
  }

  string sg10 = "";
  if(axis_cell == "b1" || axis_cell == "b2" || axis_cell == "b3" || axis_cell==setting_1 || axis_cell==setting_anrl) {
    sg10 = "Space Group 10\n (0,0,0) \n 4 o 1 (x,y,z) (-x,y,-z) (-x,-y,-z) (x,-y,z) \n 2 n m (x,.5,z) (-x,.5,-z) \n 2 m m (x,0,z) (-x,0,-z) \n 2 l 2 (.5,y,.5) (.5,-y,.5) \n 2 k 2 (0,y,.5) (0,-y,.5) \n 2 j 2 (.5,y,0) (.5,-y,0) \n 2 i 2 (0,y,0) (0,-y,0) \n 1 h 2/m (.5,.5,.5) \n 1 g 2/m (.5,0,.5) \n 1 f 2/m (0,.5,.5) \n 1 e 2/m (.5,.5,0) \n 1 d 2/m (.5,0,0) \n 1 c 2/m (0,0,.5) \n 1 b 2/m (0,.5,0) \n 1 a 2/m (0,0,0) \n ";
  } else if(axis_cell == "m1" || axis_cell == "m2" || axis_cell == "m3" || axis_cell==setting_2) { //DX20180816 - added axis-c as option 2
    sg10 = "Space Group 10\n (0,0,0) \n 4 o 1 (x,y,z) (-x,-y,z) (-x,-y,-z) (x,y,-z) \n 2 n m (x,y,.5) (-x,-y,.5) \n 2 m m (x,y,0) (-x,-y,0) \n 2 l 2 (.5,.5,z) (.5,.5,-z) \n 2 k 2 (.5,0,z) (.5,0,-z) \n 2 j 2 (0,.5,z) (0,.5,-z) \n 2 i 2 (0,0,z) (0,0,-z) \n 1 h 2/m (.5,.5,.5) \n 1 g 2/m (.5,.5,0) \n 1 f 2/m (.5,0,.5) \n 1 e 2/m (0,.5,.5) \n 1 d 2/m (0,.5,0) \n 1 c 2/m (.5,0,0) \n 1 b 2/m (0,0,.5) \n 1 a 2/m (0,0,0) \n ";
  } else {
    sg10 = "Space Group 10\n (0,0,0) \n 4 o 1 (x,y,z) (-x,y,-z) (-x,-y,-z) (x,-y,z) \n 2 n m (x,.5,z) (-x,.5,-z) \n 2 m m (x,0,z) (-x,0,-z) \n 2 l 2 (.5,y,.5) (.5,-y,.5) \n 2 k 2 (0,y,.5) (0,-y,.5) \n 2 j 2 (.5,y,0) (.5,-y,0) \n 2 i 2 (0,y,0) (0,-y,0) \n 1 h 2/m (.5,.5,.5) \n 1 g 2/m (.5,0,.5) \n 1 f 2/m (0,.5,.5) \n 1 e 2/m (.5,.5,0) \n 1 d 2/m (.5,0,0) \n 1 c 2/m (0,0,.5) \n 1 b 2/m (0,.5,0) \n 1 a 2/m (0,0,0) \n ";
  }

  string sg11 = "";
  if(axis_cell == "b1" || axis_cell == "b2" || axis_cell == "b3" || axis_cell==setting_1 || axis_cell==setting_anrl) {
    sg11 = "Space Group 11\n (0,0,0) \n 4 f 1 (x,y,z) (-x,y+.5,-z) (-x,-y,-z) (x,-y+.5,z) \n 2 e m (x,.25,z) (-x,.75,-z) \n 2 d -1 (.5,0,.5) (.5,.5,.5) \n 2 c -1 (0,0,.5) (0,.5,.5) \n 2 b -1 (.5,0,0) (.5,.5,0) \n 2 a -1 (0,0,0) (0,.5,0) \n ";
  } else if(axis_cell == "m1" || axis_cell == "m2" || axis_cell == "m3" || axis_cell==setting_2) { //DX20180816 - added axis-c as option 2
    sg11 = "Space Group 11\n (0,0,0) \n 4 f 1 (x,y,z) (-x,-y,z+.5) (-x,-y,-z) (x,y,-z+.5) \n 2 e m (x,y,.25) (-x,-y,.75) \n 2 d -1 (.5,.5,0) (.5,.5,.5) \n 2 c -1 (.5,0,0) (.5,0,.5) \n 2 b -1 (0,.5,0) (0,.5,.5) \n 2 a -1 (0,0,0) (0,0,.5) \n ";
  } else {
    sg11 = "Space Group 11\n (0,0,0) \n 4 f 1 (x,y,z) (-x,y+.5,-z) (-x,-y,-z) (x,-y+.5,z) \n 2 e m (x,.25,z) (-x,.75,-z) \n 2 d -1 (.5,0,.5) (.5,.5,.5) \n 2 c -1 (0,0,.5) (0,.5,.5) \n 2 b -1 (.5,0,0) (.5,.5,0) \n 2 a -1 (0,0,0) (0,.5,0) \n ";
  }

  string sg12 = "";
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl) {
    sg12 = "Space Group 12b1\n (0,0,0) (.5,.5,0) \n 8 j 1 (x,y,z) (-x,y,-z) (-x,-y,-z) (x,-y,z) \n 4 i m (x,0,z) (-x,0,-z) \n 4 h 2 (0,y,.5) (0,-y,.5) \n 4 g 2 (0,y,0) (0,-y,0) \n 4 f -1 (.25,.25,.5) (.75,.25,.5) \n 4 e -1 (.25,.25,0) (.75,.25,0) \n 2 d 2/m (0,.5,.5) \n 2 c 2/m (0,0,.5) \n 2 b 2/m (0,.5,0) \n 2 a 2/m (0,0,0) \n ";
  } else if(axis_cell == "b2") {
    sg12 = "Space Group 12b2\n (0,0,0) (0,.5,.5) \n 8 j 1 (x,y,z) (-x,y,-z) (-x,-y,-z) (x,-y,z) \n 4 i m (x,0,z) (-x,0,-z) \n 4 h 2 (.5,y,.5) (.5,-y,.5) \n 4 g 2 (0,y,0) (0,-y,0) \n 4 f -1 (.5,.25,.75) (.5,.25,.25) \n 4 e -1 (0,.25,.25) (0,.25,.75) \n 2 d 2/m (.5,.5,.5) \n 2 c 2/m (.5,0,.5) \n 2 b 2/m (0,.5,0) \n 2 a 2/m (0,0,0) \n ";
  } else if(axis_cell == "b3") {
    sg12 = "Space Group 12b3\n (0,0,0) (.5,.5,.5) \n 8 j 1 (x,y,z) (-x,y,-z) (-x,-y,-z) (x,-y,z) \n 4 i m (x,0,z) (-x,0,-z) \n 4 h 2 (.5,y,0) (.5,-y,0) \n 4 g 2 (0,y,0) (0,-y,0) \n 4 f -1 (.25,.25,.75) (.75,.25,.25) \n 4 e -1 (.75,.25,.75) (.25,.25,.25) \n 2 d 2/m (.5,.5,0) \n 2 c 2/m (.5,0,0) \n 2 b 2/m (0,.5,0) \n 2 a 2/m (0,0,0) \n ";
  } else if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180816 - added axis-c as option 2
    sg12 = "Space Group 12c1\n (0,0,0) (0,.5,.5) \n 8 j 1 (x,y,z) (-x,-y,z) (-x,-y,-z) (x,y,-z) \n 4 i m (x,y,0) (-x,-y,0) \n 4 h 2 (.5,0,z) (.5,0,-z) \n 4 g 2 (0,0,z) (0,0,-z) \n 4 f -1 (.5,.25,.25) (.5,.75,.25) \n 4 e -1 (0,.25,.25) (0,.75,.25) \n 2 d 2/m (.5,0,.5) \n 2 c 2/m (.5,0,0) \n 2 b 2/m (0,0,.5) \n 2 a 2/m (0,0,0) \n ";
  } else if(axis_cell == "m2") {
    sg12 = "Space Group 12c2\n (0,0,0) (.5,0,.5) \n 8 j 1 (x,y,z) (-x,-y,z) (-x,-y,-z) (x,y,-z) \n 4 i m (x,y,0) (-x,-y,0) \n 4 h 2 (.5,.5,z) (.5,.5,-z) \n 4 g 2 (0,0,z) (0,0,-z) \n 4 f -1 (.75,.5,.25) (.25,.5,.25) \n 4 e -1 (.25,0,.25) (.75,0,.25) \n 2 d 2/m (.5,.5,.5) \n 2 c 2/m (.5,.5,0) \n 2 b 2/m (0,0,.5) \n 2 a 2/m (0,0,0) \n ";
  } else if(axis_cell == "m3") {
    sg12 = "Space Group 12c3\n (0,0,0) (.5,.5,.5) \n 8 j 1 (x,y,z) (-x,-y,z) (-x,-y,-z) (x,y,-z) \n 4 i m (x,y,0) (-x,-y,0) \n 4 h 2 (0,.5,z) (0,.5,-z) \n 4 g 2 (0,0,z) (0,0,-z) \n 4 f -1 (.75,.25,.25) (.25,.75,.25) \n 4 e -1 (.75,.75,.25) (.25,.25,.25) \n 2 d 2/m (0,.5,.5) \n 2 c 2/m (0,.5,0) \n 2 b 2/m (0,0,.5) \n 2 a 2/m (0,0,0) \n ";
  } else {
    sg12 = "Space Group 12\n (0,0,0) (.5,.5,0) \n 8 j 1 (x,y,z) (-x,y,-z) (-x,-y,-z) (x,-y,z) \n 4 i m (x,0,z) (-x,0,-z) \n 4 h 2 (0,y,.5) (0,-y,.5) \n 4 g 2 (0,y,0) (0,-y,0) \n 4 f -1 (.25,.25,.5) (.75,.25,.5) \n 4 e -1 (.25,.25,0) (.75,.25,0) \n 2 d 2/m (0,.5,.5) \n 2 c 2/m (0,0,.5) \n 2 b 2/m (0,.5,0) \n 2 a 2/m (0,0,0) \n ";
  }

  string sg13 = "";
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl) {
    sg13 = "Space Group 13b1\n (0,0,0) \n 4 g 1 (x,y,z) (-x,y,-z+.5) (-x,-y,-z) (x,-y,z+.5) \n 2 f 2 (.5,y,.25) (.5,-y,.75) \n 2 e 2 (0,y,.25) (0,-y,.75) \n 2 d -1 (.5,0,0) (.5,0,.5) \n 2 c -1 (0,.5,0) (0,.5,.5) \n 2 b -1 (.5,.5,0) (.5,.5,.5) \n 2 a -1 (0,0,0) (0,0,.5) \n ";
  } else if(axis_cell == "b2") {
    sg13 = "Space Group 13b2\n (0,0,0) \n 4 g 1 (x,y,z) (-x+.5,y,-z+.5) (-x,-y,-z) (x+.5,-y,z+.5) \n 2 f 2 (.75,y,.25) (.25,-y,.75) \n 2 e 2 (.75,y,.75) (.25,-y,.25) \n 	2 d -1 (0,0,.5) (.5,0,0) \n 2 c -1 (0,.5,0) (.5,.5,.5) \n 2 b -1 (0,.5,.5) (.5,.5,0) \n 2 a -1 (0,0,0) (.5,0,.5) \n ";
  } else if(axis_cell == "b3") {
    sg13 = "Space Group 13b3\n (0,0,0) \n 4 g 1 (x,y,z) (-x+.5,y,-z) (-x,-y,-z) (x+.5,-y,z) \n 2 f 2 (.75,y,.5) (.25,-y,.5) \n 2 e 2 (.25,y,0) (.75,-y,0) \n 2 d -1 (.5,0,.5) (0,0,.5) \n 2 c -1 (0,.5,0) (.5,.5,0) \n 2 b -1 (.5,.5,.5) (0,.5,.5) \n 2 a -1 (0,0,0) (.5,0,0) \n ";
  } else if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180816 - added axis-c as option 2
    sg13 = "Space Group 13c1\n (0,0,0) \n 4 g 1 (x,y,z) (-x+.5,-y,z) (-x,-y,-z) (x+.5,y,-z) \n 2 f 2 (.25,.5,z) (.75,.5,-z) \n 2 e 2 (.25,0,z) (.75,0,-z) \n 2 d -1 (0,.5,0) (.5,.5,0) \n 2 c -1 (0,0,.5) (.5,0,.5) \n 2 b -1 (0,.5,.5) (.5,.5,.5) \n 2 a -1 (0,0,0) (.5,0,0) \n ";
  } else if(axis_cell == "m2") {
    sg13 = "Space Group 13c2\n (0,0,0) \n 4 g 1 (x,y,z) (-x+.5,-y+.5,z) (-x,-y,-z) (x+.5,y+.5,-z) \n 2 f 2 (.25,.75,z) (.75,.25,-z) \n 2 e 2 (.75,.75,z) (.25,.25,-z) \n 	2 d -1 (.5,0,0) (0,.5,0) \n 2 c -1 (0,0,.5) (.5,.5,.5) \n 2 b -1 (.5,0,.5) (0,.5,.5) \n 2 a -1 (0,0,0) (.5,.5,0) \n ";
  } else if(axis_cell == "m3") {
    sg13 = "Space Group 13c3\n (0,0,0) \n 4 g 1 (x,y,z) (-x,-y+.5,z) (-x,-y,-z) (x,y+.5,-z) \n 2 f 2 (.5,.75,z) (.5,.25,-z) \n 2 e 2 (0,.25,z) (0,.75,-z) \n 2 d -1 (.5,.5,0) (.5,0,0) \n 2 c -1 (0,0,.5) (0,.5,.5) \n 2 b -1 (.5,.5,.5) (.5,0,.5) \n 2 a -1 (0,0,0) (0,.5,0) \n ";
  } else {
    sg13 = "Space Group 13\n (0,0,0) \n 4 g 1 (x,y,z) (-x,y,-z+.5) (-x,-y,-z) (x,-y,z+.5) \n 2 f 2 (.5,y,.25) (.5,-y,.75) \n 2 e 2 (0,y,.25) (0,-y,.75) \n 2 d -1 (.5,0,0) (.5,0,.5) \n 2 c -1 (0,.5,0) (0,.5,.5) \n 2 b -1 (.5,.5,0) (.5,.5,.5) \n 2 a -1 (0,0,0) (0,0,.5) \n ";
  }

  string sg14 = "";
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl) {
    sg14 = "Space Group 14\n  (0,0,0) \n 4 e 1 (x,y,z) (-x,y+.5,-z+.5) (-x,-y,-z) (x,-y+.5,z+.5) \n 2 d -1 (.5,0,.5) (.5,.5,0) \n 2 c -1 (0,0,.5) (0,.5,0) \n 2 b -1 (.5,0,0) (.5,.5,.5) \n 2 a -1 (0,0,0) (0,.5,.5) \n ";
  } else if(axis_cell == "b2") {
    sg14 = "Space Group 14\n  (0,0,0) \n 4 e 1 (x,y,z) (-x+.5,y+.5,-z+.5) (-x,-y,-z) (x+.5,-y+.5,z+.5) \n 2 d -1 (.5,0,0) (0,.5,.5) \n 2 c -1 (.5,0,.5) (0,.5,0) \n 2 b -1 (0,0,.5) (.5,.5,0) \n 2 a -1 (0,0,0) (.5,.5,.5) \n ";
  } else if(axis_cell == "b3") {
    sg14 = "Space Group 14\n  (0,0,0) \n 4 e 1 (x,y,z) (-x+.5,y+.5,-z) (-x,-y,-z) (x+.5,-y+.5,z) \n 2 d -1 (0,0,.5) (.5,.5,.5) \n 2 c -1 (.5,0,0) (0,.5,0) \n 2 b -1 (.5,0,.5) (0,.5,.5) \n 2 a -1 (0,0,0) (.5,.5,0) \n ";
  } else if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180816 - added axis-c as option 2
    sg14 = "Space Group 14\n  (0,0,0) \n 4 e 1 (x,y,z) (-x+.5,-y,z+.5) (-x,-y,-z) (x+.5,y,-z+.5) \n 2 d -1 (.5,.5,0) (0,.5,.5) \n 2 c -1 (.5,0,0) (0,0,.5) \n 2 b -1 (0,.5,0) (.5,.5,.5) \n 2 a -1 (0,0,0) (.5,0,.5) \n ";
  } else if(axis_cell == "m2") {
    sg14 = "Space Group 14\n  (0,0,0) \n 4 e 1 (x,y,z) (-x+.5,-y+.5,z+.5) (-x,-y,-z) (x+.5,y+.5,-z+.5) \n 2 d -1 (0,.5,0) (.5,0,.5) \n 2 c -1 (.5,.5,0) (0,0,.5) \n 2 b -1 (.5,0,0) (0,.5,.5) \n 2 a -1 (0,0,0) (.5,.5,.5) \n ";
  } else if(axis_cell == "m3") {
    sg14 = "Space Group 14\n  (0,0,0) \n 4 e 1 (x,y,z) (-x,-y+.5,z+.5) (-x,-y,-z) (x,y+.5,-z+.5) \n 2 d -1 (.5,0,0) (.5,.5,.5) \n 2 c -1 (0,.5,0) (0,0,.5) \n 2 b -1 (.5,.5,0) (.5,0,.5) \n 2 a -1 (0,0,0) (0,.5,.5) \n ";
  } else {
    sg14 = "Space Group 14\n  (0,0,0) \n 4 e 1 (x,y,z) (-x,y+.5,-z+.5) (-x,-y,-z) (x,-y+.5,z+.5) \n 2 d -1 (.5,0,.5) (.5,.5,0) \n 2 c -1 (0,0,.5) (0,.5,0) \n 2 b -1 (.5,0,0) (.5,.5,.5) \n 2 a -1 (0,0,0) (0,.5,.5) \n ";
  }

  string sg15 = "";
  if(axis_cell == "b1" || axis_cell==setting_1 || axis_cell==setting_anrl) {
    sg15 = "Space Group 15b1\n (0,0,0) (.5,.5,0) \n 8 f 1 (x,y,z) (-x,y,-z+.5) (-x,-y,-z) (x,-y,z+.5) \n 4 e 2 (0,y,.25) (0,-y,.75) \n 4 d -1 (.25,.25,.5) (.75,.25,0) \n 4 c -1 (.25,.25,0) (.75,.25,.5) \n 4 b -1 (0,.5,0) (0,.5,.5) \n 4 a -1 (0,0,0) (0,0,.5) \n ";
  } else if(axis_cell == "b2") {
    sg15 = "Space Group 15b2\n (0,0,0) (0,.5,.5) \n 8 f 1 (x,y,z) (-x+.5,y,-z+.5) (-x,-y,-z) (x+.5,-y,z+.5) \n 4 e 2 (.75,y,.75) (.25,-y,.25) \n 4 d -1 (.5,.25,.75) (0,.25,.75) \n 4 c -1 (0,.25,.25) (.5,.25,.25) \n 4 b -1 (0,.5,0) (.5,.5,.5) \n 4 a -1 (0,0,0) (.5,0,.5) \n ";
  } else if(axis_cell == "b3") {
    sg15 = "Space Group 15b3\n (0,0,0) (.5,.5,.5) \n 8 f 1 (x,y,z) (-x+.5,y,-z) (-x,-y,-z) (x+.5,-y,z) \n 4 e 2 (.25,y,0) (.75,-y,0) \n 4 d -1 (.25,.25,.75) (.25,.25,.25) \n 4 c -1 (.75,.25,.75) (.75,.25,.25) \n 4 b -1 (0,.5,0) (.5,.5,0) \n 4 a -1 (0,0,0) (.5,0,0) \n ";
  } else if(axis_cell == "m1" || axis_cell==setting_2) { //DX20180816 - added axis-c as option 2
    sg15 = "Space Group 15c1\n (0,0,0) (0,.5,.5) \n 8 f 1 (x,y,z) (-x+.5,-y,z) (-x,-y,-z) (x+.5,y,-z) \n 4 e 2 (.25,0,z) (.75,0,-z) \n 4 d -1 (.5,.25,.25) (0,.75,.25) \n 4 c -1 (0,.25,.25) (.5,.75,.25) \n 4 b -1 (0,0,.5) (.5,0,.5) \n 4 a -1 (0,0,0) (.5,0,0) \n ";
  } else if(axis_cell == "m2") {
    sg15 = "Space Group 15c2\n (0,0,0) (.5,0,.5) \n 8 f 1 (x,y,z) (-x+.5,-y+.5,z) (-x,-y,-z) (x+.5,y+.5,-z) \n 4 e 2 (.75,.75,z) (.25,.25,-z) \n 4 d -1 (.75,.5,.25) (.75,0,.25) \n 4 c -1 (.25,0,.25) (.25,.5,.25) \n 4 b -1 (0,0,.5) (.5,.5,.5) \n 4 a -1 (0,0,0) (.5,.5,0) \n ";
  } else if(axis_cell == "m3") {
    sg15 = "Space Group 15c3\n (0,0,0) (.5,.5,.5) \n 8 f 1 (x,y,z) (-x,-y+.5,z) (-x,-y,-z) (x,y+.5,-z) \n 4 e 2 (0,.25,z) (0,.75,-z) \n 4 d -1 (.75,.25,.25) (.25,.25,.25) \n 4 c -1 (.75,.75,.25) (.25,.75,.25) \n 4 b -1 (0,0,.5) (0,.5,.5) \n 4 a -1 (0,0,0) (0,.5,0) \n ";
  } else {
    sg15 = "Space Group 15\n (0,0,0) (.5,.5,0) \n 8 f 1 (x,y,z) (-x,y,-z+.5) (-x,-y,-z) (x,-y,z+.5) \n 4 e 2 (0,y,.25) (0,-y,.75) \n 4 d -1 (.25,.25,.5) (.75,.25,0) \n 4 c -1 (.25,.25,0) (.75,.25,.5) \n 4 b -1 (0,.5,0) (0,.5,.5) \n 4 a -1 (0,0,0) (0,0,.5) \n ";
  }

  string sg16 = "Space Group 16\n (0,0,0) \n 4 u 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) \n 2 t ..2 (.5,.5,z) (.5,.5,-z) \n 2 s ..2 (0,.5,z) (0,.5,-z) \n 2 r ..2 (.5,0,z) (0,0,-z) \n 2 q ..2 (0,0,z) (0,0,-z) \n 2 p .2. (.5,y,.5) (.5,-y,.5) \n 2 o .2. (.5,y,0) (.5,-y,0) \n 2 m .2. (0,y,0) (0,-y,0) \n 2 l 2.. (x,.5,.5) (-x,.5,.5) \n 2 k 2.. (x,.5,0) (-x,.5,0) \n 2 j 2.. (x,0,.5) (-x,0,.5) \n 2 i 2.. (x,0,0) (-x,0,0) \n 1 h 222 (.5,.5,.5) \n 1 g 222 (0,.5,.5) \n 1 f 222 (.5,0,.5) \n 1 e 222 (.5,.5,0) \n 1 d 222 (0,0,.5) \n 1 c 222 (0,.5,0) \n 1 b 222 (.5,0,0) \n 1 a 222 (0,0,0) \n ";

  string sg17 = "Space Group 17\n (0,0,0) \n 4 e 1 (x,y,z) (-x,-y,z+.5) (-x,y,-z+.5) (x,-y,-z) \n 2 d .2. (.5,y,.25) (.5,-y,.75) \n 2 c .2. (0,y,.25) (0,-y,.75) \n 2 b 2.. (x,.5,0) (-x,.5,.5) \n 2 a 2.. (x,0,0) (-x,0,.5) \n ";

  string sg18 = "Space Group 18\n (0,0,0) \n 4 c 1 (x,y,z) (-x,-y,z) (-x+.5,y+.5,-z) (x+.5,-y+.5,-z) \n 2 b ..2 (0,.5,z) (.5,0,-z) \n 2 a ..2 (0,0,z) (.5,.5,-z) \n ";

  string sg19 = "Space Group 19\n (0,0,0) \n 4 a 1 (x,y,z) (-x+.5,-y,z+.5) (-x,y+.5,-z+.5) (x+.5,-y+.5,-z) \n ";

  string sg20 = "Space Group 20\n (0,0,0) (.5,.5,0) \n 8 c 1 (x,y,z) (-x,-y,z+.5) (-x,y,-z+.5) (x,-y,-z) \n 4 b .2. (0,y,.25) (0,-y,.75) \n 4 a 2.. (x,0,0) (-x,0,.5) \n ";

  string sg21 = "Space Group 21\n (0,0,0) (.5,.5,0) \n 8 l 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) \n 4 k ..2 (.25,.25,z) (.75,.25,-z) \n 4 j ..2 (0,.5,z) (0,.5,-z) \n 4 i ..2 (0,0,z) (0,0,-z) \n 4 h .2. (0,y,.5) (0,-y,.5) \n 4 g .2. (0,y,0) (0,-y,0) \n 4 f 2.. (x,0,.5) (-x,0,.5) \n 4 e 2.. (x,0,0) (-x,0,0) \n 	2 d 222 (0,0,.5) \n 2 c 222 (.5,0,.5) \n 2 b 222 (0,.5,0)  \n 2 a 222 (0,0,0) \n ";

  string sg22 = "Space Group 22\n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 16 k 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) \n 8 j 2.. (x,.25,.25) (-x,.75,.25) \n 8 i .2. (.25,y,.25) (.75,-y,.25) \n 8 h ..2 (.25,.25,z) (.75,.25,-z) \n 8 g ..2 (0,0,z) (0,0,-z) \n 8 f .2. (0,y,0) (0,-y,0) \n 8 e 2.. (x,0,0) (-x,0,0) \n 4 d 222 (.25,.25,.75) \n 4 c 222 (.25,.25,.25) \n 4 b 222 (0,0,.5) \n 4 a 222 (0,0,0) \n ";

  string sg23 = "Space Group 23\n (0,0,0) (.5,.5,.5) \n 8 k 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) \n 4 j ..2 (0,.5,z) (0,.5,-z) \n 4 i ..2 (0,0,z) (0,0,-z) \n 4 h .2. (.5,y,0) (.5,-y,0) \n 4 g .2. (0,y,0) (0,-y,0) \n 4 f 2.. (x,0,.5) (-x,0,.5) \n 4 e 2.. (x,0,0) (-x,0,0) \n 2 d 222 (0,.5,0) \n 2 c 222 (0,0,.5) \n 2 b 222 (.5,0,0) \n 2 a 222 (0,0,0) \n ";

  string sg24 = "Space Group 24\n (0,0,0) (.5,.5,.5) \n 8 d 1 (x,y,z) (-x+.5,-y,z+.5) (-x,y+.5,-z+.5) (x+.5,-y+.5,-z) \n 4 c ..2 (0,.25,z) (0,.75,-z+.5) \n 4 b .2. (.25,y,0) (.25,-y,.5) \n 4 a 2.. (x,0,.25) (-x+.5,0,.75) \n ";

  string sg25 = "Space Group 25\n (0,0,0) \n 4 i 1 (x,y,z) (-x,-y,z) (x,-y,z) (-x,y,z) \n 2 h m.. (.5,y,z) (.5,-y,z) \n 2 g m.. (0,y,z) (0,-y,z) \n 2 f .m. (x,.5,z) (-x,.5,z) \n 2 e .m. (x,0,z) (-x,0,z) \n 1 d mm2 (.5,.5,z) \n 1 c mm2 (.5,0,z) \n 1 b mm2 (0,.5,z) \n 1 a mm2 (0,0,z) \n ";

  string sg26 = "Space Group 26\n (0,0,0) \n 4 c 1 (x,y,z) (-x,-y,z+.5) (x,-y,z+.5) (-x,y,z) \n 2 b m.. (.5,y,z) (.5,-y,z+.5) \n 2 a m.. (0,y,z) (0,-y,z+.5) \n ";

  string sg27 = "Space Group 27\n (0,0,0) \n 4 e 1 (x,y,z) (-x,-y,z) (x,-y,z+.5) (-x,y,z+.5) \n 2 d ..2 (.5,.5,z) (.5,.5,z+.5) \n 2 c ..2 (.5,0,z) (.5,0,z+.5) \n 2 b ..2 (0,.5,z) (0,.5,z+.5) \n 2 a ..2 (0,0,z) (0,0,z+.5) \n ";

  string sg28 = "Space Group 28\n (0,0,0) \n 4 d 1 (x,y,z) (-x,-y,z) (x+.5,-y,z) (-x+.5,y,z) \n 2 c m.. (.25,y,z) (.75,-y,z) \n 2 b ..2 (0,.5,z) (.5,.5,z) \n 2 a ..2 (0,0,z) (.5,0,z) \n ";

  string sg29 = "Space Group 29\n (0,0,0) \n 4 a 1 (x,y,z) (-x,-y,z+.5) (x+.5,-y,z) (-x+.5,y,z+.5) \n ";

  string sg30 = "Space Group 30\n (0,0,0) \n 4 c 1 (x,y,z) (-x,-y,z) (x,-y+.5,z+.5) (-x,y+.5,z+.5) \n 2 b ..2 (.5,0,z) (.5,.5,z+.5) \n 2 a ..2 (0,0,z) (0,.5,z+.5) \n ";

  string sg31 = "Space Group 31\n (0,0,0) \n 4 b 1 (x,y,z) (-x+.5,-y,z+.5) (x+.5,-y,z+.5) (-x,y,z) \n 2 a m.. (0,y,z) (.5,-y,z+.5)\n ";

  string sg32 = "Space Group 32\n (0,0,0) \n 4 c 1 (x,y,z) (-x,-y,z) (x+.5,-y+.5,z) (-x+.5,y+.5,z) \n 2 b ..2 (0,.5,z) (.5,0,z) \n 2 a ..2 (0,0,z) (.5,.5,z) \n ";

  string sg33 = "Space Group 33\n (0,0,0) \n 4 a 1 (x,y,z) (-x,-y,z+.5) (x+.5,-y+.5,z) (-x+.5,y+.5,z+.5) \n ";

  string sg34 = "Space Group 34\n (0,0,0) \n 4 c 1 (x,y,z) (-x,-y,z) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) \n 2 b ..2 (0,.5,z) (.5,0,z+.5) \n 2 a ..2 (0,0,z) (.5,.5,z+.5) \n ";

  string sg35 = "Space Group 35\n (0,0,0) (.5,.5,0) \n 8 f 1 (x,y,z) (-x,-y,z) (x,-y,z) (-x,y,z) \n 4 e m.. (0,y,z) (0,-y,z) \n 4 d .m. (x,0,z) (-x,0,z) \n 4 c ..2 (.25,.25,z) (.25,.75,z) \n 2 b mm2 (0,.5,z) \n 2 a mm2 (0,0,z) \n ";

  string sg36 = "Space Group 36\n 	(0,0,0) (.5,.5,0) \n 8 b 1 (x,y,z) (-x,-y,z+.5) (x,-y,z+.5) (-x,y,z) \n 4 a m.. (0,y,z) (0,-y,z+.5) \n ";

  string sg37 = "Space Group 37\n (0,0,0) (.5,.5,0) \n 8 d 1 (x,y,z) (-x,-y,z) (x,-y,z+.5) (-x,y,z+.5) \n 4 c ..2 (.25,.25,z) (.25,.75,z+.5) \n 4 b ..2 (0,.5,z) (0,.5,z+.5) \n 4 a ..2 (0,0,z) (0,0,z+.5) \n ";

  string sg38 = "Space Group 38\n (0,0,0) (0,.5,.5) \n 8 f 1 (x,y,z) (-x,-y,z) (x,-y,z) (-x,y,z) \n 4 e m.. (.5,y,z) (.5,-y,z) \n 4 d m.. (0,y,z) (0,-y,z) \n 4 c .m. (x,0,z) (-x,0,z) \n 2 b mm2 (.5,0,z) \n 2 a mm2 (0,0,z) \n ";

  string sg39 = "Space Group 39\n (0,0,0) (0,.5,.5) \n 8 d 1 (x,y,z) (-x,-y,z) (x,-y+.5,z) (-x,y+.5,z) \n 4 c .m. (x,.25,z) (-x,.75,z) \n 4 b ..2 (.5,0,z) (.5,.5,z) \n 4 a ..2 (0,0,z) (0,.5,z) \n ";

  string sg40 = "Space Group 40\n (0,0,0) (0,.5,.5) \n 8 c 1 (x,y,z) (-x,-y,z) (x+.5,-y,z) (-x+.5,y,z) \n 4 b m.. (.25,y,z) (.75,-y,z) \n 4 a ..2 (0,0,z) (.5,0,z) \n ";

  string sg41 = "Space Group 41\n (0,0,0) (0,.5,.5) \n 8 b 1 (x,y,z) (-x,-y,z) (x+.5,-y+.5,z) (-x+.5,y+.5,z) \n 4 a ..2 (0,0,z) (.5,.5,z) \n ";

  string sg42 = "Space Group 42\n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 16 e 1 (x,y,z) (-x,-y,z) (x,-y,z) (-x,y,z) \n 8 d .m. (x,0,z) (-x,0,z) \n 8 c m.. (0,y,z) (0,-y,z) \n 8 b ..2 (.25,.25,z) (.25,.75,z) \n 4 a mm2 (0,0,z) \n ";

  string sg43 = "Space Group 43\n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 16 b 1 (x,y,z) (-x,-y,z) (x+.25,-y+.25,z+.25) (-x+.25,y+.25,z+.25) \n 8 a ..2 (0,0,z) (.25,.25,z+.25) \n ";

  string sg44 = "Space Group 44\n (0,0,0) (.5,.5,.5) \n 8 e 1 (x,y,z) (-x,-y,z) (x,-y,z) (-x,y,z) \n 4 d m.. (0,y,z) (0,-y,z) \n 4 c .m. (x,0,z) (-x,0,z) \n 2 b mm2 (0,.5,z) \n 2 a mm2 (0,0,z) \n ";

  string sg45 = "Space Group 45\n (0,0,0) (.5,.5,.5) \n 8 c 1 (x,y,z) (-x,-y,z) (x+.5,-y+.5,z) (-x+.5,y+.5,z) \n 4 b ..2 (0,.5,z) (.5,0,z) \n 4 a ..2 (0,0,z) (.5,.5,z) \n ";

  string sg46 = "Space Group 46\n (0,0,0) (.5,.5,.5) \n 8 c 1 (x,y,z) (-x,-y,z) (x+.5,-y,z) (-x+.5,y,z) \n 4 b m.. (.25,y,z) (.75,-y,z) \n 4 a ..2 (0,0,z) (.5,0,z) \n ";

  string sg47 = "Space Group 47\n (0,0,0) \n 8 A 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (-x,-y,-z) (x,y,-z) (x,-y,z) (-x,y,z) \n 4 z ..m (x,y,.5) (-x,-y,.5) (-x,y,.5) (x,-y,.5) \n 4 y ..m (x,y,0) (-x,-y,0) (-x,y,0) (x,-y,0) \n 4 x .m. (x,.5,z) (-x,.5,z) (-x,.5,-z) (x,.5,-z) \n 4 w .m. (x,0,z) (-x,0,z) (-x,0,-z) (x,0,-z) \n 4 v m.. (.5,y,z) (.5,-y,z) (.5,y,-z) (.5,-y,-z) \n 4 u m.. (0,y,z) (0,-y,z) (0,y,-z) (0,-y,-z) \n 2 t mm2 (.5,.5,z) (.5,.5,-z) \n 2 s mm2 (.5,0,z) (.5,0,-z) \n 2 r mm2 (0,.5,z) (0,.5,-z) \n 2 q mm2 (0,0,z) (0,0,-z) \n 2 p m2m (.5,y,.5) (.5,-y,.5) \n 2 o m2m (.5,y,0) (.5,-y,0) \n 2 n m2m (0,y,.5) (0,-y,.5) \n 2 m m2m (0,y,0) (0,-y,0) \n 2 l 2mm (x,.5,.5) (-x,.5,.5) \n 2 k 2mm (x,.5,0) (-x,.5,0) \n 2 j 2mm (x,0,.5) (-x,0,.5) \n 2 i 2mm (x,0,0) (-x,0,0) \n 1 h mmm (.5,.5,.5) \n 1 g mmm (0,.5,.5) \n 1 f mmm (.5,.5,0) \n 1 e mmm (0,.5,0) \n 1 d mmm (.5,0,.5) \n 1 c mmm (0,0,.5) \n 1 b mmm (.5,0,0) \n 1 a mmm (0,0,0) \n "; //DX20200530 - fixed formatting; removed extra space between "k" Wyckoff position //DX20220511 - changed "alpha" to "A" to fix prototype commands

  //DX20180723 - added origin choice 2
  string sg48 = "";
  if(axis_cell==setting_1){
    sg48 = "Space Group 48\n (0,0,0) \n 8 m 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) \n 4 l ..2 (0,.5,z) (0,.5,-z) (.5,0,-z+.5) (.5,0,z+.5) \n 4 k ..2 (0,0,z) (0,0,-z) (.5,.5,-z+.5) (.5,.5,z+.5) \n 4 j .2. (.5,y,0) (.5,-y,0) (0,-y+.5,.5) (0,y+.5,.5) \n 4 i .2. (0,y,0) (0,-y,0) (.5,-y+.5,.5) (.5,y+.5,.5) \n 4 h 2.. (x,0,.5) (-x,0,.5) (-x+.5,.5,0) (x+.5,.5,0) \n 4 g 2.. (x,0,0) (-x,0,0) (-x+.5,.5,.5) (x+.5,.5,.5) \n 4 f -1 (.75,.75,.75) (.25,.25,.75) (.25,.75,.25) (.75,.25,.25) \n 4 e -1 (.25,.25,.25) (.75,.75,.25) (.75,.25,.75) (.25,.75,.75) \n 2 d 222 (0,.5,0) (.5,0,.5) \n 2 c 222 (0,0,.5) (.5,.5,0) \n 2 b 222 (.5,0,0) (0,.5,.5) \n 2 a 222 (0,0,0) (.5,.5,.5) \n ";
  }    
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg48 = "Space Group 48\n (0,0,0) \n 8 m 1 (x,y,z) (-x+.5,-y+.5,z) (-x+.5,y,-z+.5) (x,-y+.5,-z+.5) (-x,-y,-z) (x+.5,y+.5,-z) (x+.5,-y,z+.5) (-x,y+.5,z+.5) \n 4 l ..2 (.25,.75,z) (.25,.75,-z+.5) (.75,.25,-z) (.75,.25,z+.5) \n 4 k ..2 (.25,.25,z) (.25,.25,-z+.5) (.75,.75,-z) (.75,.75,z+.5) \n 4 j .2. (.75,y,.25) (.75,-y+.5,.25) (.25,-y,.75) (.25,y+.5,.75) \n 4 i .2. (.25,y,.25) (.25,-y+.5,.25) (.75,-y,.75) (.75,y+.5,.75) \n 4 h 2.. (x,.25,.75) (-x+.5,.25,.75) (-x,.75,.25) (x+.5,.75,.25) \n 4 g 2.. (x,.25,.25) (-x+.5,.25,.25) (-x,.75,.75) (x+.5,.75,.75) \n 4 f -1 (0,0,0) (.5,.5,0) (.5,0,.5) (0,.5,.5) \n 4 e -1 (.5,.5,.5) (0,0,.5) (0,.5,0) (.5,0,0) \n 2 d 222 (.25,.75,.25) (.75,.25,.75) \n 2 c 222 (.25,.25,.75) (.75,.75,.25) \n 2 b 222 (.75,.25,.25) (.25,.75,.75) \n 2 a 222 (.25,.25,.25) (.75,.75,.75) \n ";
  }
  else {
    sg48 = "Space Group 48\n (0,0,0) \n 8 m 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) \n 4 l ..2 (0,.5,z) (0,.5,-z) (.5,0,-z+.5) (.5,0,z+.5) \n 4 k ..2 (0,0,z) (0,0,-z) (.5,.5,-z+.5) (.5,.5,z+.5) \n 4 j .2. (.5,y,0) (.5,-y,0) (0,-y+.5,.5) (0,y+.5,.5) \n 4 i .2. (0,y,0) (0,-y,0) (.5,-y+.5,.5) (.5,y+.5,.5) \n 4 h 2.. (x,0,.5) (-x,0,.5) (-x+.5,.5,0) (x+.5,.5,0) \n 4 g 2.. (x,0,0) (-x,0,0) (-x+.5,.5,.5) (x+.5,.5,.5) \n 4 f -1 (.75,.75,.75) (.25,.25,.75) (.25,.75,.25) (.75,.25,.25) \n 4 e -1 (.25,.25,.25) (.75,.75,.25) (.75,.25,.75) (.25,.75,.75) \n 2 d 222 (0,.5,0) (.5,0,.5) \n 2 c 222 (0,0,.5) (.5,.5,0) \n 2 b 222 (.5,0,0) (0,.5,.5) \n 2 a 222 (0,0,0) (.5,.5,.5) \n ";
  }    

  string sg49 = "Space Group 49\n (0,0,0) \n 8 r 1 (x,y,z) (-x,-y,z) (-x,y,-z+.5) (x,-y,-z+.5) (-x,-y,-z) (x,y,-z) (x,-y,z+.5) (-x,y,z+.5) \n 4 q ..m (x,y,0) (-x,-y,0) (-x,y,.5) (x,-y,.5) \n 4 p ..2 (.5,0,z) (.5,0,-z+.5) (.5,0,-z) (.5,0,z+.5) \n 4 o ..2 (0,.5,z) (0,.5,-z+.5) (0,.5,-z) (0,.5,z+.5) \n 4 n ..2 (.5,.5,z) (.5,.5,-z+.5) (.5,.5,-z) (.5,.5,z+.5) \n 4 m ..2 (0,0,z) (0,0,-z+.5) (0,0,-z) (0,0,z+.5) \n 4 l .2. (.5,y,.25) (.5,-y,.25) (.5,-y,.75) (.5,y,.75) \n 4 k .2. (0,y,.25) (0,-y,.25) (0,-y,.75) (0,y,.75) \n 4 j 2.. (x,.5,.25) (-x,.5,.25) (-x,.5,.75) (x,.5,.75) \n 4 i 2.. (x,0,.25) (-x,0,.25) (-x,0,.75) (x,0,.75) \n 2 h 222 (.5,.5,.25) (.5,.5,.75) \n 2 g 222 (0,.5,.25) (0,.5,.75) \n 2 f 222 (.5,.0,.25) (.5,0,.75) \n 2 e 222 (0,0,.25) (0,0,.75) \n 2 d ..2/m (.5,0,0) (.5,0,.5) \n 2 c ..2/m (0,.5,0) (0,.5,.5) \n 2 b ..2/m (.5,.5,0) (.5,.5,.5) \n 2 a ..2/m (0,0,0) (0,0,.5) \n ";

  //DX20180723 - added origin choice 2
  string sg50 = "";
  if(axis_cell==setting_1){
    sg50 = "Space Group 50\n (0,0,0) \n 8 m 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (-x+.5,-y+.5,-z) (x+.5,y+.5,-z) (x+.5,-y+.5,z) (-x+.5,y+.5,z) \n 4 l ..2 (0,.5,z) (0,.5,-z) (.5,0,-z) (.5,0,z) \n 4 k ..2 (0,0,z) (0,0,-z) (.5,.5,-z) (.5,.5,z) \n 4 j .2. (0,y,.5) (0,-y,.5) (.5,-y+.5,.5) (.5,y+.5,.5) \n 4 i .2. (0,y,0) (0,-y,0) (.5,-y+.5,0) (.5,y+.5,0) \n 4 h 2.. (x,0,.5) (-x,0,.5) (-x+.5,.5,.5) (x+.5,.5,.5) \n 4 g 2.. (x,0,0) (-x,0,0) (-x+.5,.5,0) (x+.5,.5,0) \n 4 f -1 (.25,.25,.5) (.75,.75,.5) (.75,.25,.5) (.25,.75,.5) \n 4 e -1 (.25,.25,0) (.75,.75,0) (.75,.25,0) (.25,.75,0) \n 2 d 222 (0,0,.5) (.5,.5,.5) \n 2 c 222 (.5,.0,.5) (0,.5,.5) \n 2 b 222 (.5,0,0) (0,.5,0) \n 2 a 222 (0,0,0) (.5,.5,0) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg50 = "Space Group 50\n (0,0,0) \n 8 m 1 (x,y,z) (-x+.5,-y+.5,z) (-x+.5,y,-z) (x,-y+.5,-z) (-x,-y,-z) (x+.5,y+.5,-z) (x+.5,-y,z) (-x,y+.5,z) \n 4 l ..2 (.25,.75,z) (.25,.75,-z) (.75,.25,-z) (.75,.25,z) \n 4 k ..2 (.25,.25,z) (.25,.25,-z) (.75,.75,-z) (.75,.75,z) \n 4 j .2. (.25,y,.5) (.25,-y+.5,.5) (.75,-y,.5) (.75,y+.5,.5) \n 4 i .2. (.25,y,0) (.25,-y+.5,0) (.75,-y,0) (.75,y+.5,0) \n 4 h 2.. (x,.25,.5) (-x+.5,.25,.5) (-x,.75,.5) (x+.5,.75,.5) \n 4 g 2.. (x,.25,0) (-x+.5,.25,0) (-x,.75,0) (x+.5,.75,0) \n 4 f -1 (0,0,.5) (.5,.5,.5) (.5,0,.5) (0,.5,.5) \n 4 e -1 (0,0,0) (.5,.5,0) (.5,0,0) (0,.5,0) \n 2 d 222 (.25,.25,.5) (.75,.75,.5) \n 2 c 222 (.75,.25,.5) (.25,.75,.5) \n 2 b 222 (.75,.25,0) (.25,.75,0) \n 2 a 222 (.25,.25,0) (.75,.75,0) \n "; //DX20180806 - third position in 4f is (.5,0,.5) not (.5,0,.25)
  }
  else {
    sg50 = "Space Group 50\n (0,0,0) \n 8 m 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (-x+.5,-y+.5,-z) (x+.5,y+.5,-z) (x+.5,-y+.5,z) (-x+.5,y+.5,z) \n 4 l ..2 (0,.5,z) (0,.5,-z) (.5,0,-z) (.5,0,z) \n 4 k ..2 (0,0,z) (0,0,-z) (.5,.5,-z) (.5,.5,z) \n 4 j .2. (0,y,.5) (0,-y,.5) (.5,-y+.5,.5) (.5,y+.5,.5) \n 4 i .2. (0,y,0) (0,-y,0) (.5,-y+.5,0) (.5,y+.5,0) \n 4 h 2.. (x,0,.5) (-x,0,.5) (-x+.5,.5,.5) (x+.5,.5,.5) \n 4 g 2.. (x,0,0) (-x,0,0) (-x+.5,.5,0) (x+.5,.5,0) \n 4 f -1 (.25,.25,.5) (.75,.75,.5) (.75,.25,.5) (.25,.75,.5) \n 4 e -1 (.25,.25,0) (.75,.75,0) (.75,.25,0) (.25,.75,0) \n 2 d 222 (0,0,.5) (.5,.5,.5) \n 2 c 222 (.5,.0,.5) (0,.5,.5) \n 2 b 222 (.5,0,0) (0,.5,0) \n 2 a 222 (0,0,0) (.5,.5,0) \n ";
  }

  string sg51 = "Space Group 51\n (0,0,0) \n 8 l 1 (x,y,z) (-x+.5,-y,z) (-x,y,-z) (x+.5,-y,-z) (-x,-y,-z) (x+.5,y,-z) (x,-y,z) (-x+.5,y,z) \n 4 k m.. (.25,y,z) (.25,-y,z) (.75,y,-z) (.75,-y,-z) \n 4 j .m. (x,.5,z) (-x+.5,.5,z) (-x,.5,-z) (x+.5,.5,-z) \n 4 i .m. (x,0,z) (-x+.5,0,z) (-x,0,-z) (x+.5,0,-z) \n 4 h .2. (0,y,.5) (.5,-y,.5) (0,-y,.5) (.5,y,.5) \n 4 g .2. (0,y,0) (.5,-y,0) (0,-y,0) (.5,y,0) \n 2 f mm2 (.25,.5,z) (.75,.5,-z) \n 2 e mm2 (.25,0,z) (.75,0,-z) \n 2 d .2/m. (0,.5,.5) (.5,.5,.5) \n 2 c .2/m. (0,0,.5) (.5,0,.5) \n 2 b .2/m. (0,.5,0) (.5,.5,0) \n 2 a .2/m. (0,0,0) (.5,0,0) \n ";

  string sg52 = "Space Group 52 \n (0,0,0) \n 8 e 1 (x,y,z) (-x+.5,-y,z) (-x+.5,y+.5,-z+.5) (x,-y+.5,-z+.5) (-x,-y,-z) (x+.5,y,-z) (x+.5,-y+.5,z+.5) (-x,y+.5,z+.5) \n 4 d 2.. (x,.25,.25) (-x+.5,.75,.25) (-x,.75,.75) (x+.5,.25,.75) \n 4 c ..2 (.25,0,z) (.25,.5,-z+.5) (.75,0,-z) (.75,.5,z+.5) \n 4 b -1 (0,0,.5) (.5,0,.5) (.5,.5,0) (0,.5,0) \n 4 a -1 (0,0,0) (.5,0,0) (.5,.5,.5) (0,.5,.5) \n ";

  string sg53 = "Space Group 53 \n (0,0,0) \n 8 i 1 (x,y,z) (-x+.5,-y,z+.5) (-x+.5,y,-z+.5) (x,-y,-z) (-x,-y,-z) (x+.5,y,-z+.5) (x+.5,-y,z+.5) (-x,y,z) \n 4 h m.. (0,y,z) (.5,-y,z+.5) (.5,y,-z+.5) (0,-y,-z) \n 4 g .2. (.25,y,.25) (.25,-y,.75) (.75,-y,.75) (.75,y,.25) \n 4 f 2.. (x,.5,0) (-x+.5,.5,.5) (-x,.5,0) (x+.5,.5,.5) \n 4 e 2.. (x,0,0) (-x+.5,0,.5) (-x,0,0) (x+.5,0,.5) \n 2 d 2/m.. (0,.5,0) (.5,.5,.5) \n 2 c 2/m.. (.5,.5,0) (0,.5,.5) \n 2 b 2/m.. (.5,0,0) (0,0,.5) \n 2 a 2/m.. (0,0,0) (.5,0,.5) \n ";

  string sg54 = "Space Group 54\n (0,0,0) \n 8 f 1 (x,y,z) (-x+.5,-y,z) (-x,y,-z+.5) (x+.5,-y,-z+.5) (-x,-y,-z) (x+.5,y,-z) (x,-y,z+.5) (-x+.5,y,z+.5) \n 4 e ..2 (.25,.5,z) (.75,.5,-z+.5) (.75,.5,-z) (.25,.5,z+.5) \n 4 d ..2 (.25,0,z) (.75,0,-z+.5) (.75,0,-z) (.25,0,z+.5) \n 4 c .2. (0,y,.25) (.5,-y,.25) (0,-y,.75) (.5,y,.75) \n 4 b -1 (0,.5,0) (.5,.5,0) (0,.5,.5) (.5,.5,.5) \n 4 a -1 (0,0,0) (.5,0,0) (0,0,.5) (.5,0,.5) \n ";

  string sg55 =
    "Space Group 55\n (0,0,0) \n 8 i 1 (x,y,z) (-x,-y,z) (-x+.5,y+.5,-z) (x+.5,-y+.5,-z) (-x,-y,-z) (x,y,-z) (x+.5,-y+.5,z) (-x+.5,y+.5,z) \n 4 h ..m (x,y,.5) (-x,-y,.5) (-x+.5,y+.5,.5) (x+.5,-y+.5,.5) \n 4 g ..m (x,y,0) (-x,-y,0) (-x+.5,y+.5,0) (x+.5,-y+.5,0) \n 4 f ..2 (0,.5,z) (.5,0,-z) (0,.5,-z) (.5,0,z) \n 4 e ..2 (0,0,z) (.5,.5,-z) (0,0,-z) (.5,.5,z) \n 2 d ..2/m (0,.5,.5) (.5,0,.5) \n 2 c ..2/m (0,.5,0) (.5,0,0) \n 2 b ..2/m (0,0,.5) (.5,.5,.5) \n 2 a ..2/m (0,0,0) (.5,.5,0) \n "; //DX20200422 - extraneous slash removed

  string sg56 = "Space Group 56\n (0,0,0) \n 8 e 1 (x,y,z) (-x+.5,-y+.5,z) (-x,y+.5,-z+.5) (x+.5,-y,-z+.5) (-x,-y,-z) (x+.5,y+.5,-z) (x,-y+.5,z+.5) (-x+.5,y,z+.5) \n 4 d ..2 (.25,.75,z) (.75,.25,-z+.5) (.75,.25,-z) (.25,.75,z+.5) \n 4 c ..2 (.25,.25,z) (.75,.75,-z+.5) (.75,.75,-z) (.25,.25,z+.5) \n 4 b -1 (0,0,.5) (.5,.5,.5) (0,.5,0) (.5,0,0) \n 4 a -1 (0,0,0) (.5,.5,0) (0,.5,.5) (.5,0,.5) \n ";

  string sg57 = "Space Group 57\n (0,0,0) \n 8 e 1 (x,y,z) (-x,-y,z+.5) (-x,y+.5,-z+.5) (x,-y+.5,-z) (-x,-y,-z) (x,y,-z+.5) (x,-y+.5,z+.5) (-x,y+.5,z) \n 4 d ..m (x,y,.25) (-x,-y,.75) (-x,y+.5,.25) (x,-y+.5,.75) \n 4 c 2.. (x,.25,0) (-x,.75,.5) (-x,.75,0) (x,.25,.5) \n 4 b -1 (.5,0,0) (.5,0,.5) (.5,.5,.5) (.5,.5,0) \n 4 a -1 (0,0,0) (0,0,.5) (0,.5,.5) (0,.5,0) \n ";

  string sg58 = "Space Group 58\n (0,0,0) \n 8 h 1 (x,y,z) (-x,-y,z) (-x+.5,y+.5,-z+.5) (x+.5,-y+.5,-z+.5) (-x,-y,-z) (x,y,-z) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) \n 4 g ..m (x,y,0) (-x,-y,0) (-x+.5,y+.5,.5) (x+.5,-y+.5,.5) \n 4 f ..2 (0,.5,z) (.5,0,-z+.5) (0,.5,-z) (.5,0,z+.5) \n 4 e ..2 (0,0,z) (.5,.5,-z+.5) (0,0,-z) (.5,.5,z+.5) \n 2 d ..2/m (0,.5,.5) (.5,0,0) \n 2 c ..2/m (0,.5,0) (.5,0,.5) \n 2 b ..2/m (0,0,.5) (.5,.5,0) \n 2 a ..2/m (0,0,0) (.5,.5,.5) \n ";

  //DX20180723 - added origin choice 2
  string sg59 = "";
  if(axis_cell==setting_1){
    sg59 = "Space Group 59\n (0,0,0) \n 8 g 1 (x,y,z) (-x,-y,z) (-x+.5,y+.5,-z) (x+.5,-y+.5,-z) (-x+.5,-y+.5,-z) (x+.5,y+.5,-z) (x,-y,z) (-x,y,z) \n 4 f .m. (x,0,z) (-x,0,z) (-x+.5,.5,-z) (x+.5,.5,-z) \n 4 e m.. (0,y,z) (0,-y,z) (.5,y+.5,-z) (.5,-y+.5,-z) \n 4 d -1 (.25,.25,.5) (.75,.75,.5) (.25,.75,.5) (.75,.25,.5) \n 4 c -1 (.25,.25,0) (.75,.75,0) (.25,.75,0) (.75,.25,0) \n 2 b mm2 (0,.5,z) (.5,0,-z) \n 2 a mm2 (0,0,z) (.5,.5,-z) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg59 = "Space Group 59\n (0,0,0) \n 8 g 1 (x,y,z) (-x+.5,-y+.5,z) (-x,y+.5,-z) (x+.5,-y,-z) (-x,-y,-z) (x+.5,y+.5,-z) (x,-y+.5,z) (-x+.5,y,z) \n 4 f .m. (x,.25,z) (-x+.5,.25,z) (-x,.75,-z) (x+.5,.75,-z) \n 4 e m.. (.25,y,z) (.25,-y+.5,z) (.75,y+.5,-z) (.75,-y,-z) \n 4 d -1 (0,0,.5) (.5,.5,.5) (0,.5,.5) (.5,0,.5) \n 4 c -1 (0,0,0) (.5,.5,0) (0,.5,0) (.5,0,0) \n 2 b mm2 (.25,.75,z) (.75,.25,-z) \n 2 a mm2 (.25,.25,z) (.75,.75,-z) \n ";
  }
  else {
    sg59 = "Space Group 59\n (0,0,0) \n 8 g 1 (x,y,z) (-x,-y,z) (-x+.5,y+.5,-z) (x+.5,-y+.5,-z) (-x+.5,-y+.5,-z) (x+.5,y+.5,-z) (x,-y,z) (-x,y,z) \n 4 f .m. (x,0,z) (-x,0,z) (-x+.5,.5,-z) (x+.5,.5,-z) \n 4 e m.. (0,y,z) (0,-y,z) (.5,y+.5,-z) (.5,-y+.5,-z) \n 4 d -1 (.25,.25,.5) (.75,.75,.5) (.25,.75,.5) (.75,.25,.5) \n 4 c -1 (.25,.25,0) (.75,.75,0) (.25,.75,0) (.75,.25,0) \n 2 b mm2 (0,.5,z) (.5,0,-z) \n 2 a mm2 (0,0,z) (.5,.5,-z) \n ";
  }

  string sg60 = "Space Group 60\n (0,0,0) \n 8 d 1 (x,y,z) (-x+.5,-y+.5,z+.5) (-x,y,-z+.5) (x+.5,-y+.5,-z) (-x,-y,-z) (x+.5,y+.5,-z+.5) (x,-y,z+.5) (-x+.5,y+.5,z) \n 4 c .2. (0,y,.25) (.5,-y+.5,.75) (0,-y,.75) (.5,y+.5,.25) \n 4 b -2 (0,.5,0) (.5,0,.5) (0,.5,.5) (.5,0,0) \n 4 a -1 (0,0,0) (.5,.5,.5) (0,0,.5) (.5,.5,0) \n ";

  string sg61 = "Space Group 61\n (0,0,0) \n 8 c 1 (x,y,z) (-x+.5,-y,z+.5) (-x,y+.5,-z+.5) (x+.5,-y+.5,-z) (-x,-y,-z) (x+.5,y,-z+.5) (x,-y+.5,z+.5) (-x+.5,y+.5,z) \n 4 b -1 (0,0,.5) (.5,0,0) (0,.5,0) (.5,.5,.5) \n 4 a -1 (0,0,0) (.5,0,.5) (0,.5,.5) (.5,.5,0) \n ";

  //DX20180510 [OBSOLETE] - 8 d 1 not 8 e 1 - string sg62 = "Space Group 62\n (0,0,0) \n 8 e 1 (x,y,z) (-x+.5,-y,z+.5) (-x,y+.5,-z) (x+.5,-y+.5,-z+.5) (-x,-y,-z) (x+.5,y,-z+.5) (x,-y+.5,z) (-x+.5,y+.5,z+.5) \n 4 c .m. (x,.25,z) (-x+.5,.75,z+.5) (-x,.75,-z) (x+.5,.25,-z+.5) \n 4 b -1 (0,0,.5) (.5,0,0) (0,.5,.5) (.5,.5,0) \n 4 a -1 (0,0,0) (.5,0,.5) (0,.5,0) (.5,.5,.5) \n ";
  string sg62 = "Space Group 62\n (0,0,0) \n 8 d 1 (x,y,z) (-x+.5,-y,z+.5) (-x,y+.5,-z) (x+.5,-y+.5,-z+.5) (-x,-y,-z) (x+.5,y,-z+.5) (x,-y+.5,z) (-x+.5,y+.5,z+.5) \n 4 c .m. (x,.25,z) (-x+.5,.75,z+.5) (-x,.75,-z) (x+.5,.25,-z+.5) \n 4 b -1 (0,0,.5) (.5,0,0) (0,.5,.5) (.5,.5,0) \n 4 a -1 (0,0,0) (.5,0,.5) (0,.5,0) (.5,.5,.5) \n ";

  string sg63 = "Space Group 63\n (0,0,0) (.5,.5,0)  \n 16 h 1 (x,y,z) (-x,-y,z+.5) (-x,y,-z+.5) (x,-y,-z) (-x,-y,-z) (x,y,-z+.5) (x,-y,z+.5) (-x,y,z) \n 8 g ..m (x,y,.25) (-x,-y,.75) (-x,y,.25) (x,-y,.75) \n 8 f m.. (0,y,z) (0,-y,z+.5) (0,y,-z+.5) (0,-y,-z) \n 8 e 2.. (x,0,0) (-x,0,.5) (-x,0,0) (x,0,.5) \n 8 d -1 (.25,.25,0) (.75,.75,.5) (.75,.25,.5) (.25,.75,0) \n 4 c m2m (0,y,.25) (0,-y,.75) \n 4 b 2/m.. (0,.5,0) (0,.5,.5) \n 4 a 2/m.. (0,0,0) (0,0,.5) \n ";

  string sg64 = "Space Group 64\n (0,0,0) (.5,.5,0) \n 16 g 1 (x,y,z) (-x,-y+.5,z+.5) (-x,y+.5,-z+.5) (x,-y,-z) (-x,-y,-z) (x,y+.5,-z+.5) (x,-y+.5,z+.5) (-x,y,z) \n 8 f m.. (0,y,z) (0,-y+.5,z+.5) (0,y+.5,-z+.5) (0,-y,-z) \n 8 e .2. (.25,y,.25) (.75,-y+.5,.75) (.75,-y,.75) (.25,y+.5,.25) \n 8 d 2.. (x,0,0) (-x,.5,.5) (-x,0,0) (x,.5,.5) \n 8 c -1 (.25,.25,0) (.75,.25,.5) (.75,.75,.5) (.25,.75,0) \n 4 b 2/m.. (.5,0,0) (.5,.5,.5) \n 4 a 2/m.. (0,0,0) (0,.5,.5) \n ";

  string sg65 = "Space Group 65\n (0,0,0) (.5,.5,0) \n 16 r 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (-x,-y,-z) (x,y,-z) (x,-y,z) (-x,y,z) \n 8 q ..m (x,y,.5) (-x,-y,.5) (-x,y,.5) (x,-y,.5) \n 8 p ..m (x,y,0) (-x,-y,0) (-x,y,0) (x,-y,0) \n 8 o .m. (x,0,z) (-x,0,z) (-x,0,-z) (x,0,-z) \n 8 n m.. (0,y,z) (0,-y,z) (0,y,-z) (0,-y,-z) \n 8 m ..2 (.25,.25,z) (.75,.25,-z) (.75,.75,-z) (.25,.75,z) \n 4 l mm2 (0,.5,z) (0,.5,-z) \n 4 k mm2 (0,0,z) (0,0,-z) \n 4 j m2m (0,y,.5) (0,-y,.5) \n 4 i m2m (0,y,0) (0,-y,0) \n 4 h 2mm (x,0,.5) (-x,0,.5) \n 4 g 2mm (x,0,0) (-x,0,0) \n 4 f ..2/m (.25,.25,.5) (.75,.25,.5) \n 4 e ..2/m (.25,.25,0) (.75,.25,0) \n 2 d mmm (0,0,.5) \n 2 c mmm (.5,0,.5) \n 2 b mmm (.5,0,0) \n 2 a mmm (0,0,0) \n ";

  string sg66 = "Space Group 66\n (0,0,0) (.5,.5,0) \n 16 m 1 (x,y,z) (-x,-y,z) (-x,y,-z+.5) (x,-y,-z+.5) (-x,-y,-z) (x,y,-z) (x,-y,z+.5) (-x,y,z+.5) \n 8 l ..m (x,y,0) (-x,-y,0) (-x,y,.5) (x,-y,.5) \n 8 k ..2 (.25,.25,z) (.75,.25,-z+.5) (.75,.75,-z) (.25,.75,z+.5) \n 8 j ..2 (0,.5,z) (0,.5,-z+.5) (0,.5,-z) (0,.5,z+.5) \n 8 i ..2 (0,0,z) (0,0,-z+.5) (0,0,-z) (0,0,z+.5) \n 8 h .2. (0,y,.25) (0,-y,.25) (0,-y,.75) (0,y,.75) \n 8 g 2.. (x,0,.25) (-x,0,.25) (-x,0,.75) (x,0,.75) \n 4 f ..2/m (.25,.75,0) (.75,.75,.5) \n 4 e ..2/m (.25,.25,0) (.75,.25,.5) \n 4 d ..2/m (0,.5,0) (0,.5,.5) \n 4 c ..2/m (0,0,0) (0,0,.5) \n 4 b 222 (0,.5,.25) (0,.5,.75) \n 4 a 222 (0,0,.25) (0,0,.75) \n ";

  string sg67 = "Space Group 67\n (0,0,0) (.5,.5,0) \n 16 o 1 (x,y,z) (-x,-y+.5,z) (-x,y+.5,-z) (x,-y,-z) (-x,-y,-z) (x,y+.5,-z) (x,-y+.5,z) (-x,y,z) \n 8 n .m. (x,.25,z) (-x,.25,z) (-x,.75,-z) (x,.75,-z) \n 8 m m.. (0,y,z) (0,-y+.5,z) (0,y+.5,-z) (0,-y,-z) \n 8 l ..2 (.25,0,z) (.75,.5,-z) (.75,0,-z) (.25,.5,z) \n 8 k .2. (.25,y,.5) (.75,-y+.5,.5) (.75,-y,.5) (.25,y+.5,.5) \n 8 j .2. (.25,y,0) (.75,-y+.5,0) (.75,-y,0) (.25,y+.5,0) \n 8 i 2.. (x,0,.5) (-x,.5,.5) (-x,0,.5) (x,.5,.5) \n 8 h 2.. (x,0,0) (-x,.5,0) (-x,0,0) (x,.5,0) \n 4 g mm2 (0,.25,z) (0,.75,-z) \n 4 f .2/m. (.25,.25,.5) (.75,.25,.5) \n 4 e .2/m. (.25,.25,0) (.75,.25,0) \n 4 d 2/m.. (0,0,.5) (0,.5,.5) \n 4 c 2/m.. (0,0,0) (0,.5,0) \n 4 b 222 (.25,0,.5) (.75,0,.5) \n 4 a 222 (.25,0,0) (.75,0,0) \n ";

  //DX20180723 - added origin choice 2
  string sg68 = "";
  if(axis_cell==setting_1){
    sg68 = "Space Group 68\n (0,0,0) (.5,.5,0) \n 16 i 1 (x,y,z) (-x+.5,-y+.5,z) (-x,y,-z) (x+.5,-y+.5,-z) (-x,-y+.5,-z+.5) (x+.5,y,-z+.5) (x,-y+.5,z+.5) (-x+.5,y,z+.5) \n 8 h ..2 (.25,.25,z) (.75,.25,-z) (.75,.25,-z+.5) (.25,.25,z+.5) \n 8 g ..2 (0,0,z) (0,0,-z) (0,.5,-z+.5) (0,.5,z+.5) \n 8 f .2. (0,y,0) (.5,-y+.5,0) (0,-y+.5,.5) (.5,y,.5) \n 8 e 2.. (x,0,0) (-x+.5,.5,0) (-x,.5,.5) (x+.5,0,.5) \n 8 d -1 (0,.25,.25) (.5,.25,.25) (0,.25,.75) (.5,.25,.75) \n 8 c -1 (.25,0,.25) (.25,.5,.25) (.75,0,.75) (.75,.5,.75) \n 4 b 222 (0,0,.5) (0,.5,0) \n 4 a 222 (0,0,0) (0,.5,.5) \n ";
  } 
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg68 = "Space Group 68\n (0,0,0) (.5,.5,0) \n 16 i 1 (x,y,z) (-x+.5,-y,z) (-x,y,-z+.5) (x+.5,-y,-z+.5) (-x,-y,-z) (x+.5,y,-z) (x,-y,z+.5) (-x+.5,y,z+.5) \n 8 h ..2 (.25,0,z) (.75,0,-z+.5) (.75,0,-z) (.25,0,z+.5) \n 8 g ..2 (0,.25,z) (0,.25,-z+.5) (0,.75,-z) (0,.75,z+.5) \n 8 f .2. (0,y,.25) (.5,-y,.25) (0,-y,.75) (.5,y,.75) \n 8 e 2.. (x,.25,.25) (-x+.5,.75,.25) (-x,.75,.75) (x+.5,.25,.75) \n 8 d -1 (0,0,0) (.5,0,0) (0,0,.5) (.5,0,.5) \n 8 c -1 (.25,.75,0) (.25,.25,0) (.75,.75,.5) (.75,.25,.5) \n 4 b 222 (0,.25,.75) (0,.75,.25) \n 4 a 222 (0,.25,.25) (0,.75,.75) \n ";
  }
  else {
    sg68 = "Space Group 68\n (0,0,0) (.5,.5,0) \n 16 i 1 (x,y,z) (-x+.5,-y+.5,z) (-x,y,-z) (x+.5,-y+.5,-z) (-x,-y+.5,-z+.5) (x+.5,y,-z+.5) (x,-y+.5,z+.5) (-x+.5,y,z+.5) \n 8 h ..2 (.25,.25,z) (.75,.25,-z) (.75,.25,-z+.5) (.25,.25,z+.5) \n 8 g ..2 (0,0,z) (0,0,-z) (0,.5,-z+.5) (0,.5,z+.5) \n 8 f .2. (0,y,0) (.5,-y+.5,0) (0,-y+.5,.5) (.5,y,.5) \n 8 e 2.. (x,0,0) (-x+.5,.5,0) (-x,.5,.5) (x+.5,0,.5) \n 8 d -1 (0,.25,.25) (.5,.25,.25) (0,.25,.75) (.5,.25,.75) \n 8 c -1 (.25,0,.25) (.25,.5,.25) (.75,0,.75) (.75,.5,.75) \n 4 b 222 (0,0,.5) (0,.5,0) \n 4 a 222 (0,0,0) (0,.5,.5) \n ";
  }

  string sg69 = "Space Group 69\n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 32 p 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (-x,-y,-z) (x,y,-z) (x,-y,z) (-x,y,z) \n 16 o ..m (x,y,0) (-x,-y,0) (-x,y,0) (x,-y,0) \n 16 n .m. (x,0,z) (-x,0,z) (-x,0,-z) (x,0,-z) \n 16 m m.. (0,y,z) (0,-y,z) (0,y,-z) (0,-y,-z) \n 16 l 2.. (x,.25,.25) (-x,.75,.25) (-x,.75,.75) (x,.25,.75) \n 16 k .2. (.25,y,.25) (.75,-y,.25) (.75,-y,.75) (.25,y,.75) \n 16 j ..2 (.25,.25,z) (.75,.25,-z) (.75,.75,-z) (.25,.75,z) \n 8 i mm2 (0,0,z) (0,0,-z) \n 8 h m2m (0,y,0) (0,-y,0) \n 8 g 2mm (x,0,0) (-x,0,0) \n 8 f 222 (.25,.25,.25) (.75,.75,.75) \n 8 e ..2/m (.25,.25,0) (.75,.25,0) \n 8 d .2/m. (.25,0,.25) (.75,0,.25) \n 8 c 2/m.. (0,.25,.25) (0,.75,.25) \n 4 b mmm (0,0,.5) \n 4 a mmm (0,0,0) \n ";

  //DX20180723 - added origin choice 2
  string sg70 = "";
  if(axis_cell==setting_1){
    sg70 = "Space Group 70\n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 32 h 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (-x+.25,-y+.25,-z+.25) (x+.25,y+.25,-z+.25) (x+.25,-y+.25,z+.25) (-x+.25,y+.25,z+.25) \n 16 g ..2 (0,0,z) (0,0,-z) (.25,.25,-z+.25) (.25,.25,z+.25) \n 16 f .2. (0,y,0) (0,-y,0) (.25,-y+.25,.25) (.25,y+.25,.25) \n 16 e 2.. (x,0,0) (-x,0,0) (-x+.25,.25,.25) (x+.25,.25,.25) \n 16 d -1 (.625,.625,.625) (.375,.375,.625) (.375,.625,.375) (.625,.375,.375) \n 16 c -1 (.125,.125,.125) (.875,.875,.125) (.875,.125,.875) (.125,.875,.875) \n 8 b 222 (0,0,.5) (.25,.25,.75) \n 8 a 222 (0,0,0) (.25,.25,.25) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg70 = "Space Group 70\n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 32 h 1 (x,y,z) (-x+.75,-y+.75,z) (-x+.75,y,-z+.75) (x,-y+.75,-z+.75) (-x,-y,-z) (x+.25,y+.25,-z) (x+.25,-y,z+.25) (-x,y+.25,z+.25) \n 16 g ..2 (.125,.125,z) (.625,.125,-z+.75) (.875,.875,-z) (.375,.875,z+.25) \n 16 f .2. (.125,y,.125) (.625,-y+.75,.125) (.875,-y,.875) (.375,y+.25,.875) \n 16 e 2.. (x,.125,.125) (-x+.75,.625,.125) (-x,.875,.875) (x+.25,.375,.875) \n 16 d -1 (.5,.5,.5) (.25,.25,.5) (.25,.5,.25) (.5,.25,.25) \n 16 c -1 (0,0,0) (.75,.75,0) (.75,0,.75) (0,.75,.75) \n 8 b 222 (.125,.125,.625) (.875,.875,.375) \n 8 a 222 (.125,.125,.125) (.875,.875,.875) \n ";
  }
  else {
    sg70 = "Space Group 70\n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 32 h 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (-x+.25,-y+.25,-z+.25) (x+.25,y+.25,-z+.25) (x+.25,-y+.25,z+.25) (-x+.25,y+.25,z+.25) \n 16 g ..2 (0,0,z) (0,0,-z) (.25,.25,-z+.25) (.25,.25,z+.25) \n 16 f .2. (0,y,0) (0,-y,0) (.25,-y+.25,.25) (.25,y+.25,.25) \n 16 e 2.. (x,0,0) (-x,0,0) (-x+.25,.25,.25) (x+.25,.25,.25) \n 16 d -1 (.625,.625,.625) (.375,.375,.625) (.375,.625,.375) (.625,.375,.375) \n 16 c -1 (.125,.125,.125) (.875,.875,.125) (.875,.125,.875) (.125,.875,.875) \n 8 b 222 (0,0,.5) (.25,.25,.75) \n 8 a 222 (0,0,0) (.25,.25,.25) \n ";
  }

  string sg71 = "Space Group 71\n (0,0,0) (.5,.5,.5) \n 16 o 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (-x,-y,-z) (x,y,-z) (x,-y,z) (-x,y,z) \n 8 n ..m (x,y,0) (-x,-y,0) (-x,y,0) (x,-y,0) \n 8 m .m. (x,0,z) (-x,0,z) (-x,0,-z) (x,0,-z) \n 8 l m.. (0,y,z) (0,-y,z) (0,y,-z) (0,-y,-z) \n 8 k -1 (.25,.25,.25) (.75,.75,.25) (.75,.25,.75) (.25,.75,.75) \n 4 j mm2 (.5,0,z) (.5,0,-z) \n 4 i mm2 (0,0,z) (0,0,-z) \n 4 h m2m (0,y,.5) (0,-y,.5) \n 4 g m2m (0,y,0) (0,-y,0) \n 4 f 2mm (x,.5,0) (-x,.5,0) \n 4 e 2mm (x,0,0) (-x,0,0) \n 2 d mmm (.5,0,.5) \n 2 c mmm (.5,.5,0) \n 2 b mmm (0,.5,.5) \n 2 a mmm (0,0,0) \n ";

  string sg72 = "Space Group 72\n (0,0,0) (.5,.5,.5) \n 16 k 1 (x,y,z) (-x,-y,z) (-x+.5,y+.5,-z) (x+.5,-y+.5,-z) (-x,-y,-z) (x,y,-z) (x+.5,-y+.5,z) (-x+.5,y+.5,z) \n 8 j ..m (x,y,0) (-x,-y,0) (-x+.5,y+.5,0) (x+.5,-y+.5,0) \n 8 i ..2 (0,.5,z) (.5,0,-z) (0,.5,-z) (.5,0,z) \n 8 h ..2 (0,0,z) (.5,.5,-z) (0,0,-z) (.5,.5,z) \n 8 g .2. (0,y,.25) (0,-y,.25) (0,-y,.75) (0,y,.75) \n 8 f 2.. (x,0,.25) (-x,0,.25) (-x,0,.75) (x,0,.75) \n 8 e -1 (.25,.25,.25) (.75,.75,.25) (.25,.75,.75) (.75,.25,.75) \n 4 d ..2/m (.5,0,0) (0,.5,0) \n 4 c ..2/m (0,0,0) (.5,.5,0) \n 4 b 222 (.5,0,.25) (.5,0,.75) \n 4 a 222 (0,0,.25) (0,0,.75) \n ";

  string sg73 = "Space Group 73\n (0,0,0) (.5,.5,.5) \n 16 f 1 (x,y,z) (-x+.5,-y,z+.5) (-x,y+.5,-z+.5) (x+.5,-y+.5,-z) (-x,-y,-z) (x+.5,y,-z+.5) (x,-y+.5,z+.5) (-x+.5,y+.5,z) \n 8 e ..2 (0,.25,z) (0,.75,-z+.5) (0,.75,-z) (0,.25,z+.5) \n 8 d .2. (.25,y,0) (.25,-y,.5) (.75,-y,0) (.75,y,.5) \n 8 c 2.. (x,0,.25) (-x+.5,0,.75) (-x,0,.75) (x+.5,0,.25) \n 8 b -1 (.25,.25,.25) (.25,.75,.75) (.75,.75,.25) (.75,.25,.75) \n 8 a -1 (0,0,0) (.5,0,.5) (0,.5,.5) (.5,.5,0) \n ";

  string sg74 = "Space Group 74\n (0,0,0) (.5,.5,.5) \n 16 j 1 (x,y,z) (-x,-y+.5,z) (-x,y+.5,-z) (x,-y,-z) (-x,-y,-z) (x,y+.5,-z) (x,-y+.5,z) (-x,y,z) \n 8 i .m. (x,.25,z) (-x,.25,z) (-x,.75,-z) (x,.75,-z) \n 8 h m.. (0,y,z) (0,-y+.5,z) (0,y+.5,-z) (0,-y,-z) \n 8 g .2. (.25,y,.25) (.75,-y+.5,.25) (.75,-y,.75) (.25,y+.5,.75) \n 8 f 2.. (x,0,0) (-x,.5,0) (-x,0,0) (x,.5,0) \n 4 e mm2 (0,.25,z) (0,.75,-z) \n 4 d .2/m. (.25,.25,.75) (.75,.25,.75) \n 4 c .2/m. (.25,.25,.25) (.75,.25,.25) \n 4 b 2/m.. (0,0,.5) (0,.5,.5) \n 4 a 2/m.. (0,0,0) (0,.5,0) \n ";

  string sg75 = "Space Group 75\n (0,0,0) \n 4 d 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) \n 2 c 2.. (0,.5,z) (.5,0,z) \n 1 b 4.. (.5,.5,z) \n 1 a 4.. (0,0,z) \n ";

  string sg76 = "Space Group 76\n (0,0,0) \n 4 a 1 (x,y,z) (-x,-y,z+.5) (-y,x,z+.25) (y,-x,z+.75) \n ";

  string sg77 = "Space Group 77\n (0,0,0) \n 4 d 1 (x,y,z) (-x,-y,z) (-y,x,z+.5) (y,-x,z+.5) \n 2 c 2.. (0,.5,z) (.5,0,z+.5) \n 2 b 2.. (.5,.5,z) (.5,.5,z+.5) \n 2 a 2.. (0,0,z) (0,0,z+.5) \n ";

  string sg78 = "Space Group 78\n (0,0,0) \n 4 a 1 (x,y,z) (-x,-y,z+.5) (-y,x,z+.75) (y,-x,z+.25) \n ";

  string sg79 = "Space Group 79\n (0,0,0) (.5,.5,.5) \n 8 c 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) \n 4 b 2.. (0,.5,z) (.5,0,z) \n 2 a 4.. (0,0,z) \n ";

  string sg80 = "Space Group 80\n (0,0,0) (.5,.5,.5) \n 8 b 1 (x,y,z) (-x+.5,-y+.5,z+.5) (-y,x+.5,z+.25) (y+.5,-x,z+.75) \n 4 a 2.. (0,0,z) (0,.5,z+.25) \n ";

  string sg81 = "Space Group 81\n (0,0,0) \n 4 h 1 (x,y,z) (-x,-y,z) (y,-x,-z) (-y,x,-z) \n 2 g 2.. (0,.5,z) (.5,0,-z) \n 2 f 2.. (.5,.5,z) (.5,.5,-z) \n 2 e 2.. (0,0,z) (0,0,-z) \n 1 d -4.. (.5,.5,.5) \n 1 c -4.. (.5,.5,0) \n 1 b -4.. (0,0,.5) \n 1 a -4.. (0,0,0) \n ";

  string sg82 = "Space Group 82\n (0,0,0) (.5,.5,.5) \n 8 g 1 (x,y,z) (-x,-y,z) (y,-x,-z) (-y,x,-z) \n 4 f 2.. (0,.5,z) (.5,0,-z) \n 4 e 2.. (0,0,z) (0,0,-z) \n 2 d -4.. (0,.5,.75) \n 2 c -4.. (0,.5,.25) \n 2 b -4.. (0,0,.5) \n 2 a -4.. (0,0,0) \n ";

  string sg83 = "Space Group 83\n (0,0,0) \n 8 l 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (-x,-y,-z) (x,y,-z) (y,-x,-z) (-y,x,-z) \n 4 k m.. (x,y,.5) (-x,-y,.5) (-y,x,.5) (y,-x,.5) \n 4 j m.. (x,y,0) (-x,-y,0) (-y,x,0) (y,-x,0) \n 4 i 2.. (0,.5,z) (.5,0,z) (0,.5,-z) (.5,0,-z) \n 2 h 4.. (.5,.5,z) (.5,.5,-z) \n 2 g 4.. (0,0,z) (0,0,-z) \n 2 f 2/m.. (0,.5,.5) (.5,0,.5) \n 2 e 2/m.. (0,.5,0) (.5,0,0) \n 1 d 4/m.. (.5,.5,.5) \n 1 c 4/m.. (.5,.5,0) \n 1 b 4/m.. (0,0,.5) \n 1 a 4/m.. (0,0,0) \n ";

  string sg84 = "Space Group 84\n (0,0,0) \n 8 k 1 (x,y,z) (-x,-y,z) (-y,x,z+.5) (y,-x,z+.5) (-x,-y,-z) (x,y,-z) (y,-x,-z+.5) (-y,x,-z+.5) \n 4 j m.. (x,y,0) (-x,-y,0) (-y,x,.5) (y,-x,.5) \n 4 i 2.. (0,.5,z) (.5,0,z+.5) (0,.5,-z) (.5,0,-z+.5) \n 4 h 2.. (.5,.5,z) (.5,.5,z+.5) (.5,.5,-z) (.5,.5,-z+.5) \n 4 g 2.. (0,0,z) (0,0,z+.5) (0,0,-z) (0,0,-z+.5) \n 2 f -4.. (.5,.5,.25) (.5,.5,.75) \n 2 e -4.. (0,0,.25) (0,0,.75) \n 2 d 2/m.. (0,.5,.5) (.5,0,0) \n 2 c 2/m.. (0,.5,0) (.5,0,.5) \n 2 b 2/m.. (.5,.5,0) (.5,.5,.5) \n 2 a 2/m.. (0,0,0) (0,0,.5) \n ";

  //DX20180723 - added origin choice 2
  string sg85 = "";
  if(axis_cell==setting_1){
    sg85 = "Space Group 85\n (0,0,0) \n 8 g 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z) (y+.5,-x+.5,z) (-x+.5,-y+.5,-z) (x+.5,y+.5,-z) (y,-x,-z) (-y,x,-z) \n 4 f 2.. (0,0,z) (.5,.5,z) (.5,.5,-z) (0,0,-z) \n 4 e -1 (.25,.25,.5) (.75,.75,.5) (.25,.75,.5) (.75,.25,.5) \n 4 d -1 (.25,.25,0) (.75,.75,0) (.25,.75,0) (.75,.25,0) \n 2 c 4.. (0,.5,z) (.5,0,-z) \n 2 b -4.. (0,0,.5) (.5,.5,.5) \n 2 a -4.. (0,0,0) (.5,.5,0) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg85 = "Space Group 85\n (0,0,0) \n 8 g 1 (x,y,z) (-x+.5,-y+.5,z) (-y+.5,x,z) (y,-x+.5,z) (-x,-y,-z) (x+.5,y+.5,-z) (y+.5,-x,-z) (-y,x+.5,-z) \n 4 f 2.. (.25,.75,z) (.75,.25,z) (.75,.25,-z) (.25,.75,-z) \n 4 e -1 (0,0,.5) (.5,.5,.5) (.5,0,.5) (0,.5,.5) \n 4 d -1 (0,0,0) (.5,.5,0) (.5,0,0) (0,.5,0) \n 2 c 4.. (.25,.25,z) (.75,.75,-z) \n 2 b -4.. (.25,.75,.5) (.75,.25,.5) \n 2 a -4.. (.25,.75,0) (.75,.25,0) \n "; //DX20210204 - fixed 4 f Wyckoff position (.25,75,-z) should be (.25,.75,-z)
  }
  else {
    sg85 = "Space Group 85\n (0,0,0) \n 8 g 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z) (y+.5,-x+.5,z) (-x+.5,-y+.5,-z) (x+.5,y+.5,-z) (y,-x,-z) (-y,x,-z) \n 4 f 2.. (0,0,z) (.5,.5,z) (.5,.5,-z) (0,0,-z) \n 4 e -1 (.25,.25,.5) (.75,.75,.5) (.25,.75,.5) (.75,.25,.5) \n 4 d -1 (.25,.25,0) (.75,.75,0) (.25,.75,0) (.75,.25,0) \n 2 c 4.. (0,.5,z) (.5,0,-z) \n 2 b -4.. (0,0,.5) (.5,.5,.5) \n 2 a -4.. (0,0,0) (.5,.5,0) \n ";
  }

  //DX20180723 - added origin choice 2
  string sg86 = "";
  if(axis_cell==setting_1){
    sg86 = "Space Group 86\n (0,0,0) \n 8 g 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z+.5) (y+.5,-x+.5,z+.5) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (y,-x,-z) (-y,x,-z) \n 4 f 2.. (0,0,z) (.5,.5,z+.5) (.5,.5,-z+.5) (0,0,-z) \n 4 e 2.. (0,.5,z) (0,.5,z+.5) (.5,0,-z+.5) (.5,0,-z) \n 4 d -1 (.25,.25,.75) (.75,.75,.75) (.25,.75,.25) (.75,.25,.25) \n 4 c -1 (.25,.25,.25) (.75,.75,.25) (.25,.75,.75) (.75,.25,.75) \n 2 b -4.. (0,0,.5) (.5,.5,0) \n 2 a -4.. (0,0,0) (.5,.5,.5) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg86 = "Space Group 86\n (0,0,0) \n 8 g 1 (x,y,z) (-x+.5,-y+.5,z) (-y,x+.5,z+.5) (y+.5,-x,z+.5) (-x,-y,-z) (x+.5,y+.5,-z) (y,-x+.5,-z+.5) (-y+.5,x,-z+.5) \n 4 f 2.. (.25,.25,z) (.75,.75,z+.5) (.75,.75,-z) (.25,.25,-z+.5) \n 4 e 2.. (.75,.25,z) (.75,.25,z+.5) (.25,.75,-z) (.25,.75,-z+.5) \n 4 d -1 (0,0,.5) (.5,.5,.5) (0,.5,0) (.5,0,0) \n 4 c -1 (0,0,0) (.5,.5,0) (0,.5,.5) (.5,0,.5) \n 2 b -4.. (.25,.25,.75) (.75,.75,.25) \n 2 a -4.. (.25,.25,.25) (.75,.75,.75) \n ";
  }
  else {
    sg86 = "Space Group 86\n (0,0,0) \n 8 g 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z+.5) (y+.5,-x+.5,z+.5) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (y,-x,-z) (-y,x,-z) \n 4 f 2.. (0,0,z) (.5,.5,z+.5) (.5,.5,-z+.5) (0,0,-z) \n 4 e 2.. (0,.5,z) (0,.5,z+.5) (.5,0,-z+.5) (.5,0,-z) \n 4 d -1 (.25,.25,.75) (.75,.75,.75) (.25,.75,.25) (.75,.25,.25) \n 4 c -1 (.25,.25,.25) (.75,.75,.25) (.25,.75,.75) (.75,.25,.75) \n 2 b -4.. (0,0,.5) (.5,.5,0) \n 2 a -4.. (0,0,0) (.5,.5,.5) \n ";
  }

  string sg87 = "Space Group 87\n (0,0,0) (.5,.5,.5) \n 16 i 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (-x,-y,-z) (x,y,-z) (y,-x,-z) (-y,x,-z) \n 8 h m.. (x,y,0) (-x,-y,0) (-y,x,0) (y,-x,0) \n 8 g 2.. (0,.5,z) (.5,0,z) (0,.5,-z) (.5,0,-z) \n 8 f -1 (.25,.25,.25) (.75,.75,.25) (.75,.25,.25) (.25,.75,.25) \n 4 e 4.. (0,0,z) (0,0,-z) \n 4 d -4.. (0,.5,.25) (.5,0,.25) \n 4 c 2/m.. (0,.5,0) (.5,0,0) \n 2 b 4/m.. (0,0,.5) \n 2 a 4/m.. (0,0,0) \n ";

  //DX20180723 - added origin choice 2
  string sg88 = "";
  if(axis_cell==setting_1){
    sg88 = "Space Group 88\n (0,0,0) (.5,.5,.5) \n 16 f 1 (x,y,z) (-x+.5,-y+.5,z+.5) (-y,x+.5,z+.25) (y+.5,-x,z+.75) (-x,-y+.5,-z+.25) (x+.5,y,-z+.75) (y,-x,-z) (-y+.5,x+.5,-z+.5) \n 8 e 2.. (0,0,z) (0,.5,z+.25) (0,.5,-z+.25) (0,0,-z) \n 8 d -1 (0,.25,.625) (.5,.25,.125) (.75,.5,.875) (.75,0,.375) \n 8 c -1 (0,.25,.125) (.5,.25,.625) (.75,.5,.375) (.75,0,.875) \n 4 b -4.. (0,0,.5) (0,.5,.75) \n 4 a -4.. (0,0,0) (0,.5,.25) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg88 = "Space Group 88\n (0,0,0) (.5,.5,.5) \n 16 f 1 (x,y,z) (-x+.5,-y,z+.5) (-y+.75,x+.25,z+.25) (y+.75,-x+.75,z+.75) (-x,-y,-z) (x+.5,y,-z+.5) (y+.25,-x+.75,-z+.75) (-y+.25,x+.25,-z+.25) \n 8 e 2.. (0,.25,z) (.5,.25,z+.25) (0,.75,-z) (.5,.75,-z+.75) \n 8 d -1 (0,0,.5) (.5,0,0) (.75,.25,.75) (.75,.75,.25) \n 8 c -1 (0,0,0) (.5,0,.5) (.75,.25,.25) (.75,.75,.75) \n 4 b -4.. (0,.25,.625) (.5,.25,.875) \n 4 a -4.. (0,.25,.125) (.5,.25,.375) \n "; //DX20180806 first position for 8 d is (0,0,.5) not (0,0,5)
  }
  else {
    sg88 = "Space Group 88\n (0,0,0) (.5,.5,.5) \n 16 f 1 (x,y,z) (-x+.5,-y+.5,z+.5) (-y,x+.5,z+.25) (y+.5,-x,z+.75) (-x,-y+.5,-z+.25) (x+.5,y,-z+.75) (y,-x,-z) (-y+.5,x+.5,-z+.5) \n 8 e 2.. (0,0,z) (0,.5,z+.25) (0,.5,-z+.25) (0,0,-z) \n 8 d -1 (0,.25,.625) (.5,.25,.125) (.75,.5,.875) (.75,0,.375) \n 8 c -1 (0,.25,.125) (.5,.25,.625) (.75,.5,.375) (.75,0,.875) \n 4 b -4.. (0,0,.5) (0,.5,.75) \n 4 a -4.. (0,0,0) (0,.5,.25) \n ";
  }

  string sg89 =
    "Space Group 89\n \
    (0,0,0) \n 8 p 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (-x,y,-z) (x,-y,-z) (y,x,-z) (-y,-x,-z) \n 4 o .2. (x,.5,0) (-x,.5,0) (.5,x,0) (.5,-x,0) \n 4 n .2. (x,0,.5) (-x,0,.5) (0,x,.5) (0,-x,.5) \n 4 m .2. (x,.5,.5) (-x,.5,.5) (.5,x,.5) (.5,-x,.5) \n 4 l .2. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) \n 4 k ..2 (x,x,.5) (-x,-x,.5) (-x,x,.5) (x,-x,.5) \n 4 j ..2 (x,x,0) (-x,-x,0) (-x,x,0) (x,-x,0) \n 4 i 2.. (0,.5,z) (.5,0,z) (0,.5,-z) (.5,0,-z) \n 2 h 4.. (.5,.5,z) (.5,.5,-z) \n 2 g 4.. (0,0,z) (0,0,-z) \n 2 f 222. (.5,0,.5) (0,.5,.5) \n 2 e 222. (.5,0,0) (0,.5,0) \n 1 d 422 (.5,.5,.5) \n 1 c 422 (.5,.5,0) \n 1 b 422 (0,0,.5) \n 1 a 422 (0,0,0) \n ";

  string sg90 = "Space Group 90\n (0,0,0) \n 8 g 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z) (y+.5,-x+.5,z) (-x+.5,y+.5,-z) (x+.5,-y+.5,-z) (y,x,-z) (-y,-x,-z) \n 4 f ..2 (x,x,.5) (-x,-x,.5) (-x+.5,x+.5,.5) (x+.5,-x+.5,.5) \n 4 e ..2 (x,x,0) (-x,-x,0) (-x+.5,x+.5,0) (x+.5,-x+.5,0) \n 4 d 2.. (0,0,z) (.5,.5,z) (.5,.5,-z) (0,0,-z) \n 2 c 4.. (0,.5,z) (.5,0,-z) \n 2 b 2.22 (0,0,.5) (.5,.5,.5) \n 2 a 2.22 (0,0,0) (.5,.5,0) \n ";

  string sg91 = "Space Group 91\n (0,0,0) \n 8 d 1 (x,y,z) (-x,-y,z+.5) (-y,x,z+.25) (y,-x,z+.75) (-x,y,-z) (x,-y,-z+.5) (y,x,-z+.75) (-y,-x,-z+.25) \n 4 c ..2 (x,x,.375) (-x,-x,.875) (-x,x,.625) (x,-x,.125) \n 4 b .2. (.5,y,0) (.5,-y,.5) (-y,.5,.25) (y,.5,.75) \n 4 a .2. (0,y,0) (0,-y,.5) (-y,0,.25) (y,0,.75) \n ";

  string sg92 = "Space Group 92\n (0,0,0) \n 8 b 1 (x,y,z) (-x,-y,z+.5) (-y+.5,x+.5,z+.25) (y+.5,-x+.5,z+.75) (-x+.5,y+.5,-z+.25) (x+.5,-y+.5,-z+.75) (y,x,-z) (-y,-x,-z+.5) \n 4 a ..2 (x,x,0) (-x,-x,.5) (-x+.5,x+.5,.25) (x+.5,-x+.5,.75) \n ";

  string sg93 = "Space Group 93\n (0,0,0) \n 8 p 1 (x,y,z) (-x,-y,z) (-y,x,z+.5) (y,-x,z+.5) (-x,y,-z) (x,-y,-z) (y,x,-z+.5) (-y,-x,-z+.5) \n 4 o ..2 (x,x,.75) (-x,-x,.75) (-x,x,.25) (x,-x,.25) \n 4 n ..2 (x,x,.25) (-x,-x,.25) (-x,x,.75) (x,-x,.75) \n 4 m .2. (x,.5,0) (-x,.5,0) (.5,x,.5) (.5,-x,.5) \n 4 l .2. (x,0,.5) (-x,0,.5) (0,x,0) (0,-x,0) \n 4 k .2. (x,.5,.5) (-x,.5,.5) (.5,x,0) (.5,-x,0) \n 4 j .2. (x,0,0) (-x,0,0) (0,x,.5) (0,-x,.5) \n 4 i 2.. (0,.5,z) (.5,0,z+.5) (0,.5,-z) (.5,0,-z+.5) \n 4 h 2.. (.5,.5,z) (.5,.5,z+.5) (.5,.5,-z) (.5,.5,-z+.5) \n 4 g 2.. (0,0,z) (0,0,z+.5) (0,0,-z) (0,0,-z+.5) \n 2 f 2.22 (.5,.5,.25) (.5,.5,.75) \n 2 e 2.22 (0,0,.25) (0,0,.75) \n 2 d 222. (0,.5,.5) (.5,0,0) \n 2 c 222. (0,.5,0) (.5,0,.5) \n 2 b 222. (.5,.5,0) (.5,.5,.5) \n 2 a 222. (0,0,0) (0,0,.5) \n ";

  string sg94 = "Space Group 94\n (0,0,0) \n 8 g 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z+.5) (y+.5,-x+.5,z+.5) (-x+.5,y+.5,-z+.5) (x+.5,-y+.5,-z+.5) (y,x,-z) (-y,-x,-z) \n 4 f ..2 (x,x,.5) (-x,-x,.5) (-x+.5,x+.5,0) (x+.5,-x+.5,0) \n 4 e ..2 (x,x,0) (-x,-x,0) (-x+.5,x+.5,.5) (x+.5,-x+.5,.5) \n 4 d 2.. (0,.5,z) (0,.5,z+.5) (.5,0,-z+.5) (.5,0,-z) \n 4 c 2.. (0,0,z) (.5,.5,z+.5) (.5,.5,-z+.5) (0,0,-z) \n 2 b 2.22 (0,0,.5) (.5,.5,0) \n 2 a 2.22 (0,0,0) (.5,.5,.5) \n ";

  string sg95 = "Space Group 95\n (0,0,0) \n 8 d 1 (x,y,z) (-x,-y,z+.5) (-y,x,z+.75) (y,-x,z+.25) (-x,y,-z) (x,-y,-z+.5) (y,x,-z+.25) (-y,-x,-z+.75) \n 4 c ..2 (x,x,.625) (-x,-x,.125) (-x,x,.375) (x,-x,.875) \n 4 b .2. (.5,y,0) (.5,-y,.5) (-y,.5,.75) (y,.5,.25) \n 4 a .2. (0,y,0) (0,-y,.5) (-y,0,.75) (y,0,.25) \n ";

  string sg96 = "Space Group 96\n (0,0,0) \n 8 b 1 (x,y,z) (-x,-y,z+.5) (-y+.5,x+.5,z+.75) (y+.5,-x+.5,z+.25) (-x+.5,y+.5,-z+.75) (x+.5,-y+.5,-z+.25) (y,x,-z) (-y,-x,-z+.5) \n 4 a ..2 (x,x,0) (-x,-x,.5) (-x+.5,x+.5,.75) (x+.5,-x+.5,.25) \n ";

  string sg97 = "Space Group 97\n (0,0,0) (.5,.5,.5) \n 16 k 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (-x,y,-z) (x,-y,-z) (y,x,-z) (-y,-x,-z) \n 8 j ..2 (x,x+.5,.25) (-x,-x+.5,.25) (-x+.5,x,.25) (x+.5,-x,.25) \n 8 i .2. (x,0,.5) (-x,0,.5) (0,x,.5) (0,-x,.5) \n 8 h .2. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) \n 8 g ..2 (x,x,0) (-x,-x,0) (-x,x,0) (x,-x,0) \n 8 f 2.. (0,.5,z) (.5,0,z) (0,.5,-z) (.5,0,-z) \n 4 e 4.. (0,0,z) (0,0,-z) \n 4 d 2.22 (0,.5,.25) (.5,0,.25) \n 4 c 222. (0,.5,0) (.5,0,0) \n 2 b 422 (0,0,.5) \n 2 a 422 (0,0,0) \n ";

  string sg98 = "Space Group 98\n (0,0,0) (.5,.5,.5) \n 16 g 1 (x,y,z) (-x+.5,-y+.5,z+.5) (-y,x+.5,z+.25) (y+.5,-x,z+.75) (-x+.5,y,-z+.75) (x,-y+.5,-z+.25) (y+.5,x+.5,-z+.5) (-y,-x,-z) \n 8 f .2. (x,.25,.125) (-x+.5,.25,.625) (.75,x+.5,.375) (.75,-x,.875) \n 8 e ..2 (-x,x,0) (x+.5,-x+.5,.5) (-x,-x+.5,.25) (x+.5,x,.75) \n 8 d ..2 (x,x,0) (-x+.5,-x+.5,.5) (-x,x+.5,.25) (x+.5,-x,.75) \n 8 c 2.. (0,0,z) (0,.5,z+.25) (.5,0,-z+.75) (.5,.5,-z+.5) \n 4 b 2.22 (0,0,.5) (0,.5,.75) \n 4 a 2.22 (0,0,0) (0,.5,.25) \n ";

  string sg99 = "Space Group 99\n (0,0,0) \n 8 g 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (x,-y,z) (-x,y,z) (-y,-x,z) (y,x,z) \n 4 f .m. (x,.5,z) (-x,.5,z) (.5,x,z) (.5,-x,z) \n 4 e .m. (x,0,z) (-x,0,z) (0,x,z) (0,-x,z) \n 4 d ..m (x,x,z) (-x,-x,z) (-x,x,z) (x,-x,z) \n 2 c 2mm. (.5,0,z) (0,.5,z) \n 1 b 4mm (.5,.5,z) \n 1 a 4mm (0,0,z) \n ";

  string sg100 = "Space Group 100\n (0,0,0) \n 8 d 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (x+.5,-y+.5,z) (-x+.5,y+.5,z) (-y+.5,-x+.5,z) (y+.5,x+.5,z) \n 4 c ..m (x,x+.5,z) (-x,-x+.5,z) (-x+.5,x,z) (x+.5,-x,z) \n 2 b 2.mm (.5,0,z) (0,.5,z) \n 2 a 4.. (0,0,z) (.5,.5,z) \n ";

  string sg101 = "Space Group 101\n (0,0,0) \n 8 e 1 (x,y,z) (-x,-y,z) (-y,x,z+.5) (y,-x,z+.5) (x,-y,z+.5) (-x,y,z+.5) (-y,-x,z) (y,x,z) \n 4 d ..m (x,x,z) (-x,-x,z) (-x,x,z+.5) (x,-x,z+.5) \n 4 c 2.. (0,.5,z) (.5,0,z+.5) (0,.5,z+.5) (.5,0,z) \n 2 b 2.mm (.5,.5,z) (.5,.5,z+.5) \n 2 a 2.mm (0,0,z) (0,0,z+.5) \n ";

  string sg102 = "Space Group 102\n (0,0,0) \n 8 d 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z+.5) (y+.5,-x+.5,z+.5) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) (-y,-x,z) (y,x,z) \n 4 c ..m (x,x,z) (-x,-x,z) (-x+.5,x+.5,z+.5) (x+.5,-x+.5,z+.5) \n 4 b 2.. (0,.5,z) (0,.5,z+.5) (.5,0,z+.5) (.5,0,z) \n 2 a 2.mm (0,0,z) (.5,.5,z+.5) \n ";

  string sg103 = "Space Group 103\n (0,0,0) \n 8 d 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (x,-y,z+.5) (-x,y,z+.5) (-y,-x,z+.5) (y,x,z+.5) \n 4 c 2.. (0,.5,z) (.5,0,z) (0,.5,z+.5) (.5,0,z+.5) \n 2 b 4.. (.5,.5,z) (.5,.5,z+.5) \n 2 a 4.. (0,0,z) (0,0,z+.5) \n ";

  string sg104 = "Space Group 104\n (0,0,0) \n 8 c 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) (-y+.5,-x+.5,z+.5) (y+.5,x+.5,z+.5) \n 4 b 2.. (0,.5,z) (.5,0,z) (.5,0,z+.5) (0,.5,z+.5) \n 2 a 4.. (0,0,z) (.5,.5,z+.5) \n ";

  string sg105 = "Space Group 105\n (0,0,0) \n 8 f 1 (x,y,z) (-x,-y,z) (-y,x,z+.5) (y,-x,z+.5) (x,-y,z) (-x,y,z) (-y,-x,z+.5) (y,x,z+.5) \n 4 e .m. (x,.5,z) (-x,.5,z) (.5,x,z+.5) (.5,-x,z+.5) \n 4 d .m. (x,0,z) (-x,0,z) (0,x,z+.5) (0,-x,z+.5) \n 2 c 2mm. (0,.5,z) (.5,0,z+.5) \n 2 b 2mm. (.5,.5,z) (.5,.5,z+.5) \n 2 a 2mm. (0,0,z) (0,0,z+.5) \n ";

  string sg106 = "Space Group 106\n (0,0,0) \n 8 c 1 (x,y,z) (-x,-y,z) (-y,x,z+.5) (y,-x,z+.5) (x+.5,-y+.5,z) (-x+.5,y+.5,z) (-y+.5,-x+.5,z+.5) (y+.5,x+.5,z+.5) \n 4 b 2.. (0,.5,z) (.5,0,z+.5) (.5,0,z) (0,.5,z+.5) \n 4 a 2.. (0,0,z) (0,0,z+.5) (.5,.5,z) (.5,.5,z+.5) \n ";

  string sg107 = "Space Group 107\n (0,0,0) (.5,.5,.5) \n 16 e 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (x,-y,z) (-x,y,z) (-y,-x,z) (y,x,z) \n 8 d .m. (x,0,z) (-x,0,z) (0,x,z) (0,-x,z) \n 8 c ..m (x,x,z) (-x,-x,z) (-x,x,z) (x,-x,z) \n 4 b 2mm. (0,.5,z) (.5,0,z) \n 2 a 4mm (0,0,z) \n ";

  string sg108 = "Space Group 108\n (0,0,0) (.5,.5,.5) \n 16 d 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (x,-y,z+.5) (-x,y,z+.5) (-y,-x,z+.5) (y,x,z+.5) \n 8 c ..m (x,x+.5,z) (-x,-x+.5,z) (-x+.5,x,z) (x+.5,-x,z) \n 4 b 2.mm (.5,0,z) (0,.5,z) \n 4 a 4.. (0,0,z) (0,0,z+.5) \n ";

  string sg109 = "Space Group 109\n (0,0,0) (.5,.5,.5) \n 16 c 1 (x,y,z) (-x+.5,-y+.5,z+.5) (-y,x+.5,z+.25) (y+.5,-x,z+.75) (x,-y,z) (-x+.5,y+.5,z+.5) (-y,-x+.5,z+.25) (y+.5,x,z+.75) \n 8 b .m. (0,y,z) (.5,-y+.5,z+.5) (-y,.5,z+.25) (y+.5,0,z+.75) \n 4 a 2mm. (0,0,z) (0,.5,z+.25) \n ";

  string sg110 = "Space Group 110\n (0,0,0) (.5,.5,.5) \n 16 b 1 (x,y,z) (-x+.5,-y+.5,z+.5) (-y,x+.5,z+.25) (y+.5,-x,z+.75) (x,-y,z+.5) (-x+.5,y+.5,z) (-y,-x+.5,z+.75) (y+.5,x,z+.25) \n 8 a 2.. (0,0,z) (0,.5,z+.25) (0,0,z+.5) (0,.5,z+.75) \n ";

  string sg111 = "Space Group 111\n (0,0,0) \n 8 o 1 (x,y,z) (-x,-y,z) (y,-x,-z) (-y,x,-z) (-x,y,-z) (x,-y,-z) (-y,-x,z) (y,x,z) \n 4 n ..m (x,x,z) (-x,-x,z) (x,-x,-z) (-x,x,-z) \n 4 m 2.. (0,.5,z) (.5,0,-z) (0,.5,-z) (.5,0,z) \n 4 l .2. (x,.5,0) (-x,.5,0) (.5,-x,0) (.5,x,0) \n 4 k .2. (x,0,.5) (-x,0,.5) (0,-x,.5) (0,x,.5) \n 4 j .2. (x,.5,.5) (-x,.5,.5) (.5,-x,.5) (.5,x,.5) \n 4 i .2. (x,0,0) (-x,0,0) (0,-x,0) (0,x,0) \n 2 h 2.mm (.5,.5,z) (.5,.5,-z) \n 2 g 2.mm (0,0,z) (0,0,-z) \n 2 f 222. (.5,0,.5) (0,.5,.5) \n 2 e 222. (.5,0,0) (0,.5,0) \n 1 d -42m (.5,.5,0) \n 1 c -42m (0,0,.5) \n 1 b -42m (.5,.5,.5) \n 1 a -42m (0,0,0) \n ";

  string sg112 = "Space Group 112\n (0,0,0) \n 8 n 1 (x,y,z) (-x,-y,z) (y,-x,-z) (-y,x,-z) (-x,y,-z+.5) (x,-y,-z+.5) (-y,-x,z+.5) (y,x,z+.5) \n 4 m 2.. (0,.5,z) (.5,0,-z) (0,.5,-z+.5) (.5,0,z+.5) \n 4 l 2.. (.5,.5,z) (.5,.5,-z) (.5,.5,-z+.5) (.5,.5,z+.5) \n 4 k 2.. (0,0,z) (0,0,-z) (0,0,-z+.5) (0,0,z+.5) \n 4 j .2. (0,y,.25) (0,-y,.25) (y,0,.75) (-y,0,.75) \n 4 i .2. (x,.5,.25) (-x,.5,.25) (.5,-x,.75) (.5,x,.75) \n 4 h .2. (.5,y,.25) (.5,-y,.25) (y,.5,.75) (-y,.5,.75) \n 4 g .2. (x,0,.25) (-x,0,.25) (0,-x,.75) (0,x,.75) \n 2 f -4.. (.5,.5,0) (.5,.5,.5) \n 2 e -4.. (0,0,0) (0,0,.5) \n 2 d 222. (0,.5,.25) (.5,0,.75) \n 2 c 222. (.5,.5,.25) (.5,.5,.75) \n 2 b 222. (.5,0,.25) (0,.5,.75) \n 2 a 222. (0,0,.25) (0,0,.75) \n ";

  string sg113 = "Space Group 113\n (0,0,0) \n 8 f 1 (x,y,z) (-x,-y,z) (y,-x,-z) (-y,x,-z) (-x+.5,y+.5,-z) (x+.5,-y+.5,-z) (-y+.5,-x+.5,z) (y+.5,x+.5,z) \n 4 e ..m (x,x+.5,z) (-x,-x+.5,z) (x+.5,-x,-z) (-x+.5,x,-z) \n 4 d 2.. (0,0,z) (0,0,-z) (.5,.5,-z) (.5,.5,z) \n 2 c 2.mm (0,.5,z) (.5,0,-z) \n 2 b -4.. (0,0,.5) (.5,.5,.5) \n 2 a -4.. (0,0,0) (.5,.5,0) \n ";

  string sg114 = "Space Group 114\n (0,0,0) \n 8 e 1 (x,y,z) (-x,-y,z) (y,-x,-z) (-y,x,-z) (-x+.5,y+.5,-z+.5) (x+.5,-y+.5,-z+.5) (-y+.5,-x+.5,z+.5) (y+.5,x+.5,z+.5) \n 4 d 2.. (0,.5,z) (.5,0,-z) (.5,0,-z+.5) (0,.5,z+.5) \n 4 c 2.. (0,0,z) (0,0,-z) (.5,.5,-z+.5) (.5,.5,z+.5) \n 2 b -4.. (0,0,.5) (.5,.5,0) \n 2 a -4.. (0,0,0) (.5,.5,.5) \n ";

  string sg115 = "Space Group 115\n (0,0,0) \n 8 l 1 (x,y,z) (-x,-y,z) (y,-x,-z) (-y,x,-z) (x,-y,z) (-x,y,z) (y,x,-z) (-y,-x,-z) \n 4 k .m. (x,.5,z) (-x,.5,z) (.5,-x,-z) (.5,x,-z) \n 4 j .m. (x,0,z) (-x,0,z) (0,-x,-z) (0,x,-z) \n 4 i ..2 (x,x,.5) (-x,-x,.5) (x,-x,.5) (-x,x,.5) \n 4 h ..2 (x,x,0) (-x,-x,0) (x,-x,0) (-x,x,0) \n 2 g 2mm. (0,.5,z) (.5,0,-z) \n 2 f 2mm. (.5,.5,z) (.5,.5,-z) \n 2 e 2mm. (0,0,z) (0,0,-z) \n 1 d -4m2 (0,0,.5) \n 1 c -4m2 (.5,.5,.5) \n 1 b -4m2 (.5,.5,0) \n 1 a -4m2 (0,0,0) \n ";

  string sg116 = "Space Group 116\n (0,0,0) \n 8 j 1 (x,y,z) (-x,-y,z) (y,-x,-z) (-y,x,-z) (x,-y,z+.5) (-x,y,z+.5) (y,x,-z+.5) (-y,-x,-z+.5) \n 4 i 2.. (0,.5,z) (.5,0,-z) (0,.5,z+.5) (.5,0,-z+.5) \n 4 h 2.. (.5,.5,z) (.5,.5,-z) (.5,.5,z+.5) (.5,.5,-z+.5) \n 4 g 2.. (0,0,z) (0,0,-z) (0,0,z+.5) (0,0,-z+.5) \n 4 f ..2 (x,x,.75) (-x,-x,.75) (x,-x,.25) (-x,x,.25) \n 4 e ..2 (x,x,.25) (-x,-x,.25) (x,-x,.75) (-x,x,.75) \n 2 d -4.. (.5,.5,0) (.5,.5,.5) \n 2 c -4.. (0,0,0) (0,0,.5) \n 2 b 2.22 (.5,.5,.25) (.5,.5,.75) \n 2 a 2.22 (0,0,.25) (0,0,.75) \n ";

  string sg117 = "Space Group 117\n (0,0,0) \n 8 i 1 (x,y,z) (-x,-y,z) (y,-x,-z) (-y,x,-z) (x+.5,-y+.5,z) (-x+.5,y+.5,z) (y+.5,x+.5,-z) (-y+.5,-x+.5,-z) \n 4 h ..2 (x,x+.5,.5) (-x,-x+.5,.5) (x+.5,-x,.5) (-x+.5,x,.5) \n 4 g ..2 (x,x+.5,0) (-x,-x+.5,0) (x+.5,-x,0) (-x+.5,x,0) \n 4 f 2.. (0,.5,z) (.5,0,-z) (.5,0,z) (0,.5,-z) \n 4 e 2.. (0,0,z) (0,0,-z) (.5,.5,z) (.5,.5,-z) \n 2 d 2.22 (0,.5,.5) (.5,0,.5) \n 2 c 2.22 (0,.5,0) (.5,0,0) \n 2 b -4.. (0,0,.5) (.5,.5,.5) \n 2 a -4.. (0,0,0) (.5,.5,0) \n ";

  string sg118 = "Space Group 118\n (0,0,0) \n 8 i 1 (x,y,z) (-x,-y,z) (y,-x,-z) (-y,x,-z) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) (y+.5,x+.5,-z+.5) (-y+.5,-x+.5,-z+.5) \n 4 h 2.. (0,.5,z) (.5,0,-z) (.5,0,z+.5) (0,.5,-z+.5) \n 4 g ..2 (x,x+.5,.25) (-x,-x+.5,.25) (x+.5,-x,.75) (-x+.5,x,.75) \n 4 f ..2 (x,-x+.5,.25) (-x,x+.5,.25) (-x+.5,-x,.75) (x+.5,x,.75) \n 4 e 2.. (0,0,z) (0,0,-z) (.5,.5,z+.5) (.5,.5,-z+.5) \n 2 d 2.22 (0,.5,.75) (.5,0,.25) \n 2 c 2.22 (0,.5,.25) (.5,0,.75) \n 2 b -4.. (0,0,.5) (.5,.5,0) \n 2 a -4.. (0,0,0) (.5,.5,.5) \n ";

  string sg119 = "Space Group 119\n (0,0,0) (.5,.5,.5) \n 16 j 1 (x,y,z) (-x,-y,z) (y,-x,-z) (-y,x,-z) (x,-y,z) (-x,y,z) (y,x,-z) (-y,-x,-z) \n 8 i .m. (x,0,z) (-x,0,z) (0,-x,-z) (0,x,-z) \n 8 h ..2 (x,x+.5,.25) (-x,-x+.5,.25) (x+.5,-x,.75) (-x+.5,x,.75) \n 8 g ..2 (x,x,0) (-x,-x,0) (x,-x,0) (-x,x,0) \n 4 f 2mm. (0,.5,z) (.5,0,-z) \n 4 e 2mm. (0,0,z) (0,0,-z) \n 2 d -4m2 (0,.5,.75) \n 2 c -4m2 (0,.5,.25) \n 2 b -4m2 (0,0,.5) \n 2 a -4m2 (0,0,0) \n ";

  string sg120 = "Space Group 120\n (0,0,0) (.5,.5,.5) \n 16 i 1 (x,y,z) (-x,-y,z) (y,-x,-z) (-y,x,-z) (x,-y,z+.5) (-x,y,z+.5) (y,x,-z+.5) (-y,-x,-z+.5) \n 8 h ..2 (x,x+.5,0) (-x,-x+.5,0) (x+.5,-x,0) (-x+.5,x,0) \n 8 g 2.. (0,.5,z) (.5,0,-z) (0,.5,z+.5) (.5,0,-z+.5) \n 8 f 2.. (0,0,z) (0,0,-z) (0,0,z+.5) (0,0,-z+.5) \n 8 e ..2 (x,x,.25) (-x,-x,.25) (x,-x,.75) (-x,x,.75) \n 4 d 2.22 (0,.5,0) (.5,0,0) \n 4 c -4.. (0,.5,.25) (0,.5,.75) \n 4 b -4.. (0,0,0) (0,0,.5) \n 4 a 2.22 (0,0,.25) (0,0,.75) \n ";

  string sg121 = "Space Group 121\n (0,0,0) (.5,.5,.5) \n 16 j 1 (x,y,z) (-x,-y,z) (y,-x,-z) (-y,x,-z) (-x,y,-z) (x,-y,-z) (-y,-x,z) (y,x,z) \n 8 i ..m (x,x,z) (-x,-x,z) (x,-x,-z) (-x,x,-z) \n 8 h 2.. (0,.5,z) (.5,0,-z) (0,.5,-z) (.5,0,z) \n 8 g .2. (x,0,.5) (-x,0,.5) (0,-x,.5) (0,x,.5) \n 8 f .2. (x,0,0) (-x,0,0) (0,-x,0) (0,x,0) \n 4 e 2.mm (0,0,z) (0,0,-z) \n 4 d -4.. (0,.5,.25) (0,.5,.75) \n 4 c 222. (0,.5,0) (.5,0,0) \n 2 b -42m  (0,0,.5) \n 2 a -42m (0,0,0) \n ";

  string sg122 = "Space Group 122\n (0,0,0) (.5,.5,.5) \n 16 e 1 (x,y,z) (-x,-y,z) (y,-x,-z) (-y,x,-z) (-x+.5,y,-z+.75) (x+.5,-y,-z+.75) (-y+.5,-x,z+.75) (y+.5,x,z+.75) \n 8 d .2. (x,.25,.125) (-x,.75,.125) (.25,-x,.875) (.75,x,.875) \n 8 c 2.. (0,0,z) (0,0,-z) (.5,0,-z+.75) (.5,0,z+.75) \n 4 b 2.. (0,0,.5) (.5,0,.25) \n 4 a -4.. (0,0,0) (.5,0,.75) \n ";

  string sg123 = "Space Group 123\n (0,0,0) \n 16 u 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (-x,y,-z) (x,-y,-z) (y,x,-z) (-y,-x,-z) (-x,-y,-z) (x,y,-z) (y,-x,-z) (-y,x,-z) (x,-y,z) (-x,y,z) (-y,-x,z) (y,x,z) \n 8 t .m. (x,.5,z) (-x,.5,z) (.5,x,z) (.5,-x,z) (-x,.5,-z) (x,.5,-z) (.5,x,-z) (.5,-x,-z) \n 8 s .m. (x,0,z) (-x,0,z) (0,x,z) (0,-x,z) (-x,0,-z) (x,0,-z) (0,x,-z) (0,-x,-z) \n 8 r ..m (x,x,z) (-x,-x,z) (-x,x,z) (x,-x,z) (-x,x,-z) (x,-x,-z) (x,x,-z) (-x,-x,-z) \n 8 q m.. (x,y,.5) (-x,-y,.5) (-y,x,.5) (y,-x,.5) (-x,y,.5) (x,-y,.5) (y,x,.5) (-y,-x,.5) \n 8 p m.. (x,y,0) (-x,-y,0) (-y,x,0) (y,-x,0) (-x,y,0) (x,-y,0) (y,x,0) (-y,-x,0) \n 4 o m2m. (x,.5,.5) (-x,.5,.5) (.5,x,.5) (.5,-x,.5) \n 4 n m2m. (x,.5,0) (-x,.5,0) (.5,x,0) (.5,-x,0) \n 4 m m2m. (x,0,.5) (-x,0,.5) (0,x,.5) (0,-x,.5) \n 4 l m2m. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) \n 4 k m.2m (x,x,.5) (-x,-x,.5) (-x,x,.5) (x,-x,.5) \n 4 j m.2m (x,x,0) (-x,-x,0) (-x,x,0) (x,-x,0) \n 4 i 2mm. (0,.5,z) (.5,0,z) (0,.5,-z) (.5,0,-z) \n 2 h 4mm (.5,.5,z) (.5,.5,-z) \n 2 g 4mm (0,0,z) (0,0,-z) \n 2 f mmm. (0,.5,0) (.5,0,0) \n 2 e mmm. (0,.5,.5) (.5,0,.5) \n 1 d 4/mmm (.5,.5,.5) \n 1 c 4/mmm (.5,.5,0) \n 1 b 4/mmm (0,0,.5) \n 1 a 4/mmm (0,0,0) \n ";

  string sg124 = "Space Group 124\n (0,0,0) \n 16 n 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (-x,y,-z+.5) (x,-y,-z+.5) (y,x,-z+.5) (-y,-x,-z+.5) (-x,-y,-z) (x,y,-z) (y,-x,-z) (-y,x,-z) (x,-y,z+.5) (-x,y,z+.5) (-y,-x,z+.5) (y,x,z+.5) \n 8 m m.. (x,y,0) (-x,-y,0) (-y,x,0) (y,-x,0) (-x,y,.5) (x,-y,.5) (y,x,.5) (-y,-x,.5) \n 8 l .2. (x,.5,.25) (-x,.5,.25) (.5,x,.25) (.5,-x,.25) (-x,.5,.75) (x,.5,.75) (.5,-x,.75) (.5,x,.75) \n 8 k .2. (x,0,.25) (-x,0,.25) (0,x,.25) (0,-x,.25) (-x,0,.75) (x,0,.75) (0,-x,.75) (0,x,.75) \n 8 j ..2 (x,x,.25) (-x,-x,.25) (-x,x,.25) (x,-x,.25) (-x,-x,.75) (x,x,.75) (x,-x,.75) (-x,x,.75) \n 8 i 2.. (0,.5,z) (.5,0,z) (0,.5,-z+.5) (.5,0,-z+.5) (0,.5,-z) (.5,0,-z) (0,.5,z+.5) (.5,0,z+.5) \n 4 h 4.. (.5,.5,z) (.5,.5,-z+.5) (.5,.5,-z) (.5,.5,z+.5) \n 4 g 4.. (0,0,z) (0,0,-z+.5) (0,0,-z) (0,0,z+.5) \n 4 f 222. (0,.5,.25) (.5,0,.25) (0,.5,.75) (.5,0,.75) \n 4 e 2/m.. (0,.5,0) (.5,0,0) (0,.5,.5) (.5,0,.5) \n 2 d 4/m.. (.5,.5,0) (.5,.5,.5) \n 2 c 422 (.5,.5,.25) (.5,.5,.75) \n 2 b 4/m.. (0,0,0) (0,0,.5) \n 2 a 422 (0,0,.25) (0,0,.75) \n ";

  //DX20180723 - added origin choice 2
  string sg125 = "";
  if(axis_cell==setting_1){
    sg125 = "Space Group 125\n (0,0,0) \n 16 n 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (-x,y,-z) (x,-y,-z) (y,x,-z) (-y,-x,-z) (-x+.5,-y+.5,-z) (x+.5,y+.5,-z) (y+.5,-x+.5,-z) (-y+.5,x+.5,-z) (x+.5,-y+.5,z) (-x+.5,y+.5,z) (-y+.5,-x+.5,z) (y+.5,x+.5,z) \n 8 m ..m (x,x+.5,z) (-x,-x+.5,z) (-x+.5,x,z) (x+.5,-x,z) (-x,x+.5,-z) (x,-x+.5,-z) (x+.5,x,-z) (-x+.5,-x,-z) \n 8 l .2. (x,0,.5) (-x,0,.5) (0,x,.5) (0,-x,.5) (-x+.5,.5,.5) (x+.5,.5,.5) (.5,-x+.5,.5) (.5,x+.5,.5) \n 8 k .2. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (-x+.5,.5,0) (x+.5,.5,0) (.5,-x+.5,0) (.5,x+.5,0) \n 8 j ..2 (x,x,.5) (-x,-x,.5) (-x,x,.5) (x,-x,.5) (-x+.5,-x+.5,.5) (x+.5,x+.5,.5) (x+.5,-x+.5,.5) (-x+.5,x+.5,.5) \n 8 i ..1 (x,x,0) (-x,-x,0) (-x,x,0) (x,-x,0) (-x+.5,-x+.5,0) (x+.5,x+.5,0) (x+.5,-x+.5,0) (-x+.5,x+.5,0) \n 4 h 2.mm (0,.5,z) (.5,0,z) (0,.5,-z) (.5,0,-z) \n 4 g 4.. (0,0,z) (0,0,-z) (.5,.5,-z) (.5,.5,z) \n 4 f ..2/m (.25,.25,.5) (.75,.75,.5) (.75,.25,.5) (.25,.75,.5) \n 4 e ..2/m (.25,.25,0) (.75,.75,0) (.75,.25,0) (.25,.75,0) \n 2 d -42m (0,.5,.5) (.5,0,.5) \n 2 c -42m (0,.5,.0) (.5,0,0) \n 2 b 422 (0,0,.5) (.5,.5,.5) \n 2 a 422 (0,0,0) (.5,.5,0) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg125 = "Space Group 125\n (0,0,0) \n 16 n 1 (x,y,z) (-x+.5,-y+.5,z) (-y+.5,x,z) (y,-x+.5,z) (-x+.5,y,-z) (x,-y+.5,-z) (y,x,-z) (-y+.5,-x+.5,-z) (-x,-y,-z) (x+.5,y+.5,-z) (y+.5,-x,-z) (-y,x+.5,-z) (x+.5,-y,z) (-x,y+.5,z) (-y,-x,z) (y+.5,x+.5,z) \n 8 m ..m (x,-x,z) (-x+.5,x+.5,z) (x+.5,x,z) (-x,-x+.5,z) (-x+.5,-x,-z) (x,x+.5,-z) (-x,x,-z) (x+.5,-x+.5,-z) \n 8 l .2. (x,.25,.5) (-x+.5,.25,.5) (.25,x,.5) (.25,-x+.5,.5) (-x,.75,.5) (x+.5,.75,.5) (.75,-x,.5) (.75,x+.5,.5) \n 8 k .2. (x,.25,0) (-x+.5,.25,0) (.25,x,0) (.25,-x+.5,0) (-x,.75,0) (x+.5,.75,0) (.75,-x,0) (.75,x+.5,0) \n 8 j ..2 (x,x,.5) (-x+.5,-x+.5,.5) (-x+.5,x,.5) (x,-x+.5,.5) (-x,-x,.5) (x+.5,x+.5,.5) (x+.5,-x,.5) (-x,x+.5,.5) \n 8 i ..1 (x,x,0) (-x+.5,-x+.5,0) (-x+.5,x,0) (x,-x+.5,0) (-x,-x,0) (x+.5,x+.5,0) (x+.5,-x,0) (-x,x+.5,0) \n 4 h 2.mm (.75,.25,z) (.25,.75,z) (.75,.25,-z) (.25,.75,-z) \n 4 g 4.. (.25,.25,z) (.25,.25,-z) (.75,.75,-z) (.75,.75,z) \n 4 f ..2/m (0,0,.5) (.5,.5,.5) (.5,0,.5) (0,.5,.5) \n 4 e ..2/m (0,0,0) (.5,.5,0) (.5,0,0) (0,.5,0) \n 2 d -42m (.75,.25,.5) (.25,.75,.5) \n 2 c -42m (.75,.25,.0) (.25,.75,0) \n 2 b 422 (.25,.25,.5) (.75,.75,.5) \n 2 a 422 (.25,.25,0) (.75,.75,0) \n ";
  }
  else {
    sg125 = "Space Group 125\n (0,0,0) \n 16 n 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (-x,y,-z) (x,-y,-z) (y,x,-z) (-y,-x,-z) (-x+.5,-y+.5,-z) (x+.5,y+.5,-z) (y+.5,-x+.5,-z) (-y+.5,x+.5,-z) (x+.5,-y+.5,z) (-x+.5,y+.5,z) (-y+.5,-x+.5,z) (y+.5,x+.5,z) \n 8 m ..m (x,x+.5,z) (-x,-x+.5,z) (-x+.5,x,z) (x+.5,-x,z) (-x,x+.5,-z) (x,-x+.5,-z) (x+.5,x,-z) (-x+.5,-x,-z) \n 8 l .2. (x,0,.5) (-x,0,.5) (0,x,.5) (0,-x,.5) (-x+.5,.5,.5) (x+.5,.5,.5) (.5,-x+.5,.5) (.5,x+.5,.5) \n 8 k .2. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (-x+.5,.5,0) (x+.5,.5,0) (.5,-x+.5,0) (.5,x+.5,0) \n 8 j ..2 (x,x,.5) (-x,-x,.5) (-x,x,.5) (x,-x,.5) (-x+.5,-x+.5,.5) (x+.5,x+.5,.5) (x+.5,-x+.5,.5) (-x+.5,x+.5,.5) \n 8 i ..1 (x,x,0) (-x,-x,0) (-x,x,0) (x,-x,0) (-x+.5,-x+.5,0) (x+.5,x+.5,0) (x+.5,-x+.5,0) (-x+.5,x+.5,0) \n 4 h 2.mm (0,.5,z) (.5,0,z) (0,.5,-z) (.5,0,-z) \n 4 g 4.. (0,0,z) (0,0,-z) (.5,.5,-z) (.5,.5,z) \n 4 f ..2/m (.25,.25,.5) (.75,.75,.5) (.75,.25,.5) (.25,.75,.5) \n 4 e ..2/m (.25,.25,0) (.75,.75,0) (.75,.25,0) (.25,.75,0) \n 2 d -42m (0,.5,.5) (.5,0,.5) \n 2 c -42m (0,.5,.0) (.5,0,0) \n 2 b 422 (0,0,.5) (.5,.5,.5) \n 2 a 422 (0,0,0) (.5,.5,0) \n ";
  }

  //DX20180723 - added origin choice 2
  string sg126 = "";
  if(axis_cell==setting_1){
    sg126 = "Space Group 126\n (0,0,0) \n 16 k 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (-x,y,-z) (x,-y,-z) (y,x,-z) (-y,-x,-z) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (y+.5,-x+.5,-z+.5) (-y+.5,x+.5,-z+.5) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) (-y+.5,-x+.5,z+.5) (y+.5,x+.5,z+.5) \n 8 j .2. (x,0,.5) (-x,0,.5) (0,x,.5) (0,-x,.5) (-x+.5,.5,0) (x+.5,.5,0) (.5,-x+.5,0) (.5,x+.5,0) \n 8 i .2. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (-x+.5,.5,.5) (x+.5,.5,.5) (.5,-x+.5,.5) (.5,x+.5,.5) \n 8 h ..2 (x,x,0) (-x,-x,0) (-x,x,0) (x,-x,0) (-x+.5,-x+.5,.5) (x+.5,x+.5,.5) (x+.5,-x+.5,.5) (-x+.5,x+.5,.5) \n 8 g 2.. (.5,0,z) (0,.5,z) (.5,0,-z) (0,.5,-z) (0,.5,-z+.5) (.5,0,-z+.5) (0,.5,z+.5) (.5,0,z+.5) \n 8 f -1 (.25,.25,.25) (.75,.75,.25) (.75,.25,.25) (.25,.75,.25) (.75,.25,.75) (.25,.75,.75) (.25,.25,.75) (.75,.75,.75) \n 4 e 4.. (0,0,z) (0,0,-z) (.5,.5,-z+.5) (.5,.5,z+.5) \n 4 d -4.. (.5,0,.25) (0,.5,.25) (.5,0,.75) (0,.5,.75) \n 4 c 222. (.5,0,0) (0,.5,0) (0,.5,.5) (.5,0,.5) \n 2 b 422 (0,0,.5) (.5,.5,0) \n 2 a 422 (0,0,0) (.5,.5,.5) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg126 = "Space Group 126\n (0,0,0) \n 16 k 1 (x,y,z) (-x+.5,-y+.5,z) (-y+.5,x,z) (y,-x+.5,z) (-x+.5,y,-z+.5) (x,-y+.5,-z+.5) (y,x,-z+.5) (-y+.5,-x+.5,-z+.5) (-x,-y,-z) (x+.5,y+.5,-z) (y+.5,-x,-z) (-y,x+.5,-z) (x+.5,-y,z+.5) (-x,y+.5,z+.5) (-y,-x,z+.5) (y+.5,x+.5,z+.5) \n 8 j .2. (x,.75,.25) (-x+.5,.75,.25) (.75,x,.25) (.75,-x+.5,.25) (-x,.25,.75) (x+.5,.25,.75) (.25,-x,.75) (.25,x+.5,.75) \n 8 i .2. (x,.25,.25) (-x+.5,.25,.25) (.25,x,.25) (.25,-x+.5,.25) (-x,.75,.75) (x+.5,.75,.75) (.75,-x,.75) (.75,x+.5,.75) \n 8 h ..2 (x,x,.25) (-x+.5,-x+.5,.25) (-x+.5,x,.25) (x,-x+.5,.25) (-x,-x,.75) (x+.5,x+.5,.75) (x+.5,-x,.75) (-x,x+.5,.75) \n 8 g 2.. (.25,.75,z) (.75,.25,z) (.25,.75,-z+.5) (.75,.25,-z+.5) (.75,.25,-z) (.25,.75,-z) (.75,.25,z+.5) (.25,.75,z+.5) \n 8 f -1 (0,0,0) (.5,.5,0) (.5,0,0) (0,.5,0) (.5,0,.5) (0,.5,.5) (0,0,.5) (.5,.5,.5) \n 4 e 4.. (.25,.25,z) (.25,.25,-z+.5) (.75,.75,-z) (.75,.75,z+.5) \n 4 d -4.. (.25,.75,0) (.75,.25,0) (.25,.75,.5) (.75,.25,.5) \n 4 c 222. (.25,.75,.75) (.75,.25,.75) (.75,.25,.25) (.25,.75,.25) \n 2 b 422 (.25,.25,.75) (.75,.75,.25) \n 2 a 422 (.25,.25,.25) (.75,.75,.75) \n ";
  }
  else {
    sg126 = "Space Group 126\n (0,0,0) \n 16 k 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (-x,y,-z) (x,-y,-z) (y,x,-z) (-y,-x,-z) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (y+.5,-x+.5,-z+.5) (-y+.5,x+.5,-z+.5) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) (-y+.5,-x+.5,z+.5) (y+.5,x+.5,z+.5) \n 8 j .2. (x,0,.5) (-x,0,.5) (0,x,.5) (0,-x,.5) (-x+.5,.5,0) (x+.5,.5,0) (.5,-x+.5,0) (.5,x+.5,0) \n 8 i .2. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (-x+.5,.5,.5) (x+.5,.5,.5) (.5,-x+.5,.5) (.5,x+.5,.5) \n 8 h ..2 (x,x,0) (-x,-x,0) (-x,x,0) (x,-x,0) (-x+.5,-x+.5,.5) (x+.5,x+.5,.5) (x+.5,-x+.5,.5) (-x+.5,x+.5,.5) \n 8 g 2.. (.5,0,z) (0,.5,z) (.5,0,-z) (0,.5,-z) (0,.5,-z+.5) (.5,0,-z+.5) (0,.5,z+.5) (.5,0,z+.5) \n 8 f -1 (.25,.25,.25) (.75,.75,.25) (.75,.25,.25) (.25,.75,.25) (.75,.25,.75) (.25,.75,.75) (.25,.25,.75) (.75,.75,.75) \n 4 e 4.. (0,0,z) (0,0,-z) (.5,.5,-z+.5) (.5,.5,z+.5) \n 4 d -4.. (.5,0,.25) (0,.5,.25) (.5,0,.75) (0,.5,.75) \n 4 c 222. (.5,0,0) (0,.5,0) (0,.5,.5) (.5,0,.5) \n 2 b 422 (0,0,.5) (.5,.5,0) \n 2 a 422 (0,0,0) (.5,.5,.5) \n ";
  }

  string sg127 = "Space Group 127\n (0,0,0) \n 16 l 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (-x+.5,y+.5,-z) (x+.5,-y+.5,-z) (y+.5,x+.5,-z) (-y+.5,-x+.5,-z) (-x,-y,-z) (x,y,-z) (y,-x,-z) (-y,x,-z) (x+.5,-y+.5,z) (-x+.5,y+.5,z) (-y+.5,-x+.5,z) (y+.5,x+.5,z) \n 8 k ..m (x,x+.5,z) (-x,-x+.5,z) (-x+.5,x,z) (x+.5,-x,z) (-x+.5,x,-z) (x+.5,-x,-z) (x,x+.5,-z) (-x,-x+.5,-z) \n 8 j m.. (x,y,.5) (-x,-y,.5) (-y,x,.5) (y,-x,.5) (-x+.5,y+.5,.5) (x+.5,-y+.5,.5) (y+.5,x+.5,.5) (-y+.5,-x+.5,.5) \n 8 i m.. (x,y,0) (-x,-y,0) (-y,x,0) (y,-x,0) (-x+.5,y+.5,0) (x+.5,-y+.5,0) (y+.5,x+.5,0) (-y+.5,-x+.5,0) \n 4 h m.2m (x,x+.5,.5) (-x,-x+.5,.5) (-x+.5,x,.5) (x+.5,-x,.5) \n 4 g m.2m (x,x+.5,0) (-x,-x+.5,0) (-x+.5,x,0) (x+.5,-x,0) \n 4 f 2.mm (0,.5,z) (.5,0,z) (.5,0,-z) (0,.5,-z) \n 4 e 4.. (0,0,z) (.5,.5,-z) (0,0,-z) (.5,.5,z) \n 2 d m.mm (0,.5,0) (.5,0,0) \n 2 c m.mm (0,.5,.5) (.5,0,.5) \n 2 b 4/m.. (0,0,.5) (.5,.5,.5) \n 2 a 4/m.. (0,0,0) (.5,.5,0) \n ";

  string sg128 = "Space Group 128\n (0,0,0) \n 16 i 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (-x+.5,y+.5,-z+.5) (x+.5,-y+.5,-z+.5) (y+.5,x+.5,-z+.5) (-y+.5,-x+.5,-z+.5) (-x,-y,-z) (x,y,-z) (y,-x,-z) (-y,x,-z) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) (-y+.5,-x+.5,z+.5) (y+.5,x+.5,z+.5) \n 8 h m.. (x,y,0) (-x,-y,0) (-y,x,0) (y,-x,0) (-x+.5,y+.5,.5) (x+.5,-y+.5,.5) (y+.5,x+.5,.5) (-y+.5,-x+.5,.5) \n 8 g ..2 (x,x+.5,.25) (-x,-x+.5,.25) (-x+.5,x,.25) (x+.5,-x,.25) (-x,-x+.5,.75) (x,x+.5,.75) (x+.5,-x,.75) (-x+.5,x,.75) \n 8 f 2.. (0,.5,z) (.5,0,z) (.5,0,-z+.5) (0,.5,-z+.5) (0,.5,-z) (.5,0,-z) (.5,0,z+.5) (0,.5,z+.5) \n 4 e 4.. (0,0,z) (.5,.5,-z+.5) (0,0,-z) (.5,.5,z+.5) \n 4 d 2.22 (0,.5,.25) (.5,0,.25) (0,.5,.75) (.5,0,.75) \n 4 c 2/m.. (0,.5,0) (.5,0,0) (.5,0,.5) (0,.5,.5) \n 2 b 4/m.. (0,0,.5) (.5,.5,0) \n 2 a 4/m.. (0,0,0) (.5,.5,.5) \n "; //DX20200528 - missing space between 3th and 4th position for Wyckoff letter g

  //DX20180723 - added origin choice 2
  string sg129 = "";
  if(axis_cell==setting_1){
    sg129 = "Space Group 129\n (0,0,0) \n 16 k 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z) (y+.5,-x+.5,z) (-x+.5,y+.5,-z) (x+.5,-y+.5,-z) (y,x,-z) (-y,-x,-z) (-x+.5,-y+.5,-z) (x+.5,y+.5,-z) (y,-x,-z) (-y,x,-z) (x,-y,z) (-x,y,z) (-y+.5,-x+.5,z) (y+.5,x+.5,z) \n 8 j ..m (x,x+.5,z) (-x,-x+.5,z) (-x,x+.5,z) (x,-x+.5,z) (-x+.5,x,-z) (x+.5,-x,-z) (x+.5,x,-z) (-x+.5,-x,-z) \n 8 i .m. (0,y,z) (0,-y,z) (-y+.5,.5,z) (y+.5,.5,z) (.5,y+.5,-z) (.5,-y+.5,-z) (y,0,-z) (-y,0,-z) \n 8 h ..2 (x,x,.5) (-x,-x,.5) (-x+.5,x+.5,.5) (x+.5,-x+.5,.5) (-x+.5,-x+.5,.5) (x+.5,x+.5,.5) (x,-x,.5) (-x,x,.5) \n 8 g ..2 (x,x,0) (-x,-x,0) (-x+.5,x+.5,0) (x+.5,-x+.5,0) (-x+.5,-x+.5,0) (x+.5,x+.5,0) (x,-x,0) (-x,x,0) \n 4 f 2mm. (0,0,z) (.5,.5,z) (.5,.5,-z) (0,0,-z) \n 4 e ..2/m (.25,.25,.5) (.75,.75,.5) (.25,.75,.5) (.75,.25,.5) \n 4 d ..2/m (.25,.25,0) (.75,.75,0) (.25,.75,0) (.75,.25,0) \n 2 c 4mm (0,.5,z) (.5,0,-z) \n 2 b -4m2 (0,0,.5) (.5,.5,.5) \n 2 a -42m (0,0,0) (.5,.5,0) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg129 = "Space Group 129\n (0,0,0) \n 16 k 1 (x,y,z) (-x+.5,-y+.5,z) (-y+.5,x,z) (y,-x+.5,z) (-x,y+.5,-z) (x+.5,-y,-z) (y+.5,x+.5,-z) (-y,-x,-z) (-x,-y,-z) (x+.5,y+.5,-z) (y+.5,-x,-z) (-y,x+.5,-z) (x,-y+.5,z) (-x+.5,y,z) (-y+.5,-x+.5,z) (y,x,z) \n 8 j ..m (x,x,z) (-x+.5,-x+.5,z) (-x+.5,x,z) (x,-x+.5,z) (-x,x+.5,-z) (x+.5,-x,-z) (x+.5,x+.5,-z) (-x,-x,-z) \n 8 i .m. (.25,y,z) (.25,-y+.5,z) (-y+.5,.25,z) (y,.25,z) (.75,y+.5,-z) (.75,-y,-z) (y+.5,.75,-z) (-y,.75,-z) \n 8 h ..2 (x,-x,.5) (-x+.5,x+.5,.5) (x+.5,x,.5) (-x,-x+.5,.5) (-x,x,.5) (x+.5,-x+.5,.5) (-x+.5,-x,.5) (x,x+.5,.5) \n 8 g ..2 (x,-x,0) (-x+.5,x+.5,0) (x+.5,x,0) (-x,-x+.5,0) (-x,x,0) (x+.5,-x+.5,0) (-x+.5,-x,0) (x,x+.5,0) \n 4 f 2mm. (.75,.25,z) (.25,.75,z) (.25,.75,-z) (.75,.25,-z) \n 4 e ..2/m (0,0,.5) (.5,.5,.5) (.5,0,.5) (0,.5,.5) \n 4 d ..2/m (0,0,0) (.5,.5,0) (.5,0,0) (0,.5,0) \n 2 c 4mm (.25,.25,z) (.75,.75,-z) \n 2 b -4m2 (.75,.25,.5) (.25,.75,.5) \n 2 a -42m (.75,.25,0) (.25,.75,0) \n ";
  }
  else {
    sg129 = "Space Group 129\n (0,0,0) \n 16 k 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z) (y+.5,-x+.5,z) (-x+.5,y+.5,-z) (x+.5,-y+.5,-z) (y,x,-z) (-y,-x,-z) (-x+.5,-y+.5,-z) (x+.5,y+.5,-z) (y,-x,-z) (-y,x,-z) (x,-y,z) (-x,y,z) (-y+.5,-x+.5,z) (y+.5,x+.5,z) \n 8 j ..m (x,x+.5,z) (-x,-x+.5,z) (-x,x+.5,z) (x,-x+.5,z) (-x+.5,x,-z) (x+.5,-x,-z) (x+.5,x,-z) (-x+.5,-x,-z) \n 8 i .m. (0,y,z) (0,-y,z) (-y+.5,.5,z) (y+.5,.5,z) (.5,y+.5,-z) (.5,-y+.5,-z) (y,0,-z) (-y,0,-z) \n 8 h ..2 (x,x,.5) (-x,-x,.5) (-x+.5,x+.5,.5) (x+.5,-x+.5,.5) (-x+.5,-x+.5,.5) (x+.5,x+.5,.5) (x,-x,.5) (-x,x,.5) \n 8 g ..2 (x,x,0) (-x,-x,0) (-x+.5,x+.5,0) (x+.5,-x+.5,0) (-x+.5,-x+.5,0) (x+.5,x+.5,0) (x,-x,0) (-x,x,0) \n 4 f 2mm. (0,0,z) (.5,.5,z) (.5,.5,-z) (0,0,-z) \n 4 e ..2/m (.25,.25,.5) (.75,.75,.5) (.25,.75,.5) (.75,.25,.5) \n 4 d ..2/m (.25,.25,0) (.75,.75,0) (.25,.75,0) (.75,.25,0) \n 2 c 4mm (0,.5,z) (.5,0,-z) \n 2 b -4m2 (0,0,.5) (.5,.5,.5) \n 2 a -42m (0,0,0) (.5,.5,0) \n ";
  }

  //DX20180723 - added origin choice 2
  string sg130 = "";
  if(axis_cell==setting_1){
    sg130 = "Space Group 130\n (0,0,0) \n 16 g 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z) (y+.5,-x+.5,z) (-x+.5,y+.5,-z+.5) (x+.5,-y+.5,-z+.5) (y,x,-z+.5) (-y,-x,-z+.5) (-x+.5,-y+.5,-z) (x+.5,y+.5,-z) (y,-x,-z) (-y,x,-z) (x,-y,z+.5) (-x,y,z+.5) (-y+.5,-x+.5,z+.5) (y+.5,x+.5,z+.5) \n 8 f ..2 (x,x,.25) (-x,-x,.25) (-x+.5,x+.5,.25) (x+.5,-x+.5,.25) (-x+.5,-x+.5,.75) (x+.5,x+.5,.75) (x,-x,.75) (-x,x,.75) \n 8 e 2.. (0,0,z) (.5,.5,z) (.5,.5,-z+.5) (0,0,-z+.5) (.5,.5,-z) (0,0,-z) (0,0,z+.5) (.5,.5,z+.5) \n 8 d -1 (.25,.25,0) (.75,.75,0) (.25,.75,0) (.75,.25,0) (.25,.75,.5) (.75,.25,.5) (.25,.25,.5) (.75,.75,.5) \n 4 c 4.. (0,.5,z) (.5,0,-z+.5) (.5,0,-z) (0,.5,z+.5) \n 4 b -4.. (0,0,0) (.5,.5,0) (.5,.5,.5) (0,0,.5) \n 4 a 2.22 (0,0,.25) (.5,.5,.25) (.5,.5,.75) (0,0,.75) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg130 = "Space Group 130\n (0,0,0) \n 16 g 1 (x,y,z) (-x+.5,-y+.5,z) (-y+.5,x,z) (y,-x+.5,z) (-x,y+.5,-z+.5) (x+.5,-y,-z+.5) (y+.5,x+.5,-z+.5) (-y,-x,-z+.5) (-x,-y,-z) (x+.5,y+.5,-z) (y+.5,-x,-z) (-y,x+.5,-z) (x,-y+.5,z+.5) (-x+.5,y,z+.5) (-y+.5,-x+.5,z+.5) (y,x,z+.5) \n 8 f ..2 (x,-x,.25) (-x+.5,x+.5,.25) (x+.5,x,.25) (-x,-x+.5,.25) (-x,x,.75) (x+.5,-x+.5,.75) (-x+.5,-x,.75) (x,x+.5,.75) \n 8 e 2.. (.75,.25,z) (.25,.75,z) (.25,.75,-z+.5) (.75,.25,-z+.5) (.25,.75,-z) (.75,.25,-z) (.75,.25,z+.5) (.25,.75,z+.5) \n 8 d -1 (0,0,0) (.5,.5,0) (.5,0,0) (0,.5,0) (0,.5,.5) (.5,0,.5) (.5,.5,.5) (0,0,.5) \n 4 c 4.. (.25,.25,z) (.75,.75,-z+.5) (.75,.75,-z) (.25,.25,z+.5) \n 4 b -4.. (.75,.25,0) (.25,.75,0) (.25,.75,.5) (.75,.25,.5) \n 4 a 2.22 (.75,.25,.25) (.25,.75,.25) (.25,.75,.75) (.75,.25,.75) \n ";
  }
  else {
    sg130 = "Space Group 130\n (0,0,0) \n 16 g 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z) (y+.5,-x+.5,z) (-x+.5,y+.5,-z+.5) (x+.5,-y+.5,-z+.5) (y,x,-z+.5) (-y,-x,-z+.5) (-x+.5,-y+.5,-z) (x+.5,y+.5,-z) (y,-x,-z) (-y,x,-z) (x,-y,z+.5) (-x,y,z+.5) (-y+.5,-x+.5,z+.5) (y+.5,x+.5,z+.5) \n 8 f ..2 (x,x,.25) (-x,-x,.25) (-x+.5,x+.5,.25) (x+.5,-x+.5,.25) (-x+.5,-x+.5,.75) (x+.5,x+.5,.75) (x,-x,.75) (-x,x,.75) \n 8 e 2.. (0,0,z) (.5,.5,z) (.5,.5,-z+.5) (0,0,-z+.5) (.5,.5,-z) (0,0,-z) (0,0,z+.5) (.5,.5,z+.5) \n 8 d -1 (.25,.25,0) (.75,.75,0) (.25,.75,0) (.75,.25,0) (.25,.75,.5) (.75,.25,.5) (.25,.25,.5) (.75,.75,.5) \n 4 c 4.. (0,.5,z) (.5,0,-z+.5) (.5,0,-z) (0,.5,z+.5) \n 4 b -4.. (0,0,0) (.5,.5,0) (.5,.5,.5) (0,0,.5) \n 4 a 2.22 (0,0,.25) (.5,.5,.25) (.5,.5,.75) (0,0,.75) \n ";
  }

  string sg131 = "Space Group 131\n (0,0,0) \n 16 r 1 (x,y,z) (-x,-y,z) (-y,x,z+.5) (y,-x,z+.5) (-x,y,-z) (x,-y,-z) (y,x,-z+.5) (-y,-x,-z+.5) (-x,-y,-z) (x,y,-z) (y,-x,-z+.5) (-y,x,-z+.5) (x,-y,z) (-x,y,z) (-y,-x,z+.5) (y,x,z+.5) \n 8 q m.. (x,y,0) (-x,-y,0) (-y,x,.5) (y,-x,.5) (-x,y,0) (x,-y,0) (y,x,.5) (-y,-x,.5) \n 8 p .m. (.5,y,z) (.5,-y,z) (-y,.5,z+.5) (y,.5,z+.5) (.5,y,-z) (.5,-y,-z) (y,.5,-z+.5) (-y,.5,-z+.5) \n 8 o .m. (0,y,z) (0,-y,z) (-y,0,z+.5) (y,0,z+.5) (0,y,-z) (0,-y,-z) (y,0,-z+.5) (-y,0,-z+.5) \n 8 n ..2 (x,x,.25) (-x,-x,.25) (-x,x,.75) (x,-x,.75) (-x,-x,.75) (x,x,.75) (x,-x,.25) (-x,x,.25) \n 4 m m2m. (x,.5,0) (-x,.5,0) (.5,x,.5) (.5,-x,.5) \n 4 l m2m. (x,0,.5) (-x,0,.5) (0,x,0) (0,-x,0) \n 4 k m2m. (x,.5,.5) (-x,.5,.5) (.5,x,0) (.5,-x,0) \n 4 j m2m. (x,0,0) (-x,0,0) (0,x,.5) (0,-x,.5) \n 4 i 2mm. (0,.5,z) (.5,0,z+.5) (0,.5,-z) (.5,0,-z+.5) \n 4 h 2mm. (.5,.5,z) (.5,.5,z+.5) (.5,.5,-z) (.5,.5,-z+.5) \n 4 g 2mm. (0,0,z) (0,0,z+.5) (0,0,-z) (0,0,-z+.5) \n 2 f -4m2 (.5,.5,.25) (.5,.5,.75) \n 2 e -4m2 (0,0,.25) (0,0,.75) \n 2 d mmm. (0,.5,.5) (.5,0,0) \n 2 c mmm. (0,.5,0) (.5,0,.5) \n 2 b mmm. (.5,.5,0) (.5,.5,.5) \n 2 a mmm. (0,0,0) (0,0,.5) \n ";

  string sg132 = "Space Group 132\n (0,0,0) \n 16 p 1 (x,y,z) (-x,-y,z) (-y,x,z+.5) (y,-x,z+.5) (-x,y,-z+.5) (x,-y,-z+.5) (y,x,-z) (-y,-x,-z) (-x,-y,-z) (x,y,-z) (y,-x,-z+.5) (-y,x,-z+.5) (x,-y,z+.5) (-x,y,z+.5) (-y,-x,z) (y,x,z) \n 8 o ..m (x,x,z) (-x,-x,z) (-x,x,z+.5) (x,-x,z+.5) (-x,x,-z+.5) (x,-x,-z+.5) (x,x,-z) (-x,-x,-z) \n 8 n m.. (x,y,0) (-x,-y,0) (-y,x,.5) (y,-x,.5) (-x,y,.5) (x,-y,.5) (y,x,0) (-y,-x,0) \n 8 m .2. (x,.5,.25) (-x,.5,.25) (.5,x,.75) (.5,-x,.75) (-x,.5,.75) (x,.5,.75) (.5,-x,.25) (.5,x,.25) \n 8 l .2. (x,0,.25) (-x,0,.25) (0,x,.75) (0,-x,.75) (-x,0,.75) (x,0,.75) (0,-x,.25) (0,x,.25) \n 8 k 2.. (0,.5,z) (.5,0,z+.5) (0,.5,-z+.5) (.5,0,-z) (0,.5,-z) (.5,0,-z+.5) (0,.5,z+.5) (.5,0,z) \n 4 j m.2m (x,x,.5) (-x,-x,.5) (-x,x,0) (x,-x,0) \n 4 i m.2m (x,x,0) (-x,-x,0) (-x,x,.5) (x,-x,.5) \n 4 h 2.mm (.5,.5,z) (.5,.5,z+.5) (.5,.5,-z+.5) (.5,.5,-z) \n 4 g 2.mm (0,0,z) (0,0,z+.5) (0,0,-z+.5) (0,0,-z) \n 4 f 2/m.. (0,.5,0) (.5,0,.5) (0,.5,.5) (.5,0,0) \n 4 e 222. (0,.5,.25) (.5,0,.75) (0,.5,.75) (.5,0,.25) \n 2 d -42m (.5,.5,.25) (.5,.5,.75) \n 2 c m.mm (.5,.5,0) (.5,.5,.5) \n 2 b -42m (0,0,.25) (0,0,.75) \n 2 a m.mm (0,0,0) (0,0,.5) \n ";

  //DX20180723 - added origin choice 2
  string sg133 = "";
  if(axis_cell==setting_1){
    sg133 = "Space Group 133\n (0,0,0) \n 16 k 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z+.5) (y+.5,-x+.5,z+.5) (-x,y,-z+.5) (x,-y,-z+.5) (y+.5,x+.5,-z) (-y+.5,-x+.5,-z) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (y,-x,-z) (-y,x,-z) (x+.5,-y+.5,z) (-x+.5,y+.5,z) (-y,-x,z+.5) (y,x,z+.5) \n 8 j ..2 (x,x+.5,0) (-x,-x+.5,0) (-x,x+.5,.5) (x,-x+.5,.5) (-x+.5,-x,.5) (x+.5,x,.5) (x+.5,-x,0) (-x+.5,x,0) \n 8 i .2. (x,0,.75) (-x,0,.75) (.5,x+.5,.25) (.5,-x+.5,.25) (-x+.5,.5,.75) (x+.5,.5,.75) (0,-x,.25) (0,x,.25) \n 8 h .2. (x,0,.25) (-x,0,.25) (.5,x+.5,.75) (.5,-x+.5,.75) (-x+.5,.5,.25) (x+.5,.5,.25) (0,-x,.75) (0,x,.75) \n 8 g 2.. (0,0,z) (.5,.5,z+.5) (0,0,-z+.5) (.5,.5,-z) (.5,.5,-z+.5) (0,0,-z) (.5,.5,z) (0,0,z+.5) \n 8 f 2.. (0,.5,z) (0,.5,z+.5) (0,.5,-z+.5) (0,.5,-z) (.5,0,-z+.5) (.5,0,-z) (.5,0,z) (.5,0,z+.5) \n 8 e -1 (.25,.25,.25) (.75,.75,.25) (.25,.75,.75) (.75,.25,.25) (.25,.75,.25) (.75,.75,.75) (.25,.25,.75) \n 4 d -4.. (0,0,0) (.5,.5,.5) (0,0,.5) (.5,.5,0) \n 4 c 2.22 (0,.5,0) (0,.5,.5) (.5,0,.5) (.5,0,0) \n 4 b 222. (0,0,.25) (.5,.5,.75) (.5,.5,.25) (0,0,.75) \n 4 a 222. (0,.5,.25) (0,.5,.75) (.5,0,.25) (.5,0,.75) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg133 = "Space Group 133\n (0,0,0) \n 16 k 1 (x,y,z) (-x+.5,-y+.5,z) (-y+.5,x,z+.5) (y,-x+.5,z+.5) (-x+.5,y,-z) (x,-y+.5,-z) (y,x,-z+.5) (-y+.5,-x+.5,-z+.5) (-x,-y,-z) (x+.5,y+.5,-z) (y+.5,-x,-z+.5) (-y,x+.5,-z+.5) (x+.5,-y,z) (-x,y+.5,z) (-y,-x,z+.5) (y+.5,x+.5,z+.5) \n 8 j ..2 (x,x,.25) (-x+.5,-x+.5,.25) (-x+.5,x,.75) (x,-x+.5,.75) (-x,-x,.75) (x+.5,x+.5,.75) (x+.5,-x,.25) (-x,x+.5,.25) \n 8 i .2. (x,.25,.5) (-x+.5,.25,.5) (.25,x,0) (.25,-x+.5,0) (-x,.75,.5) (x+.5,.75,.5) (.75,-x,0) (.75,x+.5,0) \n 8 h .2. (x,.25,0) (-x+.5,.25,0) (.25,x,.5) (.25,-x+.5,.5) (-x,.75,0) (x+.5,.75,0) (.75,-x,.5) (.75,x+.5,.5) \n 8 g 2.. (.75,.25,z) (.25,.75,z+.5) (.75,.25,-z) (.25,.75,-z+.5) (.25,.75,-z) (.75,.25,-z+.5) (.25,.75,z) (.75,.25,z+.5) \n 8 f 2.. (.25,.25,z) (.25,.25,z+.5) (.25,.25,-z) (.25,.25,-z+.5) (.75,.75,-z) (.75,.75,-z+.5) (.75,.75,z) (.75,.75,z+.5) \n 8 e -1 (0,0,0) (.5,.5,0) (.5,0,.5) (0,.5,.5) (.5,0,0) (0,.5,0) (0,0,.5) \n 4 d -4.. (.75,.25,.75) (.25,.75,.25) (.75,.25,.25) (.25,.75,.75) \n 4 c 2.22 (.25,.25,.25) (.25,.25,.75) (.75,.75,.75) (.75,.75,.25) \n 4 b 222. (.75,.25,0) (.25,.75,.5) (.25,.75,0) (.75,.25,.5) \n 4 a 222. (.25,.25,0) (.25,.25,.5) (.75,.75,0) (.75,.75,.5) \n ";
  }
  else {
    sg133 = "Space Group 133\n (0,0,0) \n 16 k 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z+.5) (y+.5,-x+.5,z+.5) (-x,y,-z+.5) (x,-y,-z+.5) (y+.5,x+.5,-z) (-y+.5,-x+.5,-z) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (y,-x,-z) (-y,x,-z) (x+.5,-y+.5,z) (-x+.5,y+.5,z) (-y,-x,z+.5) (y,x,z+.5) \n 8 j ..2 (x,x+.5,0) (-x,-x+.5,0) (-x,x+.5,.5) (x,-x+.5,.5) (-x+.5,-x,.5) (x+.5,x,.5) (x+.5,-x,0) (-x+.5,x,0) \n 8 i .2. (x,0,.75) (-x,0,.75) (.5,x+.5,.25) (.5,-x+.5,.25) (-x+.5,.5,.75) (x+.5,.5,.75) (0,-x,.25) (0,x,.25) \n 8 h .2. (x,0,.25) (-x,0,.25) (.5,x+.5,.75) (.5,-x+.5,.75) (-x+.5,.5,.25) (x+.5,.5,.25) (0,-x,.75) (0,x,.75) \n 8 g 2.. (0,0,z) (.5,.5,z+.5) (0,0,-z+.5) (.5,.5,-z) (.5,.5,-z+.5) (0,0,-z) (.5,.5,z) (0,0,z+.5) \n 8 f 2.. (0,.5,z) (0,.5,z+.5) (0,.5,-z+.5) (0,.5,-z) (.5,0,-z+.5) (.5,0,-z) (.5,0,z) (.5,0,z+.5) \n 8 e -1 (.25,.25,.25) (.75,.75,.25) (.25,.75,.75) (.75,.25,.25) (.25,.75,.25) (.75,.75,.75) (.25,.25,.75) \n 4 d -4.. (0,0,0) (.5,.5,.5) (0,0,.5) (.5,.5,0) \n 4 c 2.22 (0,.5,0) (0,.5,.5) (.5,0,.5) (.5,0,0) \n 4 b 222. (0,0,.25) (.5,.5,.75) (.5,.5,.25) (0,0,.75) \n 4 a 222. (0,.5,.25) (0,.5,.75) (.5,0,.25) (.5,0,.75) \n ";
  }

  //DX20180723 - added origin choice 2
  string sg134 = "";
  if(axis_cell==setting_1){
    sg134 = "Space Group 134\n (0,0,0) \n 16 n 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z+.5) (y+.5,-x+.5,z+.5) (-x,y,-z) (x,-y,-z) (y+.5,x+.5,-z+.5) (-y+.5,-x+.5,-z+.5) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (y,-x,-z) (-y,x,-z) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) (-y,-x,z) (y,x,z) \n 8 m ..m (x,x,z) (-x,-x,z) (-x+.5,x+.5,z+.5) (x+.5,-x+.5,z+.5) (-x,x,-z) (x,-x,-z) (x+.5,x+.5,-z+.5) (-x+.5,-x+.5,-z+.5) \n 8 l ..2 (x,x+.5,.75) (-x,-x+.5,.75) (-x,x+.5,.25) (x,-x+.5,.25) (-x+.5,-x,.75) (x+.5,x,.75) (x+.5,-x,.25) (-x+.5,x,.25) \n 8 k ..2 (x,x+.5,.25) (-x,-x+.5,.25) (-x,x+.5,.75) (x,-x+.5,.75) (-x+.5,-x,.25) (x+.5,x,.25) (x+.5,-x,.75) (-x+.5,x,.75) \n 8 j .2. (x,0,.5) (-x,0,.5) (.5,x+.5,0) (.5,-x+.5,0) (-x+.5,.5,0) (x+.5,.5,0) (0,-x,.5) (0,x,.5) \n 8 i .2. (x,0,0) (-x,0,0) (.5,x+.5,.5) (.5,-x+.5,.5) (-x+.5,.5,.5) (x+.5,.5,.5) (0,-x,0) (0,x,0) \n 8 h 2.. (0,.5,z) (0,.5,z+.5) (0,.5,-z) (0,.5,-z+.5) (.5,0,-z+.5) (.5,0,-z) (.5,0,z+.5) (.5,0,z) \n 8 g 2.mm (0,0,z) (.5,.5,z+.5) (0,0,-z) (.5,.5,-z+.5) \n 4 f ..2/m (.75,.75,.75) (.25,.25,.75) (.75,.25,.25) (.25,.75,.25) \n 4 e ..2/m (.25,.25,.25) (.75,.75,.25) (.25,.75,.75) (.75,.25,.75) \n 4 d 2.22 (0,.5,.25) (0,.5,.75) (.5,0,.25) (.5,0,.75) \n 4 c 222. (0,.5,0) (0,.5,.5) (.5,0,.5) (.5,0,0) \n 2 b -42m (0,0,.5) (.5,.5,0) \n 2 a -42m (0,0,0) (.5,.5,.5) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg134 = "Space Group 134\n (0,0,0) \n 16 n 1 (x,y,z) (-x+.5,-y+.5,z) (-y+.5,x,z+.5) (y,-x+.5,z+.5) (-x+.5,y,-z+.5) (x,-y+.5,-z+.5) (y,x,-z) (-y+.5,-x+.5,-z) (-x,-y,-z) (x+.5,y+.5,-z) (y+.5,-x,-z+.5) (-y,x+.5,-z+.5) (x+.5,-y,z+.5) (-x,y+.5,z+.5) (-y,-x,z) (y+.5,x+.5,z) \n 8 m ..m (x,-x,z) (-x+.5,x+.5,z) (x+.5,x,z+.5) (-x,-x+.5,z+.5) (-x+.5,-x,-z+.5) (x,x+.5,-z+.5) (-x,x,-z) (x+.5,-x+.5,-z) \n 8 l ..2 (x,x,.5) (-x+.5,-x+.5,.5) (-x+.5,x,0) (x,-x+.5,0) (-x,-x,.5) (x+.5,x+.5,.5) (x+.5,-x,0) (-x,x+.5,0) \n 8 k ..2 (x,x,0) (-x+.5,-x+.5,0) (-x+.5,x,.5) (x,-x+.5,.5) (-x,-x,0) (x+.5,x+.5,0) (x+.5,-x,.5) (-x,x+.5,.5) \n 8 j .2. (x,.25,.25) (-x+.5,.25,.25) (.25,x,.75) (.25,-x+.5,.75) (-x,.75,.75) (x+.5,.75,.75) (.75,-x,.25) (.75,x+.5,.25) \n 8 i .2. (x,.25,.75) (-x+.5,.25,.75) (.25,x,.25) (.25,-x+.5,.25) (-x,.75,.25) (x+.5,.75,.25) (.75,-x,.75) (.75,x+.5,.75) \n 8 h 2.. (.25,.25,z) (.25,.25,z+.5) (.25,.25,-z+.5) (.25,.25,-z) (.75,.75,-z) (.75,.75,-z+.5) (.75,.75,z+.5) (.75,.75,z) \n 8 g 2.mm (.75,.25,z) (.25,.75,z+.5) (.75,.25,-z+.5) (.25,.75,-z) \n 4 f ..2/m (0,0,0) (.5,.5,0) (.5,0,.5) (0,.5,.5) \n 4 e ..2/m (0,0,.5) (.5,.5,.5) (.5,0,0) (0,.5,0) \n 4 d 2.22 (.25,.25,0) (.25,.25,.5) (.75,.75,0) (.75,.75,.5) \n 4 c 222. (.25,.25,.25) (.25,.25,.75) (.75,.75,.75) (.75,.75,.25) \n 2 b -42m (.75,.25,.25) (.25,.75,.75) \n 2 a -42m (.25,.75,.25) (.75,.25,.75) \n ";
  }
  else {
    sg134 = "Space Group 134\n (0,0,0) \n 16 n 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z+.5) (y+.5,-x+.5,z+.5) (-x,y,-z) (x,-y,-z) (y+.5,x+.5,-z+.5) (-y+.5,-x+.5,-z+.5) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (y,-x,-z) (-y,x,-z) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) (-y,-x,z) (y,x,z) \n 8 m ..m (x,x,z) (-x,-x,z) (-x+.5,x+.5,z+.5) (x+.5,-x+.5,z+.5) (-x,x,-z) (x,-x,-z) (x+.5,x+.5,-z+.5) (-x+.5,-x+.5,-z+.5) \n 8 l ..2 (x,x+.5,.75) (-x,-x+.5,.75) (-x,x+.5,.25) (x,-x+.5,.25) (-x+.5,-x,.75) (x+.5,x,.75) (x+.5,-x,.25) (-x+.5,x,.25) \n 8 k ..2 (x,x+.5,.25) (-x,-x+.5,.25) (-x,x+.5,.75) (x,-x+.5,.75) (-x+.5,-x,.25) (x+.5,x,.25) (x+.5,-x,.75) (-x+.5,x,.75) \n 8 j .2. (x,0,.5) (-x,0,.5) (.5,x+.5,0) (.5,-x+.5,0) (-x+.5,.5,0) (x+.5,.5,0) (0,-x,.5) (0,x,.5) \n 8 i .2. (x,0,0) (-x,0,0) (.5,x+.5,.5) (.5,-x+.5,.5) (-x+.5,.5,.5) (x+.5,.5,.5) (0,-x,0) (0,x,0) \n 8 h 2.. (0,.5,z) (0,.5,z+.5) (0,.5,-z) (0,.5,-z+.5) (.5,0,-z+.5) (.5,0,-z) (.5,0,z+.5) (.5,0,z) \n 8 g 2.mm (0,0,z) (.5,.5,z+.5) (0,0,-z) (.5,.5,-z+.5) \n 4 f ..2/m (.75,.75,.75) (.25,.25,.75) (.75,.25,.25) (.25,.75,.25) \n 4 e ..2/m (.25,.25,.25) (.75,.75,.25) (.25,.75,.75) (.75,.25,.75) \n 4 d 2.22 (0,.5,.25) (0,.5,.75) (.5,0,.25) (.5,0,.75) \n 4 c 222. (0,.5,0) (0,.5,.5) (.5,0,.5) (.5,0,0) \n 2 b -42m (0,0,.5) (.5,.5,0) \n 2 a -42m (0,0,0) (.5,.5,.5) \n ";
  }

  string sg135 = "Space Group 135\n (0,0,0) \n 16 i 1 (x,y,z) (-x,-y,z) (-y,x,z+.5) (y,-x,z+.5) (-x+.5,y+.5,-z) (x+.5,-y+.5,-z) (y+.5,x+.5,-z+.5) (-y+.5,-x+.5,-z+.5) (-x,-y,-z) (x,y,-z) (y,-x,-z+.5) (-y,x,-z+.5) (x+.5,-y+.5,z) (-x+.5,y+.5,z) (-y+.5,-x+.5,z+.5) (y+.5,x+.5,z+.5) \n 8 h m.. (x,y,0) (-x,-y,0) (-y,x,.5) (y,-x,.5) (-x+.5,y+.5,0) (x+.5,-y+.5,0) (y+.5,x+.5,.5) (-y+.5,-x+.5,.5) \n 8 g ..2 (x,x+.5,.25) (-x,-x+.5,.25) (-x+.5,x,.75) (x+.5,-x,.75) (-x,-x+.5,.75) (x,x+.5,.75) (x+.5,-x,.25) (-x+.5,x,.25) \n 8 f 2.. (0,.5,z) (.5,0,z+.5) (.5,0,-z) (0,.5,-z+.5) (0,.5,-z) (.5,0,-z+.5) (.5,0,z) (0,.5,z+.5) \n 8 e 2.. (0,0,z) (0,0,z+.5) (.5,.5,-z) (.5,.5,-z+.5) (0,0,-z) (0,0,-z+.5) (.5,.5,z) (.5,.5,z+.5) \n 4 d 2.22 (0,.5,.25) (.5,0,.75) (0,.5,.75) (.5,0,.25) \n 4 c 2/m.. (0,.5,0) (.5,0,.5) (.5,0,0) (0,.5,.5) \n 4 b -4.. (0,0,.25) (0,0,.75) (.5,.5,.75) (.5,.5,.25) \n 4 a 2/m.. (0,0,0) (0,0,.5) (.5,.5,0) (.5,.5,.5) \n ";

  string sg136 = "Space Group 136\n (0,0,0) \n 16 k 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z+.5) (y+.5,-x+.5,z+.5) (-x+.5,y+.5,-z+.5) (x+.5,-y+.5,-z+.5) (y,x,-z) (-y,-x,-z) (-x,-y,-z) (x,y,-z) (y+.5,-x+.5,-z+.5) (-y+.5,x+.5,-z+.5) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) (-y,-x,z) (y,x,z) \n 8 j ..m (x,x,z) (-x,-x,z) (-x+.5,x+.5,z+.5) (x+.5,-x+.5,z+.5) (-x+.5,x+.5,-z+.5) (x+.5,-x+.5,-z+.5) (x,x,-z) (-x,-x,-z) \n 8 i m.. (x,y,0) (-x,-y,0) (-y+.5,x+.5,.5) (y+.5,-x+.5,.5) (-x+.5,y+.5,.5) (x+.5,-y+.5,.5) (y,x,0) (-y,-x,0) \n 8 h 2.. (0,.5,z) (0,.5,z+.5) (.5,0,-z+.5) (.5,0,-z) (0,.5,-z) (0,.5,-z+.5) (.5,0,z+.5) (.5,0,z) \n 4 g m.2m (x,-x,0) (-x,x,0) (x+.5,x+.5,.5) (-x+.5,-x+.5,.5) \n 4 f m.2m (x,x,0) (-x,-x,0) (-x+.5,x+.5,.5) (x+.5,-x+.5,.5) \n 4 e 2.mm (0,0,z) (.5,.5,z+.5) (.5,.5,-z+.5) (0,0,-z) \n 4 d -4.. (0,.5,.25) (0,.5,.75) (.5,0,.25) (.5,0,.75) \n 4 c 2/m.. (0,.5,0) (0,.5,.5) (.5,0,.5) (.5,0,0) \n 2 b m.mm (0,0,.5) (.5,.5,0) \n 2 a m.mm (0,0,0) (.5,.5,.5) \n ";

  //DX20180723 - added origin choice 2
  string sg137 = "";
  if(axis_cell==setting_1){
    sg137 = "Space Group 137\n (0,0,0) \n 16 h 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z+.5) (y+.5,-x+.5,z+.5) (-x+.5,y+.5,-z+.5) (x+.5,-y+.5,-z+.5) (y,x,-z) (-y,-x,-z) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (y,-x,-z) (-y,x,-z) (x,-y,z) (-x,y,z) (-y+.5,-x+.5,z+.5) (y+.5,x+.5,z+.5) \n 8 g .m. (0,y,z) (0,-y,z) (-y+.5,.5,z+.5) (y+.5,.5,z+.5) (.5,y+.5,-z+.5) (.5,-y+.5,-z+.5) (y,0,-z) (-y,0,-z) \n 8 f ..2 (x,x,0) (-x,-x,0) (-x+.5,x+.5,.5) (x+.5,-x+.5,.5) (-x+.5,-x+.5,.5) (x+.5,x+.5,.5) (x,-x,0) (-x,x,0) \n 8 e -1 (.25,.25,.25) (.75,.75,.25) (.25,.75,.75) (.75,.25,.75) (.25,.75,.25) (.75,.25,.25) (.25,.25,.75) (.75,.75,.75) \n 4 d 2mm. (0,.5,z) (0,.5,z+.5) (.5,0,-z+.5) (.5,0,-z) \n 4 c 2mm. (0,0,z) (.5,.5,z+.5) (.5,.5,-z+.5) (0,0,-z) \n 2 b -4m2 (0,0,.5) (.5,.5,0) \n 2 a -4m2 (0,0,0) (.5,.5,.5) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg137 = "Space Group 137\n (0,0,0) \n 16 h 1 (x,y,z) (-x+.5,-y+.5,z) (-y+.5,x,z+.5) (y,-x+.5,z+.5) (-x,y+.5,-z) (x+.5,-y,-z) (y+.5,x+.5,-z+.5) (-y,-x,-z+.5) (-x,-y,-z) (x+.5,y+.5,-z) (y+.5,-x,-z+.5) (-y,x+.5,-z+.5) (x,-y+.5,z) (-x+.5,y,z) (-y+.5,-x+.5,z+.5) (y,x,z+.5) \n 8 g .m. (.25,y,z) (.25,-y+.5,z) (-y+.5,.25,z+.5) (y,.25,z+.5) (.75,y+.5,-z) (.75,-y,-z) (y+.5,.75,-z+.5) (-y,.75,-z+.5) \n 8 f ..2 (x,-x,.25) (-x+.5,x+.5,.25) (x+.5,x,.75) (-x,-x+.5,.75) (-x,x,.75) (x+.5,-x+.5,.75) (-x+.5,-x,.25) (x,x+.5,.25) \n 8 e -1 (0,0,0) (.5,.5,0) (.5,0,.5) (0,.5,.5) (0,.5,0) (.5,0,0) (.5,.5,.5) (0,0,.5) \n 4 d 2mm. (.25,.25,z) (.25,.25,z+.5) (.75,.75,-z) (.75,.75,-z+.5) \n 4 c 2mm. (.75,.25,z) (.25,.75,z+.5) (.25,.75,-z) (.75,.25,-z+.5) \n 2 b -4m2 (.75,.25,.25) (.25,.75,.75) \n 2 a -4m2 (.75,.25,.75) (.25,.75,.25) \n ";
  }
  else {
    sg137 = "Space Group 137\n (0,0,0) \n 16 h 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z+.5) (y+.5,-x+.5,z+.5) (-x+.5,y+.5,-z+.5) (x+.5,-y+.5,-z+.5) (y,x,-z) (-y,-x,-z) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (y,-x,-z) (-y,x,-z) (x,-y,z) (-x,y,z) (-y+.5,-x+.5,z+.5) (y+.5,x+.5,z+.5) \n 8 g .m. (0,y,z) (0,-y,z) (-y+.5,.5,z+.5) (y+.5,.5,z+.5) (.5,y+.5,-z+.5) (.5,-y+.5,-z+.5) (y,0,-z) (-y,0,-z) \n 8 f ..2 (x,x,0) (-x,-x,0) (-x+.5,x+.5,.5) (x+.5,-x+.5,.5) (-x+.5,-x+.5,.5) (x+.5,x+.5,.5) (x,-x,0) (-x,x,0) \n 8 e -1 (.25,.25,.25) (.75,.75,.25) (.25,.75,.75) (.75,.25,.75) (.25,.75,.25) (.75,.25,.25) (.25,.25,.75) (.75,.75,.75) \n 4 d 2mm. (0,.5,z) (0,.5,z+.5) (.5,0,-z+.5) (.5,0,-z) \n 4 c 2mm. (0,0,z) (.5,.5,z+.5) (.5,.5,-z+.5) (0,0,-z) \n 2 b -4m2 (0,0,.5) (.5,.5,0) \n 2 a -4m2 (0,0,0) (.5,.5,.5) \n ";
  }

  //DX20180723 - added origin choice 2
  string sg138 = "";
  if(axis_cell==setting_1){
    sg138 = "Space Group 138\n (0,0,0) \n 16 j 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z+.5) (y+.5,-x+.5,z+.5) (-x+.5,y+.5,-z) (x+.5,-y+.5,-z) (y,x,-z+.5) (-y,-x,-z+.5) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (y,-x,-z) (-y,x,-z) (x,-y,z+.5) (-x,y,z+.5) (-y+.5,-x+.5,z) (y+.5,x+.5,z) \n 8 i ..m (x,x+.5,z) (-x,-x+.5,z) (-x,x+.5,z+.5) (x,-x+.5,z+.5) (-x+.5,x,-z) (x+.5,-x,-z) (x+.5,x,-z+.5) (-x+.5,-x,-z+.5) \n 8 h ..2 (x,x,.75) (-x,-x,.75) (-x+.5,x+.5,.25) (x+.5,-x+.5,.25) (-x+.5,-x+.5,.75) (x+.5,x+.5,.75) (x,-x,.25) (-x,x,.25) \n 8 g ..2 (x,x,.25) (-x,-x,.25) (-x+.5,x+.5,.75) (x+.5,-x+.5,.75) (-x+.5,-x+.5,.25) (x+.5,x+.5,.25) (x,-x,.75) (-x,x,.75) \n 8 f 2.. (0,0,z) (.5,.5,z+.5) (.5,.5,-z) (0,0,-z+.5) (.5,.5,-z+.5) (0,0,-z) (0,0,z+.5) (.5,.5,z) \n 4 e 2.mm (0,.5,z) (0,.5,z+.5) (.5,0,-z) (.5,0,-z+.5) \n 4 d ..2/m (.25,.25,.75) (.75,.75,.75) (.25,.75,.25) (.75,.25,.25)\n 4 c ..2/m (.25,.25,.25) (.75,.75,.25) (.25,.75,.75) (.75,.25,.75) \n 4 b -4.. (0,0,0) (.5,.5,.5) (.5,.5,0) (0,0,.5) \n 4 a 2.22 (0,0,.25) (.5,.5,.75) (.5,.5,.25) (0,0,.75) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg138 = "Space Group 138\n (0,0,0) \n 16 j 1 (x,y,z) (-x+.5,-y+.5,z) (-y+.5,x,z+.5) (y,-x+.5,z+.5) (-x,y+.5,-z+.5) (x+.5,-y,-z+.5) (y+.5,x+.5,-z) (-y,-x,-z) (-x,-y,-z) (x+.5,y+.5,-z) (y+.5,-x,-z+.5) (-y,x+.5,-z+.5) (x,-y+.5,z+.5) (-x+.5,y,z+.5) (-y+.5,-x+.5,z) (y,x,z) \n 8 i ..m (x,x,z) (-x+.5,-x+.5,z) (-x+.5,x,z+.5) (x,-x+.5,z+.5) (-x,x+.5,-z+.5) (x+.5,-x,-z+.5) (x+.5,x+.5,-z) (-x,-x,-z) \n 8 h ..2 (x,-x,0) (-x+.5,x+.5,0) (x+.5,x,.5) (-x,-x+.5,.5) (-x,x,0) (x+.5,-x+.5,0) (-x+.5,-x,.5) (x,x+.5,.5) \n 8 g ..2 (x,-x,.5) (-x+.5,x+.5,.5) (x+.5,x,0) (-x,-x+.5,0) (-x,x,.5) (x+.5,-x+.5,.5) (-x+.5,-x,0) (x,x+.5,0) \n 8 f 2.. (.75,.25,z) (.25,.75,z+.5) (.25,.75,-z+.5) (.75,.25,-z) (.25,.75,-z) (.75,.25,-z+.5) (.75,.25,z+.5) (.25,.75,z) \n 4 e 2.mm (.25,.25,z) (.25,.25,z+.5) (.75,.75,-z+.5) (.75,.75,-z) \n 4 d ..2/m (0,0,0) (.5,.5,0) (.5,0,.5) (0,.5,.5)\n 4 c ..2/m (0,0,.5) (.5,.5,.5) (.5,0,0) (0,.5,0) \n 4 b -4.. (.75,.25,.75) (.25,.75,.25) (.25,.75,.75) (.75,.25,.25) \n 4 a 2.22 (.75,.25,0) (.25,.75,.5) (.25,.75,0) (.75,.25,.5) \n "; //DX20180806 - second 8 h is (-x+.5,x+.5,0)  not (-x+.5,-x+.5,0)
  }
  else {
    sg138 = "Space Group 138\n (0,0,0) \n 16 j 1 (x,y,z) (-x,-y,z) (-y+.5,x+.5,z+.5) (y+.5,-x+.5,z+.5) (-x+.5,y+.5,-z) (x+.5,-y+.5,-z) (y,x,-z+.5) (-y,-x,-z+.5) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (y,-x,-z) (-y,x,-z) (x,-y,z+.5) (-x,y,z+.5) (-y+.5,-x+.5,z) (y+.5,x+.5,z) \n 8 i ..m (x,x+.5,z) (-x,-x+.5,z) (-x,x+.5,z+.5) (x,-x+.5,z+.5) (-x+.5,x,-z) (x+.5,-x,-z) (x+.5,x,-z+.5) (-x+.5,-x,-z+.5) \n 8 h ..2 (x,x,.75) (-x,-x,.75) (-x+.5,x+.5,.25) (x+.5,-x+.5,.25) (-x+.5,-x+.5,.75) (x+.5,x+.5,.75) (x,-x,.25) (-x,x,.25) \n 8 g ..2 (x,x,.25) (-x,-x,.25) (-x+.5,x+.5,.75) (x+.5,-x+.5,.75) (-x+.5,-x+.5,.25) (x+.5,x+.5,.25) (x,-x,.75) (-x,x,.75) \n 8 f 2.. (0,0,z) (.5,.5,z+.5) (.5,.5,-z) (0,0,-z+.5) (.5,.5,-z+.5) (0,0,-z) (0,0,z+.5) (.5,.5,z) \n 4 e 2.mm (0,.5,z) (0,.5,z+.5) (.5,0,-z) (.5,0,-z+.5) \n 4 d ..2/m (.25,.25,.75) (.75,.75,.75) (.25,.75,.25) (.75,.25,.25)\n 4 c ..2/m (.25,.25,.25) (.75,.75,.25) (.25,.75,.75) (.75,.25,.75) \n 4 b -4.. (0,0,0) (.5,.5,.5) (.5,.5,0) (0,0,.5) \n 4 a 2.22 (0,0,.25) (.5,.5,.75) (.5,.5,.25) (0,0,.75) \n ";
  }

  string sg139 = "Space Group 139\n (0,0,0) (.5,.5,.5) \n 32 o 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (-x,y,-z) (x,-y,-z) (y,x,-z) (-y,-x,-z) (-x,-y,-z) (x,y,-z) (y,-x,-z) (-y,x,-z) (x,-y,z) (-x,y,z) (-y,-x,z) (y,x,z) \n 16 n .m. (0,y,z) (0,-y,z) (-y,0,z) (y,0,z) (0,y,-z) (0,-y,-z) (y,0,-z) (-y,0,-z) \n 16 m ..m (x,x,z) (-x,-x,z) (-x,x,z) (x,-x,z) (-x,x,-z) (x,-x,-z) (x,x,-z) (-x,-x,-z) \n 16 l m.. (x,y,0) (-x,-y,0) (-y,x,0) (y,-x,0) (-x,y,0) (x,-y,0) (y,x,0) (-y,-x,0) \n 16 k ..2 (x,x+.5,.25) (-x,-x+.5,.25) (-x+.5,x,.25) (x+.5,-x,.25) (-x,-x+.5,.75) (x,x+.5,.75) (x+.5,-x,.75) (-x+.5,x,.75) \n 8 j m2m. (x,.5,0) (-x,.5,0) (.5,x,0) (.5,-x,0) \n 8 i m2m. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) \n 8 h m.2m (x,x,0) (-x,-x,0) (-x,x,0) (x,-x,0) \n 8 g 2mm. (0,.5,z) (.5,0,z) (0,.5,-z) (.5,0,-z) \n 8 f ..2/m (.25,.25,.25) (.75,.75,.25) (.75,.25,.25) (.25,.75,.25) \n 4 e 4mm (0,0,z) (0,0,-z) \n 4 d -4m2 (0,.5,.25) (.5,.0,.25) \n 4 c mmm. (0,.5,0) (.5,0,0) \n 2 b 4/mmm (0,0,.5) \n 2 a 4/mmm (0,0,0) \n ";

  string sg140 = "Space Group 140\n (0,0,0) (.5,.5,.5) \n 32 m 1 (x,y,z) (-x,-y,z) (-y,x,z) (y,-x,z) (-x,y,-z+.5) (x,-y,-z+.5) (y,x,-z+.5) (-y,-x,-z+.5) (-x,-y,-z) (x,y,-z) (y,-x,-z) (-y,x,-z) (x,-y,z+.5) (-x,y,z+.5) (-y,-x,z+.5) (y,x,z+.5) \n 16 l ..m (x,x+.5,z) (-x,-x+.5,z) (-x+.5,x,z) (x+.5,-x,z) (-x,x+.5,-z+.5) (x,-x+.5,-z+.5) (x+.5,x,-z+.5) (-x+.5,-x,-z+.5) \n 16 k m.. (x,y,0) (-x,-y,0) (-y,x,0) (y,-x,0) (-x,y,.5) (x,-y,.5) (y,x,.5) (-y,-x,.5) \n 16 j .2. (x,0,.25) (-x,0,.25) (0,x,.25) (0,-x,.25) (-x,0,.75) (x,0,.75) (0,-x,.75) (0,x,.75) \n 16 i ..2 (x,x,.25) (-x,-x,.25) (-x,x,.25) (x,-x,.25) (-x,-x,.75) (x,x,.75) (x,-x,.75) (-x,x,.75) \n 8 h m.2m (x,x+.5,0) (-x,-x+.5,0) (-x+.5,x,0) (x+.5,-x,0) \n 8 g 2.mm (0,.5,z) (.5,0,z) (0,.5,-z+.5) (.5,0,-z+.5) \n 8 f 4.. (0,0,z) (0,0,-z+.5) (0,0,-z) (0,0,z+.5) \n 8 e ..2/m (.25,.25,.25) (.75,.75,.25) (.75,.25,.25) (.25,.75,.25) \n 4 d m.mm (0,.5,0) (.5,0,0) \n 4 c 4/m.. (0,0,0) (0,0,.5) \n 4 b -42m (0,.5,.25) (.5,0,.25) \n 4 a 422 (0,0,.25) (0,0,.75) \n ";

  //DX20180723 - added origin choice 2
  string sg141 = "";
  if(axis_cell==setting_1){
    sg141 = "Space Group 141\n (0,0,0) (.5,.5,.5) \n 32 i 1 (x,y,z) (-x+.5,-y+.5,z+.5) (-y,x+.5,z+.25) (y+.5,-x,z+.75) (-x+.5,y,-z+.75) (x,-y+.5,-z+.25) (y+.5,x+.5,-z+.5) (-y,-x,-z) (-x,-y+.5,-z+.25) (x+.5,y,-z+.75) (y,-x,-z) (-y+.5,x+.5,-z+.5) (x+.5,-y+.5,z+.5) (-x,y,z) (-y+.5,-x,z+.75) (y,x+.5,z+.25) \n 16 h .m. (0,y,z) (.5,-y+.5,z+.5) (-y,.5,z+.25) (y+.5,0,z+.75) (.5,y,-z+.75) (0,-y+.5,-z+.25) (y+.5,.5,-z+.5) (-y,0,-z) \n 16 g ..2 (x,x,0) (-x+.5,-x+.5,.5) (-x,x+.5,.25) (x+.5,-x,.75) (-x,-x+.5,.25) (x+.5,x,.75) (x,-x,0) (-x+.5,x+.5,.5) \n 16 f .2. (x,.25,.125) (-x+.5,.25,.625) (.75,x+.5,.375) (.75,-x,.875) (-x,.25,.125) (x+.5,.25,.625) (.25,-x,.875) (.25,x+.5,.375) \n 8 e 2mm. (0,0,z) (0,.5,z+.25) (.5,0,-z+.75) (.5,.5,-z+.5) \n 8 d .2/m. (0,.25,.625) (.5,.25,.125) (.75,.5,.875) (.75,0,.375) \n 8 c .2/m. (0,.25,.125) (.5,.25,.625) (.75,.5,.375) (.75,0,.875) \n 4 b -4m2 (0,0,.5) (0,.5,.75) \n 4 a -4m2 (0,0,0) (0,.5,.25) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg141 = "Space Group 141\n (0,0,0) (.5,.5,.5) \n 32 i 1 (x,y,z) (-x+.5,-y,z+.5) (-y+.25,x+.75,z+.25) (y+.25,-x+.25,z+.75) (-x+.5,y,-z+.5) (x,-y,-z) (y+.25,x+.75,-z+.25) (-y+.25,-x+.25,-z+.75) (-x,-y,-z) (x+.5,y,-z+.5) (y+.75,-x+.25,-z+.75) (-y+.75,x+.75,-z+.25) (x+.5,-y,z+.5) (-x,y,z) (-y+.75,-x+.25,z+.75) (y+.75,x+.75,z+.25) \n 16 h .m. (0,y,z) (.5,-y,z+.5) (-y+.25,.75,z+.25) (y+.25,.25,z+.75) (.5,y,-z+.5) (0,-y,-z) (y+.25,.75,-z+.25) (-y+.25,.25,-z+.75) \n 16 g ..2 (x,x+.25,.875) (-x+.5,-x+.75,.375) (-x,x+.75,.125) (x+.5,-x+.25,.625) (-x,-x+.75,.125) (x+.5,x+.25,.625) (x,-x+.25,.875) (-x+.5,x+.75,.375) \n 16 f .2. (x,0,0) (-x+.5,0,.5) (.25,x+.75,.25) (.25,-x+.25,.75) (-x,0,0) (x+.5,0,.5) (.75,-x+.25,.75) (.75,x+.75,.25) \n 8 e 2mm. (0,.25,z) (0,.75,z+.25) (.5,.25,-z+.5) (.5,.75,-z+.25) \n 8 d .2/m. (0,0,.5) (.5,0,0) (.25,.75,.75) (.25,.25,.25) \n 8 c .2/m. (0,0,0) (.5,0,.5) (.25,.75,.25) (.25,.25,.75) \n 4 b -4m2 (0,.25,.375) (0,.75,.625) \n 4 a -4m2 (0,.75,.125) (.5,.75,.375) \n "; //DX20180806 - changed 4 b -4m2 (0,.25,.375) (0,.625,.75) to 4 b -4m2 (0,.25,.375) (0,.75,.625)
  }
  else {
    sg141 = "Space Group 141\n (0,0,0) (.5,.5,.5) \n 32 i 1 (x,y,z) (-x+.5,-y+.5,z+.5) (-y,x+.5,z+.25) (y+.5,-x,z+.75) (-x+.5,y,-z+.75) (x,-y+.5,-z+.25) (y+.5,x+.5,-z+.5) (-y,-x,-z) (-x,-y+.5,-z+.25) (x+.5,y,-z+.75) (y,-x,-z) (-y+.5,x+.5,-z+.5) (x+.5,-y+.5,z+.5) (-x,y,z) (-y+.5,-x,z+.75) (y,x+.5,z+.25) \n 16 h .m. (0,y,z) (.5,-y+.5,z+.5) (-y,.5,z+.25) (y+.5,0,z+.75) (.5,y,-z+.75) (0,-y+.5,-z+.25) (y+.5,.5,-z+.5) (-y,0,-z) \n 16 g ..2 (x,x,0) (-x+.5,-x+.5,.5) (-x,x+.5,.25) (x+.5,-x,.75) (-x,-x+.5,.25) (x+.5,x,.75) (x,-x,0) (-x+.5,x+.5,.5) \n 16 f .2. (x,.25,.125) (-x+.5,.25,.625) (.75,x+.5,.375) (.75,-x,.875) (-x,.25,.125) (x+.5,.25,.625) (.25,-x,.875) (.25,x+.5,.375) \n 8 e 2mm. (0,0,z) (0,.5,z+.25) (.5,0,-z+.75) (.5,.5,-z+.5) \n 8 d .2/m. (0,.25,.625) (.5,.25,.125) (.75,.5,.875) (.75,0,.375) \n 8 c .2/m. (0,.25,.125) (.5,.25,.625) (.75,.5,.375) (.75,0,.875) \n 4 b -4m2 (0,0,.5) (0,.5,.75) \n 4 a -4m2 (0,0,0) (0,.5,.25) \n ";
  }

  //DX20180723 - added origin choice 2
  string sg142 = "";
  if(axis_cell==setting_1){
    sg142 = "Space Group 142\n (0,0,0) (.5,.5,.5) \n 32 g 1 (x,y,z) (-x+.5,-y+.5,z+.5) (-y,x+.5,z+.25) (y+.5,-x,z+.75) (-x+.5,y,-z+.25) (x,-y+.5,-z+.75) (y+.5,x+.5,-z) (-y,-x,-z+.5) (-x,-y+.5,-z+.25) (x+.5,y,-z+.75) (y,-x,-z) (-y+.5,x+.5,-z+.5) (x+.5,-y+.5,z) (-x,y,z+.5) (-y+.5,-x,z+.25) (y,x+.5,z+.75) \n 16 f ..2 (x,x,.25) (-x+.5,-x+.5,.75) (-x,x+.5,.5) (x+.5,-x,0) (-x,-x+.5,0) (x+.5,x,.5) (x,-x,.75) (-x+.5,x+.5,.25) \n 16 e .2. (.25,y,.125) (.25,-y+.5,.625) (-y,.75,.375) (y+.5,.75,.875) (.75,-y+.5,.125) (.75,y,.625) (y,.75,.875) (-y+.5,.75,.375) \n 16 d 2.. (0,0,z) (0,.5,z+.25) (.5,0,-z+.25) (.5,.5,-z) (0,.5,-z+.25) (0,0,-z) (.5,.5,z) (.5,0,z+.25) \n 16 c -1 (0,.25,.125) (.5,.25,.625) (.75,.5,.375) (.75,0,.875) (.5,.25,.125) (0,.25,.625) (.75,.5,.875) (.75,0,.375) \n 8 b 2.22 (0,0,.25) (0,.5,.5) (0,.5,0) (0,0,.75) \n 8 a -4.. (0,0,0) (0,.5,.25) (.5,0,.25) (.5,.5,0) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg142 = "Space Group 142\n (0,0,0) (.5,.5,.5) \n 32 g 1 (x,y,z) (-x+.5,-y,z+.5) (-y+.25,x+.75,z+.25) (y+.25,-x+.25,z+.75) (-x+.5,y,-z) (x,-y,-z+.5) (y+.25,x+.75,-z+.75) (-y+.25,-x+.25,-z+.25) (-x,-y,-z) (x+.5,y,-z+.5) (y+.75,-x+.25,-z+.75) (-y+.75,x+.75,-z+.25) (x+.5,-y,z) (-x,y,z+.5) (-y+.75,-x+.25,z+.25) (y+.75,x+.75,z+.75) \n 16 f ..2 (x,x+.25,.125) (-x+.5,-x+.75,.625) (-x,x+.75,.375) (x+.5,-x+.25,.875) (-x,-x+.75,.875) (x+.5,x+.25,.375) (x,-x+.25,.625) (-x+.5,x+.75,.125) \n 16 e .2. (x,0,.25) (-x+.5,0,.75) (.25,x+.75,.5) (.25,-x+.25,0) (-x,0,.75) (x+.5,0,.25) (.75,-x+.25,.5) (.75,x+.75,0) \n 16 d 2.. (0,.25,z) (0,.75,z+.25) (.5,.25,-z) (.5,.75,-z+.75) (0,.75,-z) (0,.25,-z+.75) (.5,.75,z) (.5,.25,z+.25) \n 16 c -1 (0,0,0) (.5,0,.5) (.25,.75,.25) (.25,.25,.75) (.5,0,0) (0,0,.5) (.25,.75,.75) (.25,.25,.25) \n 8 b 2.22 (0,.25,.125) (0,.75,.375) (0,.75,.875) (0,.25,.625) \n 8 a -4.. (0,.25,.375) (0,.75,.625) (.5,.25,.625) (.5,.75,.375) \n "; //DX201806 - second from end in 16 e .2. is (.75,-x+.25,.5) not (.75,-x+.25,y,.5)
  }
  else {
    sg142 = "Space Group 142\n (0,0,0) (.5,.5,.5) \n 32 g 1 (x,y,z) (-x+.5,-y+.5,z+.5) (-y,x+.5,z+.25) (y+.5,-x,z+.75) (-x+.5,y,-z+.25) (x,-y+.5,-z+.75) (y+.5,x+.5,-z) (-y,-x,-z+.5) (-x,-y+.5,-z+.25) (x+.5,y,-z+.75) (y,-x,-z) (-y+.5,x+.5,-z+.5) (x+.5,-y+.5,z) (-x,y,z+.5) (-y+.5,-x,z+.25) (y,x+.5,z+.75) \n 16 f ..2 (x,x,.25) (-x+.5,-x+.5,.75) (-x,x+.5,.5) (x+.5,-x,0) (-x,-x+.5,0) (x+.5,x,.5) (x,-x,.75) (-x+.5,x+.5,.25) \n 16 e .2. (.25,y,.125) (.25,-y+.5,.625) (-y,.75,.375) (y+.5,.75,.875) (.75,-y+.5,.125) (.75,y,.625) (y,.75,.875) (-y+.5,.75,.375) \n 16 d 2.. (0,0,z) (0,.5,z+.25) (.5,0,-z+.25) (.5,.5,-z) (0,.5,-z+.25) (0,0,-z) (.5,.5,z) (.5,0,z+.25) \n 16 c -1 (0,.25,.125) (.5,.25,.625) (.75,.5,.375) (.75,0,.875) (.5,.25,.125) (0,.25,.625) (.75,.5,.875) (.75,0,.375) \n 8 b 2.22 (0,0,.25) (0,.5,.5) (0,.5,0) (0,0,.75) \n 8 a -4.. (0,0,0) (0,.5,.25) (.5,0,.25) (.5,.5,0) \n ";
  }

  string sg143 = "Space Group 143\n (0,0,0) \n 3 d 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) \n 1 c 3.. (2/3,1/3,z) \n 1 b 3.. (1/3,2/3,z) \n 1 a 3.. (0,0,z) \n ";

  string sg144 = "Space Group 144\n (0,0,0) \n 3 a 1 (x,y,z) (-y,x-y,z+1/3) (-x+y,-x,z+2/3) \n ";

  string sg145 = "Space Group 145\n (0,0,0) \n 3 a 1 (x,y,z) (-y,x-y,z+2/3) (-x+y,-x,z+1/3) \n ";

  //DX20180723 - added rhombohedral setting
  string sg146 = "";
  if(axis_cell==setting_2){
    sg146 = "Space Group 146\n (0,0,0) (2/3,1/3,1/3) (1/3,2/3,2/3) \n 9 b 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) \n 3 a 3. (0,0,z) \n ";
  }
  else if(axis_cell==setting_1 || axis_cell==setting_anrl){
    sg146 = "Space Group 146\n (0,0,0) \n 3 b 1 (x,y,z) (z,x,y) (y,z,x) \n 1 a 3. (x,x,x) \n ";
  }
  else {
    sg146 = "Space Group 146\n (0,0,0) (2/3,1/3,1/3) (1/3,2/3,2/3) \n 9 b 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) \n 3 a 3. (0,0,z) \n ";
  }
  //string sg146R = "Space Group 146R\n (0,0,0) \n 3 b 1 (x,y,z) (z,x,y) (y,z,x) \n 1 a 3. (x,x,x) \n ";

  string sg147 = "Space Group 147\n (0,0,0) \n 6 g 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) \n 3 f -1 (.5,0,.5) (0,.5,.5) (.5,.5,.5) \n 3 e -1 (.5,0,0) (0,.5,0) (.5,.5,0) \n 2 d 3.. (1/3,2/3,z) (2/3,1/3,-z) \n 2 c 3.. (0,0,z) (0,0,-z) \n 1 b -3.. (0,0,.5) \n 1 a -3.. (0,0,0) \n ";

  //DX20180723 - added rhombohedral setting
  string sg148 = "";
  if(axis_cell==setting_2){
    sg148 = "Space Group 148\n (0,0,0) (2/3,1/3,1/3) (1/3,2/3,2/3) \n 18 f 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) \n 9 e -1 (.5,0,0) (0,.5,0) (.5,.5,0) \n 9 d -1 (.5,0,.5) (0,.5,.5) (.5,.5,.5) \n 6 c 3. (0,0,z) (0,0,-z) \n 3 b -3. (0,0,.5) \n 3 a -3. (0,0,0) \n ";
  }
  else if(axis_cell==setting_1 || axis_cell==setting_anrl){
    sg148 = "Space Group 148\n (0,0,0) \n 6 f 1 (x,y,z) (z,x,y) (y,z,x) (-x,-y,-z) (-z,-x,-y) (-y,-z,-x) \n 3 e -1 (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 3 d -1 (.5,0,0) (0,.5,0) (0,0,.5) \n 2 c 3. (x,x,x) (-x,-x,-x) \n 1 b -3. (.5,.5,.5) \n 1 a -3. (0,0,0) \n ";
  }
  else {
    sg148 = "Space Group 148\n (0,0,0) (2/3,1/3,1/3) (1/3,2/3,2/3) \n 18 f 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) \n 9 e -1 (.5,0,0) (0,.5,0) (.5,.5,0) \n 9 d -1 (.5,0,.5) (0,.5,.5) (.5,.5,.5) \n 6 c 3. (0,0,z) (0,0,-z) \n 3 b -3. (0,0,.5) \n 3 a -3. (0,0,0) \n ";
  }
  //string sg148R = "Space Group 148R\n (0,0,0) \n 6 f 1 (x,y,z) (z,x,y) (y,z,x) (-x,-y,-z) (-z,-x,-y) (-y,-z,-x) \n 3 e -1 (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 3 d -1 (.5,0,0) (0,.5,0) (0,0,.5) \n 2 c 3. (x,x,x) (-x,-x,-x) \n 1 b -3. (.5,.5,.5) \n 1 a -3. (0,0,0) \n ";

  string sg149 = "Space Group 149\n (0,0,0) \n 6 l 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-y,-x,-z) (-x+y,y,-z) (x,x-y,-z) \n 3 k ..2 (x,-x,.5) (x,2*x,.5) (-2*x,-x,.5) \n 3 j ..2 (x,-x,0) (x,2*x,0) (-2*x,-x,0) \n 2 i 3.. (2/3,1/3,z) (2/3,1/3,-z) \n 2 h 3.. (1/3,2/3,z) (1/3,2/3,-z) \n 2 g 3.. (0,0,z) (0,0,-z) \n 1 f 3.2 (2/3,1/3,.5) \n 1 e 3.2 (2/3,1/3,0) \n 1 d 3.2 (1/3,2/3,.5) \n 1 c 3.2 (1/3,2/3,0) \n 1 b 3.2 (0,0,.5) \n 1 a 3.2 (0,0,0) \n ";

  string sg150 = "Space Group 150\n (0,0,0) \n 6 g 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (y,x,-z) (x-y,-y,-z) (-x,-x+y,-z) \n 3 f .2. (x,0,.5) (0,x,.5) (-x,-x,.5) \n 3 e .2. (x,0,0) (0,x,0) (-x,-x,0) \n 2 d 3.. (1/3,2/3,z) (2/3,1/3,-z) \n 2 c 3.. (0,0,z) (0,0,-z) \n 1 b 32. (0,0,.5) \n 1 a 32. (0,0,0) \n ";

  string sg151 = "Space Group 151\n (0,0,0) \n 6 c 1 (x,y,z) (-y,x-y,z+1/3) (-x+y,-x,z+2/3) (-y,-x,-z+2/3) (-x+y,y,-z+1/3) (x,x-y,-z) \n 3 b ..2 (x,-x,5/6) (x,2*x,1/6) (-2*x,-x,.5) \n 3 a ..2 (x,-x,1/3) (x,2*x,2/3) (-2*x,-x,0) \n ";

  string sg152 = "Space Group 152\n (0,0,0) \n 6 c 1 (x,y,z) (-y,x-y,z+1/3) (-x+y,-x,z+2/3) (y,x,-z) (x-y,-y,-z+2/3) (-x,-x+y,-z+1/3) \n 3 b .2. (x,0,5/6) (0,x,1/6) (-x,-x,.5) \n 3 a .2. (x,0,1/3) (0,x,2/3) (-x,-x,0) \n ";

  string sg153 = "Space Group 153\n (0,0,0) \n 6 c 1 (x,y,z) (-y,x-y,z+2/3) (-x+y,-x,z+1/3) (-y,-x,-z+1/3) (-x+y,y,-z+2/3) (x,x-y,-z) \n 3 b ..2 (x,-x,1/6) (x,2*x,5/6) (-2*x,-x,.5) \n 3 a ..2 (x,-x,2/3) (x,2*x,1/3) (-2*x,-x,0) \n ";

  string sg154 = "Space Group 154\n (0,0,0) \n 6 c 1 (x,y,z) (-y,x-y,z+2/3) (-x+y,-x,z+1/3) (y,x,-z) (x-y,-y,-z+1/3) (-x,-x+y,-z+2/3) \n 3 b .2. (x,0,1/6) (0,x,5/6) (-x,-x,.5) \n 3 a .2. (x,0,2/3) (0,x,1/3) (-x,-x,0) \n ";

  //DX20180723 - added rhombohedral setting
  string sg155 = "";
  if(axis_cell==setting_2){
    sg155 = "Space Group 155\n (0,0,0) (2/3,1/3,1/3) (1/3,2/3,2/3) \n 18 f 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (y,x,-z) (x-y,-y,-z) (-x,-x+y,-z) \n 9 e .2 (x,0,.5) (0,x,.5) (-x,-x,.5) \n 9 d .2 (x,0,0) (0,x,0) (-x,-x,0) \n 6 c 3. (0,0,z) (0,0,-z) \n 3 b 32 (0,0,.5) \n 3 a 32 (0,0,0) \n ";
  }
  else if(axis_cell==setting_1 || axis_cell==setting_anrl){
    sg155 = "Space Group 155\n (0,0,0) \n 6 f 1 (x,y,z) (z,x,y) (y,z,x) (-z,-y,-x) (-y,-x,-z) (-x,-z,-y) \n 3 e .2 (.5,y,-y) (-y,.5,y) (y,-y,.5) \n 3 d .2 (0,y,-y) (-y,0,y) (y,-y,0) \n 2 c 3. (x,x,x) (-x,-x,-x) \n 1 b 32 (.5,.5,.5) \n 1 a 32 (0,0,0) \n ";
  }
  else {
    sg155 = "Space Group 155\n (0,0,0) (2/3,1/3,1/3) (1/3,2/3,2/3) \n 18 f 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (y,x,-z) (x-y,-y,-z) (-x,-x+y,-z) \n 9 e .2 (x,0,.5) (0,x,.5) (-x,-x,.5) \n 9 d .2 (x,0,0) (0,x,0) (-x,-x,0) \n 6 c 3. (0,0,z) (0,0,-z) \n 3 b 32 (0,0,.5) \n 3 a 32 (0,0,0) \n ";
  }
  //string sg155R = "Space Group 155R\n (0,0,0) \n 6 f 1 (x,y,z) (z,x,y) (y,z,x) (-z,-y,-x) (-y,-x,-z) (-x,-z,-y) \n 3 e .2 (.5,y,-y) (-y,.5,y) (y,-y,.5) \n 3 d .2 (0,y,-y) (-y,0,y) (y,-y,0) \n 2 c 3. (x,x,x) (-x,-x,-x) \n 1 b 32 (.5,.5,.5) \n 1 a 32 (0,0,0) \n ";

  string sg156 = "Space Group 156\n (0,0,0) \n 6 e 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-y,-x,z) (-x+y,y,z) (x,x-y,z) \n 3 d .m. (x,-x,z) (x,2*x,z) (-2*x,-x,z) \n 1 c 3m. (2/3,1/3,z) \n 1 b 3m. (1/3,2/3,z) \n 1 a 3m. (0,0,z) \n ";

  string sg157 = "Space Group 157\n (0,0,0) \n 6 d 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (y,x,z) (x-y,-y,z) (-x,-x+y,z) \n 3 c ..m (x,0,z) (0,x,z) (-x,-x,z) \n 2 b 3.. (1/3,2/3,z) (2/3,1/3,z) \n 1 a 3.m (0,0,z) \n ";

  string sg158 = "Space Group 158\n (0,0,0) \n 6 d 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-y,-x,z+.5) (-x+y,y,z+.5) (x,x-y,z+.5) \n 2 c 3.. (2/3,1/3,z) (2/3,1/3,z+.5) \n 2 b 3.. (1/3,2/3,z) (1/3,2/3,z+.5) \n 2 a 3.. (0,0,z) (0,0,z+.5) \n ";

  string sg159 = "Space Group 159\n (0,0,0) \n 6 c 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (y,x,z+.5) (x-y,-y,z+.5) (-x,-x+y,z+.5) \n 2 b 3.. (1/3,2/3,z) (2/3,1/3,z+.5) \n 2 a 3.. (0,0,z) (0,0,z+.5) \n ";

  //DX20180723 - added rhombohedral setting
  string sg160 = "";
  if(axis_cell==setting_2){
    sg160 = "Space Group 160\n (0,0,0) (2/3,1/3,1/3) (1/3,2/3,2/3) \n 18 c 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-y,-x,z) (-x+y,y,z) (x,x-y,z) \n 9 b .m (x,-x,z) (x,2*x,z) (-2*x,-x,z) \n 3 a 3m (0,0,z) \n ";
  }
  else if(axis_cell==setting_1 || axis_cell==setting_anrl){
    sg160 = "Space Group 160\n (0,0,0) \n 6 c 1 (x,y,z) (z,x,y) (y,z,x) (z,y,x) (y,x,z) (x,z,y) \n 3 b .m (x,x,z) (z,x,x) (x,z,x) \n 1 a 3m (x,x,x) \n ";
  }
  else {
    sg160 = "Space Group 160\n (0,0,0) (2/3,1/3,1/3) (1/3,2/3,2/3) \n 18 c 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-y,-x,z) (-x+y,y,z) (x,x-y,z) \n 9 b .m (x,-x,z) (x,2*x,z) (-2*x,-x,z) \n 3 a 3m (0,0,z) \n ";
  }
  //string sg160R = "Space Group 160R\n (0,0,0) \n 6 c 1 (x,y,z) (z,x,y) (y,z,x) (z,y,x) (y,x,z) (x,z,y) \n 3 b .m (x,x,z) (z,x,x) (x,z,x) \n 1 a 3m (x,x,x) \n ";

  //DX20180723 - added rhombohedral setting
  string sg161 = "";
  if(axis_cell==setting_2){
    sg161 = "Space Group 161\n (0,0,0) (2/3,1/3,1/3) (1/3,2/3,2/3) \n 18 b 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-y,-x,z+.5) (-x+y,y,z+.5) (x,x-y,z+.5) \n 6 a 3. (0,0,z) (0,0,z+.5) \n ";
  }
  else if(axis_cell==setting_1 || axis_cell==setting_anrl){
    sg161 = "Space Group 161\n (0,0,0) \n 6 b 1 (x,y,z) (z,x,y) (y,z,x) (z+.5,y+.5,x+.5) (y+.5,x+.5,z+.5) (x+.5,z+.5,y+.5) \n 2 a 3. (x,x,x) (x+.5,x+.5,x+.5) \n ";
  }
  else {
    sg161 = "Space Group 161\n (0,0,0) (2/3,1/3,1/3) (1/3,2/3,2/3) \n 18 b 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-y,-x,z+.5) (-x+y,y,z+.5) (x,x-y,z+.5) \n 6 a 3. (0,0,z) (0,0,z+.5) \n ";
  }
  //string sg161R = "Space Group 161R\n (0,0,0) \n 6 b 1 (x,y,z) (z,x,y) (y,z,x) (z+.5,y+.5,x+.5) (y+.5,x+.5,z+.5) (x+.5,z+.5,y+.5) \n 2 a 3. (x,x,x) (x+.5,x+.5,x+.5) \n ";

  string sg162 = "Space Group 162\n (0,0,0) \n 12 l 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-y,-x,-z) (-x+y,y,-z) (x,x-y,-z) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) (y,x,z) (x-y,-y,z) (-x,-x+y,z) \n 6 k ..m (x,0,z) (0,x,z) (-x,-x,z) (0,-x,-z) (-x,0,-z) (x,x,-z) \n 6 j ..2 (x,-x,.5) (x,2*x,.5) (-2*x,-x,.5) (-x,x,.5) (-x,-2*x,.5) (2*x,x,.5) \n 6 i ..2 (x,-x,0) (x,2*x,0) (-2*x,-x,0) (-x,x,0) (-x,-2*x,0) (2*x,x,0) \n 4 h 3.. (1/3,2/3,z) (1/3,2/3,-z) (2/3,1/3,-z) (2/3,1/3,z) \n 3 g ..2/m (.5,0,.5) (0,.5,.5) (.5,.5,.5) \n 3 f ..2/m (.5,0,0) (0,.5,0) (.5,.5,0) \n 2 e 3.m (0,0,z) (0,0,-z) \n 2 d 3.2 (1/3,2/3,.5) (2/3,1/3,.5) \n 2 c 3.2 (1/3,2/3,0) (2/3,1/3,0) \n 1 b -3.m (0,0,.5) \n 1 a -3.m (0,0,0) \n ";

  string sg163 = "Space Group 163\n (0,0,0) \n 12 i 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-y,-x,-z+.5) (-x+y,y,-z+.5) (x,x-y,-z+.5) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) (y,x,z+.5) (x-y,-y,z+.5) (-x,-x+y,z+.5) \n 6 h ..2 (x,-x,.25) (x,2*x,.25) (-2*x,-x,.25) (-x,x,.75) (-x,-2*x,.75) (2*x,x,.75) \n 6 g -1 (.5,0,0) (0,.5,0) (.5,.5,0) (0,.5,.5) (.5,0,.5) (.5,.5,.5) \n 4 f 3.. (1/3,2/3,z) (1/3,2/3,-z+.5) (2/3,1/3,-z) (2/3,1/3,z+.5) \n 4 e 3.. (0,0,z) (0,0,-z+.5) (0,0,-z) (0,0,z+.5) \n 2 d 3.2 (2/3,1/3,.25) (1/3,2/3,.75) \n 2 c 3.2 (1/3,2/3,.25) (2/3,1/3,.75) \n 2 b -3.. (0,0,0) (0,0,.5) \n 2 a 3.2 (0,0,.25) (0,0,.75) \n ";

  string sg164 = "Space Group 164\n (0,0,0) \n 12 j 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (y,x,-z) (x-y,-y,-z) (-x,-x+y,-z) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) (-y,-x,z) (-x+y,y,z) (x,x-y,z) \n 6 i .m. (x,-x,z) (x,2*x,z) (-2*x,-x,z) (-x,x,-z) (2*x,x,-z) (-x,-2*x,-z) \n 6 h .2. (x,0,.5) (0,x,.5) (-x,-x,.5) (-x,0,.5) (0,-x,.5) (x,x,.5) \n 6 g .2. (x,0,0) (0,x,0) (-x,-x,0) (-x,0,0) (0,-x,0) (x,x,0) \n 3 f .2/m. (.5,0,.5) (0,.5,.5) (.5,.5,.5) \n 3 e .2/m. (.5,0,0) (0,.5,0) (.5,.5,0) \n 2 d 3m. (1/3,2/3,z) (2/3,1/3,-z) \n 2 c 3m. (0,0,z) (0,0,-z) \n 1 b -3m. (0,0,.5) \n 1 a -3m. (0,0,0) \n "; //DX20201103 - need multiplication (*) symbol

  string sg165 = "Space Group 165\n (0,0,0) \n 12 g 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (y,x,-z+.5) (x-y,-y,-z+.5) (-x,-x+y,-z+.5) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) (-y,-x,z+.5) (-x+y,y,z+.5) (x,x-y,z+.5) \n 6 f .2. (x,0,.25) (0,x,.25) (-x,-x,.25) (-x,0,.75) (0,-x,.75) (x,x,.75) \n 6 e -1 (.5,0,0) (0,.5,0) (.5,.5,0) (0,.5,.5) (.5,0,.5) (.5,.5,.5) \n 4 d 3.. (1/3,2/3,z) (2/3,1/3,-z+.5) (2/3,1/3,-z) (1/3,2/3,z+.5) \n 4 c 3.. (0,0,z) (0,0,-z+.5) (0,0,-z) (0,0,z+.5) \n 2 b -3.. (0,0,0) (0,0,.5) \n 2 a 32. (0,0,.25) (0,0,.75) \n ";

  //DX20180723 - added rhombohedral setting
  string sg166 = "";
  if(axis_cell==setting_2){
    sg166 = "Space Group 166\n (0,0,0) (2/3,1/3,1/3) (1/3,2/3,2/3) \n 36 i 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (y,x,-z) (x-y,-y,-z) (-x,-x+y,-z) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) (-y,-x,z) (-x+y,y,z) (x,x-y,z) \n 18 h .m (x,-x,z) (x,2*x,z) (-2*x,-x,z) (-x,x,-z) (2*x,x,-z) (-x,-2*x,-z) \n 18 g .2 (x,0,.5) (0,x,.5) (-x,-x,.5) (-x,0,.5) (0,-x,.5) (x,x,.5) \n 18 f .2 (x,0,0) (0,x,0) (-x,-x,0) (-x,0,0) (0,-x,0) (x,x,0) \n 9 e .2/m (.5,0,0) (0,.5,0) (.5,.5,0) \n 9 d .2/m (.5,0,.5) (0,.5,.5) (.5,.5,.5) \n 6 c 3m (0,0,z) (0,0,-z) \n 3 b -3m (0,0,.5) \n 3 a -3m (0,0,0) \n ";
  }
  else if(axis_cell==setting_1 || axis_cell==setting_anrl){
    sg166 = "Space Group 166\n (0,0,0) \n 12 i 1 (x,y,z) (z,x,y) (y,z,x) (-z,-y,-x) (-y,-x,-z) (-x,-z,-y) (-x,-y,-z) (-z,-x,-y) (-y,-z,-x) (z,y,x) (y,x,z) (x,z,y) \n 6 h .m (x,x,z) (z,x,x) (x,z,x) (-z,-x,-x) (-x,-x,-z) (-x,-z,-x) \n 6 g .2 (x,-x,.5) (.5,x,-x) (-x,.5,x) (-x,x,.5) (.5,-x,x) (x,.5,-x) \n 6 f .2 (x,-x,0) (0,x,-x) (-x,0,x) (-x,x,0) (0,-x,x) (x,0,-x) \n 3 e .2/m (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 3 d .2/m (.5,0,0) (0,.5,0) (0,0,.5) \n 2 c 3m (x,x,x) (-x,-x,-x) \n 1 b -3m (.5,.5,.5) \n 1 a -3m (0,0,0) \n ";
  }
  else {
    sg166 = "Space Group 166\n (0,0,0) (2/3,1/3,1/3) (1/3,2/3,2/3) \n 36 i 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (y,x,-z) (x-y,-y,-z) (-x,-x+y,-z) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) (-y,-x,z) (-x+y,y,z) (x,x-y,z) \n 18 h .m (x,-x,z) (x,2*x,z) (-2*x,-x,z) (-x,x,-z) (2*x,x,-z) (-x,-2*x,-z) \n 18 g .2 (x,0,.5) (0,x,.5) (-x,-x,.5) (-x,0,.5) (0,-x,.5) (x,x,.5) \n 18 f .2 (x,0,0) (0,x,0) (-x,-x,0) (-x,0,0) (0,-x,0) (x,x,0) \n 9 e .2/m (.5,0,0) (0,.5,0) (.5,.5,0) \n 9 d .2/m (.5,0,.5) (0,.5,.5) (.5,.5,.5) \n 6 c 3m (0,0,z) (0,0,-z) \n 3 b -3m (0,0,.5) \n 3 a -3m (0,0,0) \n ";
  }
  //string sg166R = "Space Group 166R\n (0,0,0) \n 12 i 1 (x,y,z) (z,x,y) (y,z,x) (-z,-y,-x) (-y,-x,-z) (-x,-z,-y) (-x,-y,-z) (-z,-x,-y) (-y,-z,-x) (z,y,x) (y,x,z) (x,z,y) \n 6 h .m (x,x,z) (z,x,x) (x,z,x) (-z,-x,-x) (-x,-x,-z) (-x,-z,-x) \n 6 g .2 (x,-x,.5) (.5,x,-x) (-x,.5,x) (-x,x,.5) (.5,-x,x) (x,.5,-x) \n 6 f .2 (x,-x,0) (0,x,-x) (-x,0,x) (-x,x,0) (0,-x,x) (x,0,-x) \n 3 e .2/m (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 3 d .2/m (.5,0,0) (0,.5,0) (0,0,.5) \n 2 c 3m (x,x,x) (-x,-x,-x) \n 1 b -3m (.5,.5,.5) \n 1 a -3m (0,0,0) \n ";

  //DX20180723 - added rhombohedral setting
  string sg167 = "";
  if(axis_cell==setting_2){
    sg167 = "Space Group 167\n (0,0,0) (2/3,1/3,1/3) (1/3,2/3,2/3) \n 36 f 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (y,x,-z+.5) (x-y,-y,-z+.5) (-x,-x+y,-z+.5) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) (-y,-x,z+.5) (-x+y,y,z+.5) (x,x-y,z+.5) \n 18 e .2 (x,0,.25) (0,x,.25) (-x,-x,.25) (-x,0,.75) (0,-x,.75) (x,x,.75) \n 18 d -1 (.5,0,0) (0,.5,0) (.5,.5,0) (0,.5,.5) (.5,0,.5) (.5,.5,.5) \n 12 c 3. (0,0,z) (0,0,-z+.5) (0,0,-z) (0,0,z+.5) \n 6 b -3. (0,0,0) (0,0,.5) \n 6 a 32 (0,0,.25) (0,0,.75) \n ";
  }
  else if(axis_cell==setting_1 || axis_cell==setting_anrl){
    sg167 = "Space Group 167\n (0,0,0) \n 12 f 1 (x,y,z) (z,x,y) (y,z,x) (-z+.5,-y+.5,-x+.5) (-y+.5,-x+.5,-z+.5) (-x+.5,-z+.5,-y+.5) (-x,-y,-z) (-z,-x,-y) (-y,-z,-x) (z+.5,y+.5,x+.5) (y+.5,x+.5,z+.5) (x+.5,z+.5,y+.5) \n 6 e .2 (x,-x+.5,.25) (.25,x,-x+.5) (-x+.5,.25,x) (-x,x+.5,.75) (.75,-x,x+.5) (x+.5,.75,-x) \n 6 d -1 (.5,0,0) (0,.5,0) (0,0,.5) (.5,.5,0) (.5,0,.5) (0,.5,.5) \n 4 c 3. (x,x,x) (-x+.5,-x+.5,-x+.5) (-x,-x,-x) (x+.5,x+.5,x+.5) \n 2 b -3. (0,0,0) (.5,.5,.5) \n 2 a 32 (.25,.25,.25) (.75,.75,.75) \n ";
  }
  else {
    sg167 = "Space Group 167\n (0,0,0) (2/3,1/3,1/3) (1/3,2/3,2/3) \n 36 f 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (y,x,-z+.5) (x-y,-y,-z+.5) (-x,-x+y,-z+.5) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) (-y,-x,z+.5) (-x+y,y,z+.5) (x,x-y,z+.5) \n 18 e .2 (x,0,.25) (0,x,.25) (-x,-x,.25) (-x,0,.75) (0,-x,.75) (x,x,.75) \n 18 d -1 (.5,0,0) (0,.5,0) (.5,.5,0) (0,.5,.5) (.5,0,.5) (.5,.5,.5) \n 12 c 3. (0,0,z) (0,0,-z+.5) (0,0,-z) (0,0,z+.5) \n 6 b -3. (0,0,0) (0,0,.5) \n 6 a 32 (0,0,.25) (0,0,.75) \n ";
  }
  //string sg167R = "Space Group 167R\n (0,0,0) \n 12 f 1 (x,y,z) (z,x,y) (y,z,x) (-z+.5,-y+.5,-x+.5) (-y+.5,-x+.5,-z+.5) (-x+.5,-z+.5,-y+.5) (-x,-y,-z) (-z,-x,-y) (-y,-z,-x) (z+.5,y+.5,x+.5) (y+.5,x+.5,z+.5) (x+.5,z+.5,y+.5) \n 6 e .2 (x,-x+.5,.25) (.25,x,-x+.5) (-x+.5,.25,x) (-x,x+.5,.75) (.75,-x,x+.5) (x+.5,.75,-x) \n 6 d -1 (.5,0,0) (0,.5,0) (0,0,.5) (.5,.5,0) (.5,0,.5) (0,.5,.5) \n 4 c 3. (x,x,x) (-x+.5,-x+.5,-x+.5) (-x,-x,-x) (x+.5,x+.5,x+.5) \n 2 b -3. (0,0,0) (.5,.5,.5) \n 2 a 32 (.25,.25,.25) (.75,.75,.75) \n ";

  string sg168 = "Space Group 168\n (0,0,0) \n 6 d 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,z) (y,-x+y,z) (x-y,x,z) \n 3 c 2.. (.5,0,z) (0,.5,z) (.5,.5,z) \n 2 b 3.. (1/3,2/3,z) (2/3,1/3,z) \n 1 a 6.. (0,0,z) \n ";

  string sg169 = "Space Group 169\n (0,0,0) \n 6 a 1 (x,y,z) (-y,x-y,z+1/3) (-x+y,-x,z+2/3) (-x,-y,z+.5) (y,-x+y,z+5/6) (x-y,x,z+1/6) \n ";

  string sg170 = "Space Group 170\n (0,0,0) \n 6 a 1 (x,y,z) (-y,x-y,z+2/3) (-x+y,-x,z+1/3) (-x,-y,z+.5) (y,-x+y,z+1/6) (x-y,x,z+5/6) \n ";

  string sg171 = "Space Group 171\n (0,0,0) \n 6 c 1 (x,y,z) (-y,x-y,z+2/3) (-x+y,-x,z+1/3) (-x,-y,z) (y,-x+y,z+2/3) (x-y,x,z+1/3) \n 3 b 2.. (.5,.5,z) (.5,0,z+2/3) (0,.5,z+1/3) \n 3 a 2.. (0,0,z) (0,0,z+2/3) (0,0,z+1/3) \n ";

  string sg172 = "Space Group 172\n (0,0,0) \n 6 c 1 (x,y,z) (-y,x-y,z+1/3) (-x+y,-x,z+2/3) (-x,-y,z) (y,-x+y,z+1/3) (x-y,x,z+2/3) \n 3 b 2.. (.5,.5,z) (.5,0,z+1/3) (0,.5,z+2/3) \n 3 a 2.. (0,0,z) (0,0,z+1/3) (0,0,z+2/3) \n ";

  string sg173 = "Space Group 173\n (0,0,0) \n 6 c 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,z+.5) (y,-x+y,z+.5) (x-y,x,z+.5) \n 2 b 3.. (1/3,2/3,z) (2/3,1/3,z+.5) \n 2 a 3.. (0,0,z) (0,0,z+.5) \n ";

  string sg174 = "Space Group 174\n (0,0,0) \n 6 l 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (x,y,-z) (-y,x-y,-z) (-x+y,-x,-z) \n 3 k m.. (x,y,.5) (-y,x-y,.5) (-x+y,-x,.5) \n 3 j m.. (x,y,0) (-y,x-y,0) (-x+y,-x,0) \n 2 i 3.. (2/3,1/3,z) (2/3,1/3,-z) \n 2 h 3.. (1/3,2/3,z) (1/3,2/3,-z) \n 2 g 3.. (0,0,z) (0,0,-z) \n 1 f -6.. (2/3,1/3,.5) \n 1 e -6.. (2/3,1/3,0) \n 1 d -6.. (1/3,2/3,.5) \n 1 c -6.. (1/3,2/3,0) \n 1 b -6.. (0,0,.5) \n 1 a -6.. (0,0,0) \n ";

  string sg175 = "Space Group 175\n (0,0,0) \n 12 l 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,z) (y,-x+y,z) (x-y,x,z) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) (x,y,-z) (-y,x-y,-z) (-x+y,-x,-z) \n 6 k m.. (x,y,.5) (-y,x-y,.5) (-x+y,-x,.5) (-x,-y,.5) (y,-x+y,.5) (x-y,x,.5) \n 6 j m.. (x,y,0) (-y,x-y,0) (-x+y,-x,0) (-x,-y,0) (y,-x+y,0) (x-y,x,0) \n 6 i 2.. (.5,0,z) (0,.5,z) (.5,.5,z) (.5,0,-z) (0,.5,-z) (.5,.5,-z) \n 4 h 3.. (1/3,2/3,z) (2/3,1/3,z) (2/3,1/3,-z) (1/3,2/3,-z) \n 3 g 2/m.. (.5,0,.5) (0,.5,.5) (.5,.5,.5) \n 3 f 2/m.. (.5,0,0) (0,.5,0) (.5,.5,0) \n 2 e 6.. (0,0,z) (0,0,-z) \n 2 d -6.. (1/3,2/3,.5) (2/3,1/3,.5) \n 2 c -6.. (1/3,2/3,0) (2/3,1/3,0) \n 1 b 6/m.. (0,0,.5) \n 1 a 6/m.. (0,0,0) \n ";

  string sg176 = "Space Group 176\n (0,0,0) \n 12 i 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,z+.5) (y,-x+y,z+.5) (x-y,x,z+.5) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) (x,y,-z+.5) (-y,x-y,-z+.5) (-x+y,-x,-z+.5) \n 6 h m.. (x,y,.25) (-y,x-y,.25) (-x+y,-x,.25) (-x,-y,.75) (y,-x+y,.75) (x-y,x,.75) \n 6 h m.. (x,y,.25) (-y,x-y,.25) (-x+y,-x,.25) (-x,-y,.75) (y,-x+y,.75) (x-y,x,.75) \n 6 g -1 (.5,0,0) (0,.5,0) (.5,.5,0) (.5,0,.5) (0,.5,.5) (.5,.5,.5) \n 4 f 3.. (1/3,2/3,z) (2/3,1/3,z+.5) (2/3,1/3,-z) (1/3,2/3,-z+.5) \n 4 e 3.. (0,0,z) (0,0,z+.5) (0,0,-z) (0,0,-z+.5) \n 2 d -6.. (2/3,1/3,.25) (1/3,2/3,.75) \n 2 c -6.. (1/3,2/3,.25) (2/3,1/3,.75) \n 2 b -3.. (0,0,0) (0,0,.5) \n 2 a -6.. (0,0,.25) (0,0,.75) \n ";

  string sg177 = "Space Group 177 \n (0,0,0) \n 12 n 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,z) (y,-x+y,z) (x-y,x,z) (y,x,-z) (x-y,-y,-z) (-x,-x+y,-z) (-y,-x,-z) (-x+y,y,-z) (x,x-y,-z) \n 6 m ..2 (x,-x,.5) (x,2*x,.5) (-2*x,-x,.5) (-x,x,.5) (-x,-2*x,.5) (2*x,x,.5) \n 6 l ..2 (x,-x,0) (x,2*x,0) (-2*x,-x,0) (-x,x,0) (-x,-2*x,0) (2*x,x,0) \n 6 k .2. (x,0,.5) (0,x,.5) (-x,-x,.5) (-x,0,.5) (0,-x,.5) (x,x,.5) \n 6 j .2. (x,0,0) (0,x,0) (-x,-x,0) (-x,0,0) (0,-x,0) (x,x,0) \n 6 i 2.. (.5,0,z) (0,.5,z) (.5,.5,z) (0,.5,-z) (.5,0,-z) (.5,.5,-z) \n 4 h 3.. (1/3,2/3,z) (2/3,1/3,z) (2/3,1/3,-z) (1/3,2/3,-z) \n 3 g 222 (.5,0,.5) (0,.5,.5) (.5,.5,.5) \n 3 f 222 (.5,0,0) (0,.5,0) (.5,.5,0) \n 2 e 6.. (0,0,z) (0,0,-z) \n 2 d 3.2 (1/3,2/3,.5) (2/3,1/3,.5) \n 2 c 3.2 (1/3,2/3,0) (2/3,1/3,0) \n 1 b 622 (0,0,.5) \n 1 a 622 (0,0,0) \n ";

  string sg178 = "Space Group 178 \n (0,0,0) \n 12 c 1 (x,y,z) (-y,x-y,z+1/3) (-x+y,-x,z+2/3) (-x,-y,z+.5) (y,-x+y,z+5/6) (x-y,x,z+1/6) (y,x,-z+1/3) (x-y,-y,-z) (-x,-x+y,-z+2/3) (-y,-x,-z+5/6) (-x+y,y,-z+.5) (x,x-y,-z+1/6) \n 6 b ..2 (x,2*x,.25) (-2*x,-x,7/12) (x,-x,11/12) (-x,-2*x,.75) (2*x,x,1/12) (-x,x,5/12) \n 6 a .2. (x,0,0) (0,x,1/3) (-x,-x,2/3) (-x,0,.5) (0,-x,5/6) (x,x,1/6) \n ";

  string sg179 = "Space Group 179 \n (0,0,0) \n 12 c 1 (x,y,z) (-y,x-y,z+2/3) (-x+y,-x,z+1/3) (-x,-y,z+.5) (y,-x+y,z+1/6) (x-y,x,z+5/6) (y,x,-z+2/3) (x-y,-y,-z) (-x,-x+y,-z+1/3) (-y,-x,-z+1/6) (-x+y,y,-z+.5) (x,x-y,-z+5/6) \n 6 b ..2 (x,2*x,.75) (-2*x,-x,5/12) (x,-x,1/12) (-x,-2*x,.25) (2*x,x,11/12) (-x,x,7/12) \n 6 a .2. (x,0,0) (0,x,2/3) (-x,-x,1/3) (-x,0,.5) (0,-x,1/6) (x,x,5/6) \n ";

  string sg180 = "Space Group 180 \n (0,0,0) \n 12 k 1 (x,y,z) (-y,x-y,z+2/3) (-x+y,-x,z+1/3) (-x,-y,z) (y,-x+y,z+2/3) (x-y,x,z+1/3) (y,x,-z+2/3) (x-y,-y,-z) (-x,-x+y,-z+1/3) (-y,-x,-z+2/3) (-x+y,y,-z) (x,x-y,-z+1/3) \n 6 j ..2 (x,2*x,.5) (-2*x,-x,1/6) (x,-x,5/6) (-x,-2*x,.5) (2*x,x,1/6) (-x,x,5/6) \n 6 i ..2 (x,2*x,0) (-2*x,-x,2/3) (x,-x,1/3) (-x,-2*x,0) (2*x,x,2/3) (-x,x,1/3) \n 6 h .2. (x,0,.5) (0,x,1/6) (-x,-x,5/6) (-x,0,.5) (0,-x,1/6) (x,x,5/6) \n 6 g .2. (x,0,0) (0,x,2/3) (-x,-x,1/3) (-x,0,0) (0,-x,2/3) (x,x,1/3) \n 6 f 2.. (.5,0,z) (0,.5,z+2/3) (.5,.5,z+1/3) (0,.5,-z+2/3) (.5,0,-z) (.5,.5,-z+1/3) \n 6 e 2.. (0,0,z) (0,0,z+2/3) (0,0,z+1/3) (0,0,-z+2/3) (0,0,-z) (0,0,-z+1/3) \n 3 d 222 (.5,0,.5) (0,.5,1/6) (.5,.5,5/6) \n 3 c 222 (.5,0,0) (0,.5,2/3) (.5,.5,1/3) \n 3 b 222 (0,0,.5) (0,0,1/6) (0,0,5/6) \n 3 a 222 (0,0,0) (0,0,2/3) (0,0,1/3) \n ";

  string sg181 = "Space Group 181 \n (0,0,0) \n 12 k 1 (x,y,z) (-y,x-y,z+1/3) (-x+y,-x,z+2/3) (-x,-y,z) (y,-x+y,z+1/3) (x-y,x,z+2/3) (y,x,-z+1/3) (x-y,-y,-z) (-x,-x+y,-z+2/3) (-y,-x,-z+1/3) (-x+y,y,-z) (x,x-y,-z+2/3) \n 6 j ..2 (x,2*x,.5) (-2*x,-x,5/6) (x,-x,1/6) (-x,-2*x,.5) (2*x,x,5/6) (-x,x,1/6) \n 6 i ..2 (x,2*x,0) (-2*x,-x,1/3) (x,-x,2/3) (-x,-2*x,0) (2*x,x,1/3) (-x,x,2/3) \n 6 h .2. (x,0,.5) (0,x,5/6) (-x,-x,1/6) (-x,0,.5) (0,-x,5/6) (x,x,1/6) \n 6 g .2. (x,0,0) (0,x,1/3) (-x,-x,2/3) (-x,0,0) (0,-x,1/3) (x,x,2/3) \n 6 f 2.. (.5,0,z) (0,.5,z+1/3) (.5,.5,z+2/3) (0,.5,-z+1/3) (.5,0,-z) (.5,.5,-z+2/3) \n 6 e 2.. (0,0,z) (0,0,z+1/3) (0,0,z+2/3) (0,0,-z+1/3) (0,0,-z) (0,0,-z+2/3) \n 3 d 222 (.5,0,.5) (0,.5,5/6) (.5,.5,1/6) \n 3 c 222 (.5,0,0) (0,.5,1/3) (.5,.5,2/3) \n 3 b 222 (0,0,.5) (0,0,5/6) (0,0,1/6) \n 3 a 222 (0,0,0) (0,0,1/3) (0,0,2/3) \n ";

  string sg182 = "Space Group 182 \n (0,0,0) \n 12 i 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,z+.5) (y,-x+y,z+.5) (x-y,x,z+.5) (y,x,-z) (x-y,-y,-z) (-x,-x+y,-z) (-y,-x,-z+.5) (-x+y,y,-z+.5) (x,x-y,-z+.5) \n 6 h ..2 (x,2*x,.25) (-2*x,-x,.25) (x,-x,.25) (-x,-2*x,.75) (2*x,x,.75) (-x,x,.75) \n 6 g .2. (x,0,0) (0,x,0) (-x,-x,0) (-x,0,.5) (0,-x,.5) (x,x,.5) \n 4 f 3.. (1/3,2/3,z) (2/3,1/3,z+.5) (2/3,1/3,-z) (1/3,2/3,-z+.5) \n 4 e 3.. (0,0,z) (0,0,z+.5) (0,0,-z) (0,0,-z+.5) \n 2 d 3.2 (1/3,2/3,.75) (2/3,1/3,.25) \n 2 c 3.2 (1/3,2/3,.25) (2/3,1/3,.75) \n 2 b 3.2 (0,0,.25) (0,0,.75) \n 2 a 32. (0,0,0) (0,0,.5) \n ";

  string sg183 = "Space Group 183 \n (0,0,0) \n 12 f 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,z) (y,-x+y,z) (x-y,x,z) (-y,-x,z) (-x+y,y,z) (x,x-y,z) (y,x,z) (x-y,-y,z) (-x,-x+y,z) \n 6 e .m. (x,-x,z) (x,2*x,z) (-2*x,-x,z) (-x,x,z) (-x,-2*x,z) (2*x,x,z) \n 6 d ..m (x,0,z) (0,x,z) (-x,-x,z) (-x,0,z) (0,-x,z) (x,x,z) \n 3 c 2mm (.5,0,z) (0,.5,z) (.5,.5,z) \n 2 b 3m. (1/3,2/3,z) (2/3,1/3,z) \n 1 a 6mm (0,0,z) \n ";

  string sg184 = "Space Group 184 \n (0,0,0) \n 12 d 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,z) (y,-x+y,z) (x-y,x,z) (-y,-x,z+.5) (-x+y,y,z+.5) (x,x-y,z+.5) (y,x,z+.5) (x-y,-y,z+.5) (-x,-x+y,z+.5) \n 6 c 2.. (.5,0,z) (0,.5,z) (.5,.5,z) (0,.5,z+.5) (.5,0,z+.5) (.5,.5,z+.5) \n 4 b 3.. (1/3,2/3,z) (2/3,1/3,z) (1/3,2/3,z+.5) (2/3,1/3,z+.5) \n 2 a 6.. (0,0,z) (0,0,z+.5) \n ";

  string sg185 = "Space Group 185 \n (0,0,0) \n 12 d 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,z+.5) (y,-x+y,z+.5) (x-y,x,z+.5) (-y,-x,z+.5) (-x+y,y,z+.5) (x,x-y,z+.5) (y,x,z) (x-y,-y,z) (-x,-x+y,z) \n 6 c ..m (x,0,z) (0,x,z) (-x,-x,z) (-x,0,z+.5) (0,-x,z+.5) (x,x,z+.5) \n 4 b 3.. (1/3,2/3,z) (2/3,1/3,z+.5) (1/3,2/3,z+.5) (2/3,1/3,z) \n 2 a 3.m (0,0,z) (0,0,z+.5) \n ";

  string sg186 = "Space Group 186 \n (0,0,0) \n 12 d 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,z+.5) (y,-x+y,z+.5) (x-y,x,z+.5) (-y,-x,z) (-x+y,y,z) (x,x-y,z) (y,x,z+.5) (x-y,-y,z+.5) (-x,-x+y,z+.5) \n 6 c .m. (x,-x,z) (x,2*x,z) (-2*x,-x,z) (-x,x,z+.5) (-x,-2*x,z+.5) (2*x,x,z+.5) \n 2 b 3m. (1/3,2/3,z) (2/3,1/3,z+.5) \n 2 a 3m. (0,0,z) (0,0,z+.5) \n ";

  string sg187 = "Space Group 187 \n (0,0,0) \n 12 o 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (x,y,-z) (-y,x-y,-z) (-x+y,-x,-z) (-y,-x,z) (-x+y,y,z) (x,x-y,z) (-y,-x,-z) (-x+y,y,-z) (x,x-y,-z) \n 6 n .m. (x,-x,z) (x,2*x,z) (-2*x,-x,z) (x,-x,-z) (x,2*x,-z) (-2*x,-x,-z) \n 6 m m.. (x,y,.5) (-y,x-y,.5) (-x+y,-x,.5) (-y,-x,.5) (-x+y,y,.5) (x,x-y,.5) \n 6 l m.. (x,y,0) (-y,x-y,0) (-x+y,-x,0) (-y,-x,0) (-x+y,y,0) (x,x-y,0) \n 3 k mm2 (x,-x,.5) (x,2*x,.5) (-2*x,-x,.5) \n 3 j mm2 (x,-x,0) (x,2*x,0) (-2*x,-x,0) \n 2 i 3m. (2/3,1/3,z) (2/3,1/3,-z) \n 2 h 3m. (1/3,2/3,z) (1/3,2/3,-z) \n 2 g 3m. (0,0,z) (0,0,-z) \n 1 f -6m2 (2/3,1/3,.5) \n 1 e -6m2 (2/3,1/3,0) \n 1 d -6m2 (1/3,2/3,.5) \n 1 c -6m2 (1/3,2/3,0) \n 1 b -6m2 (0,0,.5) \n 1 a -6m2 (0,0,0) \n ";

  string sg188 = "Space Group 188 \n (0,0,0) \n 12 l 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (x,y,-z+.5) (-y,x-y,-z+.5) (-x+y,-x,-z+.5) (-y,-x,z+.5) (-x+y,y,z+.5) (x,x-y,z+.5) (-y,-x,-z) (-x+y,y,-z) (x,x-y,-z) \n 6 k m.. (x,y,.25) (-y,x-y,.25) (-x+y,-x,.25) (-y,-x,.75) (-x+y,y,.75) (x,x-y,.75) \n 6 j ..2 (x,-x,0) (x,2*x,0) (-2*x,-x,0) (x,-x,.5) (x,2*x,.5) (-2*x,-x,.5) \n 4 i 3.. (2/3,1/3,z) (2/3,1/3,-z+.5) (2/3,1/3,z+.5) (2/3,1/3,-z) \n 4 h 3.. (1/3,2/3,z) (1/3,2/3,-z+.5) (1/3,2/3,z+.5) (1/3,2/3,-z) \n 4 g 3.. (0,0,z) (0,0,-z+.5) (0,0,z+.5) (0,0,-z) \n 2 f -6.. (2/3,1/3,.25) (2/3,1/3,.75) \n 2 e 3.2 (2/3,1/3,0) (2/3,1/3,.5) \n 2 d -6.. (1/3,2/3,.25) (1/3,2/3,.75) \n 2 c 3.2 (1/3,2/3,0) (1/3,2/3,.5) \n 2 b -6.. (0,0,.25) (0,0,.75) \n 2 a 3.2 (0,0,0) (0,0,.5) \n ";

  string sg189 = "Space Group 189 \n (0,0,0) \n 12 l 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (x,y,-z) (-y,x-y,-z) (-x+y,-x,-z) (y,x,-z) (x-y,-y,-z) (-x,-x+y,-z) (y,x,z) (x-y,-y,z) (-x,-x+y,z) \n 6 k m.. (x,y,.5) (-y,x-y,.5) (-x+y,-x,.5) (y,x,.5) (x-y,-y,.5) (-x,-x+y,.5) \n 6 j m.. (x,y,0) (-y,x-y,0) (-x+y,-x,0) (y,x,0) (x-y,-y,0) (-x,-x+y,0) \n 6 i ..m (x,0,z) (0,x,z) (-x,-x,z) (x,0,-z) (0,x,-z) (-x,-x,-z) \n 4 h 3.. (1/3,2/3,z) (1/3,2/3,-z) (2/3,1/3,-z) (2/3,1/3,z) \n 3 g m2m (x,0,.5) (0,x,.5) (-x,-x,.5) \n 3 f m2m (x,0,0) (0,x,0) (-x,-x,0) \n 2 e 3.m (0,0,z) (0,0,-z) \n 2 d -6.. (1/3,2/3,.5) (2/3,1/3,.5) \n 2 c -6.. (1/3,2/3,0) (2/3,1/3,0) \n 1 b -62m (0,0,.5) \n 1 a -62m (0,0,0) \n ";

  string sg190 = "Space Group 190 \n (0,0,0) \n 12 i 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (x,y,-z+.5) (-y,x-y,-z+.5) (-x+y,-x,-z+.5) (y,x,-z) (x-y,-y,-z) (-x,-x+y,-z) (y,x,z+.5) (x-y,-y,z+.5) (-x,-x+y,z+.5) \n 6 h m.. (x,y,.25) (-y,x-y,.25) (-x+y,-x,.25) (y,x,.75) (x-y,-y,.75) (-x,-x+y,.75) \n 6 g .2. (x,0,0) (0,x,0) (-x,-x,0) (x,0,.5) (0,x,.5) (-x,-x,.5) \n 4 f 3.. (1/3,2/3,z) (1/3,2/3,-z+.5) (2/3,1/3,-z) (2/3,1/3,z+.5) \n 4 e 3.. (0,0,z) (0,0,-z+.5) (0,0,-z) (0,0,z+.5) \n 2 d -6.. (2/3,1/3,.25) (1/3,2/3,.75) \n 2 c -6.. (1/3,2/3,.25) (2/3,1/3,.75) \n 2 b -6.. (0,0,.25) (0,0,.75) \n 2 a 32. (0,0,0) (0,0,.5) \n ";

  string sg191 = "Space Group 191 \n (0,0,0) \n 24 r 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,z) (y,-x+y,z) (x-y,x,z) (y,x,-z) (x-y,-y,-z) (-x,-x+y,-z) (-y,-x,-z) (-x+y,y,-z) (x,x-y,-z) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) (x,y,-z) (-y,x-y,-z) (-x+y,-x,-z) (-y,-x,z) (-x+y,y,z) (x,x-y,z) (y,x,z) (x-y,-y,z) (-x,-x+y,z) \n 12 q m.. (x,y,.5) (-y,x-y,.5) (-x+y,-x,.5) (-x,-y,.5) (y,-x+y,.5) (x-y,x,.5) (y,x,.5) (x-y,-y,.5) (-x,-x+y,.5) (-y,-x,.5) (-x+y,y,.5) (x,x-y,.5) \n 12 p m.. (x,y,0) (-y,x-y,0) (-x+y,-x,0) (-x,-y,0) (y,-x+y,0) (x-y,x,0) (y,x,0) (x-y,-y,0) (-x,-x+y,0) (-y,-x,0) (-x+y,y,0) (x,x-y,0) \n 12 o .m. (x,2*x,z) (-2*x,-x,z) (x,-x,z) (-x,-2*x,z) (2*x,x,z) (-x,x,z) (2*x,x,-z) (-x,-2*x,-z) (-x,x,-z) (-2*x,-x,-z) (x,2*x,-z) (x,-x,-z) \n 12 n ..m (x,0,z) (0,x,z) (-x,-x,z) (-x,0,z) (0,-x,z) (x,x,z) (0,x,-z) (x,0,-z) (-x,-x,-z) (0,-x,-z) (-x,0,-z) (x,x,-z) \n 6 m mm2 (x,2*x,.5) (-2*x,-x,.5) (x,-x,.5) (-x,-2*x,.5) (2*x,x,.5) (-x,x,.5) \n 6 l mm2 (x,2*x,0) (-2*x,-x,0) (x,-x,0) (-x,-2*x,0) (2*x,x,0) (-x,x,0) \n 6 k m2m (x,0,.5) (0,x,.5) (-x,-x,.5) (-x,0,.5) (0,-x,.5) (x,x,.5) \n 6 j m2m (x,0,0) (0,x,0) (-x,-x,0) (-x,0,0) (0,-x,0) (x,x,0) \n 6 i 2mm (.5,0,z) (0,.5,z) (.5,.5,z) (0,.5,-z) (.5,0,-z) (.5,.5,-z) \n 4 h 3m. (1/3,2/3,z) (2/3,1/3,z) (2/3,1/3,-z) (1/3,2/3,-z) \n 3 g mmm (.5,0,.5) (0,.5,.5) (.5,.5,.5) \n 3 f mmm (.5,0,0) (0,.5,0) (.5,.5,0) \n 2 e 6mm (0,0,z) (0,0,-z) \n 2 d -6m2 (1/3,2/3,.5) (2/3,1/3,.5) \n 2 c -6m2 (1/3,2/3,0) (2/3,1/3,0) \n 1 b 6/mmm (0,0,.5) \n 1 a 6/mmm (0,0,0) \n "; //DX20201103 - added multiplication (*)

  string sg192 = "Space Group 192 \n (0,0,0) \n 24 m 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,z) (y,-x+y,z) (x-y,x,z) (y,x,-z+.5) (x-y,-y,-z+.5) (-x,-x+y,-z+.5) (-y,-x,-z+.5) (-x+y,y,-z+.5) (x,x-y,-z+.5) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) (x,y,-z) (-y,x-y,-z) (-x+y,-x,-z) (-y,-x,z+.5) (-x+y,y,z+.5) (x,x-y,z+.5) (y,x,z+.5) (x-y,-y,z+.5) (-x,-x+y,z+.5) \n 12 l m.. (x,y,0) (-y,x-y,0) (-x+y,-x,0) (-x,-y,0) (y,-x+y,0) (x-y,x,0) (y,x,.5) (x-y,-y,.5) (-x,-x+y,.5) (-y,-x,.5) (-x+y,y,.5) (x,x-y,.5) \n 12 k ..2 (x,2*x,.25) (-2*x,-x,.25) (x,-x,.25) (-x,-2*x,.25) (2*x,x,.25) (-x,x,.25) (-x,-2*x,.75) (2*x,x,.75) (-x,x,.75) (x,2*x,.75) (-2*x,-x,.75) (x,-x,.75) \n 12 j .2. (x,0,.25) (0,x,.25) (-x,-x,.25) (-x,0,.25) (0,-x,.25) (x,x,.25) (-x,0,.75) (0,-x,.75) (x,x,.75) (x,0,.75) (0,x,.75) (-x,-x,.75) \n 12 i 2.. (.5,0,z) (0,.5,z) (.5,.5,z) (0,.5,-z+.5) (.5,0,-z+.5) (.5,.5,-z+.5) (.5,0,-z) (0,.5,-z) (.5,.5,-z) (0,.5,z+.5) (.5,0,z+.5) (.5,.5,z+.5) \n 8 h 3.. (1/3,2/3,z) (2/3,1/3,z) (2/3,1/3,-z+.5) (1/3,2/3,-z+.5) (2/3,1/3,-z) (1/3,2/3,-z) (1/3,2/3,z+.5) (2/3,1/3,z+.5) \n 6 g 2/m.. (.5,0,0) (0,.5,0) (.5,.5,0) (0,.5,.5) (.5,0,.5) (.5,.5,.5) \n 6 f 222 (.5,0,.25) (0,.5,.25) (.5,.5,.25) (.5,0,.75) (0,.5,.75) (.5,.5,.75) \n 4 e 6.. (0,0,z) (0,0,-z+.5) (0,0,-z) (0,0,z+.5) \n 4 d -6.. (1/3,2/3,0) (2/3,1/3,0) (2/3,1/3,.5) (1/3,2/3,.5) \n 4 c 3.2 (1/3,2/3,.25) (2/3,1/3,.25) (2/3,1/3,.75) (1/3,2/3,.75) \n 2 b 6/m.. (0,0,0) (0,0,.5) \n 2 a 622 (0,0,.25) (0,0,.75) \n ";

  string sg193 = "Space Group 193 \n (0,0,0) \n 24 l 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,z+.5) (y,-x+y,z+.5) (x-y,x,z+.5) (y,x,-z+.5) (x-y,-y,-z+.5) (-x,-x+y,-z+.5) (-y,-x,-z) (-x+y,y,-z) (x,x-y,-z) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) (x,y,-z+.5) (-y,x-y,-z+.5) (-x+y,-x,-z+.5) (-y,-x,z+.5) (-x+y,y,z+.5) (x,x-y,z+.5) (y,x,z) (x-y,-y,z) (-x,-x+y,z) \n 12 k ..m (x,0,z) (0,x,z) (-x,-x,z) (-x,0,z+.5) (0,-x,z+.5) (x,x,z+.5) (0,x,-z+.5) (x,0,-z+.5) (-x,-x,-z+.5) (0,-x,-z) (-x,0,-z) (x,x,-z) \n 12 j m.. (x,y,.25) (-y,x-y,.25) (-x+y,-x,.25) (-x,-y,.75) (y,-x+y,.75) (x-y,x,.75) (y,x,.25) (x-y,-y,.25) (-x,-x+y,.25) (-y,-x,.75) (-x+y,y,.75) (x,x-y,.75) \n 12 i ..2 (x,2*x,0) (-2*x,-x,0) (x,-x,0) (-x,-2*x,.5) (2*x,x,.5) (-x,x,.5) (-x,-2*x,0) (2*x,x,0) (-x,x,0) (x,2*x,.5) (-2*x,-x,.5) (x,-x,.5) \n 8 h 3.. (1/3,2/3,z) (2/3,1/3,z+.5) (2/3,1/3,-z+.5) (1/3,2/3,-z) (2/3,1/3,-z) (1/3,2/3,-z+.5) (1/3,2/3,z+.5) (2/3,1/3,z) \n 6 g m2m (x,0,.25) (0,x,.25) (-x,-x,.25) (-x,0,.75) (0,-x,.75) (x,x,.75) \n 6 f ..2/m (.5,0,0) (0,.5,0) (.5,.5,0) (.5,0,.5) (0,.5,.5) (.5,.5,.5) \n 4 e 3.m (0,0,z) (0,0,z+.5) (0,0,-z+.5) (0,0,-z) \n 4 d 3.2 (1/3,2/3,0) (2/3,1/3,.5) (2/3,1/3,0) (1/3,2/3,.5) \n 4 c -6.. (1/3,2/3,.25) (2/3,1/3,.75) (2/3,1/3,.25) (1/3,2/3,.75) \n 2 b -3.m (0,0,0) (0,0,.5) \n 2 a -62m (0,0,.25) (0,0,.75) \n ";

  string sg194 = "Space Group 194 \n (0,0,0) \n 24 l 1 (x,y,z) (-y,x-y,z) (-x+y,-x,z) (-x,-y,z+.5) (y,-x+y,z+.5) (x-y,x,z+.5) (y,x,-z) (x-y,-y,-z) (-x,-x+y,-z) (-y,-x,-z+.5) (-x+y,y,-z+.5) (x,x-y,-z+.5) (-x,-y,-z) (y,-x+y,-z) (x-y,x,-z) (x,y,-z+.5) (-y,x-y,-z+.5) (-x+y,-x,-z+.5) (-y,-x,z) (-x+y,y,z) (x,x-y,z) (y,x,z+.5) (x-y,-y,z+.5) (-x,-x+y,z+.5) \n 12 k .m. (x,2*x,z) (-2*x,-x,z) (x,-x,z) (-x,-2*x,z+.5) (2*x,x,z+.5) (-x,x,z+.5) (2*x,x,-z) (-x,-2*x,-z) (-x,x,-z) (-2*x,-x,-z+.5) (x,2*x,-z+.5) (x,-x,-z+.5) \n 12 j m.. (x,y,.25) (-y,x-y,.25) (-x+y,-x,.25) (-x,-y,.75) (y,-x+y,.75) (x-y,x,.75) (y,x,.75) (x-y,-y,.75) (-x,-x+y,.75) (-y,-x,.25) (-x+y,y,.25) (x,x-y,.25) \n 12 i .2. (x,0,0) (0,x,0) (-x,-x,0) (-x,0,.5) (0,-x,.5) (x,x,.5) (-x,0,0) (0,-x,0) (x,x,0) (x,0,.5) (0,x,.5) (-x,-x,.5) \n 6 h mm2 (x,2*x,.25) (-2*x,-x,.25) (x,-x,.25) (-x,-2*x,.75) (2*x,x,.75) (-x,x,.75) \n 6 g .2/m. (.5,0,0) (0,.5,0) (.5,.5,0) (.5,0,.5) (0,.5,.5) (.5,.5,.5) \n 4 f 3m. (1/3,2/3,z) (2/3,1/3,z+.5) (2/3,1/3,-z) (1/3,2/3,-z+.5) \n 4 e 3m. (0,0,z) (0,0,z+.5) (0,0,-z) (0,0,-z+.5) \n 2 d -6m2 (1/3,2/3,.75) (2/3,1/3,.25) \n 2 c -6m2 (1/3,2/3,.25) (2/3,1/3,.75) \n 2 b -6m2 (0,0,.25) (0,0,.75) \n 2 a -3m. (0,0,0) (0,0,.5) \n ";

  string sg195 = "Space Group 195 \n (0,0,0) \n 12 j 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) \n 6 i 2.. (x,.5,.5) (-x,.5,.5) (.5,x,.5) (.5,-x,.5) (.5,.5,x) (.5,.5,-x) \n 6 h 2.. (x,.5,0) (-x,.5,0) (0,x,.5) (0,-x,.5) (.5,0,x) (.5,0,-x) \n 6 g 2.. (x,0,.5) (-x,0,.5) (.5,x,0) (.5,-x,0) (0,.5,x) (0,.5,-x) \n 6 f 2.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) \n 4 e .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) \n 3 d 222.. (.5,0,0) (0,.5,0) (0,0,.5) \n 3 c 222.. (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 1 b 23. (.5,.5,.5) \n 1 a 23. (0,0,0) \n ";

  string sg196 = "Space Group 196 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 48 h 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) \n 24 g 2.. (x,.25,.25) (-x,.75,.25) (.25,x,.25) (.25,-x,.75) (.25,.25,x) (.75,.25,-x) \n 24 f 2.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) \n 16 e .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) \n 4 d 23. (.75,.75,.75) \n 4 c 23. (.25,.25,.25) \n 4 b 23. (.5,.5,.5) \n 4 a 23. (0,0,0) \n ";

  string sg197 = "Space Group 197 \n (0,0,0) (.5,.5,.5) \n 24 f 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) \n 12 e 2.. (x,.5,0) (-x,.5,0) (0,x,.5) (0,-x,.5) (.5,0,x) (.5,0,-x) \n 12 d 2.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) \n 8 c .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) \n 6 b 222.. (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 2 a 23. (0,0,0) \n ";

  string sg198 = "Space Group 198 \n (0,0,0) \n 12 b 1 (x,y,z) (-x+.5,-y,z+.5) (-x,y+.5,-z+.5) (x+.5,-y+.5,-z) (z,x,y) (z+.5,-x+.5,-y) (-z+.5,-x,y+.5) (-z,x+.5,-y+.5) (y,z,x) (-y,z+.5,-x+.5) (y+.5,-z+.5,-x) (-y+.5,-z,x+.5) \n 4 a .3. (x,x,x) (-x+.5,-x,x+.5) (-x,x+.5,-x+.5) (x+.5,-x+.5,-x) \n ";

  string sg199 = "Space Group 199 \n (0,0,0) (.5,.5,.5) \n 24 c 1 (x,y,z) (-x+.5,-y,z+.5) (-x,y+.5,-z+.5) (x+.5,-y+.5,-z) (z,x,y) (z+.5,-x+.5,-y) (-z+.5,-x,y+.5) (-z,x+.5,-y+.5) (y,z,x) (-y,z+.5,-x+.5) (y+.5,-z+.5,-x) (-y+.5,-z,x+.5) \n 12 b 2.. (x,0,.25) (-x+.5,0,.75) (.25,x,0) (.75,-x+.5,0) (0,.25,x) (0,.75,-x+.5) \n 8 a .3. (x,x,x) (-x+.5,-x,x+.5) (-x,x+.5,-x+.5) (x+.5,-x+.5,-x) \n ";

  string sg200 = "Space Group 200 \n (0,0,0) \n 24 l 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (-x,-y,-z) (x,y,-z) (x,-y,z) (-x,y,z) (-z,-x,-y) (-z,x,y) (z,x,-y) (z,-x,y) (-y,-z,-x) (y,-z,x) (-y,z,x) (y,z,-x) \n 12 k m.. (.5,y,z) (.5,-y,z) (.5,y,-z) (.5,-y,-z) (z,.5,y) (z,.5,-y) (-z,.5,y) (-z,.5,-y) (y,z,.5) (-y,z,.5) (y,-z,.5) (-y,-z,.5) \n 12 j m.. (0,y,z) (0,-y,z) (0,y,-z) (0,-y,-z) (z,0,y) (z,0,-y) (-z,0,y) (-z,0,-y) (y,z,0) (-y,z,0) (y,-z,0) (-y,-z,0) \n 8 i .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (-x,-x,-x) (x,x,-x) (x,-x,x) (-x,x,x) \n 6 h mm2.. (x,.5,.5) (-x,.5,.5) (.5,x,.5) (.5,-x,.5) (.5,.5,x) (.5,.5,-x) \n 6 g mm2.. (x,.5,0) (-x,.5,0) (0,x,.5) (0,-x,.5) (.5,0,x) (.5,0,-x) \n 6 f mm2.. (x,0,.5) (-x,0,.5) (.5,x,0) (.5,-x,0) (0,.5,x) (0,.5,-x) \n 6 e mm2.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) \n 3 d mmm.. (.5,0,0) (0,.5,0) (0,0,.5) \n 3 c mmm.. (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 1 b m-3. (.5,.5,.5) \n 1 a m-3. (0,0,0) \n ";

  //DX20180723 - added origin choice 2
  string sg201 = "";
  if(axis_cell==setting_1){
    sg201 = "Space Group 201 \n (0,0,0) \n 24 h 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) (-z+.5,-x+.5,-y+.5) (-z+.5,x+.5,y+.5) (z+.5,x+.5,-y+.5) (z+.5,-x+.5,y+.5) (-y+.5,-z+.5,-x+.5) (y+.5,-z+.5,x+.5) (-y+.5,z+.5,x+.5) (y+.5,z+.5,-x+.5) \n 12 g 2.. (x,.5,0) (-x,.5,0) (0,x,.5) (0,-x,.5) (.5,0,x) (.5,0,-x) (-x+.5,0,.5) (x+.5,0,.5) (.5,-x+.5,0) (.5,x+.5,0) (0,.5,-x+.5) (0,.5,x+.5) \n 12 f 2.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) (-x+.5,.5,.5) (x+.5,.5,.5) (.5,-x+.5,.5) (.5,x+.5,.5) (.5,.5,-x+.5) (.5,.5,x+.5) \n 8 e .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (-x+.5,-x+.5,-x+.5) (x+.5,x+.5,-x+.5) (x+.5,-x+.5,x+.5) (-x+.5,x+.5,x+.5) \n 6 d 222.. (0,.5,.5) (.5,0,.5) (.5,.5,0) (.5,0,0) (0,.5,0) (0,0,.5) \n 4 c .-3. (.75,.75,.75) (.25,.25,.75) (.25,.75,.25) (.75,.25,.25) \n 4 b .-3. (.25,.25,.25) (.75,.75,.25) (.75,.25,.75) (.25,.75,.75) \n 2 a 23. (0,0,0) (.5,.5,.5) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg201 = "Space Group 201 \n (0,0,0) \n 24 h 1 (x,y,z) (-x+.5,-y+.5,z) (-x+.5,y,-z+.5) (x,-y+.5,-z+.5) (z,x,y) (z,-x+.5,-y+.5) (-z+.5,-x+.5,y) (-z+.5,x,-y+.5) (y,z,x) (-y+.5,z,-x+.5) (y,-z+.5,-x+.5) (-y+.5,-z+.5,x) (-x,-y,-z) (x+.5,y+.5,-z) (x+.5,-y,z+.5) (-x,y+.5,z+.5) (-z,-x,-y) (-z,x+.5,y+.5) (z+.5,x+.5,-y) (z+.5,-x,y+.5) (-y,-z,-x) (y+.5,-z,x+.5) (-y,z+.5,x+.5) (y+.5,z+.5,-x) \n 12 g 2.. (x,.75,.25) (-x+.5,.75,.25) (.25,x,.75) (.25,-x+.5,.75) (.75,.25,x) (.75,.25,-x+.5) (-x,.25,.75) (x+.5,.25,.75) (.75,-x,.25) (.75,x+.5,.25) (.25,.75,-x) (.25,.75,x+.5) \n 12 f 2.. (x,.25,.25) (-x+.5,.25,.25) (.25,x,.25) (.25,-x+.5,.25) (.25,.25,x) (.25,.25,-x+.5) (-x,.75,.75) (x+.5,.75,.75) (.75,-x,.75) (.75,x+.5,.75) (.75,.75,-x) (.75,.75,x+.5) \n 8 e .3. (x,x,x) (-x+.5,-x+.5,x) (-x+.5,x,-x+.5) (x,-x+.5,-x+.5) (-x,-x,-x) (x+.5,x+.5,-x) (x+.5,-x,x+.5) (-x,x+.5,x+.5) \n 6 d 222.. (.25,.75,.75) (.75,.25,.75) (.75,.75,.25) (.75,.25,.25) (.25,.75,.25) (.25,.25,.75) \n 4 c .-3. (.5,.5,.5) (0,0,.5) (0,.5,0) (.5,0,0) \n 4 b .-3. (0,0,0) (.5,.5,0) (.5,0,.5) (0,.5,.5) \n 2 a 23. (.25,.25,.25) (.75,.75,.75) \n ";
  }
  else {
    sg201 = "Space Group 201 \n (0,0,0) \n 24 h 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) (-z+.5,-x+.5,-y+.5) (-z+.5,x+.5,y+.5) (z+.5,x+.5,-y+.5) (z+.5,-x+.5,y+.5) (-y+.5,-z+.5,-x+.5) (y+.5,-z+.5,x+.5) (-y+.5,z+.5,x+.5) (y+.5,z+.5,-x+.5) \n 12 g 2.. (x,.5,0) (-x,.5,0) (0,x,.5) (0,-x,.5) (.5,0,x) (.5,0,-x) (-x+.5,0,.5) (x+.5,0,.5) (.5,-x+.5,0) (.5,x+.5,0) (0,.5,-x+.5) (0,.5,x+.5) \n 12 f 2.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) (-x+.5,.5,.5) (x+.5,.5,.5) (.5,-x+.5,.5) (.5,x+.5,.5) (.5,.5,-x+.5) (.5,.5,x+.5) \n 8 e .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (-x+.5,-x+.5,-x+.5) (x+.5,x+.5,-x+.5) (x+.5,-x+.5,x+.5) (-x+.5,x+.5,x+.5) \n 6 d 222.. (0,.5,.5) (.5,0,.5) (.5,.5,0) (.5,0,0) (0,.5,0) (0,0,.5) \n 4 c .-3. (.75,.75,.75) (.25,.25,.75) (.25,.75,.25) (.75,.25,.25) \n 4 b .-3. (.25,.25,.25) (.75,.75,.25) (.75,.25,.75) (.25,.75,.75) \n 2 a 23. (0,0,0) (.5,.5,.5) \n ";
  }

  string sg202 = "Space Group 202 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 96 i 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (-x,-y,-z) (x,y,-z) (x,-y,z) (-x,y,z) (-z,-x,-y) (-z,x,y) (z,x,-y) (z,-x,y) (-y,-z,-x) (y,-z,x) (-y,z,x) (y,z,-x) \n 48 h m.. (0,y,z) (0,-y,z) (0,y,-z) (0,-y,-z) (z,0,y) (z,0,-y) (-z,0,y) (-z,0,-y) (y,z,0) (-y,z,0) (y,-z,0) (-y,-z,0) \n 48 g 2.. (x,.25,.25) (-x,.75,.25) (.25,x,.25) (.25,-x,.75) (.25,.25,x) (.75,.25,-x) (-x,.75,.75) (x,.25,.75) (.75,-x,.75) (.75,x,.25) (.75,.75,-x) (.25,.75,x) \n 32 f .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (-x,-x,-x) (x,x,-x) (x,-x,x) (-x,x,x) \n 24 e mm2.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) \n 24 d 2/m.. (0,.25,.25) (0,.75,.25) (.25,0,.25) (.25,0,.75) (.25,.25,0) (.75,.25,0) \n 8 c 23. (.25,.25,.25) (.75,.75,.75) \n 4 b m-3. (.5,.5,.5) \n 4 a m-3. (0,0,0) \n ";

  //DX20180723 - added origin choice 2
  string sg203 = "";
  if(axis_cell==setting_1){
    sg203 = "Space Group 203 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 96 g 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (-x+.25,-y+.25,-z+.25) (x+.25,y+.25,-z+.25) (x+.25,-y+.25,z+.25) (-x+.25,y+.25,z+.25) (-z+.25,-x+.25,-y+.25) (-z+.25,x+.25,y+.25) (z+.25,x+.25,-y+.25) (z+.25,-x+.25,y+.25) (-y+.25,-z+.25,-x+.25) (y+.25,-z+.25,x+.25) (-y+.25,z+.25,x+.25) (y+.25,z+.25,-x+.25) \n 48 f 2.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) (-x+.25,.25,.25) (x+.25,.25,.25) (.25,-x+.25,.25) (.25,x+.25,.25) (.25,.25,-x+.25) (.25,.25,x+.25) \n 32 e .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (-x+.25,-x+.25,-x+.25) (x+.25,x+.25,-x+.25) (x+.25,-x+.25,x+.25) (-x+.25,x+.25,x+.25) \n 16 d .-3. (.625,.625,.625) (.375,.375,.625) (.375,.625,.375) (.625,.375,.375) \n 16 c .-3. (.125,.125,.125) (.875,.875,.125) (.875,.125,.875) (.125,.875,.875) \n 8 b 23. (.5,.5,.5) (.75,.75,.75) \n 8 a 23. (0,0,0) (.25,.25,.25) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg203 = "Space Group 203 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 96 g 1 (x,y,z) (-x+.75,-y+.75,z) (-x+.75,y,-z+.75) (x,-y+.75,-z+.75) (z,x,y) (z,-x+.75,-y+.75) (-z+.75,-x+.75,y) (-z+.75,x,-y+.75) (y,z,x) (-y+.75,z,-x+.75) (y,-z+.75,-x+.75) (-y+.75,-z+.75,x) (-x,-y,-z) (x+.25,y+.25,-z) (x+.25,-y,z+.25) (-x,y+.25,z+.25) (-z,-x,-y) (-z,x+.25,y+.25) (z+.25,x+.25,-y) (z+.25,-x,y+.25) (-y,-z,-x) (y+.25,-z,x+.25) (-y,z+.25,x+.25) (y+.25,z+.25,-x) \n 48 f 2.. (x,.125,.125) (-x+.75,.625,.125) (.125,x,.125) (.125,-x+.75,.625) (.125,.125,x) (.625,.125,-x+.75) (-x,.875,.875) (x+.25,.375,.875) (.875,-x,.875) (.875,x+.25,.375) (.875,.875,-x) (.375,.875,x+.25) \n 32 e .3. (x,x,x) (-x+.75,-x+.75,x) (-x+.75,x,-x+.75) (x,-x+.75,-x+.75) (-x,-x,-x) (x+.25,x+.25,-x) (x+.25,-x,x+.25) (-x,x+.25,x+.25) \n 16 d .-3. (.5,.5,.5) (.25,.25,.5) (.25,.5,.25) (.5,.25,.25) \n 16 c .-3. (0,0,0) (.75,.75,0) (.75,0,.75) (0,.75,.75) \n 8 b 23. (.625,.625,.625) (.375,.375,.375) \n 8 a 23. (.125,.125,.125) (.875,.875,.875) \n ";
  }
  else {
    sg203 = "Space Group 203 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 96 g 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (-x+.25,-y+.25,-z+.25) (x+.25,y+.25,-z+.25) (x+.25,-y+.25,z+.25) (-x+.25,y+.25,z+.25) (-z+.25,-x+.25,-y+.25) (-z+.25,x+.25,y+.25) (z+.25,x+.25,-y+.25) (z+.25,-x+.25,y+.25) (-y+.25,-z+.25,-x+.25) (y+.25,-z+.25,x+.25) (-y+.25,z+.25,x+.25) (y+.25,z+.25,-x+.25) \n 48 f 2.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) (-x+.25,.25,.25) (x+.25,.25,.25) (.25,-x+.25,.25) (.25,x+.25,.25) (.25,.25,-x+.25) (.25,.25,x+.25) \n 32 e .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (-x+.25,-x+.25,-x+.25) (x+.25,x+.25,-x+.25) (x+.25,-x+.25,x+.25) (-x+.25,x+.25,x+.25) \n 16 d .-3. (.625,.625,.625) (.375,.375,.625) (.375,.625,.375) (.625,.375,.375) \n 16 c .-3. (.125,.125,.125) (.875,.875,.125) (.875,.125,.875) (.125,.875,.875) \n 8 b 23. (.5,.5,.5) (.75,.75,.75) \n 8 a 23. (0,0,0) (.25,.25,.25) \n ";
  }

  string sg204 = "Space Group 204 \n (0,0,0) (.5,.5,.5) \n 48 h 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (-x,-y,-z) (x,y,-z) (x,-y,z) (-x,y,z) (-z,-x,-y) (-z,x,y) (z,x,-y) (z,-x,y) (-y,-z,-x) (y,-z,x) (-y,z,x) (y,z,-x) \n 24 g m.. (0,y,z) (0,-y,z) (0,y,-z) (0,-y,-z) (z,0,y) (z,0,-y) (-z,0,y) (-z,0,-y) (y,z,0) (-y,z,0) (y,-z,0) (-y,-z,0) \n 16 f .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (-x,-x,-x) (x,x,-x) (x,-x,x) (-x,x,x) \n 12 e mm2.. (x,0,.5) (-x,0,.5) (.5,x,0) (.5,-x,0) (0,.5,x) (0,.5,-x) \n 12 d mm2.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) \n 8 c .-3. (.25,.25,.25) (.75,.75,.25) (.75,.25,.75) (.25,.75,.75) \n 6 b mmm.. (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 2 a m-3. (0,0,0) \n ";

  string sg205 = "Space Group 205 \n (0,0,0) \n 24 d 1 (x,y,z) (-x+.5,-y,z+.5) (-x,y+.5,-z+.5) (x+.5,-y+.5,-z) (z,x,y) (z+.5,-x+.5,-y) (-z+.5,-x,y+.5) (-z,x+.5,-y+.5) (y,z,x) (-y,z+.5,-x+.5) (y+.5,-z+.5,-x) (-y+.5,-z,x+.5) (-x,-y,-z) (x+.5,y,-z+.5) (x,-y+.5,z+.5) (-x+.5,y+.5,z) (-z,-x,-y) (-z+.5,x+.5,y) (z+.5,x,-y+.5) (z,-x+.5,y+.5) (-y,-z,-x) (y,-z+.5,x+.5) (-y+.5,z+.5,x) (y+.5,z,-x+.5) \n 8 c .3. (x,x,x) (-x+.5,-x,x+.5) (-x,x+.5,-x+.5) (x+.5,-x+.5,-x) (-x,-x,-x) (x+.5,x,-x+.5) (x,-x+.5,x+.5) (-x+.5,x+.5,x) \n 4 b .-3. (.5,.5,.5) (0,.5,0) (.5,0,0) (0,0,.5) \n 4 a .-3. (0,0,0) (.5,0,.5) (0,.5,.5) (.5,.5,0) \n ";

  string sg206 = "Space Group 206 \n (0,0,0) (.5,.5,.5) \n 48 e 1 (x,y,z) (-x+.5,-y,z+.5) (-x,y+.5,-z+.5) (x+.5,-y+.5,-z) (z,x,y) (z+.5,-x+.5,-y) (-z+.5,-x,y+.5) (-z,x+.5,-y+.5) (y,z,x) (-y,z+.5,-x+.5) (y+.5,-z+.5,-x) (-y+.5,-z,x+.5) (-x,-y,-z) (x+.5,y,-z+.5) (x,-y+.5,z+.5) (-x+.5,y+.5,z) (-z,-x,-y) (-z+.5,x+.5,y) (z+.5,x,-y+.5) (z,-x+.5,y+.5) (-y,-z,-x) (y,-z+.5,x+.5) (-y+.5,z+.5,x) (y+.5,z,-x+.5) \n 24 d 2.. (x,0,.25) (-x+.5,0,.75) (.25,x,0) (.75,-x+.5,0) (0,.25,x) (0,.75,-x+.5) (-x,0,.75) (x+.5,0,.25) (.75,-x,0) (.25,x+.5,0) (0,.75,-x) (0,.25,x+.5) \n 16 c .3. (x,x,x) (-x+.5,-x,x+.5) (-x,x+.5,-x+.5) (x+.5,-x+.5,-x) (-x,-x,-x) (x+.5,x,-x+.5) (x,-x+.5,x+.5) (-x+.5,x+.5,x) \n 8 b .-3. (.25,.25,.25) (.25,.75,.75) (.75,.75,.25) (.75,.25,.75) \n 8 a .-3. (0,0,0) (.5,0,.5) (0,.5,.5) (.5,.5,0) \n ";

  string sg207 = "Space Group 207 \n (0,0,0) \n 24 k 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y,x,-z) (-y,-x,-z) (y,-x,z) (-y,x,z) (x,z,-y) (-x,z,y) (-x,-z,-y) (x,-z,y) (z,y,-x) (z,-y,x) (-z,y,x) (-z,-y,-x) \n 12 j ..2 (.5,y,y) (.5,-y,y) (.5,y,-y) (.5,-y,-y) (y,.5,y) (y,.5,-y) (-y,.5,y) (-y,.5,-y) (y,y,.5) (-y,y,.5) (y,-y,.5) (-y,-y,.5) \n 12 i ..2 (0,y,y) (0,-y,y) (0,y,-y) (0,-y,-y) (y,0,y) (y,0,-y) (-y,0,y) (-y,0,-y) (y,y,0) (-y,y,0) (y,-y,0) (-y,-y,0) \n 12 h 2.. (x,.5,0) (-x,.5,0) (0,x,.5) (0,-x,.5) (.5,0,x) (.5,0,-x) (.5,x,0) (.5,-x,0) (x,0,.5) (-x,0,.5) (0,.5,-x) (0,.5,x) \n 8 g .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x,x,-x) (-x,-x,-x) (x,-x,x) (-x,x,x) \n 6 f 4.. (x,.5,.5) (-x,.5,.5) (.5,x,.5) (.5,-x,.5) (.5,.5,x) (.5,.5,-x) \n 6 e 4.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) \n 3 d 42.2 (.5,0,0) (0,.5,0) (0,0,.5) \n 3 c 42.2 (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 1 b 432 (.5,.5,.5) \n 1 a 432 (0,0,0) \n ";

  string sg208 = "Space Group 208 \n (0,0,0) \n 24 m 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y+.5,x+.5,-z+.5) (-y+.5,-x+.5,-z+.5) (y+.5,-x+.5,z+.5) (-y+.5,x+.5,z+.5) (x+.5,z+.5,-y+.5) (-x+.5,z+.5,y+.5) (-x+.5,-z+.5,-y+.5) (x+.5,-z+.5,y+.5) (z+.5,y+.5,-x+.5) (z+.5,-y+.5,x+.5) (-z+.5,y+.5,x+.5) (-z+.5,-y+.5,-x+.5) \n 12 l ..2 (.25,y,y+.5) (.75,-y,y+.5) (.75,y,-y+.5) (.25,-y,-y+.5) (y+.5,.25,y) (y+.5,.75,-y) (-y+.5,.75,y) (-y+.5,.25,-y) (y,y+.5,.25) (-y,y+.5,.75) (y,-y+.5,.75) (-y,-y+.5,.25) \n 12 k ..2 (.25,y,-y+.5) (.75,-y,-y+.5) (.75,y,y+.5) (.25,-y,y+.5) (-y+.5,.25,y) (-y+.5,.75,-y) (y+.5,.75,y) (y+.5,.25,-y) (y,-y+.5,.25) (-y,-y+.5,.75) (y,y+.5,.75) (-y,y+.5,.25) \n 12 j 2.. (x,.5,0) (-x,.5,0) (0,x,.5) (0,-x,.5) (.5,0,x) (.5,0,-x) (0,x+.5,.5) (0,-x+.5,.5) (x+.5,.5,0) (-x+.5,.5,0) (.5,0,-x+.5) (.5,0,x+.5) \n 12 i 2.. (x,0,.5) (-x,0,.5) (.5,x,0) (.5,-x,0) (0,.5,x) (0,.5,-x) (.5,x+.5,0) (.5,-x+.5,0) (x+.5,0,.5) (-x+.5,0,.5) (0,.5,-x+.5) (0,.5,x+.5) \n 12 h 2.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) (.5,x+.5,.5) (.5,-x+.5,.5) (x+.5,.5,.5) (-x+.5,.5,.5) (.5,.5,-x+.5) (.5,.5,x+.5) \n 8 g .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x+.5,x+.5,-x+.5) (-x+.5,-x+.5,-x+.5) (x+.5,-x+.5,x+.5) (-x+.5,x+.5,x+.5) \n 6 f 2.22 (.25,.5,0) (.75,.5,0) (0,.25,.5) (0,.75,.5) (.5,0,.25) (.5,0,.75) \n 6 e 2.22 (.25,0,.5) (.75,0,.5) (.5,.25,0) (.5,.75,0) (0,.5,.25) (0,.5,.75) \n 6 d 222.. (0,.5,.5) (.5,0,.5) (.5,.5,0) (0,.5,0) (.5,0,0) (0,0,.5) \n 4 c .32 (.75,.75,.75) (.25,.25,.75) (.25,.75,.25) (.75,.25,.25) \n 4 b .32 (.25,.25,.25) (.75,.75,.25) (.75,.25,.75) (.25,.75,.75) \n 2 a 23. (0,0,0) (.5,.5,.5) \n ";

  string sg209 = "Space Group 209 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 96 j 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y,x,-z) (-y,-x,-z) (y,-x,z) (-y,x,z) (x,z,-y) (-x,z,y) (-x,-z,-y) (x,-z,y) (z,y,-x) (z,-y,x) (-z,y,x) (-z,-y,-x) \n 48 i 2.. (x,.25,.25) (-x,.75,.25) (.25,x,.25) (.25,-x,.75) (.25,.25,x) (.75,.25,-x) (.25,x,.75) (.75,-x,.75) (x,.25,.75) (-x,.25,.25) (.25,.25,-x) (.25,.75,x) \n 48 h ..2 (.5,y,y) (.5,-y,y) (.5,y,-y) (.5,-y,-y) (y,.5,y) (y,.5,-y) (-y,.5,y) (-y,.5,-y) (y,y,.5) (-y,y,.5) (y,-y,.5) (-y,-y,.5) \n 48 g ..2 (0,y,y) (0,-y,y) (0,y,-y) (0,-y,-y) (y,0,y) (y,0,-y) (-y,0,y) (-y,0,-y) (y,y,0) (-y,y,0) (y,-y,0) (-y,-y,0) \n 32 f .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x,x,-x) (-x,-x,-x) (x,-x,x) (-x,x,x) \n 24 e 4.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) \n 24 d 2.22 (0,.25,.25) (0,.75,.25) (.25,0,.25) (.25,0,.75) (.25,.25,0) (.75,.25,0) \n 8 c 23. (.25,.25,.25) (.25,.25,.75) \n 4 b 432 (.5,.5,.5) \n 4 a 432 (0,0,0) \n "; //DX20201105 - fixed formatting

  string sg210 = "Space Group 210 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 96 h 1 (x,y,z) (-x,-y+.5,z+.5) (-x+.5,y+.5,-z) (x+.5,-y,-z+.5) (z,x,y) (z+.5,-x,-y+.5) (-z,-x+.5,y+.5) (-z+.5,x+.5,-y) (y,z,x) (-y+.5,z+.5,-x) (y+.5,-z,-x+.5) (-y,-z+.5,x+.5) (y+.75,x+.25,-z+.75) (-y+.25,-x+.25,-z+.25) (y+.25,-x+.75,z+.75) (-y+.75,x+.75,z+.25) (x+.75,z+.25,-y+.75) (-x+.75,z+.75,y+.25) (-x+.25,-z+.25,-y+.25) (x+.25,-z+.75,y+.75) (z+.75,y+.25,-x+.75) (z+.25,-y+.75,x+.75) (-z+.75,y+.75,x+.25) (-z+.25,-y+.25,-x+.25) \n 48 g ..2 (.125,y,-y+.25) (.875,-y+.5,-y+.75) (.375,y+.5,y+.75) (.625,-y,y+.25) (-y+.25,.125,y) (-y+.75,.875,-y+.5) (y+.75,.375,y+.5) (y+.25,.625,-y) (y,-y+.25,.125) (-y+.5,-y+.75,.875) (y+.5,y+.75,.375) (-y,y+.25,.625) \n 48 f 2.. (x,0,0) (-x,.5,.5) (0,x,0) (.5,-x,.5) (0,0,x) (.5,.5,-x) (.75,x+.25,.75) (.25,-x+.25,.25) (x+.75,.25,.75) (-x+.75,.75,.25) (.75,.25,-x+.75) (.25,.75,x+.75) \n 32 e .3. (x,x,x) (-x,-x+.5,x+.5) (-x+.5,x+.5,-x) (x+.5,-x,-x+.5) (x+.75,x+.25,-x+.75) (-x+.25,-x+.25,-x+.25) (x+.25,-x+.75,x+.75) (-x+.75,x+.75,x+.25) \n 16 d .32 (.625,.625,.625) (.375,.875,.125) (.875,.125,.375) (.125,.375,.875) \n 16 c .32 (.125,.125,.125) (.875,.375,.625) (.375,.625,.875) (.625,.875,.375) \n 8 b 23. (.5,.5,.5) (.25,.75,.25) \n 8 a 23. (0,0,0) (.75,.25,.75) \n ";

  string sg211 = "Space Group 211 \n (0,0,0) (.5,.5,.5) \n 48 j 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y,x,-z) (-y,-x,-z) (y,-x,z) (-y,x,z) (x,z,-y) (-x,z,y) (-x,-z,-y) (x,-z,y) (z,y,-x) (z,-y,x) (-z,y,x) (-z,-y,-x) \n 24 i ..2 (.25,y,-y+.5) (.75,-y,-y+.5) (.75,y,y+.5) (.25,-y,y+.5) (-y+.5,.25,y) (-y+.5,.75,-y) (y+.5,.75,y) (y+.5,.25,-y) (y,-y+.5,.25) (-y,-y+.5,.75) (y,y+.5,.75) (-y,y+.5,.25) \n 24 h ..2 (0,y,y) (0,-y,y) (0,y,-y) (0,-y,-y) (y,0,y) (y,0,-y) (-y,0,y) (-y,0,-y) (y,y,0) (-y,y,0) (y,-y,0) (-y,-y,0) \n 24 g 2.. (x,.5,0) (-x,.5,0) (0,x,.5) (0,-x,.5) (.5,0,x) (.5,0,-x) (.5,x,0) (.5,-x,0) (x,0,.5) (-x,0,.5) (0,.5,-x) (0,.5,x) \n 16 f .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x,x,-x) (-x,-x,-x) (x,-x,x) (-x,x,x) \n 12 e 4.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) \n 12 d 2.22 (.25,.5,0) (.75,.5,0) (0,.25,.5) (0,.75,.5) (.5,0,.25) (.5,0,.75) \n 8 c .32 (.25,.25,.25) (.75,.75,.25) (.75,.25,.75) (.25,.75,.75) \n 6 b 42.2 (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 2 a 432 (0,0,0) \n ";

  string sg212 = "Space Group 212 \n (0,0,0) \n 24 e 1 (x,y,z) (-x+.5,-y,z+.5) (-x,y+.5,-z+.5) (x+.5,-y+.5,-z) (z,x,y) (z+.5,-x+.5,-y) (-z+.5,-x,y+.5) (-z,x+.5,-y+.5) (y,z,x) (-y,z+.5,-x+.5) (y+.5,-z+.5,-x) (-y+.5,-z,x+.5) (y+.25,x+.75,-z+.75) (-y+.25,-x+.25,-z+.25) (y+.75,-x+.75,z+.25) (-y+.75,x+.25,z+.75) (x+.25,z+.75,-y+.75) (-x+.75,z+.25,y+.75) (-x+.25,-z+.25,-y+.25) (x+.75,-z+.75,y+.25) (z+.25,y+.75,-x+.75) (z+.75,-y+.75,x+.25) (-z+.75,y+.25,x+.75) (-z+.25,-y+.25,-x+.25) \n 12 d ..2 (.125,y,-y+.25) (.375,-y,-y+.75) (.875,y+.5,y+.25) (.625,-y+.5,y+.75) (-y+.25,.125,y) (-y+.75,.375,-y) (y+.25,.875,y+.5) (y+.75,.625,-y+.5) (y,-y+.25,.125) (-y,-y+.75,.375) (y+.5,y+.25,.875) (-y+.5,y+.75,.625) \n 8 c .3. (x,x,x) (-x+.5,-x,x+.5) (-x,x+.5,-x+.5) (x+.5,-x+.5,-x) (x+.25,x+.75,-x+.75) (-x+.25,-x+.25,-x+.25) (x+.75,-x+.75,x+.25) (-x+.75,x+.25,x+.75) \n 4 b .32 (.625,.625,.625) (.875,.375,.125) (.375,.125,.875) (.125,.875,.375) \n 4 a .32 (.125,.125,.125) (.375,.875,.625) (.875,.625,.375) (.625,.375,.875) \n ";

  string sg213 = "Space Group 213 \n (0,0,0) \n 24 e 1 (x,y,z) (-x+.5,-y,z+.5) (-x,y+.5,-z+.5) (x+.5,-y+.5,-z) (z,x,y) (z+.5,-x+.5,-y) (-z+.5,-x,y+.5) (-z,x+.5,-y+.5) (y,z,x) (-y,z+.5,-x+.5) (y+.5,-z+.5,-x) (-y+.5,-z,x+.5) (y+.75,x+.25,-z+.25) (-y+.75,-x+.75,-z+.75) (y+.25,-x+.25,z+.75) (-y+.25,x+.75,z+.25) (x+.75,z+.25,-y+.25) (-x+.25,z+.75,y+.25) (-x+.75,-z+.75,-y+.75) (x+.25,-z+.25,y+.75) (z+.75,y+.25,-x+.25) (z+.25,-y+.25,x+.75) (-z+.25,y+.75,x+.25) (-z+.75,-y+.75,-x+.75) \n 12 d ..2 (.125,y,y+.25) (.375,-y,y+.75) (.875,y+.5,-y+.25) (.625,-y+.5,-y+.75) (y+.25,.125,y) (y+.75,.375,-y) (-y+.25,.875,y+.5) (-y+.75,.625,-y+.5) (y,y+.25,.125) (-y,y+.75,.375) (y+.5,-y+.25,.875) (-y+.5,-y+.75,.625) \n 8 c .3. (x,x,x) (-x+.5,-x,x+.5) (-x,x+.5,-x+.5) (x+.5,-x+.5,-x) (x+.75,x+.25,-x+.25) (-x+.75,-x+.75,-x+.75) (x+.25,-x+.25,x+.75) (-x+.25,x+.75,x+.25) \n 4 b .32 (.875,.875,.875) (.625,.125,.375) (.125,.375,.625) (.375,.625,.125) \n 4 a .32 (.375,.375,.375) (.125,.625,.875) (.625,.875,.125) (.875,.125,.625) \n ";

  string sg214 = "Space Group 214 \n (0,0,0) (.5,.5,.5) \n 48 i 1 (x,y,z) (-x+.5,-y,z+.5) (-x,y+.5,-z+.5) (x+.5,-y+.5,-z) (z,x,y) (z+.5,-x+.5,-y) (-z+.5,-x,y+.5) (-z,x+.5,-y+.5) (y,z,x) (-y,z+.5,-x+.5) (y+.5,-z+.5,-x) (-y+.5,-z,x+.5) (y+.75,x+.25,-z+.25) (-y+.75,-x+.75,-z+.75) (y+.25,-x+.25,z+.75) (-y+.25,x+.75,z+.25) (x+.75,z+.25,-y+.25) (-x+.25,z+.75,y+.25) (-x+.75,-z+.75,-y+.75) (x+.25,-z+.25,y+.75) (z+.75,y+.25,-x+.25) (z+.25,-y+.25,x+.75) (-z+.25,y+.75,x+.25) (-z+.75,-y+.75,-x+.75) \n 24 h ..2 (.125,y,-y+.25) (.375,-y,-y+.75) (.875,y+.5,y+.25) (.625,-y+.5,y+.75) (-y+.25,.125,y) (-y+.75,.375,-y) (y+.25,.875,y+.5) (y+.75,.625,-y+.5) (y,-y+.25,.125) (-y,-y+.75,.375) (y+.5,y+.25,.875) (-y+.5,y+.75,.625) \n 24 g ..2 (.125,y,y+.25) (.375,-y,y+.75) (.875,y+.5,-y+.25) (.625,-y+.5,-y+.75) (y+.25,.125,y) (y+.75,.375,-y) (-y+.25,.875,y+.5) (-y+.75,.625,-y+.5) (y,y+.25,.125) (-y,y+.75,.375) (y+.5,-y+.25,.875) (-y+.5,-y+.75,.625) \n 24 f 2.. (x,0,.25) (-x+.5,0,.75) (.25,x,0) (.75,-x+.5,0) (0,.25,x) (0,.75,-x+.5) (.75,x+.25,0) (.75,-x+.75,.5) (x+.75,.5,.25) (-x+.25,0,.25) (0,.25,-x+.25) (.5,.25,x+.75) \n 16 e .3. (x,x,x) (-x+.5,-x,x+.5) (-x,x+.5,-x+.5) (x+.5,-x+.5,-x) (x+.75,x+.25,-x+.25) (-x+.75,-x+.75,-x+.75) (x+.25,-x+.25,x+.75) (-x+.25,x+.75,x+.25) \n 12 d 2.22 (.625,0,.25) (.875,0,.75) (.25,.625,0) (.75,.875,0) (0,.25,.625) (0,.75,.875) \n 12 c 2.22 (.125,0,.25) (.375,0,.75) (.25,.125,0) (.75,.375,0) (0,.25,.125) (0,.75,.375) \n 8 b .32 (.875,.875,.875) (.625,.125,.375) (.125,.375,.625) (.375,.625,.125) \n 8 a .32 (.125,.125,.125) (.375,.875,.625) (.875,.625,.375) (.625,.375,.875) \n ";

  string sg215 = "Space Group 215 \n (0,0,0) \n 24 j 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y,x,z) (-y,-x,z) (y,-x,-z) (-y,x,-z) (x,z,y) (-x,z,-y) (-x,-z,y) (x,-z,-y) (z,y,x) (z,-y,-x) (-z,y,-x) (-z,-y,x) \n 12 i ..m (x,x,z) (-x,-x,z) (-x,x,-z) (x,-x,-z) (z,x,x) (z,-x,-x) (-z,-x,x) (-z,x,-x) (x,z,x) (-x,z,-x) (x,-z,-x) (-x,-z,x) \n 12 h 2.. (x,.5,0) (-x,.5,0) (0,x,.5) (0,-x,.5) (.5,0,x) (.5,0,-x) (.5,x,0) (.5,-x,0) (x,0,.5) (-x,0,.5) (0,.5,x) (0,.5,-x) \n 6 g 2.mm (x,.5,.5) (-x,.5,.5) (.5,x,.5) (.5,-x,.5) (.5,.5,x) (.5,.5,-x) \n 6 f 2.mm (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) \n 4 e .3m (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) \n 3 d -42.m (.5,0,0) (0,.5,0) (0,0,.5) \n 3 c -42.m (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 1 b -43m (.5,.5,.5) \n 1 a -43m (0,0,0) \n ";

  string sg216 = "Space Group 216 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 96 i 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y,x,z) (-y,-x,z) (y,-x,-z) (-y,x,-z) (x,z,y) (-x,z,-y) (-x,-z,y) (x,-z,-y) (z,y,x) (z,-y,-x) (-z,y,-x) (-z,-y,x) \n 48 h ..m (x,x,z) (-x,-x,z) (-x,x,-z) (x,-x,-z) (z,x,x) (z,-x,-x) (-z,-x,x) (-z,x,-x) (x,z,x) (-x,z,-x) (x,-z,-x) (-x,-z,x) \n 24 g 2.mm (x,.25,.25) (-x,.75,.25) (.25,x,.25) (.25,-x,.75) (.25,.25,x) (.75,.25,-x) \n 24 f 2.mm (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) \n 16 e .3m (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) \n 4 d -43m (.75,.75,.75) \n 4 c -43m (.25,.25,.25) \n 4 b -43m (.5,.5,.5) \n 4 a -43m (0,0,0) \n ";

  string sg217 = "Space Group 217 \n (0,0,0) (.5,.5,.5) \n 48 h 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y,x,z) (-y,-x,z) (y,-x,-z) (-y,x,-z) (x,z,y) (-x,z,-y) (-x,-z,y) (x,-z,-y) (z,y,x) (z,-y,-x) (-z,y,-x) (-z,-y,x) \n 24 g ..m (x,x,z) (-x,-x,z) (-x,x,-z) (x,-x,-z) (z,x,x) (z,-x,-x) (-z,-x,x) (-z,x,-x) (x,z,x) (-x,z,-x) (x,-z,-x) (-x,-z,x) \n 24 f 2.. (x,.5,0) (-x,.5,0) (0,x,.5) (0,-x,.5) (.5,0,x) (.5,0,-x) (.5,x,0) (.5,-x,0) (x,0,.5) (-x,0,.5) (0,.5,x) (0,.5,-x) \n 12 e 2.mm (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) \n 12 d -4.. (.25,.5,0) (.75,.5,0) (0,.25,.5) (0,.75,.5) (.5,0,.25) (.5,0,.75) \n 8 c .3m (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) \n 6 b -42.m (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 2 a -43m (0,0,0) \n ";

  string sg218 = "Space Group 218 \n (0,0,0) \n 24 i 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y+.5,x+.5,z+.5) (-y+.5,-x+.5,z+.5) (y+.5,-x+.5,-z+.5) (-y+.5,x+.5,-z+.5) (x+.5,z+.5,y+.5) (-x+.5,z+.5,-y+.5) (-x+.5,-z+.5,y+.5) (x+.5,-z+.5,-y+.5) (z+.5,y+.5,x+.5) (z+.5,-y+.5,-x+.5) (-z+.5,y+.5,-x+.5) (-z+.5,-y+.5,x+.5) \n 12 h 2.. (x,0,.5) (-x,0,.5) (.5,x,0) (.5,-x,0) (0,.5,x) (0,.5,-x) (.5,x+.5,0) (.5,-x+.5,0) (x+.5,0,.5) (-x+.5,0,.5) (0,.5,x+.5) (0,.5,-x+.5) \n 12 g 2.. (x,.5,0) (-x,.5,0) (0,x,.5) (0,-x,.5) (.5,0,x) (.5,0,-x) (0,x+.5,.5) (0,-x+.5,.5) (x+.5,.5,0) (-x+.5,.5,0) (.5,0,x+.5) (.5,0,-x+.5) \n 12 f 2.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) (.5,x+.5,.5) (.5,-x+.5,.5) (x+.5,.5,.5) (-x+.5,.5,.5) (.5,.5,x+.5) (.5,.5,-x+.5) \n 8 e .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x+.5,x+.5,x+.5) (-x+.5,-x+.5,x+.5) (x+.5,-x+.5,-x+.5) (-x+.5,x+.5,-x+.5) \n 6 d -4.. (.25,0,.5) (.75,0,.5) (.5,.25,0) (.5,.75,0) (0,.5,.25) (0,.5,.75) \n 6 c -4.. (.25,.5,0) (.75,.5,0) (0,.25,.5) (0,.75,.5) (.5,0,.25) (.5,0,.75) \n 6 b 222.. (0,.5,.5) (.5,0,.5) (.5,.5,0) (0,.5,0) (.5,0,0) (0,0,.5) \n 2 a 23. (0,0,0) (.5,.5,.5) \n ";

  string sg219 = "Space Group 219 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 96 h 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y+.5,x+.5,z+.5) (-y+.5,-x+.5,z+.5) (y+.5,-x+.5,-z+.5) (-y+.5,x+.5,-z+.5) (x+.5,z+.5,y+.5) (-x+.5,z+.5,-y+.5) (-x+.5,-z+.5,y+.5) (x+.5,-z+.5,-y+.5) (z+.5,y+.5,x+.5) (z+.5,-y+.5,-x+.5) (-z+.5,y+.5,-x+.5) (-z+.5,-y+.5,x+.5) \n 48 g 2.. (x,.25,.25) (-x,.75,.25) (.25,x,.25) (.25,-x,.75) (.25,.25,x) (.75,.25,-x) (.75,x+.5,.75) (.25,-x+.5,.75) (x+.5,.75,.75) (-x+.5,.75,.25) (.75,.75,x+.5) (.75,.25,-x+.5) \n 48 f 2.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) (.5,x+.5,.5) (.5,-x+.5,.5) (x+.5,.5,.5) (-x+.5,.5,.5) (.5,.5,x+.5) (.5,.5,-x+.5) \n 32 e .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x+.5,x+.5,x+.5) (-x+.5,-x+.5,x+.5) (x+.5,-x+.5,-x+.5) (-x+.5,x+.5,-x+.5) \n 24 d -4.. (.25,0,0) (.75,0,0) (0,.25,0) (0,.75,0) (0,0,.25) (0,0,.75) \n 24 c -4.. (0,.25,.25) (0,.75,.25) (.25,0,.25) (.25,0,.75) (.25,.25,0) (.75,.25,0) \n 8 b 23. (.25,.25,.25) (.75,.75,.75) \n 8 a 23. (0,0,0) (.5,.5,.5) \n ";

  string sg220 = "Space Group 220 \n (0,0,0) (.5,.5,.5) \n 48 e 1 (x,y,z) (-x+.5,-y,z+.5) (-x,y+.5,-z+.5) (x+.5,-y+.5,-z) (z,x,y) (z+.5,-x+.5,-y) (-z+.5,-x,y+.5) (-z,x+.5,-y+.5) (y,z,x) (-y,z+.5,-x+.5) (y+.5,-z+.5,-x) (-y+.5,-z,x+.5) (y+.25,x+.25,z+.25) (-y+.25,-x+.75,z+.75) (y+.75,-x+.25,-z+.75) (-y+.75,x+.75,-z+.25) (x+.25,z+.25,y+.25) (-x+.75,z+.75,-y+.25) (-x+.25,-z+.75,y+.75) (x+.75,-z+.25,-y+.75) (z+.25,y+.25,x+.25) (z+.75,-y+.25,-x+.75) (-z+.75,y+.75,-x+.25) (-z+.25,-y+.75,x+.75) \n 24 d 2.. (x,0,.25) (-x+.5,0,.75) (.25,x,0) (.75,-x+.5,0) (0,.25,x) (0,.75,-x+.5) (.25,x+.25,.5) (.25,-x+.75,0) (x+.25,.5,.25) (-x+.75,0,.25) (.5,.25,x+.25) (0,.25,-x+.75) \n 16 c .3. (x,x,x) (-x+.5,-x,x+.5) (-x,x+.5,-x+.5) (x+.5,-x+.5,-x) (x+.25,x+.25,x+.25) (-x+.25,-x+.75,x+.75) (x+.75,-x+.25,-x+.75) (-x+.75,x+.75,-x+.25) \n 12 b -4.. (.875,0,.25) (.625,0,.75) (.25,.875,0) (.75,.625,0) (0,.25,.875) (0,.75,.625) \n 12 a -4.. (.375,0,.25) (.125,0,.75) (.25,.375,0) (3/4,.125,0) (0,.25,.375) (0,.75,.125) \n ";

  string sg221 = "Space Group 221 \n (0,0,0) \n 48 n 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y,x,-z) (-y,-x,-z) (y,-x,z) (-y,x,z) (x,z,-y) (-x,z,y) (-x,-z,-y) (x,-z,y) (z,y,-x) (z,-y,x) (-z,y,x) (-z,-y,-x) (-x,-y,-z) (x,y,-z) (x,-y,z) (-x,y,z) (-z,-x,-y) (-z,x,y) (z,x,-y) (z,-x,y) (-y,-z,-x) (y,-z,x) (-y,z,x) (y,z,-x) (-y,-x,z) (y,x,z) (-y,x,-z) (y,-x,-z) (-x,-z,y) (x,-z,-y) (x,z,y) (-x,z,-y) (-z,-y,x) (-z,y,-x) (z,-y,-x) (z,y,x) \n 24 m ..m (x,x,z) (-x,-x,z) (-x,x,-z) (x,-x,-z) (z,x,x) (z,-x,-x) (-z,-x,x) (-z,x,-x) (x,z,x) (-x,z,-x) (x,-z,-x) (-x,-z,x) (x,x,-z) (-x,-x,-z) (x,-x,z) (-x,x,z) (x,z,-x) (-x,z,x) (-x,-z,-x) (x,-z,x) (z,x,-x) (z,-x,x) (-z,x,x) (-z,-x,-x) \n 24 l m.. (.5,y,z) (.5,-y,z) (.5,y,-z) (.5,-y,-z) (z,.5,y) (z,.5,-y) (-z,.5,y) (-z,.5,-y) (y,z,.5) (-y,z,.5) (y,-z,.5) (-y,-z,.5) (y,.5,-z) (-y,.5,-z) (y,.5,z) (-y,.5,z) (.5,z,-y) (.5,z,y) (.5,-z,-y) (.5,-z,y) (z,y,.5) (z,-y,.5) (-z,y,.5) (-z,-y,.5) \n 24 k m.. (0,y,z) (0,-y,z) (0,y,-z) (0,-y,-z) (z,0,y) (z,0,-y) (-z,0,y) (-z,0,-y) (y,z,0) (-y,z,0) (y,-z,0) (-y,-z,0) (y,0,-z) (-y,0,-z) (y,0,z) (-y,0,z) (0,z,-y) (0,z,y) (0,-z,-y) (0,-z,y) (z,y,0) (z,-y,0) (-z,y,0) (-z,-y,0) \n 12 j m.m2 (.5,y,y) (.5,-y,y) (.5,y,-y) (.5,-y,-y) (y,.5,y) (y,.5,-y) (-y,.5,y) (-y,.5,-y) (y,y,.5) (-y,y,.5) (y,-y,.5) (-y,-y,.5) \n 12 i m.m2 (0,y,y) (0,-y,y) (0,y,-y) (0,-y,-y) (y,0,y) (y,0,-y) (-y,0,y) (-y,0,-y) (y,y,0) (-y,y,0) (y,-y,0) (-y,-y,0) \n 12 h mm2.. (x,.5,0) (-x,.5,0) (0,x,.5) (0,-x,.5) (.5,0,x) (.5,0,-x) (.5,x,0) (.5,-x,0) (x,0,.5) (-x,0,.5) (0,.5,-x) (0,.5,x) \n 8 g .3m (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x,x,-x) (-x,-x,-x) (x,-x,x) (-x,x,x) \n 6 f 4m.m (x,.5,.5) (-x,.5,.5) (.5,x,.5) (.5,-x,.5) (.5,.5,x) (.5,.5,-x) \n 6 e 4m.m (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) \n 3 d 4/mm.m (.5,0,0) (0,.5,0) (0,0,.5) \n 3 c 4/mm.m (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 1 b m-3m (.5,.5,.5) \n 1 a m-3m (0,0,0) \n ";

  //DX20180723 - added origin choice 2
  string sg222 = "";
  if(axis_cell==setting_1){
    sg222 = "Space Group 222 \n (0,0,0) \n 48 i 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y,x,-z) (-y,-x,-z) (y,-x,z) (-y,x,z) (x,z,-y) (-x,z,y) (-x,-z,-y) (x,-z,y) (z,y,-x) (z,-y,x) (-z,y,x) (-z,-y,-x) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) (-z+.5,-x+.5,-y+.5) (-z+.5,x+.5,y+.5) (z+.5,x+.5,-y+.5) (z+.5,-x+.5,y+.5) (-y+.5,-z+.5,-x+.5) (y+.5,-z+.5,x+.5) (-y+.5,z+.5,x+.5) (y+.5,z+.5,-x+.5) (-y+.5,-x+.5,z+.5) (y+.5,x+.5,z+.5) (-y+.5,x+.5,-z+.5) (y+.5,-x+.5,-z+.5) (-x+.5,-z+.5,y+.5) (x+.5,-z+.5,-y+.5) (x+.5,z+.5,y+.5) (-x+.5,z+.5,-y+.5) (-z+.5,-y+.5,x+.5) (-z+.5,y+.5,-x+.5) (z+.5,-y+.5,-x+.5) (z+.5,y+.5,x+.5) \n 24 h ..2 (0,y,y) (0,-y,y) (0,y,-y) (0,-y,-y) (y,0,y) (y,0,-y) (-y,0,y) (-y,0,-y) (y,y,0) (-y,y,0) (y,-y,0) (-y,-y,0) (.5,-y+.5,-y+.5) (.5,y+.5,-y+.5) (.5,-y+.5,y+.5) (.5,y+.5,y+.5) (-y+.5,.5,-y+.5) (-y+.5,.5,y+.5) (y+.5,.5,-y+.5) (y+.5,.5,y+.5) (-y+.5,-y+.5,.5) (y+.5,-y+.5,.5) (-y+.5,y+.5,.5) (y+.5,y+.5,.5) \n 24 g 2.. (x,0,.5) (-x,0,.5) (.5,x,0) (.5,-x,0) (0,.5,x) (0,.5,-x) (0,x,.5) (0,-x,.5) (x,.5,0) (-x,.5,0) (.5,0,-x) (.5,0,x) (-x+.5,.5,0) (x+.5,.5,0) (0,-x+.5,.5) (0,x+.5,.5) (.5,0,-x+.5) (.5,0,x+.5) (.5,-x+.5,0) (.5,x+.5,0) (-x+.5,0,.5) (x+.5,0,.5) (0,.5,x+.5) (0,.5,-x+.5) \n 16 f .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x,x,-x) (-x,-x,-x) (x,-x,x) (-x,x,x) (-x+.5,-x+.5,-x+.5) (x+.5,x+.5,-x+.5) (x+.5,-x+.5,x+.5) (-x+.5,x+.5,x+.5) (-x+.5,-x+.5,x+.5) (x+.5,x+.5,x+.5) (-x+.5,x+.5,-x+.5) (x+.5,-x+.5,-x+.5) \n 12 e 4.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) (-x+.5,.5,.5) (x+.5,.5,.5) (.5,-x+.5,.5) (.5,x+.5,.5) (.5,.5,-x+.5) (.5,.5,x+.5) \n 12 d -4.. (.25,0,.5) (.75,0,.5) (.5,.25,0) (.5,.75,0) (0,.5,.25) (0,.5,.75) (0,.25,.5) (0,.75,.5) (.25,.5,0) (.75,.5,0) (.5,0,.75) (.5,0,.25) \n 8 c .-3. (.25,.25,.25) (.75,.75,.25) (.75,.25,.75) (.25,.75,.75) (.25,.25,.75) (.75,.75,.75) (.25,.75,.25) (.75,.25,.25) \n 6 b 42.2 (0,.5,.5) (.5,0,.5) (.5,.5,0) (.5,0,0) (0,.5,0) (0,0,.5) \n 2 a 432 (0,0,0) (.5,.5,.5) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg222 = "Space Group 222 \n (0,0,0) \n 48 i 1 (x,y,z) (-x+.5,-y+.5,z) (-x+.5,y,-z+.5) (x,-y+.5,-z+.5) (z,x,y) (z,-x+.5,-y+.5) (-z+.5,-x+.5,y) (-z+.5,x,-y+.5) (y,z,x) (-y+.5,z,-x+.5) (y,-z+.5,-x+.5) (-y+.5,-z+.5,x) (y,x,-z+.5) (-y+.5,-x+.5,-z+.5) (y,-x+.5,z) (-y+.5,x,z) (x,z,-y+.5) (-x+.5,z,y) (-x+.5,-z+.5,-y+.5) (x,-z+.5,y) (z,y,-x+.5) (z,-y+.5,x) (-z+.5,y,x) (-z+.5,-y+.5,-x+.5) (-x,-y,-z) (x+.5,y+.5,-z) (x+.5,-y,z+.5) (-x,y+.5,z+.5) (-z,-x,-y) (-z,x+.5,y+.5) (z+.5,x+.5,-y) (z+.5,-x,y+.5) (-y,-z,-x) (y+.5,-z,x+.5) (-y,z+.5,x+.5) (y+.5,z+.5,-x) (-y,-x,z+.5) (y+.5,x+.5,z+.5) (-y,x+.5,-z) (y+.5,-x,-z) (-x,-z,y+.5) (x+.5,-z,-y) (x+.5,z+.5,y+.5) (-x,z+.5,-y) (-z,-y,x+.5) (-z,y+.5,-x) (z+.5,-y,-x) (z+.5,y+.5,x+.5) \n 24 h ..2 (.25,y,y) (.25,-y+.5,y) (.25,y,-y+.5) (.25,-y+.5,-y+.5) (y,.25,y) (y,.25,-y+.5) (-y+.5,.25,y) (-y+.5,.25,-y+.5) (y,y,.25) (-y+.5,y,.25) (y,-y+.5,.25) (-y+.5,-y+.5,.25) (.75,-y,-y) (.75,y+.5,-y) (.75,-y,y+.5) (.75,y+.5,y+.5) (-y,.75,-y) (-y,.75,y+.5) (y+.5,.75,-y) (y+.5,.75,y+.5) (-y,-y,.75) (y+.5,-y,.75) (-y,y+.5,.75) (y+.5,y+.5,.75) \n 24 g 2.. (x,.75,.25) (-x+.5,.75,.25) (.25,x,.75) (.25,-x+.5,.75) (.75,.25,x) (.75,.25,-x+.5) (.75,x,.25) (.75,-x+.5,.25) (x,.25,.75) (-x+.5,.25,.75) (.25,.75,-x+.5) (.25,.75,x) (-x,.25,.75) (x+.5,.25,.75) (.75,-x,.25) (.75,x+.5,.25) (.25,.75,-x) (.25,.75,x+.5) (.25,-x,.75) (.25,x+.5,.75) (-x,.75,.25) (x+.5,.75,.25) (.75,.25,x+.5) (.75,.25,-x) \n 16 f .3. (x,x,x) (-x+.5,-x+.5,x) (-x+.5,x,-x+.5) (x,-x+.5,-x+.5) (x,x,-x+.5) (-x+.5,-x+.5,-x+.5) (x,-x+.5,x) (-x+.5,x,x) (-x,-x,-x) (x+.5,x+.5,-x) (x+.5,-x,x+.5) (-x,x+.5,x+.5) (-x,-x,x+.5) (x+.5,x+.5,x+.5) (-x,x+.5,-x) (x+.5,-x,-x) \n 12 e 4.. (x,.25,.25) (-x+.5,.25,.25) (.25,x,.25) (.25,-x+.5,.25) (.25,.25,x) (.25,.25,-x+.5) (-x,.75,.75) (x+.5,.75,.75) (.75,-x,.75) (.75,x+.5,.75) (.75,.75,-x) (.75,.75,x+.5) \n 12 d -4.. (0,.75,.25) (.5,.75,.25) (.25,0,.75) (.25,.5,.75) (.75,.25,0) (.75,.25,.5) (.75,0,.25) (.75,.5,.25) (0,.25,.75) (.5,.25,.75) (.25,.75,.5) (.25,.75,0) \n 8 c .-3. (0,0,0) (.5,.5,0) (.5,0,.5) (0,.5,.5) (0,0,.5) (.5,.5,.5) (0,.5,0) (.5,0,0) \n 6 b 42.2 (.75,.25,.25) (.25,.75,.25) (.25,.25,.75) (.25,.75,.75) (.75,.25,.75) (.75,.75,.25) \n 2 a 432 (.25,.25,.25) (.75,.75,.75) \n ";
  }
  else {
    sg222 = "Space Group 222 \n (0,0,0) \n 48 i 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y,x,-z) (-y,-x,-z) (y,-x,z) (-y,x,z) (x,z,-y) (-x,z,y) (-x,-z,-y) (x,-z,y) (z,y,-x) (z,-y,x) (-z,y,x) (-z,-y,-x) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) (-z+.5,-x+.5,-y+.5) (-z+.5,x+.5,y+.5) (z+.5,x+.5,-y+.5) (z+.5,-x+.5,y+.5) (-y+.5,-z+.5,-x+.5) (y+.5,-z+.5,x+.5) (-y+.5,z+.5,x+.5) (y+.5,z+.5,-x+.5) (-y+.5,-x+.5,z+.5) (y+.5,x+.5,z+.5) (-y+.5,x+.5,-z+.5) (y+.5,-x+.5,-z+.5) (-x+.5,-z+.5,y+.5) (x+.5,-z+.5,-y+.5) (x+.5,z+.5,y+.5) (-x+.5,z+.5,-y+.5) (-z+.5,-y+.5,x+.5) (-z+.5,y+.5,-x+.5) (z+.5,-y+.5,-x+.5) (z+.5,y+.5,x+.5) \n 24 h ..2 (0,y,y) (0,-y,y) (0,y,-y) (0,-y,-y) (y,0,y) (y,0,-y) (-y,0,y) (-y,0,-y) (y,y,0) (-y,y,0) (y,-y,0) (-y,-y,0) (.5,-y+.5,-y+.5) (.5,y+.5,-y+.5) (.5,-y+.5,y+.5) (.5,y+.5,y+.5) (-y+.5,.5,-y+.5) (-y+.5,.5,y+.5) (y+.5,.5,-y+.5) (y+.5,.5,y+.5) (-y+.5,-y+.5,.5) (y+.5,-y+.5,.5) (-y+.5,y+.5,.5) (y+.5,y+.5,.5) \n 24 g 2.. (x,0,.5) (-x,0,.5) (.5,x,0) (.5,-x,0) (0,.5,x) (0,.5,-x) (0,x,.5) (0,-x,.5) (x,.5,0) (-x,.5,0) (.5,0,-x) (.5,0,x) (-x+.5,.5,0) (x+.5,.5,0) (0,-x+.5,.5) (0,x+.5,.5) (.5,0,-x+.5) (.5,0,x+.5) (.5,-x+.5,0) (.5,x+.5,0) (-x+.5,0,.5) (x+.5,0,.5) (0,.5,x+.5) (0,.5,-x+.5) \n 16 f .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x,x,-x) (-x,-x,-x) (x,-x,x) (-x,x,x) (-x+.5,-x+.5,-x+.5) (x+.5,x+.5,-x+.5) (x+.5,-x+.5,x+.5) (-x+.5,x+.5,x+.5) (-x+.5,-x+.5,x+.5) (x+.5,x+.5,x+.5) (-x+.5,x+.5,-x+.5) (x+.5,-x+.5,-x+.5) \n 12 e 4.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) (-x+.5,.5,.5) (x+.5,.5,.5) (.5,-x+.5,.5) (.5,x+.5,.5) (.5,.5,-x+.5) (.5,.5,x+.5) \n 12 d -4.. (.25,0,.5) (.75,0,.5) (.5,.25,0) (.5,.75,0) (0,.5,.25) (0,.5,.75) (0,.25,.5) (0,.75,.5) (.25,.5,0) (.75,.5,0) (.5,0,.75) (.5,0,.25) \n 8 c .-3. (.25,.25,.25) (.75,.75,.25) (.75,.25,.75) (.25,.75,.75) (.25,.25,.75) (.75,.75,.75) (.25,.75,.25) (.75,.25,.25) \n 6 b 42.2 (0,.5,.5) (.5,0,.5) (.5,.5,0) (.5,0,0) (0,.5,0) (0,0,.5) \n 2 a 432 (0,0,0) (.5,.5,.5) \n ";
  }


  string sg223 = "Space Group 223 \n (0,0,0) \n 48 l 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y+.5,x+.5,-z+.5) (-y+.5,-x+.5,-z+.5) (y+.5,-x+.5,z+.5) (-y+.5,x+.5,z+.5) (x+.5,z+.5,-y+.5) (-x+.5,z+.5,y+.5) (-x+.5,-z+.5,-y+.5) (x+.5,-z+.5,y+.5) (z+.5,y+.5,-x+.5) (z+.5,-y+.5,x+.5) (-z+.5,y+.5,x+.5) (-z+.5,-y+.5,-x+.5) (-x,-y,-z) (x,y,-z) (x,-y,z) (-x,y,z) (-z,-x,-y) (-z,x,y) (z,x,-y) (z,-x,y) (-y,-z,-x) (y,-z,x) (-y,z,x) (y,z,-x) (-y+.5,-x+.5,z+.5) (y+.5,x+.5,z+.5) (-y+.5,x+.5,-z+.5) (y+.5,-x+.5,-z+.5) (-x+.5,-z+.5,y+.5) (x+.5,-z+.5,-y+.5) (x+.5,z+.5,y+.5) (-x+.5,z+.5,-y+.5) (-z+.5,-y+.5,x+.5) (-z+.5,y+.5,-x+.5) (z+.5,-y+.5,-x+.5) (z+.5,y+.5,x+.5) \n 24 k m.. (0,y,z) (0,-y,z) (0,y,-z) (0,-y,-z) (z,0,y) (z,0,-y) (-z,0,y) (-z,0,-y) (y,z,0) (-y,z,0) (y,-z,0) (-y,-z,0) (y+.5,.5,-z+.5) (-y+.5,.5,-z+.5) (y+.5,.5,z+.5) (-y+.5,.5,z+.5) (.5,z+.5,-y+.5) (.5,z+.5,y+.5) (.5,-z+.5,-y+.5) (.5,-z+.5,y+.5) (z+.5,y+.5,.5) (z+.5,-y+.5,.5) (-z+.5,y+.5,.5) (-z+.5,-y+.5,.5) \n 24 j ..2 (.25,y,y+.5) (.75,-y,y+.5) (.75,y,-y+.5) (.25,-y,-y+.5) (y+.5,.25,y) (y+.5,.75,-y) (-y+.5,.75,y) (-y+.5,.25,-y) (y,y+.5,.25) (-y,y+.5,.75) (y,-y+.5,.75) (-y,-y+.5,.25) (.75,-y,-y+.5) (.25,y,-y+.5) (.25,-y,y+.5) (.75,y,y+.5) (-y+.5,.75,-y) (-y+.5,.25,y) (y+.5,.25,-y) (y+.5,.75,y) (-y,-y+.5,.75) (y,-y+.5,.25) (-y,y+.5,.25) (y,y+.5,.75) \n 16 i .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x+.5,x+.5,-x+.5) (-x+.5,-x+.5,-x+.5) (x+.5,-x+.5,x+.5) (-x+.5,x+.5,x+.5) (-x,-x,-x) (x,x,-x) (x,-x,x) (-x,x,x) (-x+.5,-x+.5,x+.5) (x+.5,x+.5,x+.5) (-x+.5,x+.5,-x+.5) (x+.5,-x+.5,-x+.5) \n 12 h mm2.. (x,.5,0) (-x,.5,0) (0,x,.5) (0,-x,.5) (.5,0,x) (.5,0,-x) (0,x+.5,.5) (0,-x+.5,.5) (x+.5,.5,0) (-x+.5,.5,0) (.5,0,-x+.5) (.5,0,x+.5) \n 12 g mm2.. (x,0,.5) (-x,0,.5) (.5,x,0) (.5,-x,0) (0,.5,x) (0,.5,-x) (.5,x+.5,0) (.5,-x+.5,0) (x+.5,0,.5) (-x+.5,0,.5) (0,.5,-x+.5) (0,.5,x+.5) \n 12 f mm2.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) (.5,x+.5,.5) (.5,-x+.5,.5) (x+.5,.5,.5) (-x+.5,.5,.5) (.5,.5,-x+.5) (.5,.5,x+.5) \n 8 e .32 (.25,.25,.25) (.75,.75,.25) (.75,.25,.75) (.25,.75,.75) (.75,.75,.75) (.25,.25,.75) (.25,.75,.25) (.75,.25,.25) \n 6 d -4m.2 (.25,.5,0) (.75,.5,0) (0,.25,.5) (0,.75,.5) (.5,0,.25) (.5,0,.75) \n 6 c -4m.2 (.25,0,.5) (.75,0,.5) (.5,.25,0) (.5,.75,0) (0,.5,.25) (0,.5,.75) \n 6 b mmm.. (0,.5,.5) (.5,0,.5) (.5,.5,0) (0,.5,0) (.5,0,0) (0,0,.5) \n 2 a m-3. (0,0,0) (.5,.5,.5) \n "; //DX20200528 - missing space between 15th and 16th position for Wyckoff letter k

  //DX20180723 - added origin choice 2
  string sg224 = "";
  if(axis_cell==setting_1){
    sg224 = "Space Group 224 \n (0,0,0) \n 48 l 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y+.5,x+.5,-z+.5) (-y+.5,-x+.5,-z+.5) (y+.5,-x+.5,z+.5) (-y+.5,x+.5,z+.5) (x+.5,z+.5,-y+.5) (-x+.5,z+.5,y+.5) (-x+.5,-z+.5,-y+.5) (x+.5,-z+.5,y+.5) (z+.5,y+.5,-x+.5) (z+.5,-y+.5,x+.5) (-z+.5,y+.5,x+.5) (-z+.5,-y+.5,-x+.5) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) (-z+.5,-x+.5,-y+.5) (-z+.5,x+.5,y+.5) (z+.5,x+.5,-y+.5) (z+.5,-x+.5,y+.5) (-y+.5,-z+.5,-x+.5) (y+.5,-z+.5,x+.5) (-y+.5,z+.5,x+.5) (y+.5,z+.5,-x+.5) (-y,-x,z) (y,x,z) (-y,x,-z) (y,-x,-z) (-x,-z,y) (x,-z,-y) (x,z,y) (-x,z,-y) (-z,-y,x) (-z,y,-x) (z,-y,-x) (z,y,x) \n 24 k ..m (x,x,z) (-x,-x,z) (-x,x,-z) (x,-x,-z) (z,x,x) (z,-x,-x) (-z,-x,x) (-z,x,-x) (x,z,x) (-x,z,-x) (x,-z,-x) (-x,-z,x) (x+.5,x+.5,-z+.5) (-x+.5,-x+.5,-z+.5) (x+.5,-x+.5,z+.5) (-x+.5,x+.5,z+.5) (x+.5,z+.5,-x+.5) (-x+.5,z+.5,x+.5) (-x+.5,-z+.5,-x+.5) (x+.5,-z+.5,x+.5) (z+.5,x+.5,-x+.5) (z+.5,-x+.5,x+.5) (-z+.5,x+.5,x+.5) (-z+.5,-x+.5,-x+.5) \n 24 j ..2 (.25,y,y+.5) (.75,-y,y+.5) (.75,y,-y+.5) (.25,-y,-y+.5) (y+.5,.25,y) (y+.5,.75,-y) (-y+.5,.75,y) (-y+.5,.25,-y) (y,y+.5,.25) (-y,y+.5,.75) (y,-y+.5,.75) (-y,-y+.5,.25) (.25,-y+.5,-y) (.75,y+.5,-y) (.75,-y+.5,y) (.25,y+.5,y) (-y,.25,-y+.5) (-y,.75,y+.5) (y,.75,-y+.5) (y,.25,y+.5) (-y+.5,-y,.25) (y+.5,-y,.75) (-y+.5,y,.75) (y+.5,y,.25) \n 24 i ..2 (.25,y,-y+.5) (.75,-y,-y+.5) (.75,y,y+.5) (.25,-y,y+.5) (-y+.5,.25,y) (-y+.5,.75,-y) (y+.5,.75,y) (y+.5,.25,-y) (y,-y+.5,.25) (-y,-y+.5,.75) (y,y+.5,.75) (-y,y+.5,.25) (.25,-y+.5,y) (.75,y+.5,y) (.75,-y+.5,-y) (.25,y+.5,-y) (y,.25,-y+.5) (y,.75,y+.5) (-y,.75,-y+.5) (-y,.25,y+.5) (-y+.5,y,.25) (y+.5,y,.75) (-y+.5,-y,.75) (y+.5,-y,.25) \n 24 h 2.. (x,0,.5) (-x,0,.5) (.5,x,0) (.5,-x,0) (0,.5,x) (0,.5,-x) (.5,x+.5,0) (.5,-x+.5,0) (x+.5,0,.5) (-x+.5,0,.5) (0,.5,-x+.5) (0,.5,x+.5) (-x+.5,.5,0) (x+.5,.5,0) (0,-x+.5,.5) (0,x+.5,.5) (.5,0,-x+.5) (.5,0,x+.5) (0,-x,.5) (0,x,.5) (-x,.5,0) (x,.5,0) (.5,0,x) (.5,0,-x) \n 12 g 2.mm (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) (.5,x+.5,.5) (.5,-x+.5,.5) (x+.5,.5,.5) (-x+.5,.5,.5) (.5,.5,-x+.5) (.5,.5,x+.5) \n 12 f 2.22 (.25,0,.5) (.75,0,.5) (.5,.25,0) (.5,.75,0) (0,.5,.25) (0,.5,.75) (.25,.5,0) (.75,.5,0) (0,.25,.5) (0,.75,.5) (.5,0,.25) (.5,0,.75) \n 8 e .3m (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x+.5,x+.5,-x+.5) (-x+.5,-x+.5,-x+.5) (x+.5,-x+.5,x+.5) (-x+.5,x+.5,x+.5) \n 6 d -42.m (0,.5,.5) (.5,0,.5) (.5,.5,0) (0,.5,0) (.5,0,0) (0,0,.5) \n 4 c .-3m (.75,.75,.75) (.25,.25,.75) (.25,.75,.25) (.75,.25,.25) \n 4 b .-3m (.25,.25,.25) (.75,.75,.25) (.75,.25,.75) (.25,.75,.75) \n 2 a -43m (0,0,0) (.5,.5,.5) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg224 = "Space Group 224 \n (0,0,0) \n 48 l 1 (x,y,z) (-x+.5,-y+.5,z) (-x+.5,y,-z+.5) (x,-y+.5,-z+.5) (z,x,y) (z,-x+.5,-y+.5) (-z+.5,-x+.5,y) (-z+.5,x,-y+.5) (y,z,x) (-y+.5,z,-x+.5) (y,-z+.5,-x+.5) (-y+.5,-z+.5,x) (y+.5,x+.5,-z) (-y,-x,-z) (y+.5,-x,z+.5) (-y,x+.5,z+.5) (x+.5,z+.5,-y) (-x,z+.5,y+.5) (-x,-z,-y) (x+.5,-z,y+.5) (z+.5,y+.5,-x) (z+.5,-y,x+.5) (-z,y+.5,x+.5) (-z,-y,-x) (-x,-y,-z) (x+.5,y+.5,-z) (x+.5,-y,z+.5) (-x,y+.5,z+.5) (-z,-x,-y) (-z,x+.5,y+.5) (z+.5,x+.5,-y) (z+.5,-x,y+.5) (-y,-z,-x) (y+.5,-z,x+.5) (-y,z+.5,x+.5) (y+.5,z+.5,-x) (-y+.5,-x+.5,z) (y,x,z) (-y+.5,x,-z+.5) (y,-x+.5,-z+.5) (-x+.5,-z+.5,y) (x,-z+.5,-y+.5) (x,z,y) (-x+.5,z,-y+.5) (-z+.5,-y+.5,x) (-z+.5,y,-x+.5) (z,-y+.5,-x+.5) (z,y,x) \n 24 k ..m (x,x,z) (-x+.5,-x+.5,z) (-x+.5,x,-z+.5) (x,-x+.5,-z+.5) (z,x,x) (z,-x+.5,-x+.5) (-z+.5,-x+.5,x) (-z+.5,x,-x+.5) (x,z,x) (-x+.5,z,-x+.5) (x,-z+.5,-x+.5) (-x+.5,-z+.5,x) (x+.5,x+.5,-z) (-x,-x,-z) (x+.5,-x,z+.5) (-x,x+.5,z+.5) (x+.5,z+.5,-x) (-x,z+.5,x+.5) (-x,-z,-x) (x+.5,-z,x+.5) (z+.5,x+.5,-x) (z+.5,-x,x+.5) (-z,x+.5,x+.5) (-z,-x,-x) \n 24 j ..2 (.5,y,-y) (0,-y+.5,-y) (0,y,y+.5) (.5,-y+.5,y+.5) (-y,.5,y) (-y,0,-y+.5) (y+.5,0,y) (y+.5,.5,-y+.5) (y,-y,.5) (-y+.5,-y,0) (y,y+.5,0) (-y+.5,y+.5,.5) (.5,-y,y) (0,y+.5,y) (0,-y,-y+.5) (.5,y+.5,-y+.5) (y,.5,-y) (y,0,y+.5) (-y+.5,0,-y) (-y+.5,.5,y+.5) (-y,y,.5) (y+.5,y,0) (-y,-y+.5,0) (y+.5,-y+.5,.5) \n 24 i ..2 (.5,y,y+.5) (0,-y+.5,y+.5) (0,y,-y) (.5,-y+.5,-y) (y+.5,.5,y) (y+.5,0,-y+.5) (-y,0,y) (-y,.5,-y+.5) (y,y+.5,.5) (-y+.5,y+.5,0) (y,-y,0) (-y+.5,-y,.5) (.5,-y,-y+.5) (0,y+.5,-y+.5) (0,-y,y) (.5,y+.5,y) (-y+.5,.5,-y) (-y+.5,0,y+.5) (y,0,-y) (y,.5,y+.5) (-y,-y+.5,.5) (y+.5,-y+.5,0) (-y,y,0) (y+.5,y,.5) \n 24 h 2.. (x,.25,.75) (-x+.5,.25,.75) (.75,x,.25) (.75,-x+.5,.25) (.25,.75,x) (.25,.75,-x+.5) (.75,x+.5,.25) (.75,-x,.25) (x+.5,.25,.75) (-x,.25,.75) (.25,.75,-x) (.25,.75,x+.5) (-x,.75,.25) (x+.5,.75,.25) (.25,-x,.75) (.25,x+.5,.75) (.75,.25,-x) (.75,.25,x+.5) (.25,-x+.5,.75) (.25,x,.75) (-x+.5,.75,.25) (x,.75,.25) (.75,.25,x) (.75,.25,-x+.5) \n 12 g 2.mm (x,.25,.25) (-x+.5,.25,.25) (.25,x,.25) (.25,-x+.5,.25) (.25,.25,x) (.25,.25,-x+.5) (.75,x+.5,.75) (.75,-x,.75) (x+.5,.75,.75) (-x,.75,.75) (.75,.75,-x) (.75,.75,x+.5) \n 12 f 2.22 (.5,.25,.75) (0,.25,.75) (.75,.5,.25) (.75,0,.25) (.25,.75,.5) (.25,.75,0) (.5,.75,.25) (0,.75,.25) (.25,.5,.75) (.25,0,.75) (.75,.25,.5) (.75,.25,0) \n 8 e .3m (x,x,x) (-x+.5,-x+.5,x) (-x+.5,x,-x+.5) (x,-x+.5,-x+.5) (x+.5,x+.5,-x) (-x,-x,-x) (x+.5,-x,x+.5) (-x,x+.5,x+.5) \n 6 d -42.m (.25,.75,.75) (.75,.25,.75) (.75,.75,.25) (.25,.75,.25) (.75,.25,.25) (.25,.25,.75) \n 4 c .-3m (.5,.5,.5) (0,0,.5) (0,.5,0) (.5,0,0) \n 4 b .-3m (0,0,0) (.5,.5,0) (.5,0,.5) (0,.5,.5) \n 2 a -43m (.25,.25,.25) (.75,.75,.75) \n ";
  }
  else {
    sg224 = "Space Group 224 \n (0,0,0) \n 48 l 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y+.5,x+.5,-z+.5) (-y+.5,-x+.5,-z+.5) (y+.5,-x+.5,z+.5) (-y+.5,x+.5,z+.5) (x+.5,z+.5,-y+.5) (-x+.5,z+.5,y+.5) (-x+.5,-z+.5,-y+.5) (x+.5,-z+.5,y+.5) (z+.5,y+.5,-x+.5) (z+.5,-y+.5,x+.5) (-z+.5,y+.5,x+.5) (-z+.5,-y+.5,-x+.5) (-x+.5,-y+.5,-z+.5) (x+.5,y+.5,-z+.5) (x+.5,-y+.5,z+.5) (-x+.5,y+.5,z+.5) (-z+.5,-x+.5,-y+.5) (-z+.5,x+.5,y+.5) (z+.5,x+.5,-y+.5) (z+.5,-x+.5,y+.5) (-y+.5,-z+.5,-x+.5) (y+.5,-z+.5,x+.5) (-y+.5,z+.5,x+.5) (y+.5,z+.5,-x+.5) (-y,-x,z) (y,x,z) (-y,x,-z) (y,-x,-z) (-x,-z,y) (x,-z,-y) (x,z,y) (-x,z,-y) (-z,-y,x) (-z,y,-x) (z,-y,-x) (z,y,x) \n 24 k ..m (x,x,z) (-x,-x,z) (-x,x,-z) (x,-x,-z) (z,x,x) (z,-x,-x) (-z,-x,x) (-z,x,-x) (x,z,x) (-x,z,-x) (x,-z,-x) (-x,-z,x) (x+.5,x+.5,-z+.5) (-x+.5,-x+.5,-z+.5) (x+.5,-x+.5,z+.5) (-x+.5,x+.5,z+.5) (x+.5,z+.5,-x+.5) (-x+.5,z+.5,x+.5) (-x+.5,-z+.5,-x+.5) (x+.5,-z+.5,x+.5) (z+.5,x+.5,-x+.5) (z+.5,-x+.5,x+.5) (-z+.5,x+.5,x+.5) (-z+.5,-x+.5,-x+.5) \n 24 j ..2 (.25,y,y+.5) (.75,-y,y+.5) (.75,y,-y+.5) (.25,-y,-y+.5) (y+.5,.25,y) (y+.5,.75,-y) (-y+.5,.75,y) (-y+.5,.25,-y) (y,y+.5,.25) (-y,y+.5,.75) (y,-y+.5,.75) (-y,-y+.5,.25) (.25,-y+.5,-y) (.75,y+.5,-y) (.75,-y+.5,y) (.25,y+.5,y) (-y,.25,-y+.5) (-y,.75,y+.5) (y,.75,-y+.5) (y,.25,y+.5) (-y+.5,-y,.25) (y+.5,-y,.75) (-y+.5,y,.75) (y+.5,y,.25) \n 24 i ..2 (.25,y,-y+.5) (.75,-y,-y+.5) (.75,y,y+.5) (.25,-y,y+.5) (-y+.5,.25,y) (-y+.5,.75,-y) (y+.5,.75,y) (y+.5,.25,-y) (y,-y+.5,.25) (-y,-y+.5,.75) (y,y+.5,.75) (-y,y+.5,.25) (.25,-y+.5,y) (.75,y+.5,y) (.75,-y+.5,-y) (.25,y+.5,-y) (y,.25,-y+.5) (y,.75,y+.5) (-y,.75,-y+.5) (-y,.25,y+.5) (-y+.5,y,.25) (y+.5,y,.75) (-y+.5,-y,.75) (y+.5,-y,.25) \n 24 h 2.. (x,0,.5) (-x,0,.5) (.5,x,0) (.5,-x,0) (0,.5,x) (0,.5,-x) (.5,x+.5,0) (.5,-x+.5,0) (x+.5,0,.5) (-x+.5,0,.5) (0,.5,-x+.5) (0,.5,x+.5) (-x+.5,.5,0) (x+.5,.5,0) (0,-x+.5,.5) (0,x+.5,.5) (.5,0,-x+.5) (.5,0,x+.5) (0,-x,.5) (0,x,.5) (-x,.5,0) (x,.5,0) (.5,0,x) (.5,0,-x) \n 12 g 2.mm (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) (.5,x+.5,.5) (.5,-x+.5,.5) (x+.5,.5,.5) (-x+.5,.5,.5) (.5,.5,-x+.5) (.5,.5,x+.5) \n 12 f 2.22 (.25,0,.5) (.75,0,.5) (.5,.25,0) (.5,.75,0) (0,.5,.25) (0,.5,.75) (.25,.5,0) (.75,.5,0) (0,.25,.5) (0,.75,.5) (.5,0,.25) (.5,0,.75) \n 8 e .3m (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x+.5,x+.5,-x+.5) (-x+.5,-x+.5,-x+.5) (x+.5,-x+.5,x+.5) (-x+.5,x+.5,x+.5) \n 6 d -42.m (0,.5,.5) (.5,0,.5) (.5,.5,0) (0,.5,0) (.5,0,0) (0,0,.5) \n 4 c .-3m (.75,.75,.75) (.25,.25,.75) (.25,.75,.25) (.75,.25,.25) \n 4 b .-3m (.25,.25,.25) (.75,.75,.25) (.75,.25,.75) (.25,.75,.75) \n 2 a -43m (0,0,0) (.5,.5,.5) \n ";
  }

  string sg225 = "Space Group 225 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 192 l 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y,x,-z) (-y,-x,-z) (y,-x,z) (-y,x,z)  (x,z,-y) (-x,z,y) (-x,-z,-y) (x,-z,y) (z,y,-x) (z,-y,x) (-z,y,x) (-z,-y,-x)  (-x,-y,-z) (x,y,-z) (x,-y,z) (-x,y,z) (-z,-x,-y) (-z,x,y) (z,x,-y) (z,-x,y) (-y,-z,-x) (y,-z,x) (-y,z,x) (y,z,-x) (-y,-x,z) (y,x,z) (-y,x,-z) (y,-x,-z) (-x,-z,y) (x,-z,-y) (x,z,y) (-x,z,-y) (-z,-y,x) (-z,y,-x) (z,-y,-x) (z,y,x) \n 96 k ..m (x,x,z) (-x,-x,z) (-x,x,-z) (x,-x,-z) (z,x,x) (z,-x,-x) (-z,-x,x) (-z,x,-x) (x,z,x) (-x,z,-x) (x,-z,-x) (-x,-z,x) (x,x,-z) (-x,-x,-z) (x,-x,z) (-x,x,z) (x,z,-x) (-x,z,x) (-x,-z,-x) (x,-z,x) (z,x,-x) (z,-x,x) (-z,x,x) (-z,-x,-x) \n 96 j m.. (0,y,z) (0,-y,z) (0,y,-z) (0,-y,-z) (z,0,y) (z,0,-y) (-z,0,y) (-z,0,-y) (y,z,0) (-y,z,0) (y,-z,0) (-y,-z,0) (y,0,-z) (-y,0,-z) (y,0,z) (-y,0,z) (0,z,-y) (0,z,y) (0,-z,-y) (0,-z,y) (z,y,0) (z,-y,0) (-z,y,0) (-z,-y,0) \n 48 i m.m2 (.5,y,y) (.5,-y,y) (.5,y,-y) (.5,-y,-y) (y,.5,y) (y,.5,-y) (-y,.5,y) (-y,.5,-y) (y,y,.5) (-y,y,.5) (y,-y,.5) (-y,-y,.5) \n 48 h m.m2 (0,y,y) (0,-y,y) (0,y,-y) (0,-y,-y) (y,0,y) (y,0,-y) (-y,0,y) (-y,0,-y) (y,y,0) (-y,y,0) (y,-y,0) (-y,-y,0) \n 48 g 2.mm (x,.25,.25) (-x,.75,.25) (.25,x,.25) (.25,-x,.75) (.25,.25,x) (.75,.25,-x) (.25,x,.75) (.75,-x,.75) (x,.25,.75) (-x,.25,.25) (.25,.25,-x) (.25,.75,x) \n 32 f .3m (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x,x,-x) (-x,-x,-x) (x,-x,x) (-x,x,x) \n 24 e 4m.m (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) \n 24 d m.mm (0,.25,.25) (0,.75,.25) (.25,0,.25) (.25,0,.75) (.25,.25,0) (.75,.25,0) \n 8 c -43m (.25,.25,.25) (.25,.25,.75) \n 4 b m-3m (.5,.5,.5) \n 4 a m-3m (0,0,0) \n ";

  string sg226 = "Space Group 226 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 192 j 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y+.5,x+.5,-z+.5) (-y+.5,-x+.5,-z+.5) (y+.5,-x+.5,z+.5) (-y+.5,x+.5,z+.5) (x+.5,z+.5,-y+.5) (-x+.5,z+.5,y+.5) (-x+.5,-z+.5,-y+.5) (x+.5,-z+.5,y+.5) (z+.5,y+.5,-x+.5) (z+.5,-y+.5,x+.5) (-z+.5,y+.5,x+.5) (-z+.5,-y+.5,-x+.5) (-x,-y,-z) (x,y,-z) (x,-y,z) (-x,y,z) (-z,-x,-y) (-z,x,y) (z,x,-y) (z,-x,y) (-y,-z,-x) (y,-z,x) (-y,z,x) (y,z,-x) (-y+.5,-x+.5,z+.5) (y+.5,x+.5,z+.5) (-y+.5,x+.5,-z+.5) (y+.5,-x+.5,-z+.5) (-x+.5,-z+.5,y+.5) (x+.5,-z+.5,-y+.5) (x+.5,z+.5,y+.5) (-x+.5,z+.5,-y+.5) (-z+.5,-y+.5,x+.5) (-z+.5,y+.5,-x+.5) (z+.5,-y+.5,-x+.5) (z+.5,y+.5,x+.5) \n 96 i m.. (0,y,z) (0,-y,z) (0,y,-z) (0,-y,-z) (z,0,y) (z,0,-y) (-z,0,y) (-z,0,-y) (y,z,0) (-y,z,0) (y,-z,0) (-y,-z,0) (y+.5,.5,-z+.5) (-y+.5,.5,-z+.5) (y+.5,.5,z+.5) (-y+.5,.5,z+.5) (.5,z+.5,-y+.5) (.5,z+.5,y+.5) (.5,-z+.5,-y+.5) (.5,-z+.5,y+.5) (z+.5,y+.5,.5) (z+.5,-y+.5,.5) (-z+.5,y+.5,.5) (-z+.5,-y+.5,.5) \n 96 h ..2 (.25,y,y) (.75,-y,y) (.75,y,-y) (.25,-y,-y) (y,.25,y) (y,.75,-y) (-y,.75,y) (-y,.25,-y) (y,y,.25) (-y,y,.75) (y,-y,.75) (-y,-y,.25) (.75,-y,-y) (.25,y,-y) (.25,-y,y) (.75,y,y) (-y,.75,-y) (-y,.25,y) (y,.25,-y) (y,.75,y) (-y,-y,.75) (y,-y,.25) (-y,y,.25) (y,y,.75) \n 64 g .3. (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x+.5,x+.5,-x+.5) (-x+.5,-x+.5,-x+.5) (x+.5,-x+.5,x+.5) (-x+.5,x+.5,x+.5) (-x,-x,-x) (x,x,-x) (x,-x,x) (-x,x,x) (-x+.5,-x+.5,x+.5) (x+.5,x+.5,x+.5) (-x+.5,x+.5,-x+5) (x+.5,-x+.5,-x+.5) \n 48 f 4.. (x,.25,.25) (-x,.75,.25) (.25,x,.25) (.25,-x,.75) (.25,.25,x) (.75,.25,-x) (-x,.75,.75) (x,.25,.75) (.75,-x,.75) (.75,x,.25) (.75,.75,-x) (.25,.75,x) \n 48 e mm2.. (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) (.5,x+.5,.5) (.5,-x+.5,.5) (x+.5,.5,.5) (-x+.5,.5,.5) (.5,.5,-x+.5) (.5,.5,x+.5) \n 24 d 4/m.. (0,.25,.25) (0,.75,.25) (.25,0,.25) (.25,0,.75) (.25,.25,0) (.75,.25,0) \n 24 c -4m.2 (.25,0,0) (.75,0,0) (0,.25,0) (0,.75,0) (0,0,.25) (0,0,.75) \n 8 b m-3. (0,0,0) (.5,.5,.5) \n 8 a 432 (.25,.25,.25) (.75,.75,.75) \n ";

  //DX20180723 - added origin choice 2
  string sg227 = "";
  if(axis_cell==setting_1){
    sg227 = "Space Group 227 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 192 i 1 (x,y,z) (-x,-y+.5,z+.5) (-x+.5,y+.5,-z) (x+.5,-y,-z+.5) (z,x,y) (z+.5,-x,-y+.5) (-z,-x+.5,y+.5) (-z+.5,x+.5,-y) (y,z,x) (-y+.5,z+.5,-x) (y+.5,-z,-x+.5) (-y,-z+.5,x+.5) (y+.75,x+.25,-z+.75) (-y+.25,-x+.25,-z+.25) (y+.25,-x+.75,z+.75) (-y+.75,x+.75,z+.25) (x+.75,z+.25,-y+.75) (-x+.75,z+.75,y+.25) (-x+.25,-z+.25,-y+.25) (x+.25,-z+.75,y+.75) (z+.75,y+.25,-x+.75) (z+.25,-y+.75,x+.75) (-z+.75,y+.75,x+.25) (-z+.25,-y+.25,-x+.25) (-x+.25,-y+.25,-z+.25) (x+.25,y+.75,-z+.75) (x+.75,-y+.75,z+.25) (-x+.75,y+.25,z+.75) (-z+.25,-x+.25,-y+.25) (-z+.75,x+.25,y+.75) (z+.25,x+.75,-y+.75) (z+.75,-x+.75,y+.25) (-y+.25,-z+.25,-x+.25) (y+.75,-z+.75,x+.25) (-y+.75,z+.25,x+.75) (y+.25,z+.75,-x+.75) (-y+.5,-x,z+.5) (y,x,z) (-y,x+.5,-z+.5) (y+.5,-x+.5,-z) (-x+.5,-z,y+.5) (x+.5,-z+.5,-y) (x,z,y) (-x,z+.5,-y+.5) (-z+.5,-y,x+.5) (-z,y+.5,-x+.5) (z+.5,-y+.5,-x) (z,y,x) \n 96 h ..2 (.125,y,-y+.25) (.875,-y+.5,-y+.75) (.375,y+.5,y+.75) (.625,-y,y+.25) (-y+.25,.125,y) (-y+.75,.875,-y+.5) (y+.75,.375,y+.5) (y+.25,.625,-y) (y,-y+.25,.125) (-y+.5,-y+.75,.875) (y+.5,y+.75,.375) (-y,y+.25,.625) (.125,-y+.25,y) (.375,y+.75,y+.5) (.875,-y+.75,-y+.5) (.625,y+.25,-y) (y,.125,-y+.25) (y+.5,.375,y+.75) (-y+.5,.875,-y+.75) (-y,.625,y+.25) (-y+.25,y,.125) (y+.75,y+.5,.375) (-y+.75,-y+.5,.875) (y+.25,-y,.625) \n 96 g ..m (x,x,z) (-x,-x+.5,z+.5) (-x+.5,x+.5,-z) (x+.5,-x,-z+.5) (z,x,x) (z+.5,-x,-x+.5) (-z,-x+.5,x+.5) (-z+.5,x+.5,-x) (x,z,x) (-x+.5,z+.5,-x) (x+.5,-z,-x+.5) (-x,-z+.5,x+.5) (x+.75,x+.25,-z+.75) (-x+.25,-x+.25,-z+.25) (x+.25,-x+.75,z+.75) (-x+.75,x+.75,z+.25) (x+.75,z+.25,-x+.75) (-x+.75,z+.75,x+.25) (-x+.25,-z+.25,-x+.25) (x+.25,-z+.75,x+.75) (z+.75,x+.25,-x+.75) (z+.25,-x+.75,x+.75) (-z+.75,x+.75,x+.25) (-z+.25,-x+.25,-x+.25) \n 48 f 2.mm (x,0,0) (-x,.5,.5) (0,x,0) (.5,-x,.5) (0,0,x) (.5,.5,-x) (.75,x+.25,.75) (.25,-x+.25,.25) (x+.75,.25,.75) (-x+.75,.75,.25) (.75,.25,-x+.75) (.25,.75,x+.75) \n 32 e .3m (x,x,x) (-x,-x+.5,x+.5) (-x+.5,x+.5,-x) (x+.5,-x,-x+.5) (x+.75,x+.25,-x+.75) (-x+.25,-x+.25,-x+.25) (x+.25,-x+.75,x+.75) (-x+.75,x+.75,x+.25) \n 16 d .-3m (.625,.625,.625) (.375,.875,.125) (.875,.125,.375) (.125,.375,.875) \n 16 c .-3m (.125,.125,.125) (.875,.375,.625) (.375,.625,.875) (.625,.875,.375) \n 8 b -43m (.5,.5,.5) (.25,.75,.25) \n 8 a -43m (0,0,0) (.75,.25,.75) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg227 = "Space Group 227 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 192 i 1 (x,y,z) (-x+.75,-y+.25,z+.5) (-x+.25,y+.5,-z+.75) (x+.5,-y+.75,-z+.25) (z,x,y) (z+.5,-x+.75,-y+.25) (-z+.75,-x+.25,y+.5) (-z+.25,x+.5,-y+.75) (y,z,x) (-y+.25,z+.5,-x+.75) (y+.5,-z+.75,-x+.25) (-y+.75,-z+.25,x+.5) (y+.75,x+.25,-z+.5) (-y,-x,-z) (y+.25,-x+.5,z+.75) (-y+.5,x+.75,z+.25) (x+.75,z+.25,-y+.5) (-x+.5,z+.75,y+.25) (-x,-z,-y) (x+.25,-z+.5,y+.75) (z+.75,y+.25,-x+.5) (z+.25,-y+.5,x+.75) (-z+.5,y+.75,x+.25) (-z,-y,-x) (-x,-y,-z) (x+.25,y+.75,-z+.5) (x+.75,-y+.5,z+.25) (-x+.5,y+.25,z+.75) (-z,-x,-y) (-z+.5,x+.25,y+.75) (z+.25,x+.75,-y+.5) (z+.75,-x+.5,y+.25) (-y,-z,-x) (y+.75,-z+.5,x+.25) (-y+.5,z+.25,x+.75) (y+.25,z+.75,-x+.5) (-y+.25,-x+.75,z+.5) (y,x,z) (-y+.75,x+.5,-z+.25) (y+.5,-x+.25,-z+.75) (-x+.25,-z+.75,y+.5) (x+.5,-z+.25,-y+.75) (x,z,y) (-x+.75,z+.5,-y+.25) (-z+.25,-y+.75,x+.5) (-z+.75,y+.5,-x+.25) (z+.5,-y+.25,-x+.75) (z,y,x) \n 96 h ..2 (0,y,-y) (.75,-y+.25,-y+.5) (.25,y+.5,y+.75) (.5,-y+.75,y+.25) (-y,0,y) (-y+.5,.75,-y+.25) (y+.75,.25,y+.5) (y+.25,.5,-y+.75) (y,-y,0) (-y+.25,-y+.5,.75) (y+.5,y+.75,.25) (-y+.75,y+.25,.5) (0,-y,y) (.25,y+.75,y+.5) (.75,-y+.5,-y+.25) (.5,y+.25,-y+.75) (y,0,-y) (y+.5,.25,y+.75) (-y+.25,.75,-y+.5) (-y+.75,.5,y+.25) (-y,y,0) (y+.75,y+.5,.25) (-y+.5,-y+.25,.75) (y+.25,-y+.75,.5) \n 96 g ..m (x,x,z) (-x+.75,-x+.25,z+.5) (-x+.25,x+.5,-z+.75) (x+.5,-x+.75,-z+.25) (z,x,x) (z+.5,-x+.75,-x+.25) (-z+.75,-x+.25,x+.5) (-z+.25,x+.5,-x+.75) (x,z,x) (-x+.25,z+.5,-x+.75) (x+.5,-z+.75,-x+.25) (-x+.75,-z+.25,x+.5) (x+.75,x+.25,-z+.5) (-x,-x,-z) (x+.25,-x+.5,z+.75) (-x+.5,x+.75,z+.25) (x+.75,z+.25,-x+.5) (-x+.5,z+.75,x+.25) (-x,-z,-x) (x+.25,-z+.5,x+.75) (z+.75,x+.25,-x+.5) (z+.25,-x+.5,x+.75) (-z+.5,x+.75,x+.25) (-z,-x,-x) \n 48 f 2.mm (x,.125,.125) (-x+.75,.125,.625) (.125,x,.125) (.625,-x+.75,.125) (.125,.125,x) (.125,.625,-x+.75) (.875,x+.25,.375) (.875,-x,.875) (x+.75,.375,.375) (-x+.5,.875,.375) (.875,.375,-x+.5) (.375,.375,x+.75) \n 32 e .3m (x,x,x) (-x+.75,-x+.25,x+.5) (-x+.25,x+.5,-x+.75) (x+.5,-x+.75,-x+.25) (x+.75,x+.25,-x+.5) (-x,-x,-x) (x+.25,-x+.5,x+.75) (-x+.5,x+.75,x+.25) \n 16 d .-3m (.5,.5,.5) (.25,.75,0) (.75,0,.25) (0,.25,.75) \n 16 c .-3m (0,0,0) (.75,.25,.5) (.25,.5,.75) (.5,.75,.25) \n 8 b -43m (.375,.375,.375) (.125,.625,.125) \n 8 a -43m (.125,.125,.125) (.875,.375,.375) \n ";
  }
  else {
    sg227 = "Space Group 227 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 192 i 1 (x,y,z) (-x,-y+.5,z+.5) (-x+.5,y+.5,-z) (x+.5,-y,-z+.5) (z,x,y) (z+.5,-x,-y+.5) (-z,-x+.5,y+.5) (-z+.5,x+.5,-y) (y,z,x) (-y+.5,z+.5,-x) (y+.5,-z,-x+.5) (-y,-z+.5,x+.5) (y+.75,x+.25,-z+.75) (-y+.25,-x+.25,-z+.25) (y+.25,-x+.75,z+.75) (-y+.75,x+.75,z+.25) (x+.75,z+.25,-y+.75) (-x+.75,z+.75,y+.25) (-x+.25,-z+.25,-y+.25) (x+.25,-z+.75,y+.75) (z+.75,y+.25,-x+.75) (z+.25,-y+.75,x+.75) (-z+.75,y+.75,x+.25) (-z+.25,-y+.25,-x+.25) (-x+.25,-y+.25,-z+.25) (x+.25,y+.75,-z+.75) (x+.75,-y+.75,z+.25) (-x+.75,y+.25,z+.75) (-z+.25,-x+.25,-y+.25) (-z+.75,x+.25,y+.75) (z+.25,x+.75,-y+.75) (z+.75,-x+.75,y+.25) (-y+.25,-z+.25,-x+.25) (y+.75,-z+.75,x+.25) (-y+.75,z+.25,x+.75) (y+.25,z+.75,-x+.75) (-y+.5,-x,z+.5) (y,x,z) (-y,x+.5,-z+.5) (y+.5,-x+.5,-z) (-x+.5,-z,y+.5) (x+.5,-z+.5,-y) (x,z,y) (-x,z+.5,-y+.5) (-z+.5,-y,x+.5) (-z,y+.5,-x+.5) (z+.5,-y+.5,-x) (z,y,x) \n 96 h ..2 (.125,y,-y+.25) (.875,-y+.5,-y+.75) (.375,y+.5,y+.75) (.625,-y,y+.25) (-y+.25,.125,y) (-y+.75,.875,-y+.5) (y+.75,.375,y+.5) (y+.25,.625,-y) (y,-y+.25,.125) (-y+.5,-y+.75,.875) (y+.5,y+.75,.375) (-y,y+.25,.625) (.125,-y+.25,y) (.375,y+.75,y+.5) (.875,-y+.75,-y+.5) (.625,y+.25,-y) (y,.125,-y+.25) (y+.5,.375,y+.75) (-y+.5,.875,-y+.75) (-y,.625,y+.25) (-y+.25,y,.125) (y+.75,y+.5,.375) (-y+.75,-y+.5,.875) (y+.25,-y,.625) \n 96 g ..m (x,x,z) (-x,-x+.5,z+.5) (-x+.5,x+.5,-z) (x+.5,-x,-z+.5) (z,x,x) (z+.5,-x,-x+.5) (-z,-x+.5,x+.5) (-z+.5,x+.5,-x) (x,z,x) (-x+.5,z+.5,-x) (x+.5,-z,-x+.5) (-x,-z+.5,x+.5) (x+.75,x+.25,-z+.75) (-x+.25,-x+.25,-z+.25) (x+.25,-x+.75,z+.75) (-x+.75,x+.75,z+.25) (x+.75,z+.25,-x+.75) (-x+.75,z+.75,x+.25) (-x+.25,-z+.25,-x+.25) (x+.25,-z+.75,x+.75) (z+.75,x+.25,-x+.75) (z+.25,-x+.75,x+.75) (-z+.75,x+.75,x+.25) (-z+.25,-x+.25,-x+.25) \n 48 f 2.mm (x,0,0) (-x,.5,.5) (0,x,0) (.5,-x,.5) (0,0,x) (.5,.5,-x) (.75,x+.25,.75) (.25,-x+.25,.25) (x+.75,.25,.75) (-x+.75,.75,.25) (.75,.25,-x+.75) (.25,.75,x+.75) \n 32 e .3m (x,x,x) (-x,-x+.5,x+.5) (-x+.5,x+.5,-x) (x+.5,-x,-x+.5) (x+.75,x+.25,-x+.75) (-x+.25,-x+.25,-x+.25) (x+.25,-x+.75,x+.75) (-x+.75,x+.75,x+.25) \n 16 d .-3m (.625,.625,.625) (.375,.875,.125) (.875,.125,.375) (.125,.375,.875) \n 16 c .-3m (.125,.125,.125) (.875,.375,.625) (.375,.625,.875) (.625,.875,.375) \n 8 b -43m (.5,.5,.5) (.25,.75,.25) \n 8 a -43m (0,0,0) (.75,.25,.75) \n ";
  }

  //DX20180723 - added origin choice 2
  string sg228 = "";
  if(axis_cell==setting_1){
    sg228 = "Space Group 228 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 192 h 1 (x,y,z) (-x,-y+.5,z+.5) (-x+.5,y+.5,-z) (x+.5,-y,-z+.5) (z,x,y) (z+.5,-x,-y+.5) (-z,-x+.5,y+.5) (-z+.5,x+.5,-y) (y,z,x) (-y+.5,z+.5,-x) (y+.5,-z,-x+.5) (-y,-z+.5,x+.5) (y+.75,x+.25,-z+.75) (-y+.25,-x+.25,-z+.25) (y+.25,-x+.75,z+.75) (-y+.75,x+.75,z+.25) (x+.75,z+.25,-y+.75) (-x+.75,z+.75,y+.25) (-x+.25,-z+.25,-y+.25) (x+.25,-z+.75,y+.75) (z+.75,y+.25,-x+.75) (z+.25,-y+.75,x+.75) (-z+.75,y+.75,x+.25) (-z+.25,-y+.25,-x+.25) (-x+.75,-y+.75,-z+.75) (x+.75,y+.25,-z+.25) (x+.25,-y+.25,z+.75) (-x+.25,y+.75,z+.25) (-z+.75,-x+.75,-y+.75) (-z+.25,x+.75,y+.25) (z+.75,x+.25,-y+.25) (z+.25,-x+.25,y+.75) (-y+.75,-z+.75,-x+.75) (y+.25,-z+.25,x+.75) (-y+.25,z+.75,x+.25) (y+.75,z+.25,-x+.25) (-y,-x+.5,z) (y+.5,x+.5,z+.5) (-y+.5,x,-z) (y,-x,-z+.5) (-x,-z+.5,y) (x,-z,-y+.5) (x+.5,z+.5,y+.5) (-x+.5,z,-y) (-z,-y+.5,x) (-z+.5,y,-x) (z,-y,-x+.5) (z+.5,y+.5,x+.5) \n 96 g ..2 (.125,y,-y+.25) (.875,-y+.5,-y+.75) (.375,y+.5,y+.75) (.625,-y,y+.25) (-y+.25,.125,y) (-y+.75,.875,-y+.5) (y+.75,.375,y+.5) (y+.25,.625,-y) (y,-y+.25,.125) (-y+.5,-y+.75,.875) (y+.5,y+.75,.375) (-y,y+.25,.625) (.625,-y+.75,y+.5) (.875,y+.25,y) (.375,-y+.25,-y) (.125,y+.75,-y+.5) (y+.5,.625,-y+.75) (y,.875,y+.25) (-y,.375,-y+.25) (-y+.5,.125,y+.75) (-y+.75,y+.5,.625) (y+.25,y,.875) (-y+.25,-y,.375) (y+.75,-y+.5,.125) \n 96 f 2.. (x,0,0) (-x,.5,.5) (0,x,0) (.5,-x,.5) (0,0,x) (.5,.5,-x) (.75,x+.25,.75) (.25,-x+.25,.25) (x+.75,.25,.75) (-x+.75,.75,.25) (.75,.25,-x+.75) (.25,.75,x+.75) (-x+.75,.75,.75) (x+.75,.25,.25) (.75,-x+.75,.75) (.25,x+.75,.25) (.75,.75,-x+.75) (.25,.25,x+.75) (0,-x+.5,0) (.5,x+.5,.5) (-x,.5,0) (x,0,.5) (0,.5,x) (.5,0,-x) \n 64 e .3. (x,x,x) (-x,-x+.5,x+.5) (-x+.5,x+.5,-x) (x+.5,-x,-x+.5) (x+.75,x+.25,-x+.75) (-x+.25,-x+.25,-x+.25) (x+.25,-x+.75,x+.75) (-x+.75,x+.75,x+.25) (-x+.75,-x+.75,-x+.75) (x+.75,x+.25,-x+.25) (x+.25,-x+.25,x+.75) (-x+.25,x+.75,x+.25) (-x,-x+.5,x) (x+.5,x+.5,x+.5) (-x+.5,x,-x) (x,-x,-x+.5) \n 48 d -4.. (.25,0,0) (.75,.5,.5) (0,.25,0) (.5,.75,.5) (0,0,.25) (.5,.5,.75) (.75,.5,.75) (.25,0,.25) (0,.25,.75) (.5,.75,.5) (.75,.25,.5) (.25,.75,0) \n 32 c .-3. (.375,.375,.375) (.625,.125,.875) (.125,.875,.625) (.875,.625,.125) (.125,.625,.375) (.875,.875,.875) (.625,.375,.125) (.375,.125,.625) \n 32 b .32 (.125,.125,.125) (.875,.375,.625) (.375,.625,.875) (.625,.875,.375) (.625,.625,.625) (.875,.375,.125) (.375,.125,.875) (.125,.875,.375) \n 16 a 23. (0,0,0) (.75,.25,.75) (.75,.75,.75) (0,.5,0) \n ";
  }
  else if(axis_cell==setting_2 || axis_cell==setting_anrl){
    sg228 = "Space Group 228 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 192 h 1 (x,y,z) (-x+.25,-y+.75,z+.5) (-x+.75,y+.5,-z+.25) (x+.5,-y+.25,-z+.75) (z,x,y) (z+.5,-x+.25,-y+.75) (-z+.25,-x+.75,y+.5) (-z+.75,x+.5,-y+.25) (y,z,x) (-y+.75,z+.5,-x+.25) (y+.5,-z+.25,-x+.75) (-y+.25,-z+.75,x+.5) (y+.75,x+.25,-z) (-y+.5,-x+.5,-z+.5) (y+.25,-x,z+.75) (-y,x+.75,z+.25) (x+.75,z+.25,-y) (-x,z+.75,y+.25) (-x+.5,-z+.5,-y+.5) (x+.25,-z,y+.75) (z+.75,y+.25,-x) (z+.25,-y,x+.75) (-z,y+.75,x+.25) (-z+.5,-y+.5,-x+.5) (-x,-y,-z) (x+.75,y+.25,-z+.5) (x+.25,-y+.5,z+.75) (-x+.5,y+.75,z+.25) (-z,-x,-y) (-z+.5,x+.75,y+.25) (z+.75,x+.25,-y+.5) (z+.25,-x+.5,y+.75) (-y,-z,-x) (y+.25,-z+.5,x+.75) (-y+.5,z+.75,x+.25) (y+.75,z+.25,-x+.5) (-y+.25,-x+.75,z) (y+.5,x+.5,z+.5) (-y+.75,x,-z+.25) (y,-x+.25,-z+.75) (-x+.25,-z+.75,y) (x,-z+.25,-y+.75) (x+.5,z+.5,y+.5) (-x+.75,z,-y+.25) (-z+.25,-y+.75,x) (-z+.75,y,-x+.25) (z,-y+.25,-x+.75) (z+.5,y+.5,x+.5) \n 96 g ..2 (.25,y,-y) (0,-y+.75,-y+.5) (.5,y+.5,y+.25) (.75,-y+.25,y+.75) (-y,.25,y) (-y+.5,0,-y+.75) (y+.25,.5,y+.5) (y+.75,.75,-y+.25) (y,-y,.25) (-y+.75,-y+.5,0) (y+.5,y+.25,.5) (-y+.25,y+.75,.75) (.75,-y,y) (0,y+.25,y+.5) (.5,-y+.5,-y+.75) (.25,y+.75,-y+.25) (y,.75,-y) (y+.5,0,y+.25) (-y+.75,.5,-y+.5) (-y+.25,.25,y+.75) (-y,y,.75) (y+.25,y+.5,0) (-y+.5,-y+.75,.5) (y+.75,-y+.25,.25) \n 96 f 2.. (x,.125,.125) (-x+.25,.625,.625) (.125,x,.125) (.625,-x+.25,.625) (.125,.125,x) (.625,.625,-x+.25) (.875,x+.25,.875) (.375,-x+.5,.375) (x+.75,.375,.875) (-x,.875,.375) (.875,.375,-x) (.375,.875,x+.75) (-x,.875,.875) (x+.75,.375,.375) (.875,-x,.875) (.375,x+.75,.375) (.875,.875,-x) (.375,.375,x+.75) (.125,-x+.75,.125) (.625,x+.5,.625) (-x+.25,.625,.125) (x,.125,.625) (.125,.625,x) (.625,.125,-x+.25) \n 64 e .3. (x,x,x) (-x+.25,-x+.75,x+.5) (-x+.75,x+.5,-x+.25) (x+.5,-x+.25,-x+.75) (x+.75,x+.25,-x) (-x+.5,-x+.5,-x+.5) (x+.25,-x,x+.75) (-x,x+.75,x+.25) (-x,-x,-x) (x+.75,x+.25,-x+.5) (x+.25,-x+.5,x+.75) (-x+.5,x+.75,x+.25) (-x+.25,-x+.75,x) (x+.5,x+.5,x+.5) (-x+.75,x,-x+.25) (x,-x+.25,-x+.75) \n 48 d -4.. (.875,.125,.125) (.375,.625,.625) (.125,.875,.125) (.625,.375,.625) (.125,.125,.875) (.625,.625,.375) (.875,.125,.875) (.625,.625,.375) (.625,.375,.875) (.125,.875,.375) (.875,.375,.125) (.375,.875,.625) \n 32 c .-3. (0,0,0) (.25,.75,.5) (.75,.5,.25) (.5,.25,.75) (.75,.25,0) (.5,.5,.5) (.25,0,.75) (0,.75,.25) \n 32 b .32 (.25,.25,.25) (0,.5,.75) (.5,.75,0) (.75,0,.5) (.75,.75,.75) (0,.5,.25) (.5,.25,0) (.25,0,.5) \n 16 a 23. (.125,.125,.125) (.875,.375,.875) (.875,.875,.875) (.125,.625,.125) \n ";
  }
  else {
    sg228 = "Space Group 228 \n (0,0,0) (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 192 h 1 (x,y,z) (-x,-y+.5,z+.5) (-x+.5,y+.5,-z) (x+.5,-y,-z+.5) (z,x,y) (z+.5,-x,-y+.5) (-z,-x+.5,y+.5) (-z+.5,x+.5,-y) (y,z,x) (-y+.5,z+.5,-x) (y+.5,-z,-x+.5) (-y,-z+.5,x+.5) (y+.75,x+.25,-z+.75) (-y+.25,-x+.25,-z+.25) (y+.25,-x+.75,z+.75) (-y+.75,x+.75,z+.25) (x+.75,z+.25,-y+.75) (-x+.75,z+.75,y+.25) (-x+.25,-z+.25,-y+.25) (x+.25,-z+.75,y+.75) (z+.75,y+.25,-x+.75) (z+.25,-y+.75,x+.75) (-z+.75,y+.75,x+.25) (-z+.25,-y+.25,-x+.25) (-x+.75,-y+.75,-z+.75) (x+.75,y+.25,-z+.25) (x+.25,-y+.25,z+.75) (-x+.25,y+.75,z+.25) (-z+.75,-x+.75,-y+.75) (-z+.25,x+.75,y+.25) (z+.75,x+.25,-y+.25) (z+.25,-x+.25,y+.75) (-y+.75,-z+.75,-x+.75) (y+.25,-z+.25,x+.75) (-y+.25,z+.75,x+.25) (y+.75,z+.25,-x+.25) (-y,-x+.5,z) (y+.5,x+.5,z+.5) (-y+.5,x,-z) (y,-x,-z+.5) (-x,-z+.5,y) (x,-z,-y+.5) (x+.5,z+.5,y+.5) (-x+.5,z,-y) (-z,-y+.5,x) (-z+.5,y,-x) (z,-y,-x+.5) (z+.5,y+.5,x+.5) \n 96 g ..2 (.125,y,-y+.25) (.875,-y+.5,-y+.75) (.375,y+.5,y+.75) (.625,-y,y+.25) (-y+.25,.125,y) (-y+.75,.875,-y+.5) (y+.75,.375,y+.5) (y+.25,.625,-y) (y,-y+.25,.125) (-y+.5,-y+.75,.875) (y+.5,y+.75,.375) (-y,y+.25,.625) (.625,-y+.75,y+.5) (.875,y+.25,y) (.375,-y+.25,-y) (.125,y+.75,-y+.5) (y+.5,.625,-y+.75) (y,.875,y+.25) (-y,.375,-y+.25) (-y+.5,.125,y+.75) (-y+.75,y+.5,.625) (y+.25,y,.875) (-y+.25,-y,.375) (y+.75,-y+.5,.125) \n 96 f 2.. (x,0,0) (-x,.5,.5) (0,x,0) (.5,-x,.5) (0,0,x) (.5,.5,-x) (.75,x+.25,.75) (.25,-x+.25,.25) (x+.75,.25,.75) (-x+.75,.75,.25) (.75,.25,-x+.75) (.25,.75,x+.75) (-x+.75,.75,.75) (x+.75,.25,.25) (.75,-x+.75,.75) (.25,x+.75,.25) (.75,.75,-x+.75) (.25,.25,x+.75) (0,-x+.5,0) (.5,x+.5,.5) (-x,.5,0) (x,0,.5) (0,.5,x) (.5,0,-x) \n 64 e .3. (x,x,x) (-x,-x+.5,x+.5) (-x+.5,x+.5,-x) (x+.5,-x,-x+.5) (x+.75,x+.25,-x+.75) (-x+.25,-x+.25,-x+.25) (x+.25,-x+.75,x+.75) (-x+.75,x+.75,x+.25) (-x+.75,-x+.75,-x+.75) (x+.75,x+.25,-x+.25) (x+.25,-x+.25,x+.75) (-x+.25,x+.75,x+.25) (-x,-x+.5,x) (x+.5,x+.5,x+.5) (-x+.5,x,-x) (x,-x,-x+.5) \n 48 d -4.. (.25,0,0) (.75,.5,.5) (0,.25,0) (.5,.75,.5) (0,0,.25) (.5,.5,.75) (.75,.5,.75) (.25,0,.25) (0,.25,.75) (.5,.75,.5) (.75,.25,.5) (.25,.75,0) \n 32 c .-3. (.375,.375,.375) (.625,.125,.875) (.125,.875,.625) (.875,.625,.125) (.125,.625,.375) (.875,.875,.875) (.625,.375,.125) (.375,.125,.625) \n 32 b .32 (.125,.125,.125) (.875,.375,.625) (.375,.625,.875) (.625,.875,.375) (.625,.625,.625) (.875,.375,.125) (.375,.125,.875) (.125,.875,.375) \n 16 a 23. (0,0,0) (.75,.25,.75) (.75,.75,.75) (0,.5,0) \n ";
  }

  string sg229 = "Space Group 229 \n (0,0,0) (.5,.5,.5) \n 96 l 1 (x,y,z) (-x,-y,z) (-x,y,-z) (x,-y,-z) (z,x,y) (z,-x,-y) (-z,-x,y) (-z,x,-y) (y,z,x) (-y,z,-x) (y,-z,-x) (-y,-z,x) (y,x,-z) (-y,-x,-z) (y,-x,z) (-y,x,z) (x,z,-y) (-x,z,y) (-x,-z,-y) (x,-z,y) (z,y,-x) (z,-y,x) (-z,y,x) (-z,-y,-x) (-x,-y,-z) (x,y,-z) (x,-y,z) (-x,y,z) (-z,-x,-y) (-z,x,y) (z,x,-y) (z,-x,y) (-y,-z,-x) (y,-z,x) (-y,z,x) (y,z,-x) (-y,-x,z) (y,x,z) (-y,x,-z) (y,-x,-z) (-x,-z,y) (x,-z,-y) (x,z,y) (-x,z,-y) (-z,-y,x) (-z,y,-x) (z,-y,-x) (z,y,x) \n 48 k ..m (x,x,z) (-x,-x,z) (-x,x,-z) (x,-x,-z) (z,x,x) (z,-x,-x) (-z,-x,x) (-z,x,-x) (x,z,x) (-x,z,-x) (x,-z,-x) (-x,-z,x) (x,x,-z) (-x,-x,-z) (x,-x,z) (-x,x,z) (x,z,-x) (-x,z,x) (-x,-z,-x) (x,-z,x) (z,x,-x) (z,-x,x) (-z,x,x) (-z,-x,-x) \n 48 j m.. (0,y,z) (0,-y,z) (0,y,-z) (0,-y,-z) (z,0,y) (z,0,-y) (-z,0,y) (-z,0,-y) (y,z,0) (-y,z,0) (y,-z,0) (-y,-z,0) (y,0,-z) (-y,0,-z) (y,0,z) (-y,0,z) (0,z,-y) (0,z,y) (0,-z,-y) (0,-z,y) (z,y,0) (z,-y,0) (-z,y,0) (-z,-y,0) \n 48 i ..2 (.25,y,-y+.25) (.75,-y,-y+.5) (.75,y,y+.5) (.25,-y,y+.5) (-y+.5,.25,y) (-y+.5,.75,-y) (y+.5,.75,y) (y+.5,.25,-y) (y,-y+.5,.25) (-y,-y+.5,.75) (y,y+.5,.75) (-y,y+.5,.25) (.75,-y,y+.5) (.25,y,y+.5) (.25,-y,-y+.5) (.75,y,-y+.5) (y+.5,.75,-y) (y+.5,.25,y) (-y+.5,.25,-y) (-y+.5,.75,y) (-y,y+.5,.75) (y,y+.5,.25) (-y,-y+.5,.25) (y,-y+.5,.75) \n 24 h m.m2 (0,y,y) (0,-y,y) (0,y,-y) (0,-y,-y) (y,0,y) (y,0,-y) (-y,0,y) (-y,0,-y) (y,y,0) (-y,y,0) (y,-y,0) (-y,-y,0) \n 24 g mm2.. (x,0,.5) (-x,0,.5) (.5,x,0) (.5,-x,0) (0,.5,x) (0,.5,-x) (0,x,.5) (0,-x,.5) (x,.5,0) (-x,.5,0) (.5,0,-x) (.5,0,x) \n 16 f .3m (x,x,x) (-x,-x,x) (-x,x,-x) (x,-x,-x) (x,x,-x) (-x,-x,-x) (x,-x,x) (-x,x,x) \n 12 e 4m.m (x,0,0) (-x,0,0) (0,x,0) (0,-x,0) (0,0,x) (0,0,-x) \n 12 d -4m.2 (.25,0,.5) (.75,0,.5) (.5,.25,0) (.5,.75,0) (0,.5,.25) (0,.5,.75) \n 8 c .-3m (.25,.25,.25) (.75,.75,.25) (.75,.25,.75) (.25,.75,.75) \n 6 b 4/mm.m (0,.5,.5) (.5,0,.5) (.5,.5,0) \n 2 a m-3m (0,0,0) \n ";

  string sg230 = "Space Group 230 \n (0,0,0) (.5,.5,.5) \n 96 h 1 (x,y,z) (-x+.5,-y,z+.5) (-x,y+.5,-z+.5) (x+.5,-y+.5,-z) (z,x,y) (z+.5,-x+.5,-y) (-z+.5,-x,y+.5) (-z,x+.5,-y+.5) (y,z,x) (-y,z+.5,-x+.5) (y+.5,-z+.5,-x) (-y+.5,-z,x+.5) (y+.75,x+.25,-z+.25) (-y+.75,-x+.75,-z+.75) (y+.25,-x+.25,z+.75) (-y+.25,x+.75,z+.25) (x+.75,z+.25,-y+.25) (-x+.25,z+.75,y+.25) (-x+.75,-z+.75,-y+.75) (x+.25,-z+.25,y+.75) (z+.75,y+.25,-x+.25) (z+.25,-y+.25,x+.75) (-z+.25,y+.75,x+.25) (-z+.75,-y+.75,-x+.75) (-x,-y,-z) (x+.5,y,-z+.5) (x,-y+.5,z+.5) (-x+.5,y+.5,z) (-z,-x,-y) (-z+.5,x+.5,y) (z+.5,x,-y+.5) (z,-x+.5,y+.5) (-y,-z,-x) (y,-z+.5,x+.5) (-y+.5,z+.5,x) (y+.5,z,-x+.5) (-y+.25,-x+.75,z+.75) (y+.25,x+.25,z+.25) (-y+.75,x+.75,-z+.25) (y+.75,-x+.25,-z+.75) (-x+.25,-z+.75,y+.75) (x+.75,-z+.25,-y+.75) (x+.25,z+.25,y+.25) (-x+.75,z+.75,-y+.25) (-z+.25,-y+.75,x+.75) (-z+.75,y+.75,-x+.25) (z+.75,-y+.25,-x+.75) (z+.25,y+.25,x+.25) \n 48 g ..2 (.125,y,-y+.25) (.375,-y,-y+.75) (.875,y+.5,y+.25) (.625,-y+.5,y+.75) (-y+.25,.125,y) (-y+.75,.375,-y) (y+.25,.875,y+.5) (y+.75,.625,-y+.5) (y,-y+.25,.125) (-y,-y+.75,.375) (y+.5,y+.25,.875) (-y+.5,y+.75,.625) (.875,-y,y+.75) (.625,y,y+.25) (.125,-y+.5,-y+.75) (.375,y+.5,-y+.25) (y+.75,.875,-y) (y+.25,.625,y) (-y+.75,.125,-y+.5) (-y+.25,.375,y+.5) (-y,y+.75,.875) (y,y+.25,.625) (-y+.5,-y+.75,.125) (y+.5,-y+.25,.375) \n 48 f 2.. (x,0,.25) (-x+.5,0,.75) (.25,x,0) (.75,-x+.5,0) (0,.25,x) (0,.75,-x+.5) (.75,x+.25,0) (.75,-x+.75,.5) (x+.75,.5,.25) (-x+.25,0,.25) (0,.25,-x+.25) (.5,.25,x+.75) (-x,0,.75) (x+.5,0,.25) (.75,-x,0) (.25,x+.5,0) (0,.75,-x) (0,.25,x+.5) (.25,-x+.75,0) (.25,x+.25,.5) (-x+.25,.5,.75) (x+.75,0,.75) (0,.75,x+.75) (.5,.75,-x+.25) \n 32 e .3. (x,x,x) (-x+.5,-x,x+.5) (-x,x+.5,-x+.5) (x+.5,-x+.5,-x) (x+.75,x+.25,-x+.25) (-x+.75,-x+.75,-x+.75) (x+.25,-x+.25,x+.75) (-x+.25,x+.75,x+.25) (-x,-x,-x) (x+.5,x,-x+.5) (x,-x+.5,x+.5) (-x+.5,x+.5,x) (-x+.25,-x+.75,x+.75) (x+.25,x+.25,x+.25) (-x+.75,x+.75,-x+.25) (x+.75,-x+.25,-x+.75) \n 24 d -4.. (.375,0,.25) (.125,0,.75) (.25,.375,0) (.75,.125,0) (0,.25,.375) (0,.75,.125) (.75,.625,0) (.75,.375,.5) (.125,.5,.25) (.875,0,.25) (0,.25,.875) (.5,.25,.125) \n 24 c 2.22 (.125,0,.25) (.375,0,.75) (.25,.125,0) (.75,.375,0) (0,.25,.125) (0,.75,.375) (.875,0,.75) (.625,0,.25) (.75,.875,0) (.25,.625,0) (0,.75,.875) (0,.25,.625) \n 16 b .32 (.125,.125,.125) (.375,.875,.625) (.875,.625,.375) (.625,.375,.875) (.875,.875,.875) (.625,.125,.375) (.125,.375,.625) (.375,.625,.125) \n 16 a .-3. (0,0,0) (.5,0,.5) (0,.5,.5) (.5,.5,0) (.75,.25,.25) (.75,.75,.75) (.25,.25,.75) (.25,.75,.25) \n "; //DX20200520 - missing space between 12th and 13th position for Wyckoff letter h

  sgs.push_back(sg1);
  sgs.push_back(sg2);
  sgs.push_back(sg3);
  sgs.push_back(sg4);
  sgs.push_back(sg5);
  sgs.push_back(sg6);
  sgs.push_back(sg7);
  sgs.push_back(sg8);
  sgs.push_back(sg9);
  sgs.push_back(sg10);
  sgs.push_back(sg11);
  sgs.push_back(sg12);
  sgs.push_back(sg13);
  sgs.push_back(sg14);
  sgs.push_back(sg15);
  sgs.push_back(sg16);
  sgs.push_back(sg17);
  sgs.push_back(sg18);
  sgs.push_back(sg19);
  sgs.push_back(sg20);
  sgs.push_back(sg21);
  sgs.push_back(sg22);
  sgs.push_back(sg23);
  sgs.push_back(sg24);
  sgs.push_back(sg25);
  sgs.push_back(sg26);
  sgs.push_back(sg27);
  sgs.push_back(sg28);
  sgs.push_back(sg29);
  sgs.push_back(sg30);
  sgs.push_back(sg31);
  sgs.push_back(sg32);
  sgs.push_back(sg33);
  sgs.push_back(sg34);
  sgs.push_back(sg35);
  sgs.push_back(sg36);
  sgs.push_back(sg37);
  sgs.push_back(sg38);
  sgs.push_back(sg39);
  sgs.push_back(sg40);
  sgs.push_back(sg41);
  sgs.push_back(sg42);
  sgs.push_back(sg43);
  sgs.push_back(sg44);
  sgs.push_back(sg45);
  sgs.push_back(sg46);
  sgs.push_back(sg47);
  sgs.push_back(sg48);
  sgs.push_back(sg49);
  sgs.push_back(sg50);
  sgs.push_back(sg51);
  sgs.push_back(sg52);
  sgs.push_back(sg53);
  sgs.push_back(sg54);
  sgs.push_back(sg55);
  sgs.push_back(sg56);
  sgs.push_back(sg57);
  sgs.push_back(sg58);
  sgs.push_back(sg59);
  sgs.push_back(sg60);
  sgs.push_back(sg61);
  sgs.push_back(sg62);
  sgs.push_back(sg63);
  sgs.push_back(sg64);
  sgs.push_back(sg65);
  sgs.push_back(sg66);
  sgs.push_back(sg67);
  sgs.push_back(sg68);
  sgs.push_back(sg69);
  sgs.push_back(sg70);
  sgs.push_back(sg71);
  sgs.push_back(sg72);
  sgs.push_back(sg73);
  sgs.push_back(sg74);
  sgs.push_back(sg75);
  sgs.push_back(sg76);
  sgs.push_back(sg77);
  sgs.push_back(sg78);
  sgs.push_back(sg79);
  sgs.push_back(sg80);
  sgs.push_back(sg81);
  sgs.push_back(sg82);
  sgs.push_back(sg83);
  sgs.push_back(sg84);
  sgs.push_back(sg85);
  sgs.push_back(sg86);
  sgs.push_back(sg87);
  sgs.push_back(sg88);
  sgs.push_back(sg89);
  sgs.push_back(sg90);
  sgs.push_back(sg91);
  sgs.push_back(sg92);
  sgs.push_back(sg93);
  sgs.push_back(sg94);
  sgs.push_back(sg95);
  sgs.push_back(sg96);
  sgs.push_back(sg97);
  sgs.push_back(sg98);
  sgs.push_back(sg99);
  sgs.push_back(sg100);
  sgs.push_back(sg101);
  sgs.push_back(sg102);
  sgs.push_back(sg103);
  sgs.push_back(sg104);
  sgs.push_back(sg105);
  sgs.push_back(sg106);
  sgs.push_back(sg107);
  sgs.push_back(sg108);
  sgs.push_back(sg109);
  sgs.push_back(sg110);
  sgs.push_back(sg111);
  sgs.push_back(sg112);
  sgs.push_back(sg113);
  sgs.push_back(sg114);
  sgs.push_back(sg115);
  sgs.push_back(sg116);
  sgs.push_back(sg117);
  sgs.push_back(sg118);
  sgs.push_back(sg119);
  sgs.push_back(sg120);
  sgs.push_back(sg121);
  sgs.push_back(sg122);
  sgs.push_back(sg123);
  sgs.push_back(sg124);
  sgs.push_back(sg125);
  sgs.push_back(sg126);
  sgs.push_back(sg127);
  sgs.push_back(sg128);
  sgs.push_back(sg129);
  sgs.push_back(sg130);
  sgs.push_back(sg131);
  sgs.push_back(sg132);
  sgs.push_back(sg133);
  sgs.push_back(sg134);
  sgs.push_back(sg135);
  sgs.push_back(sg136);
  sgs.push_back(sg137);
  sgs.push_back(sg138);
  sgs.push_back(sg139);
  sgs.push_back(sg140);
  sgs.push_back(sg141);
  sgs.push_back(sg142);
  sgs.push_back(sg143);
  sgs.push_back(sg144);
  sgs.push_back(sg145);
  sgs.push_back(sg146);
  sgs.push_back(sg147);
  sgs.push_back(sg148);
  sgs.push_back(sg149);
  sgs.push_back(sg150);
  sgs.push_back(sg151);
  sgs.push_back(sg152);
  sgs.push_back(sg153);
  sgs.push_back(sg154);
  sgs.push_back(sg155);
  sgs.push_back(sg156);
  sgs.push_back(sg157);
  sgs.push_back(sg158);
  sgs.push_back(sg159);
  sgs.push_back(sg160);
  sgs.push_back(sg161);
  sgs.push_back(sg162);
  sgs.push_back(sg163);
  sgs.push_back(sg164);
  sgs.push_back(sg165);
  sgs.push_back(sg166);
  sgs.push_back(sg167);
  sgs.push_back(sg168);
  sgs.push_back(sg169);
  sgs.push_back(sg170);
  sgs.push_back(sg171);
  sgs.push_back(sg172);
  sgs.push_back(sg173);
  sgs.push_back(sg174);
  sgs.push_back(sg175);
  sgs.push_back(sg176);
  sgs.push_back(sg177);
  sgs.push_back(sg178);
  sgs.push_back(sg179);
  sgs.push_back(sg180);
  sgs.push_back(sg181);
  sgs.push_back(sg182);
  sgs.push_back(sg183);
  sgs.push_back(sg184);
  sgs.push_back(sg185);
  sgs.push_back(sg186);
  sgs.push_back(sg187);
  sgs.push_back(sg188);
  sgs.push_back(sg189);
  sgs.push_back(sg190);
  sgs.push_back(sg191);
  sgs.push_back(sg192);
  sgs.push_back(sg193);
  sgs.push_back(sg194);
  sgs.push_back(sg195);
  sgs.push_back(sg196);
  sgs.push_back(sg197);
  sgs.push_back(sg198);
  sgs.push_back(sg199);
  sgs.push_back(sg200);
  sgs.push_back(sg201);
  sgs.push_back(sg202);
  sgs.push_back(sg203);
  sgs.push_back(sg204);
  sgs.push_back(sg205);
  sgs.push_back(sg206);
  sgs.push_back(sg207);
  sgs.push_back(sg208);
  sgs.push_back(sg209);
  sgs.push_back(sg210);
  sgs.push_back(sg211);
  sgs.push_back(sg212);
  sgs.push_back(sg213);
  sgs.push_back(sg214);
  sgs.push_back(sg215);
  sgs.push_back(sg216);
  sgs.push_back(sg217);
  sgs.push_back(sg218);
  sgs.push_back(sg219);
  sgs.push_back(sg220);
  sgs.push_back(sg221);
  sgs.push_back(sg222);
  sgs.push_back(sg223);
  sgs.push_back(sg224);
  sgs.push_back(sg225);
  sgs.push_back(sg226);
  sgs.push_back(sg227);
  sgs.push_back(sg228);
  sgs.push_back(sg229);
  sgs.push_back(sg230);

  gl_sgs = sgs;

  return true;
}
//} //namespace SYM

#endif
