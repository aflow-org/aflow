// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Richard H. Taylor
// Updated by David Hicks
// d.hicks@duke.edu

#ifndef _AFLOW_SYMMETRY_SPACEGROUP_CPP_
#define _AFLOW_SYMMETRY_SPACEGROUP_CPP_
#include <string>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm>
#include <pthread.h>
#include "aflow_symmetry_spacegroup.h"
#include "aflow.h"

//======================= Declare global variables ======================= //
//DX20190215 [OBSOLETE]#ifdef AFLOW_SYMMETRY_MULTITHREADS_ENABLE
//DX20190215 [OBSOLETE] thread_local double _SYM_TOL_;
//DX20190215 [OBSOLETE]#else
//DX20190215 [OBSOLETE]double _SYM_TOL_;
//DX20190215 [OBSOLETE]#endif

// ***************************************************************************
// Set Tolerance Function
// ***************************************************************************
//DX20190215 [OBSOLETE] namespace SYM {
//DX20190215 [OBSOLETE]   void SetTolerance(const double& starting_tol) {
//DX20190215 [OBSOLETE]    _SYM_TOL_ = starting_tol;
//DX20190215 [OBSOLETE]  }
//DX20190215 [OBSOLETE]} //namespace SYM

// ========== Using namespace standard ========== //
using namespace std;

// **************************************************************************************************************************************************
//  Orient Conventional Cell
// **************************************************************************************************************************************************
// For use after finding conventional form of crystal. Standardizes the unit cell, correcting
// underlying Cartesian coordinate system and orientation of lattice basis vectors.
namespace SYM {
  void orient(xstructure& xstr, bool update_atom_positions) {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    stringstream message;
    double sym_tol = xstr.sym_eps; //DX20190215 - added sym_tol 
    double tol = _ZERO_TOL_;
    xvector<double> a;
    a(1) = 1;
    xvector<double> b;
    b(2) = 1;
    xvector<double> c;
    c(3) = 1;
    xvector<double> tmpvec;
    xvector<double> zero;
    zero(1) = 0;
    zero(2) = 0;
    zero(3) = 0;

    char lattice_label = xstr.lattice_label_ITC;

    //print(lattice_basis)
    xmatrix<double> lattice_basis_xmat = xstr.lattice;
    vector<xvector<double> > lattice_basis;
    lattice_basis.push_back(xstr.lattice(1)); //DX20190905 - SYM::extract_row -> aurostd::xmatrix::operator(int)
    lattice_basis.push_back(xstr.lattice(2)); //DX20190905 - SYM::extract_row -> aurostd::xmatrix::operator(int)
    lattice_basis.push_back(xstr.lattice(3)); //DX20190905 - SYM::extract_row -> aurostd::xmatrix::operator(int)

    xmatrix<double> Linv = aurostd::inverse(xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]));

    deque<_atom> atomic_basis_ = xstr.atoms;

    // ========================================== ORIENTATION ROUTINE: ========================================== //
    if(true == true) {
      //MORE ORIENTATION (CHANGING UNDERLYING CARTESIAN COORD. SYSTEM)

      // ==================== MONOCLINIC ==================== //
      // ===== UNIQUE AXIS: b ===== //
      if(lattice_label == 'b') {
        //GET ANGLE BETWEEN lattice_basis[1] AND (0 0 1)
        double angle = aurostd::angle(lattice_basis[1], b);
        if(aurostd::abs(angle - Pi_r) > tol && aurostd::abs(angle) > tol) {
          xvector<double> rotaxis = aurostd::vector_product(lattice_basis[1], b); //DX20190905 - SYM::CrossPro -> aurostd::vector_product 
          rotaxis/=aurostd::modulus(rotaxis); //DX20190905 - SYM::normalize -> divide by aurostd::modulus
          Screw S;  //use screw class for rotation
          S.get_screw_direct(rotaxis, zero, 2 * Pi_r / angle);
          lattice_basis[0] = S * lattice_basis[0];
          lattice_basis[1] = S * lattice_basis[1];
          lattice_basis[2] = S * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S, new_lattice); //DX20190905 - return changed to void
          }
        }
        //DX: ORIG but it breaks for LIB2/LIB/AgAl/f342
        //for (int i = 0; i < 3; i++) {
        //  if(lattice_basis[1](i + 1) < -tol) {
        //    lattice_basis[1] = -1 * lattice_basis[1];
        //    minus_one(atomic_basis_, i + 1);
        //  }
        //}
        //DX: Orig above, new below
        if(lattice_basis[1](2) < -tol) {
          lattice_basis[1] = -1 * lattice_basis[1];
          minus_one(atomic_basis_, 2);
        }
        if(det(xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2])) < -tol) {
          xvector<double> lattice_tmp = lattice_basis[0];
          lattice_basis[0] = lattice_basis[2];
          lattice_basis[2] = lattice_tmp;
          swap_columns(atomic_basis_, 1, 3);
        }
      }

      // [BELOW IS OBSOLETE]

      // ====== UNIQUE AXIS: c ===== //
      if(lattice_label == 'm'){
        //GET ANGLE BETWEEN lattice_basis[2] AND (0 0 1) 
        double angle = aurostd::angle(lattice_basis[2],c);
        if(aurostd::abs( angle - Pi_r ) > tol && aurostd::abs ( angle ) > tol){
          xvector<double> rotaxis = aurostd::vector_product(lattice_basis[2],c); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          rotaxis/=aurostd::modulus(rotaxis); //DX20190905 - SYM::normalize -> divide by aurostd::modulus
          Screw S; //use screw class for rotation
          S.get_screw_direct(rotaxis,zero,2*Pi_r/angle);
          lattice_basis[0] =  S*lattice_basis[0];
          lattice_basis[1] =  S*lattice_basis[1];
          lattice_basis[2] =  S*lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0],lattice_basis[1],lattice_basis[2]);
          updateAtomPositions(atomic_basis_,S,new_lattice); //DX20190905 - return changed to void
        }
        //TESTif(lattice_basis[2](1) < -tol || lattice_basis[2](2) < -tol || lattice_basis[2](3) < -tol){
        //TEST  lattice_basis[2] = -1*lattice_basis[2];
        //TEST  minus_one(atomic_basis_,3);
        //TEST}
        if(lattice_basis[2](3) < -tol) {
          lattice_basis[2] = -1*lattice_basis[2];
          minus_one(atomic_basis_,3);
        }
        if(det(xvec2xmat(lattice_basis[0],lattice_basis[1],lattice_basis[2])) <0.0){
          xvector<double> lattice_tmp=lattice_basis[0];
          lattice_basis[0]=lattice_basis[1];
          lattice_basis[1]=lattice_tmp;
          swap_columns(atomic_basis_, 1, 2);
        }
      } 


      // ==================== ORTHORHOMBIC ==================== //
      if(lattice_label == 'o') {
        //GET ANGLE BETWEEN lattice_basis[0] AND (1 0 0)
        double angle = aurostd::angle(lattice_basis[0], a);
        if(aurostd::abs(angle - Pi_r) > tol && aurostd::abs(angle) > tol) {
          xvector<double> rotaxis = aurostd::vector_product(lattice_basis[0], a); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          rotaxis/=aurostd::modulus(rotaxis); //DX20190905 - SYM::normalize -> divide by aurostd::modulus
          Screw S;  //use screw class for rotation
          S.get_screw_direct(rotaxis, zero, 2 * Pi_r / angle);
          lattice_basis[0] = S * lattice_basis[0];
          lattice_basis[1] = S * lattice_basis[1];
          lattice_basis[2] = S * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S, new_lattice); //DX20190905 - return changed to void
          }
        }
        //GET ANGLE BETWEEN lattice_basis[2] AND (0 0 1)
        double angle1 = aurostd::angle(lattice_basis[2], c);
        if(aurostd::abs(angle1 - Pi_r) > tol && aurostd::abs(angle1) > tol) {
          xvector<double> rotaxis1 = aurostd::vector_product(lattice_basis[2], c); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          rotaxis1/=aurostd::modulus(rotaxis1); //DX20190905 - SYM::normalize -> divide by aurostd::modulus;
          Screw S;  //use screw class for rotation
          S.get_screw_direct(rotaxis1, zero, 2 * Pi_r / angle1);
          lattice_basis[0] = S * lattice_basis[0];
          lattice_basis[1] = S * lattice_basis[1];
          lattice_basis[2] = S * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S, new_lattice); //DX20190905 - return changed to void
          }
        }
        //GET ANGLE BETWEEN lattice_basis[1] AND (0 1 0)
        double angle2 = aurostd::angle(lattice_basis[1], b);
        if(aurostd::abs(angle2 - Pi_r) > tol && aurostd::abs(angle2) > tol) {
          xvector<double> rotaxis2 = aurostd::vector_product(lattice_basis[1], b); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          rotaxis2/=aurostd::modulus(rotaxis2); //DX20190905 - SYM::normalize -> divide by aurostd::modulus
          Screw S2;  //use screw class for rotation
          S2.get_screw_direct(rotaxis2, zero, 2 * Pi_r / angle2);
          lattice_basis[0] = S2 * lattice_basis[0];
          lattice_basis[1] = S2 * lattice_basis[1];
          lattice_basis[2] = S2 * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S2, new_lattice); //DX20190905 - return changed to void
          }
        }
        for (int i = 0; i < 3; i++) {
          if(lattice_basis[i](i + 1) < -tol) {
            lattice_basis[i] = -1 * lattice_basis[i];
            minus_one(atomic_basis_, i + 1);
          }
        }
      }

      // ==================== TETRAGONAL ==================== //
      if(lattice_label == 't') {
        //GET ANGLE BETWEEN lattice_basis[2] AND (0 0 1)
        double angle = aurostd::angle(lattice_basis[2], c);
        if(aurostd::abs(angle - Pi_r) > tol && aurostd::abs(angle) > tol) {
          xvector<double> rotaxis = aurostd::vector_product(lattice_basis[2], c); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          rotaxis/=aurostd::modulus(rotaxis); //DX20190905 - SYM::normalize -> divide by aurostd::modulus
          Screw S;  //use screw class for rotation
          S.get_screw_direct(rotaxis, zero, 2 * Pi_r / angle);
          lattice_basis[0] = S * lattice_basis[0];
          lattice_basis[1] = S * lattice_basis[1];
          lattice_basis[2] = S * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S, new_lattice); //DX20190905 - return changed to void
          }
        }
        //GET ANGLE BETWEEN lattice_basis[1] AND (0 1 0)
        double angle2 = aurostd::angle(lattice_basis[1], b);
        if(aurostd::abs(angle2 - Pi_r) > tol && aurostd::abs(angle2) > tol) {
          xvector<double> rotaxis2 = aurostd::vector_product(lattice_basis[1], b); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          rotaxis2/=aurostd::modulus(rotaxis2); //DX20190905 - SYM::normalize -> divide by aurostd::modulus
          Screw S2;  //use screw class for rotation
          S2.get_screw_direct(rotaxis2, zero, 2 * Pi_r / angle2);
          lattice_basis[0] = S2 * lattice_basis[0];
          lattice_basis[1] = S2 * lattice_basis[1];
          lattice_basis[2] = S2 * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S2, new_lattice); //DX20190905 - return changed to void
          }
        }
        for (int i = 0; i < 3; i++) {
          if(lattice_basis[i](i + 1) < -tol) {
            lattice_basis[i] = -1 * lattice_basis[i];
            minus_one(atomic_basis_, i + 1);
          }
        }
      }

      // ==================== CUBIC ==================== //
      if(lattice_label == 'c') {
        //GET ANGLE BETWEEN lattice_basis[2] AND (0 0 1)
        double angle = aurostd::angle(lattice_basis[2], c);
        if(aurostd::abs(angle - Pi_r) > tol && aurostd::abs(angle) > tol) {
          xvector<double> rotaxis = aurostd::vector_product(lattice_basis[2], c); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          rotaxis/=aurostd::modulus(rotaxis); //DX20190905 - SYM::normalize -> divide by aurostd::modulus
          Screw S;  //use screw class for rotation
          S.get_screw_direct(rotaxis, zero, 2 * Pi_r / angle);
          lattice_basis[0] = S * lattice_basis[0];
          lattice_basis[1] = S * lattice_basis[1];
          lattice_basis[2] = S * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S, new_lattice); //DX20190905 - return changed to void
          }
        }
        //GET ANGLE BETWEEN lattice_basis[1] AND (0 1 0)
        double angle2 = aurostd::angle(lattice_basis[1], b);
        if(aurostd::abs(angle2 - Pi_r) > tol && aurostd::abs(angle2) > tol) {
          xvector<double> rotaxis2 = aurostd::vector_product(lattice_basis[1], b); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          rotaxis2/=aurostd::modulus(rotaxis2); //DX20190905 - SYM::normalize -> divide by aurostd::modulus
          Screw S2;  //use screw class for rotation
          S2.get_screw_direct(rotaxis2, zero, 2 * Pi_r / angle2);
          lattice_basis[0] = S2 * lattice_basis[0];
          lattice_basis[1] = S2 * lattice_basis[1];
          lattice_basis[2] = S2 * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S2, new_lattice); //DX20190905 - return changed to void
          }
        }
        for (int i = 0; i < 3; i++) {
          if(lattice_basis[i](i + 1) < -tol) {
            lattice_basis[i] = -1 * lattice_basis[i];
            minus_one(atomic_basis_, i + 1);
          }
        }
      }

      // ==================== HEXAGONAL ==================== //
      if(lattice_label == 'h') {
        double angle;
        xvector<double> rotaxis;
        Screw S;
        if(lattice_basis[2](3) < -tol) {
          lattice_basis[2] = -1 * lattice_basis[2];
          minus_one(atomic_basis_, 3);
        }
        xvector<double> tmpc = c * aurostd::modulus(lattice_basis[2]);
        //If lattice_basis[2] is already parallel to 001, this is not necessary
        if(aurostd::isdifferent(tmpc, lattice_basis[2],_ZERO_TOL_)) { //DX20190805 - !SYM::vec_compare -> aurostd::isdifferent
          angle = aurostd::angle(lattice_basis[2], c);
          rotaxis = aurostd::vector_product(lattice_basis[2], c); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          S.get_screw_direct(rotaxis, zero, 2 * Pi_r / angle);
          lattice_basis[0] = S * lattice_basis[0];
          lattice_basis[1] = S * lattice_basis[1];
          lattice_basis[2] = S * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S, new_lattice); //DX20190905 - return changed to void
          }
        }
        //GET ANGLE BETWEEN lattice_basis[0] AND (Sqrt(3)/2 1/2 0)
        xvector<double> fromx;  // 30 degrees above x axis
        fromx(1) = sqrt(3.0) / 2.0;
        fromx(2) = 1.0 / 2.0;
        angle = aurostd::angle(lattice_basis[0], fromx);
        if(aurostd::abs(angle - Pi_r) > tol && aurostd::abs(angle) > tol) {
          rotaxis = aurostd::vector_product(lattice_basis[0], fromx); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          rotaxis/=aurostd::modulus(rotaxis); //DX20190905 - SYM::normalize -> divide by aurostd::modulus
          S.get_screw_direct(rotaxis, zero, 2 * Pi_r / angle);
          lattice_basis[0] = S * lattice_basis[0];
          lattice_basis[1] = S * lattice_basis[1];
          lattice_basis[2] = S * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S, new_lattice); //DX20190905 - return changed to void
          }
        }
        xvector<double> tmpb = -b * aurostd::modulus(lattice_basis[1]);
        if(aurostd::identical(tmpb, lattice_basis[1], _ZERO_TOL_)) { //DX20190805 - SYM::vec_compare -> aurostd::identical
          rotaxis = aurostd::vector_product(lattice_basis[0], lattice_basis[1]); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          angle = -120.0 * Pi_r / 180.0;
          S.get_screw_direct(rotaxis, zero, 2 * Pi_r / angle);
          lattice_basis[0] = S * lattice_basis[0];
          lattice_basis[1] = S * lattice_basis[1];
          lattice_basis[2] = S * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S, new_lattice); //DX20190905 - return changed to void
          }
        }

        if(aurostd::determinant(xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2])) < -tol) {
          tmpvec = lattice_basis[0];
          lattice_basis[0] = lattice_basis[1];
          lattice_basis[1] = tmpvec;
          swap_columns(atomic_basis_, 1, 2);
        }
      }

      // ==================== RHOMBOHEDRAL ==================== //
      if(lattice_label == 'R') {
        double angle;
        xvector<double> rotaxis;
        Screw S;
        if(lattice_basis[2](3) < -tol) {
          lattice_basis[2] = -1 * lattice_basis[2];
          minus_one(atomic_basis_, 3);
        }
        xvector<double> tmpc = c * aurostd::modulus(lattice_basis[2]);
        //If lattice_basis[2] is already parallel to 001, this is not necessary
        if(aurostd::isdifferent(tmpc, lattice_basis[2], _ZERO_TOL_)) { //DX20190805 - !SYM::vec_compare -> aurostd::isdifferent
          angle = aurostd::angle(lattice_basis[2], c);
          rotaxis = aurostd::vector_product(lattice_basis[2], c); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          S.get_screw_direct(rotaxis, zero, 2 * Pi_r / angle);
          lattice_basis[0] = S * lattice_basis[0];
          lattice_basis[1] = S * lattice_basis[1];
          lattice_basis[2] = S * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S, new_lattice); //DX20190905 - return changed to void
          }
        }
        //GET ANGLE BETWEEN lattice_basis[0] AND (Sqrt(3)/2 1/2 0)
        xvector<double> fromx;  // 30 degrees above x axis
        fromx(1) = sqrt(3.0) / 2.0;
        fromx(2) = 1.0 / 2.0;
        angle = aurostd::angle(lattice_basis[0], fromx);
        if(aurostd::abs(angle) > tol && aurostd::abs(angle - Pi_r) > tol) {
          rotaxis = aurostd::vector_product(lattice_basis[0], fromx); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          rotaxis/=aurostd::modulus(rotaxis); //DX20190905 - SYM::normalize -> divide by aurostd::modulus
          S.get_screw_direct(rotaxis, zero, 2 * Pi_r / angle);
          lattice_basis[0] = S * lattice_basis[0];
          lattice_basis[1] = S * lattice_basis[1];
          lattice_basis[2] = S * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S, new_lattice); //DX20190905 - return changed to void
          }
        }
        xvector<double> tmpb = -b * aurostd::modulus(lattice_basis[1]);
        if(aurostd::identical(tmpb, lattice_basis[1], _ZERO_TOL_)) { //DX20190805 - SYM::vec_compare -> aurostd::identical
          rotaxis = aurostd::vector_product(lattice_basis[0], lattice_basis[1]); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          angle = -120.0 * Pi_r / 180.0;
          S.get_screw_direct(rotaxis, zero, 2 * Pi_r / angle);
          lattice_basis[0] = S * lattice_basis[0];
          lattice_basis[1] = S * lattice_basis[1];
          lattice_basis[2] = S * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S, new_lattice); //DX20190905 - return changed to void
          }
        }

        if(aurostd::determinant(xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2])) < -tol) {
          tmpvec = lattice_basis[0];
          lattice_basis[0] = lattice_basis[1];
          lattice_basis[1] = tmpvec;
          swap_columns(atomic_basis_, 1, 2);
        }
        lattice_basis_xmat = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
        // Check for reverse/obverse setting
        // [OBSOLETE] DX- This is done later in the code: isObverseSetting(lattice_basis_xmat,atomic_basis_,xstr.dist_nn_min,_SYM_TOL_);
        lattice_basis[0] = lattice_basis_xmat(1); //DX20190905 - SYM::extract_row -> aurostd::xmatrix::operator(int)
        lattice_basis[1] = lattice_basis_xmat(2); //DX20190905 - SYM::extract_row -> aurostd::xmatrix::operator(int)
        lattice_basis[2] = lattice_basis_xmat(3); //DX20190905 - SYM::extract_row -> aurostd::xmatrix::operator(int)
      }
      if(lattice_label == 'r') {
        //cerr << "do not NEED to orient rhombohedral cell, rhombohedral axes" << endl;
      }

      // ==================== RHOMBOHEDRAL IN HEXAGONAL SETTING ==================== //
      if(lattice_label == 'X') {  //
        if(LDEBUG) { cerr << __AFLOW_FUNC__ << " RHOMBOHEDRAL IN HEXAGONAL SETTING" << endl; }
        double angle;
        xvector<double> rotaxis;
        Screw S;
        if(lattice_basis[2](3) < -tol) {
          lattice_basis[2] = -1 * lattice_basis[2];
          minus_one(atomic_basis_, 3);
        }
        //GET ANGLE BETWEEN lattice_basis[0] AND (Sqrt(3)/2 1/2 0)
        xvector<double> fromx;  // 30 degrees above x axis
        fromx(1) = sqrt(3.0) / 2.0;
        fromx(2) = 1.0 / 2.0;
        angle = aurostd::angle(lattice_basis[0], fromx);
        if(aurostd::abs(angle - Pi_r) > tol && aurostd::abs(angle) > tol) {
          rotaxis = aurostd::vector_product(lattice_basis[0], fromx); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          rotaxis/=aurostd::modulus(rotaxis); //DX20190905 - SYM::normalize -> divide by aurostd::modulus
          S.get_screw_direct(rotaxis, zero, 2 * Pi_r / angle);
          lattice_basis[0] = S * lattice_basis[0];
          lattice_basis[1] = S * lattice_basis[1];
          lattice_basis[2] = S * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S, new_lattice); //DX20190905 - return changed to void
          }
        }
        xvector<double> tmpb = -b * aurostd::modulus(lattice_basis[1]);
        if(aurostd::identical(tmpb, lattice_basis[1],_ZERO_TOL_)) { //DX20190805 - SYM::vec_compare -> aurostd::identical
          rotaxis = aurostd::vector_product(lattice_basis[0], lattice_basis[1]); //DX20190905 - SYM::CrossPro -> aurostd::vector_product
          angle = -120.0 * Pi_r / 180.0;
          S.get_screw_direct(rotaxis, zero, 2 * Pi_r / angle);
          lattice_basis[0] = S * lattice_basis[0];
          lattice_basis[1] = S * lattice_basis[1];
          lattice_basis[2] = S * lattice_basis[2];
          xmatrix<double> new_lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
          if(update_atom_positions){
            updateAtomPositions(atomic_basis_, S, new_lattice); //DX20190905 - return changed to void
          }
        }
        //print(lattice_basis);
        lattice_basis_xmat = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
        xmatrix<double> f2c = trasp(lattice_basis_xmat);
        xmatrix<double> c2f = inverse(trasp(lattice_basis_xmat));
        bool skew = SYM::isLatticeSkewed(lattice_basis_xmat, xstr.dist_nn_min, sym_tol); //DX20190215 - _SYM_TOL_ to sym_tol

        //xb();
        //print(atomic_basis_);
        _atom obverse1;
        _atom obverse2;
        _atom obverse3;
        obverse1.fpos(1) = 2.0 / 3.0;
        obverse1.fpos(2) = 1.0 / 3.0;
        obverse1.fpos(3) = 1.0 / 3.0;
        obverse1.name = "0";
        obverse2.fpos(1) = 1.0 / 3.0;
        obverse2.fpos(2) = 2.0 / 3.0;
        obverse2.fpos(3) = 2.0 / 3.0;
        obverse2.name = "0";
        obverse3.name = "0";
        uint match_type1 = 0;
        uint match_type2 = 0;
        uint match_type3 = 0;

        if(aurostd::determinant(xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2])) < -tol) {
          tmpvec = lattice_basis[0];
          lattice_basis[0] = lattice_basis[1];
          lattice_basis[1] = tmpvec;
          swap_columns(atomic_basis_, 1, 2);
        }
        if((!SYM::FPOSMatch(atomic_basis_, obverse1, match_type1, lattice_basis_xmat, f2c, skew, sym_tol) || //DX20190215 - _SYM_TOL_ to sym_tol //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
              !SYM::FPOSMatch(atomic_basis_, obverse2, match_type2, lattice_basis_xmat, f2c, skew, sym_tol) || //DX20190215 - _SYM_TOL_ to sym_tol //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name 
              !SYM::FPOSMatch(atomic_basis_, obverse3, match_type3, lattice_basis_xmat, f2c, skew, sym_tol)) && //DX20190215 - _SYM_TOL_ to sym_tol //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
            (match_type1 != match_type2 || match_type2 != match_type3 || match_type1 != match_type3)) {
          lattice_basis[0] = -1 * lattice_basis[0];
          minus_one(atomic_basis_, 1);
          lattice_basis[1] = -1 * lattice_basis[1];
          minus_one(atomic_basis_, 2);
        }
        if(!(SYM::FPOSMatch(atomic_basis_, obverse1, match_type1, lattice_basis_xmat, f2c, skew, sym_tol) && //DX20190215 - _SYM_TOL_ to sym_tol //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
              SYM::FPOSMatch(atomic_basis_, obverse2, match_type2, lattice_basis_xmat, f2c, skew, sym_tol) && //DX20190215 - _SYM_TOL_ to sym_tol //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
              SYM::FPOSMatch(atomic_basis_, obverse3, match_type3, lattice_basis_xmat, f2c, skew, sym_tol)) && //DX20190215 - _SYM_TOL_ to sym_tol //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
            (match_type1 != match_type2 || match_type2 != match_type3 || match_type1 != match_type3)) {
          lattice_basis[0] = -1 * lattice_basis[0];
          minus_one(atomic_basis_, 1);
          lattice_basis[2] = -1 * lattice_basis[2];
          minus_one(atomic_basis_, 3);
        }
        if(!(SYM::FPOSMatch(atomic_basis_, obverse1, match_type1, lattice_basis_xmat, f2c, skew, sym_tol) && //DX20190215 - _SYM_TOL_ to sym_tol //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
              SYM::FPOSMatch(atomic_basis_, obverse2, match_type2, lattice_basis_xmat, f2c, skew, sym_tol) && //DX20190215 - _SYM_TOL_ to sym_tol //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
              SYM::FPOSMatch(atomic_basis_, obverse3, match_type3, lattice_basis_xmat, f2c, skew, sym_tol)) && //DX20190215 - _SYM_TOL_ to sym_tol //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
            (match_type1 != match_type2 || match_type2 != match_type3 || match_type1 != match_type3)) {
          lattice_basis[1] = -1 * lattice_basis[1];
          minus_one(atomic_basis_, 2);
          lattice_basis[2] = -1 * lattice_basis[2];
          minus_one(atomic_basis_, 3);
        }
        if(!(SYM::FPOSMatch(atomic_basis_, obverse1, match_type1, lattice_basis_xmat, f2c, skew, sym_tol) && //DX20190215 - _SYM_TOL_ to sym_tol //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
              SYM::FPOSMatch(atomic_basis_, obverse2, match_type2, lattice_basis_xmat, f2c, skew, sym_tol) && //DX20190215 - _SYM_TOL_ to sym_tol //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
              SYM::FPOSMatch(atomic_basis_, obverse3, match_type3, lattice_basis_xmat, f2c, skew, sym_tol)) && //DX20190215 - _SYM_TOL_ to sym_tol //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
            (match_type1 != match_type2 || match_type2 != match_type3 || match_type1 != match_type3)) {
          message << "PROBLEM TRANFORMING TO OBVERSE SETTING. QUITING PROGRAM [dir = " << xstr.directory << "].";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
        }
      }
      xstr.lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);

      xstr.atoms = atomic_basis_;
    }
  }
} //namespace SYM

// **************************************************************************************************************************************************
// Is obverse setting ?
// **************************************************************************************************************************************************
// Check and transform a rhombohedral system into the obverse setting if it is not already in this setting
namespace SYM {
  bool isObverseSetting(xstructure& xstr, double& tolerance) {
    return isObverseSetting(xstr.lattice, xstr.atoms, xstr.dist_nn_min, tolerance);
  }
} //namespace SYM

// ===== deque<_atom> input ===== //
namespace SYM {
  bool isObverseSetting(xmatrix<double>& lattice, deque<_atom>& atomic_basis_, double& dist_nn_min, double& tolerance) {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    vector<xvector<double> > lattice_basis;
    lattice_basis.push_back(lattice(1));
    lattice_basis.push_back(lattice(2));
    lattice_basis.push_back(lattice(3));
    xmatrix<double> f2c = trasp(lattice);
    xmatrix<double> c2f = inverse(trasp(lattice));
    bool skew = SYM::isLatticeSkewed(lattice, dist_nn_min, tolerance);
    //CHECK IF SETTING IS REVERSE OR OBVERSE. TRANSFORM TO OBVERSE AS NECESSARY
    _atom obverse1;
    _atom obverse2;
    _atom obverse3;
    obverse1.fpos(1) = 2.0 / 3.0;
    obverse1.fpos(2) = 1.0 / 3.0;
    obverse1.fpos(3) = 1.0 / 3.0;
    obverse1.name = "0";
    obverse2.fpos(1) = 1.0 / 3.0;
    obverse2.fpos(2) = 2.0 / 3.0;
    obverse2.fpos(3) = 2.0 / 3.0;
    obverse2.name = "0";
    obverse3.name = "0";
    _atom reverse1;
    _atom reverse2;
    _atom reverse3;
    reverse1.fpos(1) = 1.0 / 3.0;
    reverse1.fpos(2) = 2.0 / 3.0;
    reverse1.fpos(3) = 1.0 / 3.0;
    reverse1.name = "0";
    reverse2.fpos(1) = 2.0 / 3.0;
    reverse2.fpos(2) = 1.0 / 3.0;
    reverse2.fpos(3) = 2.0 / 3.0;
    reverse2.name = "0";
    reverse3.name = "0";

    // Also, check (1/2,0,0) and (1/2,0,1/2)
    _atom half_x;
    half_x.fpos(1) = 0.5;
    _atom half_x_shift1;
    half_x_shift1.fpos = half_x.fpos + reverse1.fpos;
    _atom half_x_shift2;
    half_x_shift2.fpos = half_x.fpos + reverse2.fpos;
    _atom half_xz;
    half_xz.fpos(1) = 0.5;
    half_xz.fpos(3) = 0.5;
    _atom half_xz_shift1;
    half_xz_shift1.fpos = half_xz.fpos + reverse1.fpos;
    _atom half_xz_shift2;
    half_xz_shift2.fpos = half_xz.fpos + reverse2.fpos;

    uint match_type1 = 0;
    uint match_type2 = 0;
    uint match_type3 = 0;
    bool obverse = true;
    if(SYM::FPOSMatch(atomic_basis_, reverse1, match_type1, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
        SYM::FPOSMatch(atomic_basis_, reverse2, match_type2, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
        SYM::FPOSMatch(atomic_basis_, reverse3, match_type3, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
        match_type1 == match_type2 && match_type2 == match_type3 && match_type1 == match_type3) {
      obverse = false;
    } 
    else if(SYM::FPOSMatch(atomic_basis_, half_x, match_type1, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
        SYM::FPOSMatch(atomic_basis_, half_x_shift1, match_type2, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
        SYM::FPOSMatch(atomic_basis_, half_x_shift2, match_type3, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
        match_type1 == match_type2 && match_type2 == match_type3 && match_type1 == match_type3) {
      obverse = false;
    } 
    else if(SYM::FPOSMatch(atomic_basis_, half_xz, match_type1, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name 
        SYM::FPOSMatch(atomic_basis_, half_xz_shift1, match_type2, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
        SYM::FPOSMatch(atomic_basis_, half_xz_shift2, match_type3, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
        match_type1 == match_type2 && match_type2 == match_type3 && match_type1 == match_type3) {
      obverse = false;
    }
    if(!obverse) {
      if(LDEBUG) { cerr << "SYM::isObverseSetting: RHL: In reversed setting, rotating to obverse setting" << endl; }
      //cerr << 1 << endl;
      transformToObverse(lattice, atomic_basis_);
      return false;
    } else {
      if(LDEBUG) { cerr << "SYM::isObverseSetting: TRIGONAL? Yes, but origin choice may prevent from finding orientation" << endl; }
      return true;
    }
  }
} //namespace SYM

// ===== deque<deque<_atom> > input ===== //
namespace SYM {
  bool isObverseSetting(xmatrix<double>& lattice, deque<deque<_atom> >& equivalent_atoms, double& dist_nn_min, double& tolerance) {
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    vector<xvector<double> > lattice_basis;
    lattice_basis.push_back(lattice(1));
    lattice_basis.push_back(lattice(2));
    lattice_basis.push_back(lattice(3));
    xmatrix<double> f2c = trasp(lattice);
    xmatrix<double> c2f = inverse(trasp(lattice));
    bool skew = SYM::isLatticeSkewed(lattice, dist_nn_min, tolerance);
    //CHECK IF SETTING IS REVERSE OR OBVERSE. TRANSFORM TO OBVERSE AS NECESSARY
    _atom obverse1;
    _atom obverse2;
    _atom obverse3;
    obverse1.fpos(1) = 2.0 / 3.0;
    obverse1.fpos(2) = 1.0 / 3.0;
    obverse1.fpos(3) = 1.0 / 3.0;
    obverse1.name = "0";
    obverse2.fpos(1) = 1.0 / 3.0;
    obverse2.fpos(2) = 2.0 / 3.0;
    obverse2.fpos(3) = 2.0 / 3.0;
    obverse2.name = "0";
    obverse3.name = "0";
    _atom reverse1;
    _atom reverse2;
    _atom reverse3;
    reverse1.fpos(1) = 1.0 / 3.0;
    reverse1.fpos(2) = 2.0 / 3.0;
    reverse1.fpos(3) = 1.0 / 3.0;
    reverse1.name = "0";
    reverse2.fpos(1) = 2.0 / 3.0;
    reverse2.fpos(2) = 1.0 / 3.0;
    reverse2.fpos(3) = 2.0 / 3.0;
    reverse2.name = "0";
    reverse3.name = "0";
    // Also, check (1/2,0,0) and (1/2,0,1/2)
    _atom half_x;
    half_x.fpos(1) = 0.5;
    _atom half_x_shift1;
    half_x_shift1.fpos = half_x.fpos + reverse1.fpos;
    _atom half_x_shift2;
    half_x_shift2.fpos = half_x.fpos + reverse2.fpos;
    _atom half_xz;
    half_xz.fpos(1) = 0.5;
    half_xz.fpos(3) = 0.5;
    _atom half_xz_shift1;
    half_xz_shift1.fpos = half_xz.fpos + reverse1.fpos;
    _atom half_xz_shift2;
    half_xz_shift2.fpos = half_xz.fpos + reverse2.fpos;
    uint match_type1 = 0;
    uint match_type2 = 0;
    uint match_type3 = 0;
    bool obverse = true;
    for (uint i = 0; i < equivalent_atoms.size(); i++) {
      if(SYM::FPOSMatch(equivalent_atoms[i], reverse1, match_type1, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name 
          SYM::FPOSMatch(equivalent_atoms[i], reverse2, match_type2, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
          SYM::FPOSMatch(equivalent_atoms[i], reverse3, match_type3, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
          match_type1 == match_type2 && match_type2 == match_type3 && match_type1 == match_type3) {
        obverse = false;
        break;
      } 
      else if(SYM::FPOSMatch(equivalent_atoms[i], half_x, match_type1, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
          SYM::FPOSMatch(equivalent_atoms[i], half_x_shift1, match_type2, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
          SYM::FPOSMatch(equivalent_atoms[i], half_x_shift2, match_type3, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
          match_type1 == match_type2 && match_type2 == match_type3 && match_type1 == match_type3) {
        obverse = false;
        break;
      } 
      else if(SYM::FPOSMatch(equivalent_atoms[i], half_xz, match_type1, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name 
          SYM::FPOSMatch(equivalent_atoms[i], half_xz_shift1, match_type2, lattice, f2c, skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
          SYM::FPOSMatch(equivalent_atoms[i], half_xz_shift2, match_type3, lattice, f2c,  skew, tolerance) && //DX20190619 - lattice_basis_xmat and f2c as input, remove "Atom" prefix from name
          match_type1 == match_type2 && match_type2 == match_type3 && match_type1 == match_type3) {
        obverse = false;
      }
    }
    if(!obverse) {
      if(LDEBUG) { cerr << "SYM::isObverseSetting: RHL: In reversed setting, rotating to obverse setting" << endl; }
      //cerr << 2 << endl;
      transformToObverse(lattice, equivalent_atoms);
      return false;
    } else {
      return true;
    }
  }
} //namespace SYM

// **************************************************************************************************************************************************
// Transform to obverse
// **************************************************************************************************************************************************
// Transform a rhombohedral system into the obverse setting
// ===== deque<_atom> input ===== //
namespace SYM {
  bool transformToObverse(xmatrix<double>& lattice, deque<_atom>& atoms) {
    vector<xvector<double> > lattice_basis;
    lattice_basis.push_back(lattice(1));
    lattice_basis.push_back(lattice(2));
    lattice_basis.push_back(lattice(3));
    //Rotate by 60 degrees
    xmatrix<double> R_60;
    R_60(1, 1) = cos(60.0 * Pi_r / 180);
    R_60(1, 2) = -sin(60.0 * Pi_r / 180);
    R_60(1, 3) = 0.0;
    R_60(2, 1) = sin(60.0 * Pi_r / 180);
    R_60(2, 2) = cos(60.0 * Pi_r / 180);
    R_60(2, 3) = 0.0;
    R_60(3, 1) = 0.0;
    R_60(3, 2) = 0.0;
    R_60(3, 3) = 1.0;
    xvector<double> lb1 = R_60 * lattice_basis[0];
    xvector<double> lb2 = R_60 * lattice_basis[1];
    xvector<double> lb3 = R_60 * lattice_basis[2];
    xmatrix<double> lb = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
    xmatrix<double> lb_r60_inv = aurostd::inverse(xvec2xmat(lb1, lb2, lb3));

    for (uint i = 0; i < atoms.size(); i++) {
      //DX20190905 [OBSOLETE-no more mod_one_xvec] atoms[i].fpos = SYM::mod_one_xvec(atoms[i].fpos * lb * lb_r60_inv);
      atoms[i].fpos = atoms[i].fpos * lb * lb_r60_inv; //DX20190905 
      BringInCellInPlace(atoms[i].fpos); //DX20190905 
    }
    lattice_basis[0] = lb1;
    lattice_basis[1] = lb2;
    lattice_basis[2] = lb3;
    lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
    return true;
  }
} // namespace SYM

// ===== deque<deque<_atom> > input ===== //
namespace SYM {
  bool transformToObverse(xmatrix<double>& lattice, deque<deque<_atom> >& equivalent_atoms) {
    vector<xvector<double> > lattice_basis;
    lattice_basis.push_back(lattice(1));
    lattice_basis.push_back(lattice(2));
    lattice_basis.push_back(lattice(3));
    xvector<double> tmpvec = lattice_basis[0];
    lattice_basis[0] = lattice_basis[1];
    lattice_basis[1] = tmpvec;
    for (uint i = 0; i < equivalent_atoms.size(); i++) {
      swap_columns(equivalent_atoms[i], 1, 2);
    }
    lattice = xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
    return true;
  }
} //namespace SYM

//DX20210611 [OBSOLETE - integrated into structure printing function in aflow_xatom] // **************************************************************************************************************************************************
//DX20210611 [OBSOLETE - integrated into structure printing function in aflow_xatom] // Prints Wyccar from xstructure
//DX20210611 [OBSOLETE - integrated into structure printing function in aflow_xatom] // **************************************************************************************************************************************************
//DX20210611 [OBSOLETE - integrated into structure printing function in aflow_xatom] namespace SYM {
//DX20210611 [OBSOLETE - integrated into structure printing function in aflow_xatom]   void printWyccar(ofstream& FileMESSAGE, xstructure& str, const bool& osswrite, ostream& oss) {
//DX20210611 [OBSOLETE - integrated into structure printing function in aflow_xatom]     ostringstream aus;
//DX20210611 [OBSOLETE - integrated into structure printing function in aflow_xatom]     for (uint i = 0; i < str.wyccar_ITC.size(); i++) {
//DX20210611 [OBSOLETE - integrated into structure printing function in aflow_xatom]       aus << str.wyccar_ITC[i];
//DX20210611 [OBSOLETE - integrated into structure printing function in aflow_xatom]       if(i < str.wyccar_ITC.size() - 1) {
//DX20210611 [OBSOLETE - integrated into structure printing function in aflow_xatom]         aus << endl;
//DX20210611 [OBSOLETE - integrated into structure printing function in aflow_xatom]       }
//DX20210611 [OBSOLETE - integrated into structure printing function in aflow_xatom]     }
//DX20210611 [OBSOLETE - integrated into structure printing function in aflow_xatom]     aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, osswrite, oss);
//DX20210611 [OBSOLETE - integrated into structure printing function in aflow_xatom]   }
//DX20210611 [OBSOLETE - integrated into structure printing function in aflow_xatom] } // namespace SYM

// **************************************************************************************************************************************************
// Primitive Cell Routine (Created by RHT, slightly different from SC version)
// **************************************************************************************************************************************************
uint xstructure::GetPrimitiveCell(void) {
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  if(LDEBUG) {
    cerr << "xstructure::GetPrimitiveCell(): Original structure (before primitivization): " << endl;
    cerr << (*this) << endl;
  }
  (*this).ReScale(1.0);//DX20180526
  xmatrix<double> prim_lattice = (*this).lattice;
  double prim_det = aurostd::det(prim_lattice);

  vector<xvector<double> > lattice_basis;
  lattice_basis.push_back(SYM::extract_row(prim_lattice, 1));
  lattice_basis.push_back(SYM::extract_row(prim_lattice, 2));
  lattice_basis.push_back(SYM::extract_row(prim_lattice, 3));
  xmatrix<double> lattice_basis_xmat = SYM::xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);

  xmatrix<double> f2c = trasp(lattice_basis_xmat);
  xmatrix<double> c2f = inverse(trasp(lattice_basis_xmat));
  //DX20210315 [OBSOLETE - now using transformation method for cell conversion] bool skew = SYM::isLatticeSkewed(lattice_basis_xmat, (*this).dist_nn_min, (*this).sym_eps); //DX20190215 - _SYM_TOL to (*this).sym_eps
  vector<xvector<double> > diff_vectors;

  // Determine the Greatest Common Denominator between the atom types
  int atoms_GCD = 0;
  deque<deque<_atom> > orig_split_atom_types;
  SYM::getAtomGCD((*this).atoms, orig_split_atom_types, atoms_GCD);
  (*this).BringInCell();

  int atom_count = (*this).atoms.size();
  deque<_atom> atomic_basis_;
  //DX20210129 [OBSOLETE] deque<int> basistypes;

  // Copy deque<_atoms> atoms (IN XSTRUCTURE) to atomic_basis_ to manipulate free of xstructure
  for (int i = 0; i < atom_count; i++) {
    _atom tmp_atom = (*this).atoms[i]; //DX20191011 - initialized here instead of explicitly setting individual attributes below
    //DX20191011 [OBSOLETE] tmp_atom.fpos = (*this).atoms[i].fpos;
    //DX20191011 [OBSOLETE] tmp_atom.cpos = (*this).atoms[i].cpos;
    //DX20191011 [OBSOLETE] tmp_atom.type = (*this).atoms[i].type;
    //DX20191011 [OBSOLETE] tmp_atom.spin = (*this).atoms[i].spin; //DX20170921 - magnetic sym
    //DX20191011 [OBSOLETE] tmp_atom.spin_is_given = (*this).atoms[i].spin_is_given; //DX20170921 - magnetic sym
    //DX20191011 [OBSOLETE] tmp_atom.noncoll_spin = (*this).atoms[i].noncoll_spin; //DX20171205 - magnetic sym (non-collinear)
    //DX20191011 [OBSOLETE] tmp_atom.noncoll_spin_is_given = (*this).atoms[i].noncoll_spin_is_given; //DX20171205 - magnetic sym (non-collinear)
    string tmpname = (*this).atoms[i].name;
    SYM::cleanupstring(tmpname);
    tmp_atom.name = tmpname;
    atomic_basis_.push_back(tmp_atom);
  }

  // ===== Expand the lattice and atomic basis ===== //
  // Expand crystal structure (expanded_basis) to check which atoms are inside reduced.
  // You must do this because reducing the cell often changes the orientation (not simply shrinking)
  // so atoms may be added that were not in the original representation.
  deque<_atom> expanded_basis;
  vector<xvector<double> > big_expanded;  //just lattice points
  int q = 3;
  big_expanded = SYM::expand_lattice(q, q, q, prim_lattice);
  expanded_basis = SYM::add_basis(big_expanded, prim_lattice, (*this));

  ////print(expanded_basis);
  //double prim_det=aurostd::det(prim_lattice);

  // Search only over smallest subset of atoms (no need to search over any other for periodic structures)
  deque<deque<_atom> > atoms_by_type = SYM::break_up_by_type(atomic_basis_);
  deque<int> numofatoms;
  uint smallest_group = atoms_by_type[0].size();
  uint index_for_smallest_group = 0;
  for (uint i = 1; i < atoms_by_type.size(); i++) {
    if(atoms_by_type[i].size() < smallest_group) {
      smallest_group = atoms_by_type[i].size();
      index_for_smallest_group = i;
    }
  }

  bool foundbasis = false;
  // If there is only one atom or one least frequently occuring atom, then the cell cannot be further reduced
  if(atom_count == 1 || atoms_by_type[index_for_smallest_group].size() == 1) {
    prim_lattice = SYM::xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
    foundbasis = true;
    sort(atomic_basis_.begin(), atomic_basis_.end(), sortAtomsTypes);
    //DX20210129 [OBSOLETE] numofatoms = ::GetNumEachType(atomic_basis_);
    //DX20210129 [OBSOLETE] basistypes = numofatoms;
  }

  deque<_atom> newbasis;
  lattice_basis_xmat = SYM::xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
  f2c = trasp(lattice_basis_xmat);
  c2f = inverse(trasp(lattice_basis_xmat));
  //DX20210315 [OBSOLETE - now using transformation method for cell conversion] skew = SYM::isLatticeSkewed(lattice_basis_xmat, (*this).dist_nn_min, (*this).sym_eps); //DX20190215 - _SYM_TOL_ to (*this).sym_eps
  vector<_atom> full_basis;  //basis that has been reduced modulo 1
  vector<xvector<double> > lattice_vector_candidates;

  // ===== Find all possible translation vectors in crystal ===== //
  if(foundbasis == false) {
    for (uint i = 0; i < atomic_basis_.size(); i++) {
      _atom tmp = atomic_basis_[i]; //DX20191011 - initialized instead of setting individually below
      //DX20190905 [OBSOLETE-no more mod_one_xvec] tmp.fpos = SYM::mod_one_xvec(atomic_basis_[i].fpos);
      tmp.fpos = ::BringInCell(atomic_basis_[i].fpos); //DX20190905 - new BringInCell function and go outside xstructure scope
      //DX20191011 [OBSOLETE] tmp.name = atomic_basis_[i].name;
      //DX20191011 [OBSOLETE] tmp.type = atomic_basis_[i].type;
      //DX20191011 [OBSOLETE] tmp.spin = atomic_basis_[i].spin; //DX20170921 - magnetic sym
      //DX20191011 [OBSOLETE] tmp.spin_is_given = atomic_basis_[i].spin_is_given; //DX20170921 - magnetic sym
      //DX20191011 [OBSOLETE] tmp.noncoll_spin = atomic_basis_[i].noncoll_spin; //DX20171205 - magnetic sym (non-collinear)
      //DX20191011 [OBSOLETE] tmp.noncoll_spin_is_given = atomic_basis_[i].noncoll_spin_is_given; //DX20171205 - magnetic sym (non-collinear)
      full_basis.push_back(tmp);
    }
    vector<xvector<double> > diff_vectors;
    vector<uint> atom_pos;
    for (uint i = 0; i < atoms_by_type[index_for_smallest_group].size(); i++) {
      //DX20190905 [OBSOLETE-no more mod_one_xvec] xvector<double> tmp = SYM::mod_one_xvec(atoms_by_type[index_for_smallest_group][i].fpos - atoms_by_type[index_for_smallest_group][0].fpos); //need to bring in cell, otherwise, checking larger cells not smaller
      xvector<double> tmp = atoms_by_type[index_for_smallest_group][i].fpos - atoms_by_type[index_for_smallest_group][0].fpos; //need to bring in cell, otherwise, checking larger cells not smaller //DX20190905
      BringInCellInPlace(tmp); //DX20190905
      diff_vectors.push_back(tmp);
      atom_pos.push_back(i);
    }

    // ---------------------------------------------------------------------------
    // translate by difference vectors and check if equivalent
    // (i.e., check if diff vector is a lattice vector) 
    bool is_frac = true;
    for (uint d = 1; d < diff_vectors.size(); d++) {  //use d=0 if we do double for loop
      if(isTranslationVector((*this), diff_vectors[d], (*this).sym_eps, is_frac)){
        lattice_vector_candidates.push_back(diff_vectors[d]);
      }
    }
    //add the original lattice vectors' translations to pool of potential translations:
    xvector<double> a1;
    a1(1) = 1;
    a1(2) = 0;
    a1(3) = 0;
    xvector<double> a2;
    a2(1) = 0;
    a2(2) = 1;
    a2(3) = 0;
    xvector<double> a3;
    a3(1) = 0;
    a3(2) = 0;
    a3(3) = 1;
    lattice_vector_candidates.push_back(a1);
    lattice_vector_candidates.push_back(a2);
    lattice_vector_candidates.push_back(a3);
    vector<double> moduli;
    uint num_of_candidates = lattice_vector_candidates.size();
    //Choose smallest L.I. potential lattice vectors. If there are only three, then these three (the original lattice three) must be the reduced basis.
    if(num_of_candidates == 3) {
      prim_lattice = SYM::xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);
      foundbasis = true;
      sort(atomic_basis_.begin(), atomic_basis_.end(), sortAtomsTypes);
      //DX20210129 [OBSOLETE] numofatoms = ::GetNumEachType(atomic_basis_);
      //DX20210129 [OBSOLETE] basistypes = numofatoms;

    } else {
      //Here because the original poscar is not primitive
      //Get vector of moduli of lattice vector candidates:
      xmatrix<double> nonsimp_reduced_basis;
      nonsimp_reduced_basis = prim_lattice;
      for (uint i = 0; i < num_of_candidates; i++) {
        moduli.push_back(aurostd::modulus(lattice_vector_candidates[i] * lattice_basis_xmat));
      }
      //Define max value inorder to find minimum in loop below:
      double sum_mods = 1e9;
      double amin = aurostd::min(aurostd::modulus(lattice_basis[0]), aurostd::min(aurostd::modulus(lattice_basis[1]), aurostd::modulus(lattice_basis[2])));
      amin = (*this).dist_nn_min;  //DX TEST
      double tol_vol = (*this).sym_eps * amin * amin * amin; //DX20190215 - _SYM_TOL_ to (*this).sym_eps

      //DX20180503 - FASTER MIN CART DISTANCE CALCULATOR - START
      //DX20180503 - only calculate multiplication once (time-saver)
      vector<xvector<double> > candidates_cpos;
      vector<double> candidates_mods;
      for(uint a=0;a<lattice_vector_candidates.size();a++){
        xvector<double> tmp = lattice_vector_candidates[a] * lattice_basis_xmat;
        candidates_cpos.push_back(tmp); candidates_mods.push_back(aurostd::modulus(tmp));
      }
      //DX20180503 - FASTER MIN CART DISTANCE CALCULATOR - END

      for (uint i = 0; i < candidates_cpos.size(); i++) {
        for (uint j = i+1; j < candidates_cpos.size(); j++) {
          for (uint k = j+1; k < candidates_cpos.size(); k++) {
            xmatrix<double> candidate_lattice = SYM::xvec2xmat(candidates_cpos[i], candidates_cpos[j], candidates_cpos[k]);
            double det = aurostd::determinant(candidate_lattice);
            if(aurostd::abs(det) > tol_vol){  //DX20180613 - can be left handed; Minkowski/Niggli will fix; necessary since we aren't checking permutations of vectors
              candidate_lattice=::MinkowskiBasisReduction(candidate_lattice);      // Minkowski first  //DX20180526
              candidate_lattice=::NiggliUnitCellForm(candidate_lattice);           // Niggli Second //DX20180526
              det = aurostd::determinant(candidate_lattice);
              //In order to find minimum, determinant should be positive,less than the first, and an integer multiple of less than first
              double atom_number_ratio = (double)smallest_group/aurostd::nint(prim_det/det); // atom ratio before and after reducing //DX20180723
              //DX [OBSOLETE] 20180723 - if(det > tol_vol && det < prim_det && aurostd::isinteger(prim_det / det, 0.05)) //DX20180614 - now it should be right-handed
              if(det > tol_vol && det < prim_det && aurostd::isinteger(prim_det / det, 0.05) && aurostd::isinteger(atom_number_ratio, 0.05) && (atom_number_ratio-1.0)>-_ZERO_TOL_) //DX20180614 - now it should be right-handed // check atom ratio after reducing //DX20180723
              { //CO20200106 - patching for auto-indenting
                //double sum_candidate_mods = candidates_mods[i] + candidates_mods[j] + candidates_mods[k];
                double sum_candidate_mods = aurostd::modulus(candidate_lattice(1)) + aurostd::modulus(candidate_lattice(2)) + aurostd::modulus(candidate_lattice(3)); //DX20180613 - need to recalculate mods, they may have changed during Minkowski/Niggli 
                //candidate_lattice=::MinkowskiBasisReduction(candidate_lattice);      // Minkowski first  //DX20180526
                //candidate_lattice=::NiggliUnitCellForm(candidate_lattice);           // Niggli Second //DX20180526
                if(sum_candidate_mods < sum_mods) {
                  sum_mods = sum_candidate_mods; 
                  nonsimp_reduced_basis = candidate_lattice;
                }
              }
            }
          }
        }
      }
      prim_lattice = nonsimp_reduced_basis;
      lattice_basis[0] = SYM::extract_row(prim_lattice, 1);
      lattice_basis[1] = SYM::extract_row(prim_lattice, 2);
      lattice_basis[2] = SYM::extract_row(prim_lattice, 3);

      lattice_basis_xmat = SYM::xvec2xmat(lattice_basis[0], lattice_basis[1], lattice_basis[2]);

      // ---------------------------------------------------------------------------
      //DX20210316 - used transformation method (more efficient)
      xmatrix<double> transformation_matrix = GetBasisTransformation((*this).lattice, lattice_basis_xmat);
      xmatrix<double> rotation_matrix = aurostd::eye<double>(3,3);
      xstructure prim = (*this);
      // try to primitivize, it may fail, so return original structure
      try{
        prim.TransformStructure(transformation_matrix, rotation_matrix);
        (*this) = prim;
      }
      catch(aurostd::xerror& re){}
      return 0;

      //DX20210315 - foldInCell method (too slow for large, skewed cells)
      //DX20210315 [OBSOLETE]xmatrix<double> f2c = trasp(lattice_basis_xmat);
      //DX20210315 [OBSOLETE]xmatrix<double> c2f = inverse(trasp(lattice_basis_xmat));
      //DX20210315 [OBSOLETE]bool skew = SYM::isLatticeSkewed(lattice_basis_xmat, (*this).dist_nn_min, (*this).sym_eps); //DX20190215 - _SYM_TOL_ to (*this).sym_eps
      //DX20210315 [OBSOLETE]xmatrix<double> PL_inv = aurostd::inverse(prim_lattice);
      //DX20210315 [OBSOLETE]double same_atom_tol = (*this).dist_nn_min - 0.1;  //min_dist itself will consider nn atom to be the same, needs to be slightly smaller.`
      //DX20210315 [OBSOLETE]//Now get atoms inside of new, reduced basis:
      //DX20210315 [OBSOLETE]//[CO20190520]newbasis = foldAtomsInCell(expanded_basis, c2f, f2c, skew, same_atom_tol); //CO20180409
      //DX20210315 [OBSOLETE]newbasis = ::foldAtomsInCell(expanded_basis, (*this).lattice, lattice_basis_xmat, skew, same_atom_tol,false); //CO20180409 //DX20190619 - false->do not check min dists (expensive)

      //DEBUG
      //cerr << "newbasis.size(): " << newbasis.size() << endl;
      //for(uint n=0;n<newbasis.size();n++){
      //  cerr << newbasis[n].fpos << endl;
      //}
      //DEBUG

      //DX20210315 [OBSOLETE]bool consistent_ratio = SYM::GCD_conventional_atomic_basis(newbasis, orig_split_atom_types, atoms_GCD);
      //DX20210315 [OBSOLETE]if(consistent_ratio == false) {
      //DX20210315 [OBSOLETE]  foundbasis = false;
      //DX20210315 [OBSOLETE]} else {
      //DX20210315 [OBSOLETE]  //Overwrite basis types (update it for primitive cell):
      //DX20210315 [OBSOLETE]  sort(newbasis.begin(), newbasis.end(), sortAtomsTypes);
      //DX20210315 [OBSOLETE]  //DX20210129 [OBSOLETE] numofatoms = ::GetNumEachType(newbasis);
      //DX20210315 [OBSOLETE]  //DX20210129 [OBSOLETE] basistypes = numofatoms;
      //DX20210315 [OBSOLETE]  atomic_basis_.clear();
      //DX20210315 [OBSOLETE]  atomic_basis_ = newbasis;
      //DX20210315 [OBSOLETE]  if(atomic_basis_.size() == 0) {
      //DX20210315 [OBSOLETE]    moduli.clear();
      //DX20210315 [OBSOLETE]  }
      //DX20210315 [OBSOLETE]  foundbasis = true;
      //DX20210315 [OBSOLETE]}
      //DX20210315 [OBSOLETE]cerr << "(*this): " << (*this) << endl;
    }
  }

  if(foundbasis == false) { return 1; }

  // ===== Put contents back in xstructure ===== //
  //DX20210129 [OBSOLETE] deque<int> numtypes;
  //DX20210129 [OBSOLETE] for (uint i = 0; i < basistypes.size(); i++) {
  //DX20210129 [OBSOLETE]   numtypes.push_back(basistypes[i]);
  //DX20210129 [OBSOLETE] }

  (*this).lattice = prim_lattice;
  if(newbasis.size() > 0) {  //IF A NEW BASIS EXISTS (CELL WAS NOT PRIMITIVE) THEN OVERWRITE a.atoms
    sort(newbasis.begin(), newbasis.end(), sortAtomsTypes);
    (*this).ReplaceAtoms(newbasis, false);
    //DX20210129 [OBSOLETE] // cerr << "USING NEW BASIS" << endl;
    //DX20210129 [OBSOLETE] // Remove old basis
    //DX20210129 [OBSOLETE] while ((*this).atoms.size() > 0) {
    //DX20210129 [OBSOLETE]   (*this).RemoveAtom((uint)0);
    //DX20210129 [OBSOLETE] }
    //DX20210129 [OBSOLETE] (*this).species.clear(); //DX
    //DX20210129 [OBSOLETE] for (uint i = 0; i < newbasis.size(); i++) {
    //DX20210129 [OBSOLETE]   _atom tmp_atom = newbasis[i]; //DX20191011 - initialized instead of setting individually below
    //DX20210129 [OBSOLETE]   //DX20191011 [OBSOLETE] tmp_atom.fpos = newbasis[i].fpos;
    //DX20210129 [OBSOLETE]   //DX20191011 [OBSOLETE] tmp_atom.cpos = newbasis[i].cpos;
    //DX20210129 [OBSOLETE]   //DX20191011 [OBSOLETE] tmp_atom.name = newbasis[i].name;
    //DX20210129 [OBSOLETE]   //DX20191011 [OBSOLETE] tmp_atom.type = newbasis[i].type;
    //DX20210129 [OBSOLETE]   //DX20191011 [OBSOLETE] tmp_atom.name_is_given = true;
    //DX20210129 [OBSOLETE]   //DX20191011 [OBSOLETE] tmp_atom.spin = newbasis[i].spin; //DX20170921 - magnetic sym
    //DX20210129 [OBSOLETE]   tmp_atom.spin_is_given = false; //DX20170921 - magnetic sym
    //DX20210129 [OBSOLETE]   if(aurostd::abs(tmp_atom.spin)>_ZERO_TOL_){
    //DX20210129 [OBSOLETE]     tmp_atom.spin_is_given = true; //DX20170921 - magnetic sym
    //DX20210129 [OBSOLETE]   }
    //DX20210129 [OBSOLETE]   tmp_atom.noncoll_spin = newbasis[i].noncoll_spin; //DX20171205 - magnetic sym (non-collinear)
    //DX20210129 [OBSOLETE]   tmp_atom.noncoll_spin_is_given = false; //DX20171205 - magnetic sym (non-collinear)
    //DX20210129 [OBSOLETE]   if(aurostd::abs(tmp_atom.noncoll_spin(1))>_ZERO_TOL_ || aurostd::abs(tmp_atom.noncoll_spin(2))>_ZERO_TOL_ || aurostd::abs(tmp_atom.noncoll_spin(3))>_ZERO_TOL_){
    //DX20210129 [OBSOLETE]     tmp_atom.noncoll_spin_is_given = true; //DX20171205 - magnetic sym (non-collinear) //DX20191108 - fixed typo, should be noncoll_spin_is_given not spin_is_given
    //DX20210129 [OBSOLETE]   }
    //DX20210129 [OBSOLETE]   (*this).AddAtom(tmp_atom,false);  //CO20230319 - add by type
    //DX20210129 [OBSOLETE] }
  }
  //DX20210129 [OBSOLETE] (*this).num_each_type = numtypes;
  //DX20210129 [moved earlier] (*this).lattice = prim_lattice;
  if(LDEBUG) {
    cerr << "xstructure::GetPrimitiveCell(): New primitivized cell:" << endl;
    cerr << (*this) << endl;
  }
  return 0;
}

// **************************************************************************************************************************************************
// Shift Origin (with L input)
// **************************************************************************************************************************************************
namespace SYM {
  vector<xvector<double> > expand_cell(xmatrix<double>& L) {
    // This expands the CELL (lattice vertices), so the expansion is in terms of unit cells.
    vector<xvector<double> > out;
    xvector<double> tmp;
    out.push_back(tmp);
    tmp(1) = (L(1, 1));
    tmp(2) = (L(1, 2));
    tmp(3) = (L(1, 3));
    out.push_back(tmp);
    tmp.clear();
    tmp(1) = (L(2, 1));
    tmp(2) = (L(2, 2));
    tmp(3) = (L(2, 3));
    out.push_back(tmp);
    tmp.clear();
    tmp(1) = (L(3, 1));
    tmp(2) = (L(3, 2));
    tmp(3) = (L(3, 3));
    out.push_back(tmp);
    tmp.clear();
    tmp(1) = (L(1, 1) + L(2, 1));
    tmp(2) = (L(1, 2) + L(2, 2));
    tmp(3) = (L(1, 3) + L(2, 3));
    out.push_back(tmp);
    tmp.clear();
    tmp(1) = (L(1, 1) + L(3, 1));
    tmp(2) = (L(1, 2) + L(3, 2));
    tmp(3) = (L(1, 3) + L(3, 3));
    out.push_back(tmp);
    tmp.clear();
    tmp(1) = (L(1, 1) + L(2, 1) + L(3, 1));
    tmp(2) = (L(1, 2) + L(2, 2) + L(3, 2));
    tmp(3) = (L(1, 3) + L(2, 3) + L(3, 3));
    out.push_back(tmp);
    tmp.clear();
    tmp(1) = (L(2, 1) + L(3, 1));
    tmp(2) = (L(2, 2) + L(3, 2));
    tmp(3) = (L(2, 3) + L(3, 3));
    out.push_back(tmp);
    tmp.clear();
    return out;
  }
} //namespace SYM

// **************************************************************************************************************************************************
// Add Basis
// **************************************************************************************************************************************************
namespace SYM {
  deque<_atom> add_basis(vector<xvector<double> >& expanded_lattice_points, xmatrix<double>& L, xstructure& xstr) {
    stringstream message;
    deque<_atom> out;
    if(xstr.atoms.size() == 0) {
      message << __AFLOW_FUNC__ << " atoms deque is empty! [dir=" << xstr.directory << "]";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_MISSING_);
    }
    for (uint i = 0; i < expanded_lattice_points.size(); i++) {
      //cerr << "lattice point: " << expanded_lattice_points[i] << endl; //cartesian
      for (uint j = 0; j < xstr.atoms.size(); j++) {  // (OBSOLETE) go from 1 because
        // the origin at 0 0 0 is taken care of by the expanded lattice
        _atom tmp = xstr.atoms[j]; //DX20191011 initialized instead of setting explictly (below)
        //tmp.print_RHT = true; //CO20190405 - get rid of this!
        tmp.cpos = expanded_lattice_points[i] + F2C(L, xstr.atoms[j].fpos);
        //DX20191011 [OBSOLETE] tmp.name = xstr.atoms[j].name;
        //DX20191011 [OBSOLETE] tmp.type = xstr.atoms[j].type;
        //DX20191011 [OBSOLETE] tmp.spin = xstr.atoms[j].spin; //DX20170921 - magnetic sym
        //DX20191011 [OBSOLETE] tmp.spin_is_given = xstr.atoms[j].spin_is_given; //DX20170921 - magnetic sym
        //DX20191011 [OBSOLETE] tmp.noncoll_spin = xstr.atoms[j].noncoll_spin; //DX20171205 - magnetic sym (non-collinear)
        //DX20191011 [OBSOLETE] tmp.noncoll_spin_is_given = xstr.atoms[j].noncoll_spin_is_given; //DX20171205 - magnetic sym (non-collinear)
        out.push_back(tmp);
      }
    }
    return out;
  }
}

// **************************************************************************************************************************************************
// Expand Lattice (added L as input)
// **************************************************************************************************************************************************
namespace SYM {
  vector<xvector<double> > expand_lattice_positive_only(int& a, int& b, int& c, xmatrix<double>& L) {
    // We are choosing a lattice point as origin and then expanding from that point.
    // The expanded points are in cartesian coordinates. This is a generic "stupid" expansion.
    vector<xvector<double> > out;
    xvector<double> tmp;
    for (int x = 0; x <= a; x++) {
      for (int y = 0; y <= b; y++) {
        for (int z = 0; z <= c; z++) {
          tmp(1) = L(1, 1) * x + L(2, 1) * y + L(3, 1) * z;
          tmp(2) = L(1, 2) * x + L(2, 2) * y + L(3, 2) * z;
          tmp(3) = L(1, 3) * x + L(2, 3) * y + L(3, 3) * z;
          out.push_back(tmp);
        }
      }
    }
    return out;
  }
} //namespace SYM

// **************************************************************************************************************************************************
// Expand Lattice (added L as input)
// **************************************************************************************************************************************************
namespace SYM {
  vector<xvector<double> > expand_lattice(int& a, int& b, int& c, xmatrix<double>& L) {
    // We are choosing a lattice point as origin and then expanding from that point.
    // The expanded points are in cartesian coordinates. This is a generic "stupid" expansion.
    vector<xvector<double> > out;
    xvector<double> tmp;
    for (int x = -a; x <= a; x++) {
      for (int y = -b; y <= b; y++) {
        for (int z = -c; z <= c; z++) {
          tmp(1) = L(1, 1) * x + L(2, 1) * y + L(3, 1) * z;
          tmp(2) = L(1, 2) * x + L(2, 2) * y + L(3, 2) * z;
          tmp(3) = L(1, 3) * x + L(2, 3) * y + L(3, 3) * z;
          out.push_back(tmp);
        }
      }
    }
    return out;
  }
} //namespace SYM

//**********************************************************************************************************************************************
// Conventional Cell Function -- Determines the conventional cell consistent with symmetry
//**********************************************************************************************************************************************
namespace SYM {
  xstructure ConventionalCell(xstructure& xstr, int& IT, int& cell_choice, bool& last_orientation, string& crystalsystem_prev,
      xstructure& CrystOut_prev, vector<xmatrix<double> >& candidate_lattice_vectors_prev,
      vector<char>& candidate_lattice_chars_prev, symfolder& checkops, SymmetryInformationITC& ITC_sym_info, bool& lattice_reformed,
      vector<int>& lattice_pgroups, vector<xmatrix<double> >& lattice_sym_mats,
      vector<xmatrix<double> >& crystal_sym_mats, bool& symmetryfound, bool foundspacegroup_global) {
    bool LDEBUG = (FALSE || XHOST.DEBUG);

    double sym_tol = xstr.sym_eps; //DX20190215 - get sym_tol

    // ========== Objects containing crystal stucture information ========== //
    xstructure xstr_out;
    //DX20210401 [OBSOLETE] xstr.GetPrimitiveCell();
    xstr.GetPrimitive(); //DX20210401
    xstructure prim = xstr; //DX20210316 - store primitive

    //vector<xmatrix<double> > lattice_vector_choices; //DX20220827
    xmatrix<double> L = xstr.lattice;
    deque<_atom> primitive_atomic_basis = xstr.atoms;
    for (uint i = 0; i < primitive_atomic_basis.size(); i++) {
      primitive_atomic_basis[i].cpos = F2C(xstr.lattice, primitive_atomic_basis[i].fpos);
    }
    //print(primitive_atomic_basis);
    xmatrix<double> LPRIM = xstr.lattice;
    xmatrix<double> f2c = trasp(LPRIM);
    xmatrix<double> c2f = inverse(trasp(LPRIM));
    bool skew = SYM::isLatticeSkewed(LPRIM, xstr.dist_nn_min, sym_tol); //DX20190215 - _SYM_TOL_ to sym_tol
    xmatrix<double> Linv = aurostd::inverse(L);


    bool atoms_overlapping = false;
    // === Check if atoms occupying the same space === //
    for (uint i = 0; i < primitive_atomic_basis.size(); i++) {
      for (uint j = i + 1; j < primitive_atomic_basis.size(); j++) {
        if(SYM::MapAtom(primitive_atomic_basis[i], primitive_atomic_basis[j], FALSE, LPRIM, f2c, skew, sym_tol)) { //DX20190215 - added sym_tol //DX20190619 - LPRIM and f2c as input
          if(LDEBUG) { cerr << "SYM::ConventionalCell: WARNING: Atoms occupying the same space (given tolerance " << sym_tol << "). Check structure [dir=" << xstr.directory << "]." << endl; } //DX20190215 - _SYM_TOL_ to sym_tol
          if(LDEBUG) { cerr << "SYM::ConventionalCell: atom #" << i << " symbol: " << primitive_atomic_basis[i].name << " coord: " << primitive_atomic_basis[i].fpos << endl; }
          if(LDEBUG) { cerr << "SYM::ConventionalCell: atom #" << j << " symbol: " << primitive_atomic_basis[j].name << " coord: " << primitive_atomic_basis[j].fpos << endl; }
          atoms_overlapping = true;
        }
      }
    }

    //========== Greatest common denominator (GCD) in the primitive atomic basis ==========//
    // This is used to check the consistency between the ratio of atoms when we find a conventional cell
    deque<deque<_atom> > prim_split_atom_types;
    int prim_GCD = 0;
    getAtomGCD(primitive_atomic_basis, prim_split_atom_types, prim_GCD);

    // ===== CRYSTAL SYSTEM AND LATTICE SYSTEM VARIABLES===== //
    string crystalsystem = "";
    string latticesystem = "";
    vector<xmatrix<double> > candidate_lattice_vectors;
    vector<char> candidate_lattice_chars;
    bool lattice_vectors_found = false;
    
    uint sym_count = 0;
    if(IT > 0) {
        crystalsystem = crystalsystem_prev;
        if((uint)(IT + sym_count) < candidate_lattice_vectors_prev.size()) {
          candidate_lattice_vectors = candidate_lattice_vectors_prev;
          candidate_lattice_chars = candidate_lattice_chars_prev;
          lattice_vectors_found = true;
        }
      }

    // ===== Cell/Lattice Expansion coefficients ===== //
    int b = 1;
    int q = 3;

    //// ===== Declare expanded lattice and expanded cell variables ===== //
    //vector<xvector<double> > expanded_lattice_1 = expand_lattice(b, b, b, L);
    //vector<xvector<double> > expanded_cell_1 = expand_cell(L);
    
    // ===== Declare expanded lattice and expanded cell variables ===== //
    vector<xvector<double> > expanded_lattice_1;
    vector<xvector<double> > expanded_cell_1;
    
    // ===== Declare symmetry operation variables ===== //
    vector<Screw> twofold_ops_vec;
    vector<Screw> rot_ops_vec;
    vector<Glide> mirror_ops_vec;
    vector<xvector<double> > mirror_lattice_vectors;
    vector<xvector<double> > twofold_lattice_vectors;
    vector<xvector<double> > rot_lattice_vectors;

    if(!lattice_vectors_found) {
      expanded_lattice_1 = expand_lattice(b, b, b, L);
      expanded_cell_1 = expand_cell(L);


    //double radius = RadiusSphereLattice(L);
    double radius = 2.0 * RadiusSphereLattice(L);

    // ========== Determine symmetry operations consistent with the lattice ==========//
    //DX20170926 [OBSOLETE] mirror_ops_vec = mirror_operations(expanded_lattice_1, expanded_cell_1, L, Linv, mirror_lattice_vectors, radius, skew);
    //DX20170926 [OBSOLETE] rot_ops_vec = triplet_operations(expanded_lattice_1, expanded_cell_1, L, Linv, rot_lattice_vectors, radius, skew);
    //DX20170926 [OBSOLETE] twofold_ops_vec = twofold_operations(expanded_lattice_1, expanded_cell_1, L, Linv, twofold_lattice_vectors, radius, skew);
    mirror_ops_vec = mirror_operations(expanded_lattice_1, expanded_cell_1, L, Linv, mirror_lattice_vectors, radius, skew, sym_tol); //DX20190215 - added sym_tol
    for(uint m=0;m<mirror_ops_vec.size();m++){
      Screw candidate_rotation;
      xvector<double> axis_direction = mirror_ops_vec[m].return_direction();
      xvector<double> point_on_axis; //Null by default
      candidate_rotation.get_screw_direct(axis_direction, point_on_axis, 2);
      twofold_ops_vec.push_back(candidate_rotation);
      twofold_lattice_vectors.push_back(mirror_lattice_vectors[m]);
    }
    if(mirror_ops_vec.size()>2){
      rot_ops_vec = triplet_operations(expanded_lattice_1, expanded_cell_1, L, Linv, rot_lattice_vectors, radius, skew, sym_tol); //DX20190215 - added sym_tol
    }
    if(LDEBUG) {
      cerr << "SYM::ConventionalCell: Number of mirror operations: " << mirror_ops_vec.size() << " [dir=" << xstr.directory << "]." << endl;
      for(uint m=0;m<mirror_ops_vec.size();m++){
        cerr << "SYM::ConventionalCell: mirror #" << m+1 << " : " << mirror_ops_vec[m] << endl;
      }
      cerr << "SYM::ConventionalCell: Number of twofold operations: " << twofold_ops_vec.size() << " [dir=" << xstr.directory << "]." << endl;
      for(uint m=0;m<twofold_ops_vec.size();m++){
        cerr << "SYM::ConventionalCell: twofold operator #" << m+1 << " : " << twofold_ops_vec[m] << endl;
      }
      cerr << "SYM::ConventionalCell: Number of n-fold (n>2) operations: " << rot_ops_vec.size() << " [dir=" << xstr.directory << "]." << endl;
      for(uint m=0;m<rot_ops_vec.size();m++){
        cerr << "SYM::ConventionalCell: n-fold (n>2) operator #" << m+1 << " : " << rot_ops_vec[m] << endl;
      }
    }
    // Note: It will not always work to check the lattice symmetry to get to correct crystal system.
    // For example, you may have a trigonal structure with a cubic lattice. This would result in
    // having nine mirror planes in the lattice, leading to the incorrect assumption that the lattice is
    // not hexagonal. You cannot go from a cubic lattice to the hexagonal crystal system by the addition
    // of internal atoms because a subgroup cannot be formed with a sixfold symmetry. However, threefold
    // symmetry axes exists in the original lattice allowing for the transition to the trigonal crystal system.
    }

    // ========== Stores the lattice vector possibilities (different orientations, etc.) for a given lattice type ========== //
    //DX20220914 vector<xmatrix<double> > candidate_lattice_vectors;
    //DX20220914 vector<char> candidate_lattice_chars;

    // ===== Variables to keep track of status of symmetry search ===== //
    //DX20220914 [MOVED] uint sym_count = 0;
    bool reform = true;

    // ****************************************************************************************************************************************** //
    //                                                            Conventional Cell Loop
    // ****************************************************************************************************************************************** //
    // This section determines the conventional cell. Note: Lattice symmetry >= crystal symmetry.  The first loop determines the conventional cell
    // based on the lattice symmetry.  If the crystal symmetry is inconsistent with the conventional cell chosen by the lattice symmetry (i.e. the
    // crystal has lower symmetry), then the conventional cell is reformed.

    while (reform == true) {
      // ===== If not the first iteration, use the symmetry elements from the first iteration ===== //
      if(IT >= 1) {
        xstr = CrystOut_prev;
        //DX20210401 [OBSOLETE] xstr.GetPrimitiveCell();
        xstr.GetPrimitive(); //DX20210401
        prim = xstr; //DX20210316 - store primitive
        L = xstr.lattice;
        Linv = aurostd::inverse(L);
      }
      xstr_out = xstr;
      CrystOut_prev = xstr;

      // ===== Declare expanded lattice and expanded cell variables ===== //
      //DX20220914 [OBSOLETE] vector<xvector<double> > expanded_lattice = expand_lattice(b, b, b, L);
      //DX20220914 [OBSOLETE] vector<xvector<double> > expanded_cell = expand_cell(L);
      vector<xvector<double> > big_expanded = expand_lattice(q, q, q, L);
      deque<_atom> expanded_basis = add_basis(big_expanded, L, xstr);

      // ===== Conventional cell variables  ===== //
      vector<xvector<double> > bravais_basis;

      // ========== Determine the number of particular symmetry elements to screen/determine the possible conventional cell ========== //
      int mcount = mirror_ops_vec.size();
      int fourcount = 0;
      int three_six_count = 0;
      for (uint i = 0; i < rot_ops_vec.size(); i++) {
        if(aurostd::abs(rot_ops_vec[i].return_order() - 4.0) < _ZERO_TOL_) {
          fourcount++;
        }
        if(aurostd::abs(rot_ops_vec[i].return_order() - 3.0) < _ZERO_TOL_ || aurostd::abs(rot_ops_vec[i].return_order() - 6.0) < _ZERO_TOL_) {
          three_six_count++;
        }
      }

      // ===== If atoms overlapping, then set to TRICLINIC (i.e. spacegroup=1) ===== //
      if(atoms_overlapping == true) {
        mcount = 0;
        fourcount = 0;
        three_six_count = 0;
      }

      //When the lattice symetry leads to an incorrect crystalsystem, reformation will occur.
      //When this happens, mcount should be set to zero so that only the variables
      //crystalsystem and latticesystem are used in logic.

      // ===== Determine if symmetry is found ===== //
      symmetryfound = false;
      //DX20220914 bool lattice_vectors_found = false;
      //When on second iteration, use only crystal system to designate cell (having both mcount and crystal system can lead to going into hexagonal
      //and orthorhombic sections which adds redundant atoms) 20130821
      if(IT == 0) {
        if(sym_count > 0) {
          mcount = -1;
          if(sym_count < candidate_lattice_vectors_prev.size()) {
            //cerr << "WARNING: Should never be here..." << endl;
            candidate_lattice_vectors = candidate_lattice_vectors_prev;
            candidate_lattice_chars = candidate_lattice_chars_prev;
            lattice_vectors_found = true;
          }
        }
      }
      if(IT > 0) {
        crystalsystem = crystalsystem_prev;
        if((uint)(IT + sym_count) < candidate_lattice_vectors_prev.size()) {
          candidate_lattice_vectors = candidate_lattice_vectors_prev;
          candidate_lattice_chars = candidate_lattice_chars_prev;
          lattice_vectors_found = true;
        }
      }

      bool continue_ortho = false;
      // === Check orthogonality defect ===== //
      uint totalcount = mcount;
      int totalnum = CombinationNr(3, totalcount);
      vector<int> combiset;
      combiset = AllCombination41(3, totalcount, 1);
      for (int i = 0; i < totalnum; i++) {
        if(aurostd::abs(orthogonality_defect(xvec2xmat(mirror_ops_vec[combiset[0]].return_direction(), mirror_ops_vec[combiset[1]].return_direction(), mirror_ops_vec[combiset[2]].return_direction())) - 1.0) < sym_tol) { //DX20190215 - _SYM_TOL_ to sym_tol
          continue_ortho = true;
          break;
        }
        combiset = AllCombination42(3, totalcount, combiset);
      }
      // ********************************************************************************************************************************
      // Symmetry/Conventional Cell Routine
      // ********************************************************************************************************************************
      //    1) Use the number of mirror/fourfold symmetry elements to determine the possible conventional cell
      //    2) Find the lattice vectors
      //    3) Find number of lattice points in cell (lattice centering)
      //	  4) Populate the new conventional cell with the conventional atomic basis (check if ratio is consistent)
      //    5) Check if the lattice symmetry is consistent with the crystal symmetry
      //	     ==> IF NOT
      //	          a) Check all possible orientations, if multiple are possible, and recheck crystal symmetry
      //		  b) If still not commensurate, reform the lattice to be consistent with the crystal symmetry and repeat 1-5 using
      //		     the crystalsystem variable as the indicator of the correct lattice type.
      //       ==> IF CONSISTENT
      //	          a) Store crystallographic information and proceed to check generators/origin shift/Wyckoff positions
      // ********************************************************************************************************************************

      // DEBUG
      // cerr << "SYM::ConventionalCell: mcount: " << mcount << endl;
      // cerr << "SYM::ConventionalCell: twofold: " << twofold_ops_vec.size() << endl;
      // cerr << "SYM::ConventionalCell: fourcount: " << fourcount << endl;
      // cerr << "SYM::ConventionalCell: crystalsystem: " << crystalsystem << endl;
      // DEBUG
      if(lattice_vectors_found == false) {
        // *****************************************************************************************************************************
        // Find the lattice vectors.
        // *****************************************************************************************************************************

        // =============== CUBIC =============== //
        if((mcount == 9 && fourcount == 3) || crystalsystem == "CUBIC") {
          if(LDEBUG) { cerr << "SYM::ConventionalCell: CUBIC [dir=" << xstr.directory << "]." << endl; }
          symmetryfound = true;
          crystalsystem = "CUBIC";

          if(!findCubicLattice(rot_lattice_vectors, rot_ops_vec, candidate_lattice_vectors, candidate_lattice_chars, sym_tol)) { //DX20190215 - added sym_tol
            symmetryfound = false;
            return xstr_out;
          }
        }

        // =============== HEXAGONAL =============== //
        else if(((mcount == 7 || mcount == 8) && three_six_count != 0)|| crystalsystem == "TRIGONAL" || crystalsystem == "HEXAGONAL") {
          symmetryfound = true;
          crystalsystem = "HEXAGONAL";
          if(LDEBUG) { cerr << "SYM::ConventionalCell: HEX [dir=" << xstr.directory << "]." << endl; }
          if(!findTrigonalLattice(rot_lattice_vectors, rot_ops_vec, big_expanded, candidate_lattice_vectors, candidate_lattice_chars, sym_tol)) { //DX20190215 - added sym_tol
            symmetryfound = false;
            return xstr_out;
          }
          //DX20180808 - check rhl - START
          // some ambiguity whether hex or rhl if using crystalsystem as conditional; 
          // if using cell choice 2, the cell is always hex 
          // if using cell choice 1, the cell is hex for hexagonal space groups and rhl for rhombohedral space groups 
          // need to check number of lattice points in the cell
          if(cell_choice==SG_SETTING_1 || cell_choice==SG_SETTING_ANRL){ //DX20190131 - add ANRL setting
            f2c = trasp(candidate_lattice_vectors[0]);
            c2f = inverse(trasp(candidate_lattice_vectors[0]));
            skew = SYM::isLatticeSkewed(candidate_lattice_vectors[0], xstr_out.dist_nn_min, sym_tol); //DX20190215 - _SYM_TOL_ to sym_tol
            // ===== Find number of lattice points in cell (lattice centering) ===== //
            vector<xvector<double> > bravais_basis;
            int bravais_count = 0;
            determineLatticeCentering(bravais_basis, bravais_count, c2f, f2c, skew, big_expanded, crystalsystem, candidate_lattice_chars, sym_tol); //DX20190215 - added sym_tol
            if(bravais_count == 3){ // if three lattice points in the cell, then the system is rhl
              //use rhombohedral setting
              if(LDEBUG) { cerr << "SYM::ConventionalCell: Crystal system is trigonal with three lattice points and setting=1, using rhombohedral setting [dir=" << xstr.directory << "]." << endl; }
              if(!findRhombohedralSetting(big_expanded, candidate_lattice_vectors, candidate_lattice_chars, sym_tol)) { //DX20190215 - added sym_tol
                symmetryfound = false;
                return xstr_out;
              }
            } 
          }
          //DX20180808 - check rhl - END
        }

        // =============== TETRAGONAL =============== //
        //else if(((mcount == 4 || mcount == 5 || mcount == 6 || mcount == 7 || mcount == 8) && fourcount != 0) || crystalsystem == "TETRAGONAL")
        //DX20210126 [OBSOLETE - go up to mcount==9] else if(((mcount >= 3 && mcount <= 8) && fourcount != 0) || crystalsystem == "TETRAGONAL")
        else if(((mcount >= 3 && mcount <= 9) && fourcount != 0) || crystalsystem == "TETRAGONAL")
        { //CO20200106 - patching for auto-indenting
          symmetryfound = true;
          crystalsystem = "TETRAGONAL";
          if(LDEBUG) { cerr << "SYM::ConventionalCell: TET [dir=" << xstr.directory << "]." << endl; }

          if(!findTetragonalLattice(rot_lattice_vectors, twofold_lattice_vectors, rot_ops_vec, twofold_ops_vec, big_expanded,
                candidate_lattice_vectors, candidate_lattice_chars, sym_tol)) { //DX20190215 - added sym_tol
            symmetryfound = false;
            return xstr_out;
          }
        }

        // =============== ORTHORHOMBIC/RHOMBOHEDRAL =============== //
        else if(((mcount >= 3 && mcount <= 5) && continue_ortho && fourcount == 0) || crystalsystem == "ORTHORHOMBIC") {
          symmetryfound = true;
          if(LDEBUG) { cerr << "SYM::ConventionalCell: ORTHO/RHOMB [dir=" << xstr.directory << "]." << endl; }
          // ===== Check between orthorhombic and rhombohedral ===== //
          // To discern between rhombohedral and orthorhombic lattices by mirror planes,
          // check the orthogonality defect. (Rhombohedral mirror planes are not mutually orthogonal).
          //bool continue_ortho = false;
          //// === Check orthogonality defect ===== //
          //uint totalcount = mcount;
          //int totalnum = CombinationNr(3, totalcount);
          //vector<int> combiset;
          //combiset = AllCombination41(3, totalcount, 1);
          //for (int i = 0; i < totalnum; i++) {
          //  if(aurostd::abs(orthogonality_defect(xvec2xmat(mirror_ops_vec[combiset[0]].return_direction(), mirror_ops_vec[combiset[1]].return_direction(), mirror_ops_vec[combiset[2]].return_direction())) - 1.0) < sym_tol) { //DX20190215 - _SYM_TOL_ to sym_tol
          //    continue_ortho = true;
          //    break;
          //  }
          //  combiset = AllCombination42(3, totalcount, combiset);
          //}
          //if(crystalsystem == "ORTHORHOMBIC") {  //IF MCOUNT == 30 CRYSTAL CLASS IS ALREADY KNOWN
          //  continue_ortho = true;
          //}

          // =============== ORTHORHOMBIC =============== //
          //  if(continue_ortho == true) {  //[CO20200106 - close bracket for indenting]}
          if(LDEBUG) { cerr << "SYM::ConventionalCell: ORTHO true: [dir=" << xstr.directory << "]." << endl; }
          crystalsystem = "ORTHORHOMBIC";

          if(!findOrthorhombicLattice(twofold_lattice_vectors, mirror_lattice_vectors, candidate_lattice_vectors, candidate_lattice_chars, sym_tol)) { //DX20190215 - added sym_tol
            symmetryfound = false;
            return xstr_out;
          }
        }
        // =============== RHOMBOHEDRAL =============== //
        //  else {  //[CO20200106 - close bracket for indenting]}
        else if(((mcount >= 3 && mcount <= 6) && !continue_ortho && three_six_count!=0 && fourcount == 0)) { //DX20210107 - added mcount==6
          symmetryfound = true;
          // On the second iteration, this section is now equivalent to the hexagonal section above.
          // Need a slightly larger expansion to get lattice basis in HEX setting
          // Using 'r' to denote rhombohedral in primitive lattice cell and 'R' to denote rhombodehral with hexagonal axes
          crystalsystem = "TRIGONAL";
          if(LDEBUG) { cerr << "SYM::ConventionalCell: RHOMB [dir=" << xstr.directory << "]." << endl; }

          if(!findRhombohedralLattice(rot_lattice_vectors, rot_ops_vec, big_expanded, candidate_lattice_vectors, candidate_lattice_chars, sym_tol)) { //DX20190215 - added sym_tol
            symmetryfound = false;
            return xstr_out;
          }
          //DX20180808 - check rhl - START
          // need to check number of lattice points in the cell; it may be a hexagonal space group and not rhombohedral
          if(cell_choice==SG_SETTING_1 || cell_choice==SG_SETTING_ANRL){ //DX20190131 - add ANRL setting
            f2c = trasp(candidate_lattice_vectors[0]);
            c2f = inverse(trasp(candidate_lattice_vectors[0]));
            skew = SYM::isLatticeSkewed(candidate_lattice_vectors[0], xstr_out.dist_nn_min, sym_tol); //DX20190215 - _SYM_TOL_ to sym_tol
            // ===== Find number of lattice points in cell (lattice centering) ===== //
            vector<xvector<double> > bravais_basis;
            int bravais_count = 0;
            determineLatticeCentering(bravais_basis, bravais_count, c2f, f2c, skew, big_expanded, crystalsystem, candidate_lattice_chars, sym_tol); //DX20190215 - added sym_tol
            if(bravais_count == 3){ // if three lattice points in the cell, then the system is rhl
              //use rhombohedral setting
              if(!findRhombohedralSetting(big_expanded, candidate_lattice_vectors, candidate_lattice_chars, sym_tol)) { //DX20190215 - added sym_tol
                symmetryfound = false;
                return xstr_out;
              }
            }
          }
          //DX20180808 - check rhl - END
        }

        // =============== MONOCLINIC =============== //
        else if(((mcount >= 1 && mcount <= 3) && !continue_ortho && three_six_count==0) || crystalsystem == "MONOCLINIC")
          //else if(((mcount == 1 || mcount == 2 || mcount == 3) && three_six_count==0) || crystalsystem == "MONOCLINIC") //DX20180613 - can be orthogonal and still be monoclinic 
        { //CO20200106 - patching for auto-indenting
          symmetryfound = true;
          crystalsystem = "MONOCLINIC";
          if(LDEBUG) { cerr << "SYM::ConventionalCell: MONO [dir=" << xstr.directory << "]." << endl; }

          if(!findMonoclinicLattice(mirror_lattice_vectors, twofold_lattice_vectors, big_expanded,
                candidate_lattice_vectors, candidate_lattice_chars, cell_choice, sym_tol)) { //DX20180816 - added cell choice //DX20190215 - added sym_tol
            symmetryfound = false;
            return xstr_out;
          }
        }

        // =============== TRICLINIC =============== //
        else if(mcount == 0 || crystalsystem == "TRICLINIC") {
          symmetryfound = true;
          crystalsystem = "TRICLINIC";
          if(LDEBUG) { cerr << "SYM::ConventionalCell: TRI [dir=" << xstr.directory << "]." << endl; }

          if(!findTriclinicLattice(L, candidate_lattice_vectors,
                candidate_lattice_chars)) {
            symmetryfound = false;
            return xstr_out;
          }
        }
        else {
          if(LDEBUG) { cerr << "SYM::ConventionalCell: Could not identify crystal system, i.e., symmetry operation counts do not match with any of the crystal systems [dir=" << xstr.directory << "]." << endl;}
        }
      }  //end of lattice_vectors_found if statement
      // *****************************************************************************************************************************
      // END: Find the lattice vectors.
      // *****************************************************************************************************************************

      // *****************************************************************************************************************************
      // Search for conventional cell orientation.
      // *****************************************************************************************************************************
      uint start_index = 0;
      if(candidate_lattice_vectors.size() > (uint)IT) {
        start_index = IT;
      } else {
        start_index = 0;
      }
      //lattice_vector_choices.clear(); //DX20220827
      for (uint i = start_index; i < candidate_lattice_vectors.size(); i++) {
        xmatrix<double> CL = candidate_lattice_vectors[i];
        if(LDEBUG) {
          cerr << "SYM::ConventionalCell: Possible conventional lattice vectors (choice " << (i+1) << " of " << candidate_lattice_vectors.size() << "): [dir=" << xstr.directory << "]." << endl;
          cerr << CL << endl;
        }
        f2c = trasp(CL);
        c2f = inverse(trasp(CL));
        skew = SYM::isLatticeSkewed(CL, xstr_out.dist_nn_min, sym_tol); //DX20190215 - _SYM_TOL_ to sym_tol
        // ===== Find number of lattice points in cell (lattice centering) ===== //
        vector<xvector<double> > bravais_basis;
        int bravais_count = 0;
        determineLatticeCentering(bravais_basis, bravais_count, c2f, f2c, skew, big_expanded, crystalsystem, candidate_lattice_chars, sym_tol); //DX20190215 - added sym_tol
        char lattice_char = candidate_lattice_chars[i];

        // *****************************************************************************************************************************
        // Transform into conventional cell //DX20210316 - transformation method based on primitive cell instead of foldAtomsInCell 
        // *****************************************************************************************************************************
        xmatrix<double> transformation_matrix = GetBasisTransformation(prim.lattice, CL);
        xmatrix<double> rotation_matrix = aurostd::eye<double>(3,3);
        // try to convert to conventional cell, if it fails, change tolerance
        try{
          xstr_out = TransformStructure(prim, transformation_matrix, rotation_matrix);
          xstr_out.lattice_label_ITC = lattice_char;
        }
        catch(aurostd::xerror& re){
          symmetryfound = false;
          return xstr_out;
        }

        //DX20210315 [OBSOLETE]// *****************************************************************************************************************************
        //DX20210315 [OBSOLETE]// Fill in cell with conventional atomic basis
        //DX20210315 [OBSOLETE]// *****************************************************************************************************************************
        //DX20210315 [OBSOLETE]// ========== Conventional Atomic Basis ========== //
        //DX20210315 [OBSOLETE]//Now add those basis atoms that fit in conventional cell:
        //DX20210315 [OBSOLETE]//[CO20190520]deque<_atom> conventional_basis_atoms = foldAtomsInCell(expanded_basis, c2f, f2c, skew, sym_tol); //DX20190215 - _SYM_TOL_ to sym_tol
        //DX20210315 [OBSOLETE]deque<_atom> conventional_basis_atoms = foldAtomsInCell(expanded_basis, xstr.lattice, CL, skew, sym_tol, false); //DX20190215 - _SYM_TOL_ to sym_tol //DX20190619 - false->do not check min dists (expensive)

        //DX20210315 [OBSOLETE]bool consistent_ratio = true;  // Change tolerances if the ratio between atomic basis atoms is inconsistent with the primitive cell.
        //DX20210315 [OBSOLETE]// ===== Check if ratio of atoms is consistent with primitive atomic basis ===== //
        //DX20210315 [OBSOLETE]// Note: If the original POSCAR has atoms overlapping, skip this part (SG will be equal to 1)
        //DX20210315 [OBSOLETE]if(atoms_overlapping == false) {
        //DX20210315 [OBSOLETE]  consistent_ratio = GCD_conventional_atomic_basis(conventional_basis_atoms, prim_split_atom_types, prim_GCD);
        //DX20210315 [OBSOLETE]}
        //DX20210315 [OBSOLETE]if(consistent_ratio == false) {
        //DX20210315 [OBSOLETE]  if(LDEBUG) { cerr << "SYM::ConventionalCell: WARNING: Conventional cell atomic basis ratio is inconsistent with the primitive cell! Recalculating. [dir=" << xstr.directory << "]." << endl; }
        //DX20210315 [OBSOLETE]  symmetryfound = false;
        //DX20210315 [OBSOLETE]  continue;
        //DX20210315 [OBSOLETE]}

        //DX20210315 [OBSOLETE]// ========== Rearrange atoms/Store ========== //
        //DX20210315 [OBSOLETE]sort(conventional_basis_atoms.begin(), conventional_basis_atoms.end(), sortAtomsTypes);

        //DX20210315 [OBSOLETE]xstr_out.lattice_label_ITC = lattice_char;
        //DX20210315 [OBSOLETE]xstr_out.lattice = CL;

        //DX20210315 [OBSOLETE]vector<string> in_names;
        //DX20210315 [OBSOLETE]for (uint c = 0; c < conventional_basis_atoms.size(); c++) {
        //DX20210315 [OBSOLETE]  conventional_basis_atoms[c].fpos = C2F(xstr_out.lattice, conventional_basis_atoms[c].cpos);
        //DX20210315 [OBSOLETE]  in_names.push_back(conventional_basis_atoms[c].name);
        //DX20210315 [OBSOLETE]}

        //DX20210315 [OBSOLETE]xstr_out.ReplaceAtoms(conventional_basis_atoms, false); //DX20210129
        //DX20210315 [OBSOLETE]//DX20210129 [OBSOLETE] //Remove old atomic basis
        //DX20210315 [OBSOLETE]//DX20210129 [OBSOLETE] while (xstr_out.atoms.size() > 0) {
        //DX20210315 [OBSOLETE]//DX20210129 [OBSOLETE]   xstr_out.RemoveAtom((uint)0);
        //DX20210315 [OBSOLETE]//DX20210129 [OBSOLETE] }
        //DX20210315 [OBSOLETE]//DX20210129 [OBSOLETE] xstr_out.species.clear();
        //DX20210315 [OBSOLETE]//DX20210129 [OBSOLETE] // Add new atomic basis
        //DX20210315 [OBSOLETE]//DX20210129 [OBSOLETE] for (uint c = 0; c < conventional_basis_atoms.size(); c++) {
        //DX20210315 [OBSOLETE]//DX20210129 [OBSOLETE]   xstr_out.AddAtom(conventional_basis_atoms[c],false); //CO20230319 - add by type
        //DX20210315 [OBSOLETE]//DX20210129 [OBSOLETE] }

        //DX20210315 [OBSOLETE]//DX20190410 START
        //DX20210315 [OBSOLETE]// Check if AddAtom removed atoms
        //DX20210315 [OBSOLETE]if(xstr_out.atoms.size() != conventional_basis_atoms.size()){
        //DX20210315 [OBSOLETE]  if(LDEBUG) {
        //DX20210315 [OBSOLETE]    cerr << "SYM::ConventionalCell: AddAtom() removed conventional basis atoms due to the tolerance value. Change the tolerance." << endl;
        //DX20210315 [OBSOLETE]  }
        //DX20210315 [OBSOLETE]  symmetryfound = false;
        //DX20210315 [OBSOLETE]  return xstr_out;
        //DX20210315 [OBSOLETE]}
        //DX20210315 [OBSOLETE]//DX20190410 END

        //DX20210315 [OBSOLETE]// Set number of each type
        //DX20210315 [OBSOLETE]//DX20210129 [OBSOLETE - ReplaceAtoms takes care of this] xstr_out.SetNumEachType();
        //DX20210315 [OBSOLETE]if(xstr_out.num_each_type.size() != in_names.size()) {
        //DX20210315 [OBSOLETE]  xstr_out = pflow::SetAllAtomNames(xstr_out, in_names);
        //DX20210315 [OBSOLETE]}
        //DX20210315 [OBSOLETE]xstr_out.FixLattices();

        // ************************************************************************************************
        // Check if lattice and crystal symmetry are commensurate
        // ************************************************************************************************
        xstr_out.crystal_system_ITC = crystalsystem;
        orient(xstr_out);
        if(LDEBUG) {
          cerr << "SYM::ConventionalCell: Oriented conventional cell: " << endl;
          cerr << xstr_out << endl;
        }
        xstr = xstr_out;

        // ========== PEARSON SYMBOL: Record lattice type and centering ========== //
        map<int, char> centering;
        centering[1] = 'P';
        centering[2] = 'I';
        centering[3] = 'R';
        centering[4] = 'F';
        centering[5] = 'C';
        xstr_out.bravais_label_ITC = centering[bravais_count];
        //DX20210316 [OBSOLETE] string pearson_symbol = getPearsonSymbol(xstr_out.bravais_label_ITC, lattice_char, conventional_basis_atoms);
        string pearson_symbol = getPearsonSymbol(xstr_out.bravais_label_ITC, lattice_char, xstr_out.atoms);

        if(LDEBUG) { cerr << "SYM::ConventionalCell: Pearson: " << pearson_symbol << " [dir=" << xstr.directory << "]." << endl; }
        checkops = check_ccell(xstr_out,ITC_sym_info); //DX20190215 - add ITC sym info

        // Store lattice symmetry operators only if not reformed (lattice symmetry should not know anything about crystal symmetry)
        // So we use the first iteration (unreformed)
        if(lattice_reformed == false) {
          lattice_reformed = true;
          lattice_pgroups = checkops.lattice_pgroups;
          lattice_sym_mats = checkops.lattice_sym_mats;
        }
        // Store crystal symmetry operators
        crystal_sym_mats = checkops.crystal_sym_mats;

        if(checkops.latticesystem == "REDO") {
          if(LDEBUG) { cerr << "SYM::ConventionalCell: Symmetry check failed... [dir=" << xstr.directory << "]." << endl; }
        }
        if(checkops.commensurate == true) {
          //cerr <<"reform unnecessary" << endl;
          reform = false;
          xstr_out.crystal_system_ITC = checkops.crystalsystem;
          crystalsystem_prev = crystalsystem;
          candidate_lattice_vectors_prev = candidate_lattice_vectors;
          candidate_lattice_chars_prev = candidate_lattice_chars;
          symmetryfound = true;
          if(i == candidate_lattice_vectors.size() - 1) {
            last_orientation = true;
          }
          //lattice_vector_choices.push_back(CL);
          //DX20220827 - break;
          //DX2022088 - continue;
          break;
        }

        // ========== If all combinations of cell orientations are exhausted, reform the lattice ========== //
        if(i == candidate_lattice_vectors.size() - 1 && (checkops.commensurate == false || checkops.latticesystem == "REDO")) {  //added redo as test
          if(foundspacegroup_global) { last_orientation = true; return xstr_out; } //DX20220924
          if(LDEBUG) { cerr << "SYM::ConventionalCell: All permutations of cell checked. Need to reform the lattice based on the crystal symmetry [dir=" << xstr.directory << "]." << endl; }
          //cerr << "i: " << i << endl;
          //cerr << "checkops.latticesystem: " << checkops.latticesystem << endl;
          //cerr << "checkops.crystalsystem: " << checkops.crystalsystem << endl;

          // Reduce conventional cell into a primitive form
          //DX20210401 [OBSOLETE] xstr.GetPrimitiveCell();
          xstr.GetPrimitive(); //DX20210401
          prim = xstr; //DX20210316 - store primitive

          L = xstr.lattice;
          Linv = aurostd::inverse(L);
          bool skew = SYM::isLatticeSkewed(L, xstr.dist_nn_min, sym_tol); //DX20190215 - _SYM_TOL_ to sym_tol
          crystalsystem = checkops.crystalsystem;
          latticesystem = checkops.latticesystem;

          // Need to find the mirror operations corresponding the the reduced lattice
          // Need to do this since the reduced cell has similar operation directions (preserves some orientation info)
          expanded_lattice_1 = expand_lattice(b, b, b, L);
          expanded_cell_1 = expand_cell(L);
          mirror_lattice_vectors.clear();
          rot_lattice_vectors.clear();
          twofold_lattice_vectors.clear();

          double radius = RadiusSphereLattice(L);
          //DX20170926 [OBSOLETE] vector<Glide> mirror_ops_vec_new = mirror_operations(expanded_lattice_1, expanded_cell_1, L, Linv, mirror_lattice_vectors, radius, skew);
          //DX20170926 [OBSOLETE] vector<Screw> rot_ops_vec_new = triplet_operations(expanded_lattice_1, expanded_cell_1, L, Linv, rot_lattice_vectors, radius, skew);
          //DX20170926 [OBSOLETE] vector<Screw> twofold_ops_vec_new = twofold_operations(expanded_lattice_1, expanded_cell_1, L, Linv, twofold_lattice_vectors, radius, skew);
          vector<Glide> mirror_ops_vec_new = mirror_operations(expanded_lattice_1, expanded_cell_1, L, Linv, mirror_lattice_vectors, radius, skew, sym_tol); //DX20190215 - added sym_tol
          vector<Screw> twofold_ops_vec_new;
          for(uint m=0;m<mirror_ops_vec_new.size();m++){
            Screw candidate_rotation;
            xvector<double> axis_direction = mirror_ops_vec_new[m].return_direction();
            xvector<double> point_on_axis; //Null by default
            candidate_rotation.get_screw_direct(axis_direction, point_on_axis, 2);      
            twofold_ops_vec_new.push_back(candidate_rotation);
            twofold_lattice_vectors.push_back(mirror_lattice_vectors[m]);
          }
          vector<Screw> rot_ops_vec_new;
          if(mirror_ops_vec_new.size()>2){
            rot_ops_vec_new = triplet_operations(expanded_lattice_1, expanded_cell_1, L, Linv, rot_lattice_vectors, radius, skew, sym_tol); //DX20190215 - added sym_tol
          }
          bool all_matched_mirror = true;
          bool all_matched_twofold = true;
          bool all_matched_rot = true;

          // Now find the corresponding lattice vectors for the checkops operations with respect to the reduced representation
          vector<xvector<double> > tmp_mirror_lattice_vectors = mirror_lattice_vectors;
          vector<xvector<double> > tmp_twofold_lattice_vectors = twofold_lattice_vectors;
          vector<xvector<double> > tmp_rot_lattice_vectors = rot_lattice_vectors;
          mirror_lattice_vectors = getLatticeVectorsFromOriginalMirrorOperations(mirror_ops_vec_new, checkops.mirror_ops,
              tmp_mirror_lattice_vectors, all_matched_mirror);
          twofold_lattice_vectors = getLatticeVectorsFromOriginalRotationOperations(twofold_ops_vec_new, rot_ops_vec_new,
              checkops.twofold_ops, tmp_twofold_lattice_vectors,
              tmp_rot_lattice_vectors, all_matched_twofold);
          if(mirror_ops_vec_new.size()>2){
            rot_lattice_vectors = getLatticeVectorsFromOriginalRotationOperations(twofold_ops_vec_new, rot_ops_vec_new,
                checkops.rot_ops, tmp_twofold_lattice_vectors,
                tmp_rot_lattice_vectors, all_matched_rot);
          }
          // If we could not match all operations from the old representation to the new, then we need to change the tolerance
          if(!all_matched_mirror || !all_matched_twofold || !all_matched_rot) {
            if(LDEBUG) {
              cerr << "SYM::ConventionalCell: WARNING: Couldn't match new mirror operations to old vectors. Changing tolerance." << endl;
            }
            symmetryfound = false;
            return xstr_out;
          }

          // Need to set the operators to the ones in checkops. These are the new ones which contain
          // the correct symmetry for the crystal
          twofold_ops_vec = checkops.twofold_ops;
          rot_ops_vec = checkops.rot_ops;
          mirror_ops_vec = checkops.mirror_ops;

          // ===== Reset variables ===== //
          lattice_vectors_found = false;
          candidate_lattice_vectors.clear();
          candidate_lattice_chars.clear();
          candidate_lattice_vectors_prev.clear();
          candidate_lattice_chars_prev.clear();
        }
      }  //END CLs
      if(latticesystem == "REDO") {
        //cerr << "NEED TO REDO" << endl;
        if(LDEBUG) { cerr << "SYM::ConventionalCell: Could not find lattice consistent with crystal symmetry [dir=" << xstr.directory << "]." << endl; }
        symmetryfound = false;
        return xstr_out;
      }

      if(symmetryfound == true) {
        sym_count++;
      }
      // **********************************************************************************************************************
      // ============================== If symmetryfound = false, need to readjust tolerances ============================== //
      // **********************************************************************************************************************
      else if(symmetryfound == false) {
        if(LDEBUG) { cerr << "SYM::ConventionalCell: Conventional cell consistent with crystal symmetry is found [dir=" << xstr.directory << "]." << endl; }
        return xstr_out;
      }
      // if we made it here, we've reformed the lattice and this is not the last orientation
      //last_orientation = false; //DX20220924
    }
    //cerr << "lattice_vector_choices= " << lattice_vector_choices.size() << endl;
    //for(uint l=0;l<lattice_vector_choices.size(); l++){
    //  cerr << "l=" << lattice_vector_choices[l] << endl;
    //}
    xstr_out.standard_lattice_ITC = LPRIM;
    xstr_out.standard_basis_ITC = primitive_atomic_basis;

    return xstr_out;
  }
} //namespace SYM

// ***********************************************************************************************************************
// Find cubic lattice
// ***********************************************************************************************************************
namespace SYM {
  bool findCubicLattice(vector<xvector<double> >& rot_lattice_vectors, vector<Screw>& rot_ops_vec,
      vector<xmatrix<double> >& candidate_lattice_vectors, vector<char>& candidate_lattice_chars, double& tol) { //DX20190215 - added tol
    // CUBIC lattice vector criteria: Parallel to 3 equivalent fourfold axes.
    bool LDEBUG = (FALSE || XHOST.DEBUG);

    xmatrix<double> CL;
    xvector<double> conv_lattice_vec_a;
    xvector<double> conv_lattice_vec_b;
    xvector<double> conv_lattice_vec_c;
    vector<xvector<double> > conven;
    for (uint i = 0; i < rot_ops_vec.size(); i++) {
      if(aurostd::abs(rot_ops_vec[i].return_order() - 4.0) < _ZERO_TOL_) {
        conven.push_back(rot_lattice_vectors[i]);
      }
    }

    // Check if any vectors are the same
    if(latticeVectorsSame(conven[0], conven[1], conven[2], tol)) { //DX20190215 - _SYM_TOL_ to tol
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Two or more of conventional basis vectors are the same..." << endl; }
      return false;
    }

    // Ensure lengths of vectors are the same
    if(aurostd::abs(aurostd::modulus(conven[0]) - aurostd::modulus(conven[1])) > tol && //DX20190215 - _SYM_TOL_ to tol
        aurostd::abs(aurostd::modulus(conven[0]) - aurostd::modulus(conven[2])) > tol && //DX20190215 - _SYM_TOL_ to tol
        aurostd::abs(aurostd::modulus(conven[2]) - aurostd::modulus(conven[1])) > tol) { //DX20190215 - _SYM_TOL_ to tol
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Lattice vectors are not the same length..." << endl; }
      return false;
    }

    if(conven.size() != 3) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"There should be three and only three 4-fold axes",_INDEX_BOUNDS_);
    }

    // ==== Orient into positive quadrant ==== //
    orientVectorsPositiveQuadrant(conven[0], tol); //DX20190215 - _SYM_TOL_ to tol
    orientVectorsPositiveQuadrant(conven[1], tol); //DX20190215 - _SYM_TOL_ to tol
    orientVectorsPositiveQuadrant(conven[2], tol); //DX20190215 - _SYM_TOL_ to tol

    alignLatticeWithXYZ(conven[0], conven[1], conven[2], tol); //DX20190215 - _SYM_TOL_ to tol

    // ===== Store possible lattice vectors ===== //
    // Need to check all possible orientation for cubic lattice to find the correct origin shift
    CL = xvec2xmat(conven[0], conven[1], conven[2]);
    candidate_lattice_vectors.push_back(CL);
    candidate_lattice_chars.push_back('c');
    CL = xvec2xmat(conven[2], conven[0], conven[1]);
    candidate_lattice_vectors.push_back(CL);
    candidate_lattice_chars.push_back('c');
    CL = xvec2xmat(conven[2], conven[1], conven[0]);
    candidate_lattice_vectors.push_back(CL);
    candidate_lattice_chars.push_back('c');
    CL = xvec2xmat(conven[0], conven[2], conven[1]);
    candidate_lattice_vectors.push_back(CL);
    candidate_lattice_chars.push_back('c');
    CL = xvec2xmat(conven[1], conven[0], conven[2]);
    candidate_lattice_vectors.push_back(CL);
    candidate_lattice_chars.push_back('c');
    CL = xvec2xmat(conven[1], conven[2], conven[0]);
    candidate_lattice_vectors.push_back(CL);
    candidate_lattice_chars.push_back('c');

    return true;
  }
} //namespace SYM

// ***********************************************************************************************************************
// Find trigonal/hexagonal lattice
// ***********************************************************************************************************************
namespace SYM {
  bool findTrigonalLattice(vector<xvector<double> >& rot_lattice_vectors, vector<Screw>& rot_ops_vec, vector<xvector<double> >& big_expanded,
      vector<xmatrix<double> >& candidate_lattice_vectors, vector<char>& candidate_lattice_chars, double& tol) { //DX20190215 - added tol
    // HEXAGONAL lattice vector criteria:
    //       lattice vector c        : Parallel to 6-fold axis
    //       lattice vectors a, b    : Parallel to 2-fold axes (angle=120 degrees); a=b
    bool LDEBUG = (FALSE || XHOST.DEBUG);

    xmatrix<double> CL;
    xvector<double> conv_lattice_vec_a;
    xvector<double> conv_lattice_vec_b;
    xvector<double> conv_lattice_vec_c;
    vector<xvector<double> > conven;
    xvector<double> null;
    null[1] = null[2] = null[3] = 0.0;

    // ===== Conventional lattice c ===== //
    vector<xvector<double> > six_or_threefold_lattice_vectors;
    vector<xvector<double> > twofold;
    for (uint i = 0; i < rot_ops_vec.size(); i++) {
      if(aurostd::abs(rot_ops_vec[i].return_order() - 6.0) < _ZERO_TOL_) {
        six_or_threefold_lattice_vectors.push_back(rot_lattice_vectors[i]);
      }
    }
    if(six_or_threefold_lattice_vectors.size() == 0) {
      //cerr << "TRIGONAL HEXAGONAL CELL (ie, no 6fold axis)" << endl;
      for (uint i = 0; i < rot_ops_vec.size(); i++) {
        if(aurostd::abs(rot_ops_vec[i].return_order() - 3.0) < _ZERO_TOL_) {
          six_or_threefold_lattice_vectors.push_back(rot_lattice_vectors[i]);
        }
      }
    }

    if(six_or_threefold_lattice_vectors.size() == 0) {
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << " c-axis not found. Tolerances may be too tight." << endl; }
      return false;
    }

    conv_lattice_vec_c = six_or_threefold_lattice_vectors[0];

    if(conv_lattice_vec_c == null) {
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << " c-axis is null." << endl; }
      return false;
    }
    if(conv_lattice_vec_c[2] < 0.0) {
      conv_lattice_vec_c = -1.0 * conv_lattice_vec_c;
    }

    // ===== Conventional lattices a and b ===== //
    xvector<double> tmpa, tmpb;
    vector<xvector<double> > possible_latt_a_b;
    double min_norm = 1e9; //DX20180514 - added initialization
    int count = 0;

    possible_latt_a_b = find_vectors_inplane(big_expanded, conv_lattice_vec_c, tol); //DX20190215 - added tol
    for (uint i = 0; i < possible_latt_a_b.size(); i++) {
      for (uint j = 0; j < possible_latt_a_b.size(); j++) {
        tmpa = possible_latt_a_b[i];
        tmpb = possible_latt_a_b[j];
        if(aurostd::abs(aurostd::modulus(tmpa) - aurostd::modulus(tmpb)) < tol && //DX20190215 - _SYM_TOL_ to tol
            SYM::checkAngle(tmpa, tmpb, (2.0 * Pi_r / 3.0), tol) && //DX20190215 - _SYM_TOL_ to tol
            SYM::checkAngle(tmpa, conv_lattice_vec_c, (Pi_r / 2.0), tol) && //DX20190215 - _SYM_TOL_ to tol
            SYM::checkAngle(tmpb, conv_lattice_vec_c, (Pi_r / 2.0), tol) && //DX20190215 - _SYM_TOL_ to tol
            tmpa != null && tmpb != null) {
          count++;
          if(count > 1) {
            if(aurostd::modulus(tmpa) <= min_norm) {
              min_norm = aurostd::modulus(tmpa);
              conv_lattice_vec_a = tmpa;
              conv_lattice_vec_b = tmpb;
            }
          } else {
            min_norm = aurostd::modulus(tmpa);
            conv_lattice_vec_a = tmpa;
            conv_lattice_vec_b = tmpb;
          }
        }
      }
    }

    // ===== Check if any null vectors ===== //
    vector<xvector<double> > abc_vecs;
    abc_vecs.push_back(conv_lattice_vec_a);
    abc_vecs.push_back(conv_lattice_vec_b);
    abc_vecs.push_back(conv_lattice_vec_c);
    if(anyNullVectors(abc_vecs, tol)) { //DX20190215 - _SYM_TOL_ to tol
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << " One or more of the lattice vectors is null" << endl; }
      return false;
    }

    // ===== Store possible lattice vectors ===== //
    // Only one possible orientation for trigonal/hexagonal lattice
    CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
    if(aurostd::determinant(CL) < tol) { //DX20190215 - _SYM_TOL_ to tol
      CL = xvec2xmat(conv_lattice_vec_b, conv_lattice_vec_a, conv_lattice_vec_c);
    }
    candidate_lattice_vectors.push_back(CL);
    candidate_lattice_chars.push_back('h');
    //DX20221110 - TRY OBVERSE SETTING - START
    CL = xvec2xmat(-conv_lattice_vec_b, -conv_lattice_vec_a, conv_lattice_vec_c);
    candidate_lattice_vectors.push_back(CL);
    candidate_lattice_chars.push_back('h');
    //DX20221110 - TRY OBVERSE SETTING - STOP

    return true;
  }
} //namespace SYM

// ***********************************************************************************************************************
// Find tetragonal lattice
// ***********************************************************************************************************************
namespace SYM {
  bool findTetragonalLattice(vector<xvector<double> >& rot_lattice_vectors, vector<xvector<double> >& twofold_lattice_vectors,
      vector<Screw>& rot_ops_vec, vector<Screw>& twofold_ops_vec, vector<xvector<double> >& big_expanded,
      vector<xmatrix<double> >& candidate_lattice_vectors, vector<char>& candidate_lattice_chars, double& tol) { //DX20190215 - added tol
    // TETRAGONAL lattice vectors criteria:
    //        lattice vector c        : Parallel to 4-fold axis
    //        lattice vectors a, b    : Parallel to 2-fold axes (angle=90 degrees); a=b
    bool LDEBUG = (FALSE || XHOST.DEBUG);

    xmatrix<double> CL;
    xvector<double> conv_lattice_vec_a;
    xvector<double> conv_lattice_vec_b;
    xvector<double> conv_lattice_vec_c;
    vector<xvector<double> > conven;
    vector<xvector<double> > fourfold_lattice_vectors;
    vector<xvector<double> > twofold;
    xvector<double> null;
    null[1] = null[2] = null[3] = 0.0;

    // ===== Conventional lattice c ===== //
    for (uint i = 0; i < rot_ops_vec.size(); i++) {
      if(aurostd::abs(rot_ops_vec[i].return_order() - 4.0) < _ZERO_TOL_) {
        fourfold_lattice_vectors.push_back(rot_lattice_vectors[i]);
      }
    }
    for (uint f = 0; f < fourfold_lattice_vectors.size(); f++) {
      conv_lattice_vec_c = fourfold_lattice_vectors[f];
      if(conv_lattice_vec_c == null) {
        if(LDEBUG) { cerr << __AFLOW_FUNC__ << " c-axis is null." << endl; }
        return false;
      }

      // ===== Conventional lattices a and b ===== //
      xvector<double> tmpa, tmpb;
      vector<xvector<double> > possible_latt_a_b;
      double min_norm = 1e9; //DX20180514 - added initialization
      int count = 0;
      bool all_dir_equal = true;

      // ===== Check to see if all two-fold operator directions are the same ===== //
      // After the first reformation, the checkops section often neglects to find
      // different operators.  If this occurs find all the vectors in the plane perpendicular
      // to the c axis
      for (uint i = 0; i < twofold_ops_vec.size(); i++) {
        xvector<double> rot_dir;
        if(i == 0) {
          rot_dir = twofold_ops_vec[i].return_direction();
        } else {
          xvector<double> rot_tmp = twofold_ops_vec[i].return_direction();
          if((rot_dir[1] - rot_tmp[1]) > tol && (rot_dir[2] - rot_tmp[2]) > tol && (rot_dir[3] - rot_tmp[3]) > tol) { //DX20190215 - _SYM_TOL_ to tol
            all_dir_equal = false;
            break;
          }
        }
      }

      // ===== All directions the same or there is only one 2-fold operator; find vectors perpendicular to c axis
      if(twofold_ops_vec.size() <= 1 || all_dir_equal == true) {
        possible_latt_a_b = find_vectors_inplane(big_expanded, conv_lattice_vec_c, tol); //DX20190215 - added tol
      } else {
        // ===== If there are different 2-fold operators ===== //
        for (uint i = 0; i < twofold_ops_vec.size(); i++) {
          possible_latt_a_b.push_back(twofold_lattice_vectors[i]);
        }
      }
      for (uint i = 0; i < possible_latt_a_b.size(); i++) {
        for (uint j = 0; j < possible_latt_a_b.size(); j++) {
          tmpa = possible_latt_a_b[i];
          tmpb = possible_latt_a_b[j];
          if(aurostd::abs(aurostd::modulus(tmpa) - aurostd::modulus(tmpb)) < tol && //DX20190215 - _SYM_TOL_ to tol
              SYM::checkAngle(tmpa, tmpb, (Pi_r / 2.0), tol) && //DX20190215 - _SYM_TOL_ to tol
              SYM::checkAngle(tmpa, conv_lattice_vec_c, (Pi_r / 2.0), tol) && //DX20190215 - _SYM_TOL_ to tol
              SYM::checkAngle(tmpb, conv_lattice_vec_c, (Pi_r / 2.0), tol) && //DX20190215 - _SYM_TOL_ to tol
              tmpa != null && tmpb != null) {
            count++;
            if(count > 1) {
              if(aurostd::modulus(tmpa) <= min_norm) {
                min_norm = aurostd::modulus(tmpa);
                conv_lattice_vec_a = tmpa;
                conv_lattice_vec_b = tmpb;
              }
            } else {
              min_norm = aurostd::modulus(tmpa);
              conv_lattice_vec_a = tmpa;
              conv_lattice_vec_b = tmpb;
            }
          }
        }
      }

      // ===== Check if any null vectors ===== //
      vector<xvector<double> > ab_vecs;
      ab_vecs.push_back(conv_lattice_vec_a);
      ab_vecs.push_back(conv_lattice_vec_b);
      if(anyNullVectors(ab_vecs, tol)) { //DX20190215 - _SYM_TOL_ to tol
        if(f != fourfold_lattice_vectors.size() - 1) {
          continue;
        } else {
          if(LDEBUG) { cerr << __AFLOW_FUNC__ << " One or more of the lattice vectors is null" << endl; }
          return false;
        }
      }
      else {
        //DX20210107 [OBSOLETE - do not break, store all possible conventional lattice possibilites]

        // ==== Orient into positive quadrant ==== //
        orientVectorsPositiveQuadrant(conv_lattice_vec_a, tol); //DX20190215 - _SYM_TOL_ to tol
        orientVectorsPositiveQuadrant(conv_lattice_vec_b, tol); //DX20190215 - _SYM_TOL_ to tol

        // Order lattice vectors so that the a vector has a positive x coordinate
        if(conv_lattice_vec_b(1) > tol && conv_lattice_vec_b(2) < tol && conv_lattice_vec_b(3) < tol) { //DX20190215 - _SYM_TOL_ to tol
          xvector<double> tmp = conv_lattice_vec_a;
          conv_lattice_vec_a = conv_lattice_vec_b;
          conv_lattice_vec_b = tmp;
        }

        // === Proceed if lattice vectors are not null === //
        CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
        candidate_lattice_vectors.push_back(CL);
        candidate_lattice_chars.push_back('t');
        CL = xvec2xmat(conv_lattice_vec_b, conv_lattice_vec_a, conv_lattice_vec_c); //DX20231010 - check alternative lattice vector order to minimize Wyckoff enumeration
        candidate_lattice_vectors.push_back(CL); //DX20231010 - check alternative lattice vector order to minimize Wyckoff enumeration
        candidate_lattice_chars.push_back('t'); //DX20231010 - check alternative lattice vector order to minimize Wyckoff enumeration
      } //DX20210107
    } //DX20210107

    return true;
  }
} //namespace SYM

// ***********************************************************************************************************************
// Find monoclinic lattice
// ***********************************************************************************************************************
namespace SYM {
  bool findMonoclinicLattice(vector<xvector<double> >& mirror_lattice_vectors, vector<xvector<double> >& twofold_lattice_vectors,
      vector<xvector<double> >& big_expanded, vector<xmatrix<double> >& candidate_lattice_vectors,
      vector<char>& candidate_lattice_chars, int& cell_choice, double& tol) { //DX20180816 - added setting_choice //DX20190215 - added tol
    // MONOCLINIC lattice vectors criteria:
    //        lattice vector unqiue axis (b or c)        : Parallel to mirror axis or 2-fold axis
    //        lattice vectors a, (c or b)                : Perpendicular to lattice vector c; angle != 90 degree; a!=b

    // ===== MONOCLINIC has different cell choices and different unique axes in the ITC; need to check all 6 possibilities ===== //
    // num_reinitials : variable to scan through combinations of vectors in monoclinic plane (e,f, or g: shortest vectors in plane)
    // IT            : scan through unique axis b or c
    bool LDEBUG = (FALSE || XHOST.DEBUG);

    xmatrix<double> CL;
    xvector<double> conv_lattice_vec_a;
    xvector<double> conv_lattice_vec_b;
    xvector<double> conv_lattice_vec_c;
    vector<xvector<double> > conven;
    xvector<double> null;
    null[1] = null[2] = null[3] = 0.0;

    // ===== Conventional lattice unique axis ===== //
    xvector<double> h_lat;
    vector<xvector<double> > possible_h_lats;
    if(mirror_lattice_vectors.size() > 0) {
      if(mirror_lattice_vectors.size() > 1) {
        if(aurostd::modulus(mirror_lattice_vectors[1]) < aurostd::modulus(mirror_lattice_vectors[0])) {
          possible_h_lats.push_back(mirror_lattice_vectors[1]);
          possible_h_lats.push_back(-mirror_lattice_vectors[1]); //DX20221129
          possible_h_lats.push_back(mirror_lattice_vectors[0]);
          possible_h_lats.push_back(-mirror_lattice_vectors[0]); //DX20221129
        } else {
          possible_h_lats.push_back(mirror_lattice_vectors[0]);
          possible_h_lats.push_back(-mirror_lattice_vectors[0]); //DX20221129
          possible_h_lats.push_back(mirror_lattice_vectors[1]);
          possible_h_lats.push_back(-mirror_lattice_vectors[1]); //DX20221129
        }
        if(LDEBUG) {
          cerr << __AFLOW_FUNC__ << " Two choices for unique axis. " << endl;
          cerr << __AFLOW_FUNC__ << " choice 1: " << possible_h_lats[0] << " (mod=" << aurostd::modulus(possible_h_lats[0]) << ")" << endl;
          cerr << __AFLOW_FUNC__ << " choice 2: " << possible_h_lats[1] << " (mod=" << aurostd::modulus(possible_h_lats[1]) << ")" << endl;
        }
      } 
      else {
        possible_h_lats.push_back(mirror_lattice_vectors[0]);
        possible_h_lats.push_back(-mirror_lattice_vectors[0]); //DX20221129
        if(LDEBUG) {
          cerr << __AFLOW_FUNC__ << " One choice for unique axis. " << endl;
          cerr << __AFLOW_FUNC__ << " choice 1: " << possible_h_lats[0] << " (mod=" << aurostd::modulus(possible_h_lats[0]) << ")" << endl;
        }
      }
    } else if(twofold_lattice_vectors.size() > 0) {
      possible_h_lats.push_back(twofold_lattice_vectors[0]);
      possible_h_lats.push_back(-twofold_lattice_vectors[0]); //DX20221129
      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " One choice for unique axis. " << endl;
        cerr << __AFLOW_FUNC__ << " choice 1: " << possible_h_lats[0] << " (mod=" << aurostd::modulus(possible_h_lats[0]) << ")" << endl;
      }
    }
    if(possible_h_lats.size() == 0) {
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << " No unique axis found" << endl; }
      return false;
    }
    if(possible_h_lats[0] == null) {
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Unique axis is null" << endl; }
      return false;
    }
    // ===== Determine direction of unique axis ===== //

    // =============== Find vectors in monoclinic plane (e, f, and g) =============== //
    // === Get vectors e and f (shortest 2 vectors in monoclinic plane) === //
    // Criteria on e and f:
    //          1) modulus(e)+modulus(f) is as small as possible
    //          2) Angle between e and f < 120 degrees
    //          3) Angle between e and f > 90 degrees
    xvector<double> tranvec_e;
    xvector<double> tranvec_f;
    xvector<double> tranvec_g;

    for (uint h = 0; h < possible_h_lats.size(); h++) {
      vector<xvector<double> > possible_efg_vectors; //DX20180110 - scoping issue; place here, not outside for loop
      vector<vector<int> > index;  //DX20180110 - scoping issue; place here, not outside for loop
      vector<double> total_mod;  //DX20180110 - scoping issue; place here, not outside for loop
      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " Testing choice " << h+1  << ": " << possible_h_lats[h] << endl;
      }
      h_lat = possible_h_lats[h];
      possible_efg_vectors = find_vectors_inplane(big_expanded, h_lat, tol); //DX20190215 - added tol
      for (uint i = 0; i < possible_efg_vectors.size(); i++) {
        for (uint j = i + 1; j < possible_efg_vectors.size(); j++) {
          for (uint k = j + 1; k < possible_efg_vectors.size(); k++) {
            int e = -1;
            int f = -1;
            int g = -1;
            xvector<double> tranvec_i = possible_efg_vectors[i];
            xvector<double> tranvec_j = possible_efg_vectors[j];
            xvector<double> tranvec_k = possible_efg_vectors[k];

            // ===== Ensure none of the possible vectors is a null vector ===== //
            if((aurostd::abs(tranvec_i[1]) >= _ZERO_TOL_ || aurostd::abs(tranvec_i[2]) >= _ZERO_TOL_ || aurostd::abs(tranvec_i[3]) >= _ZERO_TOL_) &&
                (aurostd::abs(tranvec_j[1]) >= _ZERO_TOL_ || aurostd::abs(tranvec_j[2]) >= _ZERO_TOL_ || aurostd::abs(tranvec_j[3]) >= _ZERO_TOL_) &&
                (aurostd::abs(tranvec_k[1]) >= _ZERO_TOL_ || aurostd::abs(tranvec_k[2]) >= _ZERO_TOL_ || aurostd::abs(tranvec_k[3]) >= _ZERO_TOL_)) {
              vector<int> int_ijk;
              int_ijk.push_back(i);
              int_ijk.push_back(j);
              int_ijk.push_back(k);
              vector<xvector<double> > tranvec_ijk;
              tranvec_ijk.push_back(tranvec_i);
              tranvec_ijk.push_back(tranvec_j);
              tranvec_ijk.push_back(tranvec_k);
              vector<double> mod_tranvec_ijk;
              mod_tranvec_ijk.push_back(aurostd::modulus(tranvec_i));
              mod_tranvec_ijk.push_back(aurostd::modulus(tranvec_j));
              mod_tranvec_ijk.push_back(aurostd::modulus(tranvec_k));
              int tmp_int_e = -1;
              int tmp_int_f = -1;
              int tmp_int_g = -1;

              // ===== Determine the 1st smallest vector : "e" ===== //
              tmp_int_e = smallest_gt_min_index(0, -1, -1, mod_tranvec_ijk);
              tranvec_e = tranvec_ijk[tmp_int_e];
              e = int_ijk[tmp_int_e];

              // ===== Determine the 2nd smallest vector : "f" ===== //
              tmp_int_f = smallest_gt_min_index(aurostd::modulus(tranvec_e), tmp_int_e, -1, mod_tranvec_ijk);
              tranvec_f = tranvec_ijk[tmp_int_f];
              f = int_ijk[tmp_int_f];

              // ===== Determine the 2nd smallest vector : "g" ===== //
              tmp_int_g = smallest_gt_min_index(aurostd::modulus(tranvec_f), tmp_int_e, tmp_int_f, mod_tranvec_ijk);
              tranvec_g = tranvec_ijk[tmp_int_g];
              g = int_ijk[tmp_int_g];

              xvector<double> closed_triangle = tranvec_e + tranvec_f + tranvec_g;
              if(nullVector(closed_triangle, tol)) { //DX20190215 - _SYM_TOL_ to tol
                //cerr << "closed angle satisfied" << endl;
                // ===== Neglect anti-parallel vector combinations ===== //
                if(!SYM::checkAngle(tranvec_e, tranvec_f, Pi_r, tol) && //DX20190215 - _SYM_TOL_ to tol
                    !SYM::checkAngle(tranvec_f, tranvec_g, Pi_r, tol) && //DX20190215 - _SYM_TOL_ to tol
                    !SYM::checkAngle(tranvec_g, tranvec_e, Pi_r, tol)) { //DX20190215 - _SYM_TOL_ to tol
                  //cerr << "not anti-parallel" << endl;
                  // ===== Determine the determinant of e,f,g with h_lat (unique axis_b) ===== //
                  int consistent = 0;
                  if((cell_choice == 0 || cell_choice == SG_SETTING_1 || cell_choice == SG_SETTING_ANRL) && aurostd::det(xvec2xmat(tranvec_f, h_lat, tranvec_e)) > 0.0 && aurostd::det(xvec2xmat(tranvec_g, h_lat, tranvec_f)) > 0.0 && //DX20190131 - add ANRL setting choice
                      aurostd::det(xvec2xmat(tranvec_e, h_lat, tranvec_g)) > 0.0) { //DX20180816 - added cell choice; check axis b
                    consistent = 3;
                  }
                  else if(cell_choice == SG_SETTING_2 && aurostd::det(xvec2xmat(tranvec_f, h_lat, tranvec_e)) > 0.0 && aurostd::det(xvec2xmat(tranvec_g, h_lat, tranvec_f)) > 0.0 &&
                      aurostd::det(xvec2xmat(tranvec_e, h_lat, tranvec_g)) > 0.0) { //DX20180816 - added cell choice; check axis c
                    consistent = 3;
                  }
                  if(consistent != 3 && (aurostd::abs(aurostd::modulus(tranvec_e) - aurostd::modulus(tranvec_f)) < tol || //DX20190215 - _SYM_TOL_ to tol
                        aurostd::abs(aurostd::modulus(tranvec_f) - aurostd::modulus(tranvec_g)) < tol || //DX20190215 - _SYM_TOL_ to tol
                        aurostd::abs(aurostd::modulus(tranvec_g) - aurostd::modulus(tranvec_e)) < tol)) { //DX20190215 - _SYM_TOL_ to tol
                    int tmp_int = -1;
                    if(aurostd::abs(aurostd::modulus(tranvec_e) - aurostd::modulus(tranvec_f)) < tol) { //DX20190215 - _SYM_TOL_ to tol
                      tmp_int = e;
                      e = f;
                      f = tmp_int;
                    } else if(aurostd::abs(aurostd::modulus(tranvec_f) - aurostd::modulus(tranvec_g)) < tol) { //DX20190215 - _SYM_TOL_ to tol
                      tmp_int = f;
                      f = g;
                      g = tmp_int;
                    } else if(aurostd::abs(aurostd::modulus(tranvec_g) - aurostd::modulus(tranvec_e)) < tol) { //DX20190215 - _SYM_TOL_ to tol
                      tmp_int = g;
                      g = e;
                      e = tmp_int;
                    }
                    tranvec_e = possible_efg_vectors[e];
                    tranvec_f = possible_efg_vectors[f];
                    tranvec_g = possible_efg_vectors[g];
                    if((cell_choice == 0 || cell_choice == SG_SETTING_1 || cell_choice == SG_SETTING_ANRL) && aurostd::det(xvec2xmat(tranvec_f, h_lat, tranvec_e)) > 0.0 && aurostd::det(xvec2xmat(tranvec_g, h_lat, tranvec_f)) > 0.0 && //DX20190131 - add ANRL setting choice
                        aurostd::det(xvec2xmat(tranvec_e, h_lat, tranvec_g)) > 0.0) { //DX20180816 - added cell choice; check axis b
                      consistent = 3;
                    }
                    else if(cell_choice == SG_SETTING_2 && aurostd::det(xvec2xmat(tranvec_f, h_lat, tranvec_e)) > 0.0 && aurostd::det(xvec2xmat(tranvec_g, h_lat, tranvec_f)) > 0.0 &&
                        aurostd::det(xvec2xmat(tranvec_e, h_lat, tranvec_g)) > 0.0) { //DX20180816 - added cell choice; check axis c
                      consistent = 3;
                    }
                  }
                  if(consistent == 3) {
                    vector<int> tmp;
                    tmp.push_back(e);
                    tmp.push_back(f);
                    tmp.push_back(g);
                    index.push_back(tmp);
                    total_mod.push_back(aurostd::modulus(tranvec_e) + aurostd::modulus(tranvec_f) + aurostd::modulus(tranvec_g));
                  }
                }
              }
            }
          }
        }
      }
      // === Check if 2 or more vectors were chosen to comprise the monoclinic plane === //
      if(index.size() <= 1 && h == possible_h_lats.size() - 1) {
        if(LDEBUG) { cerr << __AFLOW_FUNC__ << " No lattice vectors found in the monoclinic plane" << endl; }
        return false;
      } else if(index.size() <= 1 && h < possible_h_lats.size() - 1) {
        if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Not enought lattice vectors found in monoclinic plane. Try the other unique axis." << endl; }
        possible_efg_vectors.clear();
        index.clear();
        total_mod.clear();
        continue;
      }
      //DX20171212 - Consider all possible unique axis choices - else {
      //DX20171212 - Consider all possible unique axis choices -   break;
      //DX20171212 - Consider all possible unique axis choices - }

      int e = index[smallest_gt_min_index(0, -1, -1, total_mod)][0];
      tranvec_e = possible_efg_vectors[e];
      int f = index[smallest_gt_min_index(0, -1, -1, total_mod)][1];
      tranvec_f = possible_efg_vectors[f];
      int g = index[smallest_gt_min_index(0, -1, -1, total_mod)][2];
      tranvec_g = possible_efg_vectors[g];
      double mod_e = aurostd::modulus(tranvec_e);
      double mod_f = aurostd::modulus(tranvec_f);
      double mod_g = aurostd::modulus(tranvec_g);

      if(LDEBUG) { 
        cerr << __AFLOW_FUNC__ << " 1st shortest vector in monoclinic plane (e): " << tranvec_e << " (mod=" << mod_e << ")." << endl;
        cerr << __AFLOW_FUNC__ << " 2nd shortest vector in monoclinic plane (f): " << tranvec_f << " (mod=" << mod_f << ")." << endl;
        cerr << __AFLOW_FUNC__ << " 3rd shortest vector in monoclinic plane (g): " << tranvec_g << " (mod=" << mod_g << ")." << endl;
        cerr << __AFLOW_FUNC__ << " Angle between e and f: " << aurostd::angle(tranvec_e,tranvec_f)*180.0/Pi_r << endl;
        cerr << __AFLOW_FUNC__ << " Angle between f and g: " << aurostd::angle(tranvec_f,tranvec_g)*180.0/Pi_r << endl;
        cerr << __AFLOW_FUNC__ << " Angle between g and e: " << aurostd::angle(tranvec_g,tranvec_e)*180.0/Pi_r << endl;
      }

      if(mod_e <= _ZERO_TOL_ || mod_f <= _ZERO_TOL_ || mod_g <= _ZERO_TOL_) {
        if(LDEBUG) { cerr << __AFLOW_FUNC__ << " One or more lattice vectors in the monoclinic plane is null" << endl; }
        return false;
      }

      // ========== Pick unique axis-b (See ITC p.38-40. Use option (i) for finding the monoclinic cell discussed on p.40) ========== //
      if(cell_choice == 0 || cell_choice == SG_SETTING_1 || cell_choice == SG_SETTING_ANRL){ //DX20190131 - add ANRL setting
        conv_lattice_vec_b = h_lat;
        // ===== e and f ===== //
        conv_lattice_vec_a = tranvec_f;  // 2nd shortest (order determined by RH coordinate system requirement)
        conv_lattice_vec_c = tranvec_e;  // 1st shortest (order determined by RH coordinate system requirement)
        CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
        if(det(CL) < 0.0) {
          //cerr << "[0] IN MONO: Negative det" << endl;
          conv_lattice_vec_a = tranvec_e;  // 2nd shortest (order determined by RH coordinate system requirement)
          conv_lattice_vec_c = tranvec_f;  // 1st shortest (order determined by RH coordinate system requirement)
          CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
        }
        candidate_lattice_vectors.push_back(CL);
        candidate_lattice_chars.push_back('b');

        // ===== e and g ===== //
        conv_lattice_vec_a = tranvec_e;  // 3rd shortest (order determined by RH coordinate system requirement)
        conv_lattice_vec_c = tranvec_g;  // 1sd shortest (order determined by RH coordinate system requirement)
        CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
        if(det(CL) < 0.0) {
          //cerr << "[1] IN MONO: Negative det" << endl;
          conv_lattice_vec_a = tranvec_g;  // 3rd shortest (order determined by RH coordinate system requirement)
          conv_lattice_vec_c = tranvec_e;  // 1sd shortest (order determined by RH coordinate system requirement)
          CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
        }
        candidate_lattice_vectors.push_back(CL);
        candidate_lattice_chars.push_back('b');

        // ===== g and f ===== //
        conv_lattice_vec_a = tranvec_g;  // 3rd shortest (order determined by RH coordinate system requirement)
        conv_lattice_vec_c = tranvec_f;  // 1sd shortest (order determined by RH coordinate system requirement)
        CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
        if(det(CL) < 0.0) {
          //cerr << "[2] IN MONO: Negative det" << endl;
          conv_lattice_vec_a = tranvec_f;  // 3rd shortest (order determined by RH coordinate system requirement)
          conv_lattice_vec_c = tranvec_g;  // 1sd shortest (order determined by RH coordinate system requirement)
          CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
        }
        candidate_lattice_vectors.push_back(CL);
        candidate_lattice_chars.push_back('b');
      }
      // ========== Pick unique axis-c (See ITC p.38-40. Use option (i) for finding the monoclinic cell discussed on p.40) ========== //
      else if(cell_choice == SG_SETTING_2){
        conv_lattice_vec_c = h_lat;
        // ===== e and f ===== //
        conv_lattice_vec_a = tranvec_e;  // 1st shortest (order determined by RH coordinate system requirement)
        conv_lattice_vec_b = tranvec_f;  // 2nd shortest (order determined by RH coordinate system requirement)
        CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
        if(det(CL) < 0.0) {
          //cerr << "[0] IN MONO: Negative det" << endl;
          conv_lattice_vec_a = tranvec_f;  // 2nd shortest (order determined by RH coordinate system requirement)
          conv_lattice_vec_b = tranvec_e;  // 1st shortest (order determined by RH coordinate system requirement)
          CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
        }
        candidate_lattice_vectors.push_back(CL);
        candidate_lattice_chars.push_back('m'); //unique axis c; use 'm' since 'c' is reserved for cubic

        // ===== e and g ===== //
        conv_lattice_vec_a = tranvec_g;  // 3rd shortest (order determined by RH coordinate system requirement)
        conv_lattice_vec_b = tranvec_e;  // 1st shortest (order determined by RH coordinate system requirement)
        CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
        if(det(CL) < 0.0) {
          //cerr << "[1] IN MONO: Negative det" << endl;
          conv_lattice_vec_a = tranvec_e;  // 1st shortest (order determined by RH coordinate system requirement)
          conv_lattice_vec_b = tranvec_g;  // 3rd shortest (order determined by RH coordinate system requirement)
          CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
        }
        candidate_lattice_vectors.push_back(CL);
        candidate_lattice_chars.push_back('m'); //unique axis c; use 'm' since 'c' is reserved for cubic

        // ===== g and f ===== //
        conv_lattice_vec_a = tranvec_f;  // 2nd shortest (order determined by RH coordinate system requirement)
        conv_lattice_vec_b = tranvec_g;  // 3rd shortest (order determined by RH coordinate system requirement)
        CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
        if(det(CL) < 0.0) {
          //cerr << "[2] IN MONO: Negative det" << endl;
          conv_lattice_vec_a = tranvec_g;  // 3rd shortest (order determined by RH coordinate system requirement)
          conv_lattice_vec_b = tranvec_f;  // 2nd shortest (order determined by RH coordinate system requirement)
          CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
        }
        candidate_lattice_vectors.push_back(CL);
        candidate_lattice_chars.push_back('m'); //unique axis c; use 'm' since 'c' is reserved for cubic
      }
    } //DX20171212 - Consider all possible unique axis choices 
    return true;
  }
} //namespace SYM

// ***********************************************************************************************************************
// Find triclinic lattice
// ***********************************************************************************************************************
namespace SYM {
  bool findTriclinicLattice(xmatrix<double>& lattice, vector<xmatrix<double> >& candidate_lattice_vectors,
      vector<char>& candidate_lattice_chars) {
    // TRICLINIC lattice vectors criteria: a!=b!=c; alpha,beta,gamma!=90
    candidate_lattice_vectors.push_back(lattice);
    candidate_lattice_chars.push_back('a');

    return true;
  }
}

// ***********************************************************************************************************************
// Find orthorhombic lattice
// ***********************************************************************************************************************
namespace SYM {
  bool findOrthorhombicLattice(vector<xvector<double> >& twofold_lattice_vectors, vector<xvector<double> >& mirror_lattice_vectors,
      vector<xmatrix<double> >& candidate_lattice_vectors, vector<char>& candidate_lattice_chars, double& tol) { //DX20190215 - added tol
    // ORTHORHOMBIC lattice vector criteria:
    //        lattice vectors a,b,c   : Parallel to 2-fold axis (or two 2-fold and one mirror)
    bool LDEBUG = (FALSE || XHOST.DEBUG);

    xmatrix<double> CL;
    vector<xvector<double> > conven;

    // ===== Search all combinations of twofold and mirror lattice vectors
    uint totalcount = twofold_lattice_vectors.size() + mirror_lattice_vectors.size();
    int totalnum = CombinationNr(3, totalcount);
    vector<int> combiset;
    combiset = AllCombination41(3, totalcount, 1);
    for (int i = 0; i < totalnum; i++) {
      conven.clear();
      for (int j = 0; j < 3; j++) {
        if(combiset[j] < (int)(twofold_lattice_vectors.size())) {
          conven.push_back(twofold_lattice_vectors[combiset[j]]);
        } else {
          conven.push_back(mirror_lattice_vectors[combiset[j] - twofold_lattice_vectors.size()]);
        }
      }
      if(SYM::checkAngle(conven[0], conven[1], (Pi_r / 2.0), tol) && //DX20190215 - _SYM_TOL_ to tol
          SYM::checkAngle(conven[0], conven[2], (Pi_r / 2.0), tol) && //DX20190215 - _SYM_TOL_ to tol
          SYM::checkAngle(conven[1], conven[2], (Pi_r / 2.0), tol)) { //DX20190215 - _SYM_TOL_ to tol
        //print(conven);
        if(latticeVectorsSame(conven[0], conven[1], conven[2], tol)) { //DX20190215 - _SYM_TOL_ to tol
          if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Two or more of conventional basis vectors are the same..." << endl; }
          conven.clear();
        } else if(SYM::checkAngle(conven[0], conven[1], Pi_r, tol) && //DX20190215 - _SYM_TOL_ to tol
            SYM::checkAngle(conven[0], conven[2], Pi_r, tol) && //DX20190215 - _SYM_TOL_ to tol
            SYM::checkAngle(conven[1], conven[2], Pi_r, tol)) { //DX20190215 - _SYM_TOL_ to tol
          if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Two or more of conventional basis vectors are the same..." << endl; }
          conven.clear();
        } else {
          break;
        }
      } else {
        conven.clear();
      }
      combiset = AllCombination42(3, totalcount, combiset);
    }
    if(conven.size() == 0) {
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << " No conventional basis vectors found...." << endl; }
      return false;
    }

    // ===== Possible orientations of lattice vectors ===== //
    CL = xvec2xmat(conven[0], conven[1], conven[2]);
    candidate_lattice_vectors.push_back(CL);
    candidate_lattice_chars.push_back('o');
    CL = xvec2xmat(conven[0], conven[2], conven[1]);
    candidate_lattice_vectors.push_back(CL);
    candidate_lattice_chars.push_back('o');
    CL = xvec2xmat(conven[1], conven[0], conven[2]);
    candidate_lattice_vectors.push_back(CL);
    candidate_lattice_chars.push_back('o');
    CL = xvec2xmat(conven[1], conven[2], conven[0]);
    candidate_lattice_vectors.push_back(CL);
    candidate_lattice_chars.push_back('o');
    CL = xvec2xmat(conven[2], conven[1], conven[0]);
    candidate_lattice_vectors.push_back(CL);
    candidate_lattice_chars.push_back('o');
    CL = xvec2xmat(conven[2], conven[0], conven[1]);
    candidate_lattice_vectors.push_back(CL);
    candidate_lattice_chars.push_back('o');

    return true;
  }
} //namespace SYM

// ***********************************************************************************************************************
// Find rhombohedral lattice
// ***********************************************************************************************************************
namespace SYM {
  bool findRhombohedralLattice(vector<xvector<double> >& rot_lattice_vectors, vector<Screw>& rot_ops_vec,
      vector<xvector<double> >& big_expanded, vector<xmatrix<double> >& candidate_lattice_vectors,
      vector<char>& candidate_lattice_chars, double& tol) { //DX20190215 - added tol
    // RHOMBOHEDRAL lattice vectors criteria:
    //        lattice vector c       : Parallel to 6- or 3-fold axis
    //        lattice vectors a,b    : Parallel to 2-fold axes; Angle ab = 120 degrees
    bool LDEBUG = (FALSE || XHOST.DEBUG);

    xmatrix<double> CL;
    xvector<double> conv_lattice_vec_a;
    xvector<double> conv_lattice_vec_b;
    xvector<double> conv_lattice_vec_c;
    vector<xvector<double> > conven;
    xvector<double> null;
    null[1] = null[2] = null[3] = 0.0;

    vector<xvector<double> > six_or_threefold_lattice_vectors;
    vector<xvector<double> > twofold;
    // ===== Lattice vector 'c' ===== //
    for (uint i = 0; i < rot_ops_vec.size(); i++) {
      if(aurostd::abs(rot_ops_vec[i].return_order() - 6.0) < _ZERO_TOL_) {
        six_or_threefold_lattice_vectors.push_back(rot_lattice_vectors[i]);
      }
    }
    if(six_or_threefold_lattice_vectors.size() == 0) {
      //cerr << "TRIGONAL HEXAGONAL CELL (ie, no 6fold axis)" << endl;
      for (uint i = 0; i < rot_ops_vec.size(); i++) {
        if(aurostd::abs(rot_ops_vec[i].return_order() - 3.0) < _ZERO_TOL_) {
          six_or_threefold_lattice_vectors.push_back(rot_lattice_vectors[i]);
        }
      }
    }
    if(six_or_threefold_lattice_vectors.size() == 0) {
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << " No 3- or 6-fold axis found." << endl; }
      return false;
    }

    conv_lattice_vec_c = six_or_threefold_lattice_vectors[0];
    if(conv_lattice_vec_c == null) {
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << " No 3- or 6-fold axis found." << endl; }
      return false;
    }

    vector<xvector<double> > possible_latt_a_b = find_vectors_inplane(big_expanded, conv_lattice_vec_c, tol); //DX20190215 - added tol
    // ===== Lattice vectors 'a' and 'b' ===== //
    xvector<double> tmpa, tmpb;
    double min_norm = 1e9; //DX20180514 - added initialization
    int count = 0;
    for (uint i = 0; i < possible_latt_a_b.size(); i++) {
      for (uint j = i; j < possible_latt_a_b.size(); j++) {
        tmpa = possible_latt_a_b[i];
        tmpb = possible_latt_a_b[j];
        if(aurostd::abs(aurostd::modulus(tmpa) - aurostd::modulus(tmpb)) < tol && //DX20190215 - _SYM_TOL_ to tol
            SYM::checkAngle(tmpa, tmpb, ((2.0 * Pi_r) / 3.0), tol)) { //DX20190215 - _SYM_TOL_ to tol
          count++;
          if(count > 1) {
            if(aurostd::modulus(tmpa) <= min_norm) {
              min_norm = aurostd::modulus(tmpa);
              conv_lattice_vec_a = tmpa;
              conv_lattice_vec_b = tmpb;
            }
          } else {
            min_norm = aurostd::modulus(tmpa);
            conv_lattice_vec_a = tmpa;
            conv_lattice_vec_b = tmpb;
          }
        }
      }
    }

    vector<xvector<double> > abc_vecs;
    abc_vecs.push_back(conv_lattice_vec_a);
    abc_vecs.push_back(conv_lattice_vec_b);
    abc_vecs.push_back(conv_lattice_vec_c);
    if(anyNullVectors(abc_vecs, tol)) { //DX20190215 - _SYM_TOL_ to tol
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << " One or more of the lattice vectors is null." << endl; }
      return false;
    }

    CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
    if(aurostd::determinant(CL) < tol) { //DX20190215 - _SYM_TOL_ to tol
      CL = xvec2xmat(conv_lattice_vec_b, conv_lattice_vec_a, conv_lattice_vec_c);
    }
    candidate_lattice_vectors.push_back(CL);
    candidate_lattice_chars.push_back('R');

    return true;
  }
}

// ***********************************************************************************************************************
// Find rhombohedral setting
// ***********************************************************************************************************************
namespace SYM {
  bool findRhombohedralSetting(vector<xvector<double> >& big_expanded, vector<xmatrix<double> >& candidate_lattice_vectors,
      vector<char>& candidate_lattice_chars, double& tol) { //DX20190215 - added tol
    // RHOMBOHEDRAL setting lattice vectors criteria:
    //        lattice vector a,b,c   : a=b=c=sqrt((hex_a^2/3) + (hex_c^2/9))
    //                                 Angle: alpha=beta=gamma=acos((2.0c^2 - 3a^2)/(2(c^2 + 3.0a^2)))

    // This function requires that findRhombohedralLattice() was run beforehand, since it takes 
    // in candidate_lattice_vectors for the hexagonal setting
    // In particular, the hexagonal c-axis is used to orient the rhl lattice, and the moduli of a and c
    // help determine the length of and angles between the rhombohedral lattice vectors

    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Finding rhombohedral setting of structure ..." << endl; }

    xmatrix<double> CL;
    xvector<double> conv_lattice_vec_a;
    xvector<double> conv_lattice_vec_b;
    xvector<double> conv_lattice_vec_c;
    vector<xvector<double> > conven;
    xvector<double> null;
    null[1] = null[2] = null[3] = 0.0;

    // store hexagonal c-axis
    xvector<double> hex_c = candidate_lattice_vectors[0](3);

    // use moduli of hex a and c to determine rhl a and alpha
    double mod_hex_a=0.0; double mod_hex_c=0.0; double mod_aprime=0.0;
    mod_hex_a = aurostd::modulus(candidate_lattice_vectors[0](1));
    mod_hex_c = aurostd::modulus(candidate_lattice_vectors[0](3));
    mod_aprime = std::sqrt((mod_hex_a*mod_hex_a/3.0) + (mod_hex_c*mod_hex_c/9.0));

    double alpha_prime = std::acos(((2.0*mod_hex_c*mod_hex_c) - (3.0*mod_hex_a*mod_hex_a))/(2.0*(mod_hex_c*mod_hex_c + 3.0*mod_hex_a*mod_hex_a)));

    if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Modulus of rhombohedral lattice vectors (a prime) is " << mod_aprime << "." << endl; }
    if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Angle between all rhombohedral lattice vectors (alpha prime) is " << alpha_prime << "." << endl; }

    // find lattice vectors that are of size mod_aprime
    vector<xvector<double> > candidate_rhl_setting_vectors;
    for (uint i = 0; i < big_expanded.size(); i++) {
      if(aurostd::abs(aurostd::modulus(big_expanded[i]) - mod_aprime) < tol){ //DX20190215 - _SYM_TOL_ to tol
        if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Found candidate vector: " << big_expanded[i]  << "." << endl; }
        candidate_rhl_setting_vectors.push_back(big_expanded[i]);
      }
    }

    xvector<double> tmpa, tmpb, tmpc;
    double moda, modb, modc, tmpalpha, tmpbeta, tmpgamma = 0.0;
    double min_norm = 1e9; //DX20180514 - added initialization
    int count = 0;
    for (uint i = 0; i < candidate_rhl_setting_vectors.size(); i++) {
      tmpa = candidate_rhl_setting_vectors[i]; moda = aurostd::modulus(tmpa);
      for (uint j = 0; j < candidate_rhl_setting_vectors.size(); j++) {
        if(j!=i){
          tmpb = candidate_rhl_setting_vectors[j]; modb = aurostd::modulus(tmpb);
          for (uint k = 0; k < candidate_rhl_setting_vectors.size(); k++) {
            if(k!=j && k!=i){
              tmpc = candidate_rhl_setting_vectors[k]; modc = aurostd::modulus(tmpc);
              // check if all vectors are coplanar
              if(aurostd::abs(aurostd::scalar_product(tmpa,aurostd::vector_product(tmpb,tmpc)))>tol){ //DX20190215 - _SYM_TOL_ to tol
                tmpalpha = aurostd::angle(tmpb,tmpc);
                tmpbeta = aurostd::angle(tmpa,tmpc);
                tmpgamma = aurostd::angle(tmpa,tmpb);
                // check if alpha=beta=gamma and != 0
                if(SYM::checkAngle(modb, modc, tmpalpha, alpha_prime, tol) && //DX20190215 - _SYM_TOL_ to tol
                    SYM::checkAngle(moda, modc, tmpbeta, alpha_prime, tol) && //DX20190215 - _SYM_TOL_ to tol
                    SYM::checkAngle(moda, modb, tmpgamma, alpha_prime, tol) && //DX20190215 - _SYM_TOL_ to tol
                    tmpalpha>_ZERO_TOL_ && tmpbeta>_ZERO_TOL_ && tmpgamma>_ZERO_TOL_ &&
                    aurostd::scalar_product(tmpa,hex_c)>tol &&                    // we want to orient the rhl cell to be in the positive direction of hex-c //DX20190215 - _SYM_TOL_ to tol
                    aurostd::scalar_product(tmpb,hex_c)>tol &&                    // we want to orient the rhl cell to be in the positive direction of hex-c //DX20190215 - _SYM_TOL_ to tol
                    aurostd::scalar_product(tmpc,hex_c)>tol) {                    // we want to orient the rhl cell to be in the positive direction of hex-c //DX20190215 - _SYM_TOL_ to tol 
                  count++;
                  if(count > 1) {
                    if(aurostd::modulus(tmpa) <= min_norm) {
                      min_norm = aurostd::modulus(tmpa);
                      conv_lattice_vec_a = tmpa;
                      conv_lattice_vec_b = tmpb;
                      conv_lattice_vec_c = tmpc;
                    }
                  } 
                  else {
                    min_norm = aurostd::modulus(tmpa);
                    conv_lattice_vec_a = tmpa;
                    conv_lattice_vec_b = tmpb;
                    conv_lattice_vec_c = tmpc;
                  }
                }
              }
            }
          }
        }
      }
    }
    CL = xvec2xmat(conv_lattice_vec_a, conv_lattice_vec_b, conv_lattice_vec_c);
    if(aurostd::determinant(CL) < tol) { //DX20190215 - _SYM_TOL_ to tol
      CL = xvec2xmat(conv_lattice_vec_b, conv_lattice_vec_a, conv_lattice_vec_c);
    }
    candidate_lattice_vectors[0]=CL;
    candidate_lattice_chars[0]='r';
    return true; 
  }
}



// ***********************************************************************************************************************
// Determine lattice centering
// ***********************************************************************************************************************
namespace SYM {
  bool determineLatticeCentering(vector<xvector<double> >& bravais_basis, int& bravais_count, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, vector<xvector<double> >& big_expanded, string& crystalsystem, vector<char>& candidate_lattice_chars, double& tol) { //DX20190215 - _SYM_TOL_ to tol
    // Find number of lattice points in cell (lattice centering)

    xmatrix<double> lattice = trasp(f2c); //DX20190619 
    for (uint j = 0; j < big_expanded.size(); j++) {
      if(bravais_basis.size() == 0) {
        bravais_basis.push_back(big_expanded[j]);
      } else {
        bool duplicate_lattice_point = false;
        for (uint a = 0; a < bravais_basis.size(); a++) {
          xvector<double> bravais_fpos = BringInCell(c2f * bravais_basis[a]);
          xvector<double> tmp = BringInCell(c2f * big_expanded[j]);
          //if(MapAtoms(bravais_fpos, tmp, f2c, skew, tol)) //DX20190215 - _SYM_TOL_ to tol
          if(SYM::MapAtom(bravais_fpos,tmp,lattice,f2c,skew,tol)) //DX20190215 - _SYM_TOL_ to tol //DX20190619 - lattice and f2c as input
          { //CO20200106 - patching for auto-indenting
            duplicate_lattice_point = true;
            break;
          }
        }
        if(duplicate_lattice_point == false) {
          bravais_basis.push_back(big_expanded[j]);
        }
      }
    }
    bravais_count = bravais_basis.size();

    // ===== If MONOCLINIC, and two lattice points, then C-centered ===== //
    if(crystalsystem == "MONOCLINIC" && bravais_count == 2) {
      bravais_count = 5;
    }
    // ===== If ORTHORHOMBIC, could be P-,F-,C-,or I-centered ===== //
    else if(crystalsystem == "ORTHORHOMBIC") {
      vector<double> zs;
      vector<double> ys;
      vector<double> xs;
      for (uint a = 0; a < bravais_basis.size(); a++) {
        xvector<double> bravais_direct = c2f * bravais_basis[a];
        zs.push_back(bravais_direct[3]);
        ys.push_back(bravais_direct[2]);
        xs.push_back(bravais_direct[1]);
      }
      // ===== Differentiate between C and I centering ===== //
      //Use Bravais Count = 5 to identify Orthorhombic C:
      if(bravais_count == 2 &&
          (aurostd::abs(zs[0] - zs[1]) < tol || //DX20190215 - _SYM_TOL_ to tol
           aurostd::abs(ys[0] - ys[1]) < tol || //DX20190215 - _SYM_TOL_ to tol
           aurostd::abs(xs[0] - xs[1]) < tol)) { //DX20190215 - _SYM_TOL_ to tol
        // === If two lattice points are on the same plane --> //DX C centering === //
        bravais_count = 5;
      }
    }
    // ===== If RHOMBOHEDRAL and the bravais_count is only 1, then it is actually HEXAGONAL ===== //
    else if(crystalsystem == "RHOMBOHEDRAL" && bravais_count == 1) {
      crystalsystem = "HEXAGONAL";
      for (uint i = 0; i < candidate_lattice_chars.size(); i++) {
        candidate_lattice_chars[i] = 'h';
      }
    } else if(crystalsystem == "HEXAGONAL") {
      if(bravais_count == 3) {
        // Check if HEX has obverse/reverse to know if this is actually rhombohedral
        xvector<double> centering1;
        centering1(1) = 2.0 / 3.0;
        centering1(2) = 1.0 / 3.0;
        centering1(3) = 1.0 / 3.0;
        xvector<double> centering2;
        centering2(1) = 1.0 / 3.0;
        centering2(2) = 2.0 / 3.0;
        centering2(3) = 2.0 / 3.0;
        xvector<double> centering1_swap;
        centering1_swap(1) = 1.0 / 3.0;
        centering1_swap(2) = 2.0 / 3.0;
        centering1_swap(3) = 1.0 / 3.0;
        xvector<double> centering2_swap;
        centering2_swap(1) = 2.0 / 3.0;
        centering2_swap(2) = 1.0 / 3.0;
        centering2_swap(3) = 2.0 / 3.0;
        bool centering1_found = false;
        bool centering2_found = false;
        for (uint a = 0; a < bravais_basis.size(); a++) {
          xvector<double> bravais_fpos = BringInCell(c2f * bravais_basis[a]);
          if(SYM::MapAtom(bravais_fpos, centering1, lattice, f2c, skew, tol) || SYM::MapAtom(bravais_fpos, centering1_swap, lattice, f2c, skew, tol)) { //DX20190215 - _SYM_TOL_ to tol //DX20190619 - lattice and f2c as input
            centering1_found = true;
          } else if(SYM::MapAtom(bravais_fpos, centering2, lattice, f2c, skew, tol) || SYM::MapAtom(bravais_fpos, centering2_swap, lattice, f2c, skew, tol)) { //DX20190215 - _SYM_TOL_ to tol //DX20190619 - lattice and f2c as input
            centering2_found = true;
          }
        }
        if(centering1_found && centering2_found) {
          for (uint i = 0; i < candidate_lattice_chars.size(); i++) {
            candidate_lattice_chars[i] = 'R';
          }
        }
      }
    }
    return true;
  }
} //namespace SYM

// ***********************************************************************************************************************
// Check conventional cell symmetry
// ***********************************************************************************************************************
namespace SYM {
  symfolder check_ccell(xstructure& xstr, SymmetryInformationITC& ITC_sym_info) {
    bool LDEBUG = (FALSE || XHOST.DEBUG);

    symfolder SOps;
    xstructure xstr_CCell = xstr;
    double tol = xstr.sym_eps; //DX20190215 - added tol

    // ===== Symmetry elements variables ===== //
    ITC_sym_info.initsymmats();  //LOOK IN MAIN (CONVENTIONALCELL)
    ITC_sym_info.initglides();
    ITC_sym_info.initsymops();
    xvector<double> tmp;
    tmp(1) = 1;
    //extern vector<xvector<double> > glideplanes;
    //extern vector<xvector<double> > glidetrans;
    //extern vector<string> glidesymbols;
    //extern vector<xmatrix<double> > sym_mats;
    //extern vector<string> symbol;
    //extern vector<string> dirparam;
    //extern hash sym_mats_direction;
    //HEX DATA
    //extern hash sym_mats_direction_hex;
    //extern vector<xmatrix<double> > sym_mats_hex;
    //extern vector<string> symbol_hex;
    //extern vector<string> dirparam_hex;
    //extern vector<xvector<double> > glideplanes_hex;
    //extern vector<xvector<double> > glidetrans_hex;
    //extern vector<string> glidesymbols_hex;
    ////INDEX FOR THE LATTICES
    //extern vector<int> index_cubic;
    //extern vector<int> index_hex;
    //extern vector<int> index_rhom;
    //extern vector<int> index_tetr;
    //extern vector<int> index_ortho;
    //extern vector<int> index_mono_b;
    //extern vector<int> index_mono_c;
    //extern vector<int> index_tric;
    vector<Screw> twofold_;
    vector<Screw> other_;
    vector<Glide> mirrors_;
    vector<xvector<double> > twofold_lattice_vectors;
    vector<xvector<double> > rot_lattice_vectors;
    vector<xvector<double> > mirror_lattice_vectors;

    // ===== Number of symmetry elements ===== //
    int mcount = 0;  //count number of mirror operations to determine crystal lattice
    int twocount = 0;
    int threecount = 0;
    int fourcount = 0;
    int sixcount = 0;
    int ident_trans = 0;  //count number of centering operations (inludes (000))
    vector<string> pointgroupops;
    vector<xvector<double> > mcount_dirs;
    vector<xvector<double> > latticecell;
    xvector<double> a;
    a(1) = 1;
    xvector<double> b;
    b(2) = 1;
    xvector<double> c;
    c(3) = 1;
    latticecell.push_back(a);
    latticecell.push_back(b);
    latticecell.push_back(c);

    deque<_atom> atomicbasis = xstr_CCell.atoms;
    xmatrix<double> L = xstr_CCell.lattice;
    xmatrix<double> f2c = trasp(L);
    xmatrix<double> c2f = inverse(trasp(L));
    bool skew = SYM::isLatticeSkewed(L, xstr_CCell.dist_nn_min, tol); //DX20190215 - _SYM_TOL_ to tol
    xmatrix<double> Linv = aurostd::inverse(L);
    char latticetypechar = xstr_CCell.lattice_label_ITC;
    vector<int> SYMINDEX;
    string latticesystem = "";

    // ===== Use lattice type to determine symmetry elements to check ===== //
    if(latticetypechar == 'c') {  //CUBIC
      SYMINDEX = ITC_sym_info.index_cubic;
    }
    if(latticetypechar == 'r') {  //RHOMBOHEDRAL
      SYMINDEX = ITC_sym_info.index_rhom;
    }
    if(latticetypechar == 't') {  //TETRAGONAL
      SYMINDEX = ITC_sym_info.index_tetr;
    }
    if(latticetypechar == 'o') {  //ORTHORHOMBIC
      SYMINDEX = ITC_sym_info.index_ortho;
    }
    if(latticetypechar == 'b') {  //MONOCLINIC (Unique axis b)
      SYMINDEX = ITC_sym_info.index_mono_b;
    }
    if(latticetypechar == 'm') {  //MONOCLINIC (Unique axis c)
      SYMINDEX = ITC_sym_info.index_mono_c;
    }
    if(latticetypechar == 'a') {  //TRICLINIC
      SYMINDEX = ITC_sym_info.index_tric;
    }
    if(latticetypechar == 'h' || latticetypechar == 'R') {  //HEXAGONAL (RHOMB)
      SYMINDEX = ITC_sym_info.index_hex;
      ITC_sym_info.sym_mats = ITC_sym_info.sym_mats_hex;
      ITC_sym_info.symbol = ITC_sym_info.symbol_hex;
      ITC_sym_info.dirparam = ITC_sym_info.dirparam_hex;
      ITC_sym_info.sym_mats_direction = ITC_sym_info.sym_mats_direction_hex;
      ITC_sym_info.glideplanes = ITC_sym_info.glideplanes_hex;
      ITC_sym_info.glidetrans = ITC_sym_info.glidetrans_hex;
      ITC_sym_info.glidesymbols = ITC_sym_info.glidesymbols_hex;
    }

    // SYMINDEX contains the point groups of the LATTICE.  Below, we check for the point group of the CRYSTAL.
    // We do not want the crystal symmetry to influence the lattice symmetry, thus if the crystal needed to be
    // reformed, we do not store the new point groups.
    SOps.lattice_pgroups = SYMINDEX;
    SOps.lattice_sym_mats = ITC_sym_info.sym_mats;
    SOps.crystal_sym_mats = ITC_sym_info.sym_mats;

    // ===== Determine the smallest group of an atom type to search for internal translations ===== //
    deque<deque<_atom> > atoms_by_type = SYM::break_up_by_type(atomicbasis);
    uint smallest_group = atoms_by_type[0].size();
    uint index_for_smallest_group = 0;
    for (uint i = 1; i < atoms_by_type.size(); i++) {
      if(atoms_by_type[i].size() < smallest_group) {
        smallest_group = atoms_by_type[i].size();
        index_for_smallest_group = i;
      }
    }
    vector<xvector<double> > translations;
    vector<xvector<double> > centeringops;
    vector<int> insym;  //vector of index of symmetry operations found inside cell
    //LOOP OVER POSSIBLE SYMMETRY OPERATIONS //Select smallest group of atom types to search for internal translations

    bool pointgroupmap_success = false;
    vector<string> pointgroup_crystalsystem;
    vector<vector<int> > all_atom_maps;
    vector<vector<int> > all_type_maps;
    // ****************************************************************************************************************************************
    // Check if the symmetry found by the lattice is consistent with the crystal
    // ****************************************************************************************************************************************
    // Applies the symmetry elements consistent with the lattice to each atom in the cell.
    // If a symmetry element is able to map an atom (one-to-one mapping) onto an equivalent atom,
    // then the symmetry element is stored.  Each type of symmetry element is stored and categorized.
    // This is then used to determine the point group symmetry and the lattice type using the
    // "pointgroupmap" variable.  If the lattice type suggested by this section is different than
    // that found in the conventional cell routine, the cell must be reformed.

    // ===== Loop over tolerance to ensure a pointgroupmap is found (otherwise segmentation fault) ===== //
    //CO START
    _sym_op symOp;        //CO
    symOp.is_fgroup = 1;  //CO
    for (uint i = 0; i < atomicbasis.size(); i++) {
      symOp.basis_atoms_map.push_back(0);
      symOp.basis_types_map.push_back(0);
    }
    //CO END
    while (pointgroupmap_success == false) {
      ident_trans = 0;
      // ===== There must be an identity operator otherwise a tolerance problem ===== //
      while (ident_trans == 0) {
        int symcount = 0;

        // ===== Loop over symmetry elements ===== //
        for (uint i = 0; i < SYMINDEX.size(); i++) {
          symcount++;
          xmatrix<double> Rf = ITC_sym_info.sym_mats[SYMINDEX[i]];  //Uf  //CO
          string sym_symbol = ITC_sym_info.symbol[SYMINDEX[i]];
          //cerr << i << " testing " << ITC_sym_info.symbol[SYMINDEX[i]] << endl;
          xvector<double> T;
          vector<int> atom_map;
          vector<int> type_map;
          //CO START
          symOp.Uf = Rf;
          symOp.Uc = f2c * Rf * c2f;
          symOp.str_Hermann_Mauguin = sym_symbol;
          //CO END
          // ===== Use the smallest group of an atom type to find the possible translations ===== //
          for (uint j = 0; j < atoms_by_type[index_for_smallest_group].size(); j++) {
            //DX20190905 [OBSOLETE-no more mod_one_xvec] T = SYM::mod_one_xvec(atoms_by_type[index_for_smallest_group][0].fpos - Rf * atoms_by_type[index_for_smallest_group][j].fpos);  //ftau
            T = atoms_by_type[index_for_smallest_group][0].fpos - Rf * atoms_by_type[index_for_smallest_group][j].fpos;  //ftau //DX20190905
            BringInCellInPlace(T); //DX20190905
            //CO START
            symOp.ftau = T;
            symOp.ctau = f2c * T;  //atoms_by_type[index_for_smallest_group][0].cpos - Rf*atoms_by_type[index_for_smallest_group][j].cpos;//f2c*T;
            //CO END
            //if(SYM::getFullSymBasis(atomicbasis,Rf,c2f,f2c,sym_symbol,T,skew,tol,atom_map,type_map)) //DX20190215 - _SYM_TOL_ to tol
            //DX if(SYM::getFullSymBasis(atomicbasis,L,c2f,f2c,symOp,skew,tol,false,atom_map,type_map)) //DX20190215 - _SYM_TOL_ to tol
            if(SYM::getFullSymBasis(atomicbasis, L, c2f, f2c, symOp, TRUE, skew, tol, atom_map, type_map)) //DX20190215 - _SYM_TOL_ to tol
            { //CO20200106 - patching for auto-indenting
              //cerr << "---------" << endl;
              //cerr << "Uf " << endl << symOp.Uf << endl;
              //cerr << "FTAU: " << symOp.ftau << endl;
              //cerr << "---------" << endl;
              // ===== If one-to-one, store the symmetry operator ===== //
              //cerr << " ===> storing " << ITC_sym_info.symbol[SYMINDEX[i]] << endl;
              //print(atom_map);
              insym.push_back(SYMINDEX[i]);
              translations.push_back(T);
              all_atom_maps.push_back(atom_map);
              all_type_maps.push_back(type_map);
              atom_map.clear();
              type_map.clear();

              // ===== Categorize by type of symmetry element ===== //
              if(ITC_sym_info.symbol[SYMINDEX[i]] == "1") {
                ident_trans++;
                centeringops.push_back(T);
                pointgroupops.push_back("1");
              } else if(ITC_sym_info.symbol[SYMINDEX[i]] == "m") {
                mcount++;
                mcount_dirs.push_back(ITC_sym_info.sym_mats_direction[SYMINDEX[i] + 1]);
                pointgroupops.push_back("m");
              } else if(ITC_sym_info.symbol[SYMINDEX[i]] == "2") {
                twocount++;
                pointgroupops.push_back(ITC_sym_info.symbol[SYMINDEX[i]]);
              } else if(havechar(ITC_sym_info.symbol[SYMINDEX[i]], '3')) {
                threecount++;
                pointgroupops.push_back(ITC_sym_info.symbol[SYMINDEX[i]]);
              } else if(havechar(ITC_sym_info.symbol[SYMINDEX[i]], '4')) {
                fourcount++;
                pointgroupops.push_back(ITC_sym_info.symbol[SYMINDEX[i]]);
              } else if(havechar(ITC_sym_info.symbol[SYMINDEX[i]], '6')) {
                sixcount++;
                pointgroupops.push_back(ITC_sym_info.symbol[SYMINDEX[i]]);
              } else {
                pointgroupops.push_back(ITC_sym_info.symbol[SYMINDEX[i]]);
              }
            } else {
              atom_map.clear();
              type_map.clear();
            }
          }
        }

        if(ident_trans == 0) {
          if(LDEBUG) { cerr << __AFLOW_FUNC__ << " ERROR: Could not find identity." << endl; }
          SOps.latticesystem = "REDO";
          SOps.commensurate = false;
          return SOps;
        }
      }

      // ===== Scale number of operators based on number of identity operators ===== //
      twocount = twocount / ident_trans;
      threecount = threecount / ident_trans;
      fourcount = fourcount / ident_trans;
      sixcount = sixcount / ident_trans;
      mcount = mcount / ident_trans;  //divide by centering operation number to get appropriate mirror count

      vector<string> tmpsecondpgops;
      for (uint i = 0; i < insym.size(); i++) {
        if(i % ident_trans == 0) {
          tmpsecondpgops.push_back(ITC_sym_info.symbol[insym[i]]);
        }
      }
      //print(tmpsecondpgops);

      // ===== Determine the pointgroup and crystal system based on the symmetry operations found ===== //
      pointgroup_crystalsystem = splitstringbyspaces(crystal_system_library_search(tmpsecondpgops));

      // === If not found, change tolerances and redo symmetry search (otherwise, segmentation fault) === //
      if(pointgroup_crystalsystem[0] == "REDO") {
        if(LDEBUG) { cerr << __AFLOW_FUNC__ << " WARNING: Point group mapping failed (number of symmetry elements does not match with a given point group) [dir=" << xstr.directory << "]." << endl; }
        SOps.latticesystem = "REDO";
        SOps.commensurate = false;
        return SOps;
      } else {
        pointgroupmap_success = true;
      }
    }  // End of the pointgroupmap_success loop

    //cout << "Point group & crystal system: " << pointgroup_crystalsystem << endl;
    cleanupstring(pointgroup_crystalsystem[0]);
    cleanupstring(pointgroup_crystalsystem[1]);

    // Given the information above, can get the Crystal System and then assign an appropriate
    // mcount for the conventionalization routine: (must account for centerings)
    ostringstream crystalsystem;

    crystalsystem << pointgroup_crystalsystem[1];

    //NOTE: no need to reform cell if difference is trigonal <-> hexagonal, but should change crystal system accordingly
    if(crystalsystem.str() == "TRIGONAL" && xstr_CCell.crystal_system_ITC == "HEXAGONAL") {  // ADDITION: 20130515
      xstr_CCell.crystal_system_ITC = "TRIGONAL";
    }

    SOps.crystalsystem = pointgroup_crystalsystem[1];
    //print(insym);

    // === DIVIDE SYMMETRY BY CENTERING OPERATIONS === //
    uint centeringgroups = 1;
    if(xstr_CCell.bravais_label_ITC == 'I' || xstr_CCell.bravais_label_ITC == 'C') {
      centeringgroups = 2;
    }
    if(xstr_CCell.bravais_label_ITC == 'F') {
      centeringgroups = 4;
    }
    if(xstr_CCell.bravais_label_ITC == 'R') {
      centeringgroups = 3;
    }

    //cerr << "centeringops.size(): " << centeringops.size() << endl;
    //cerr << "centeringgroups: " << centeringgroups << endl;
    // === Ensure the centering operators (translations) are consistent with the lattice centering (# lattice points) === //
    if(centeringops.size() != centeringgroups) {
      // === If not restart from the beginning, changing the tolerance === //
      if(LDEBUG) { 
        cerr << __AFLOW_FUNC__ << " WARNING: Centering operations incommensurate... [dir=" << xstr.directory << "]." << endl;
      }
      SOps.latticesystem = "REDO";
      SOps.commensurate = false;
      return SOps;
    }

    // ===== If the crystal lattice system is not commensurate with the original lattice symmetry, prepare to reform ===== //
    if(xstr_CCell.crystal_system_ITC != crystalsystem.str()) {
      //cerr << "CCell not equal to crystalsystem.str(): " << CCell.CrystalSystem << " != " << crystalsystem.str() << endl;
      for (uint i = 0; i < insym.size(); i++) {
        xvector<double> origin;
        if(ITC_sym_info.symbol[insym[i]] == "m") {
          Glide mirror_tmp;
          mirror_tmp.get_glide_direct(ITC_sym_info.sym_mats_direction[insym[i] + 1] * L, origin);
          mirror_lattice_vectors.push_back(ITC_sym_info.sym_mats_direction[insym[i] + 1] * L);
          mirrors_.push_back(mirror_tmp);
        } else {
          Screw screw_tmp;
          if(havechar(ITC_sym_info.symbol[insym[i]], '2')) {
            screw_tmp.get_screw_direct(ITC_sym_info.sym_mats_direction[insym[i] + 1] * L, origin, 2);
            twofold_lattice_vectors.push_back(ITC_sym_info.sym_mats_direction[insym[i] + 1] * L);
            twofold_.push_back(screw_tmp);
          }
          if(havechar(ITC_sym_info.symbol[insym[i]], '3')) {
            screw_tmp.get_screw_direct(ITC_sym_info.sym_mats_direction[insym[i] + 1] * L, origin, 3);
            rot_lattice_vectors.push_back(ITC_sym_info.sym_mats_direction[insym[i] + 1] * L);
            other_.push_back(screw_tmp);
          }
          if(havechar(ITC_sym_info.symbol[insym[i]], '4')) {
            screw_tmp.get_screw_direct(ITC_sym_info.sym_mats_direction[insym[i] + 1] * L, origin, 4);
            rot_lattice_vectors.push_back(ITC_sym_info.sym_mats_direction[insym[i] + 1] * L);
            other_.push_back(screw_tmp);
          }
          if(havechar(ITC_sym_info.symbol[insym[i]], '6')) {
            screw_tmp.get_screw_direct(ITC_sym_info.sym_mats_direction[insym[i] + 1] * L, origin, 6);
            rot_lattice_vectors.push_back(ITC_sym_info.sym_mats_direction[insym[i] + 1] * L);
            other_.push_back(screw_tmp);
          }
        }
      }
      SOps.latticesystem = latticesystem;
      SOps.crystalsystem = crystalsystem.str();
      SOps.commensurate = false;
      SOps.twofold_ops = twofold_;
      SOps.rot_ops = other_;
      SOps.mirror_ops = mirrors_;
      SOps.twofold_lattice_vectors = twofold_lattice_vectors;
      SOps.rot_lattice_vectors = rot_lattice_vectors;
      SOps.mirror_lattice_vectors = mirror_lattice_vectors;
    } else {
      SOps.insym = insym;
      SOps.translations = translations;
      SOps.all_atom_maps = all_atom_maps;
      SOps.all_type_maps = all_type_maps;
      SOps.centeringops = centeringops;
      SOps.pointgroupops = pointgroupops;
      SOps.commensurate = true;
      //cerr << "Original lattice symmetry is commensurate with crystal lattice symmetry" << endl;
    }
    return SOps;
  }
} //namespace SYM

// ********************************************************************************************************************************
// calculateSpaceGroups()
// ********************************************************************************************************************************
namespace SYM {
  void calculateSpaceGroups(vector<xstructure>& vxstrs, uint start_index, uint end_index, uint setting){ //DX20191230 add setting option

    bool no_scan = false;

    // if end index is default (i.e., AUROSTD_MAX_UINT), then compute symmetry analysis for all structures
    if(end_index == AUROSTD_MAX_UINT){ end_index=vxstrs.size(); }

    for(uint i=start_index;i<end_index;i++){ //DX20191108 - [switching from compare::splitTaskIntoThreads to getThreadDistribution, end index conventions differ, <= vs <]
      double use_tol = SYM::defaultTolerance(vxstrs[i]);
      vxstrs[i].SpaceGroup_ITC(use_tol, -1, setting, no_scan);
    }
  }
}

// ********************************************************************************************************************************
// AFLOW2SG:: Prints space group in string format (Void Input)
// ********************************************************************************************************************************
string xstructure::aflow2sg(void) {
  //DX20170905 - [OBSOLETE] int change_sym_count = 1;
  double use_tol = SYM::defaultTolerance((*this));
  //DX20170905 - [OBSOLETE] double orig_tolerance = use_tol;
  bool no_scan = false;
  //DX20170905 - [OBSOLETE] uint sgroup = (*this).SpaceGroup_ITC(use_tol, orig_tolerance, -1, change_sym_count, no_scan);
  uint sgroup = (*this).SpaceGroup_ITC(use_tol, -1, 0, no_scan); //DX20180806 - added setting; 0 means default
  return GetSpaceGroupName(sgroup,(*this).directory) + " #" + aurostd::utype2string(sgroup); //DX20180526 - add directory
}

// ********************************************************************************************************************************
// AFLOW2SG:: Prints space group in string format (Tolerance Input)
// ********************************************************************************************************************************
string xstructure::aflow2sg(double& use_tol) {
  //DX20170905 - [OBSOLETE] double orig_tolerance = use_tol;
  //DX20170905 - [OBSOLETE] int change_sym_count = 1;
  bool no_scan = false;
  //DX20170905 - [OBSOLETE] uint sgroup = (*this).SpaceGroup_ITC(use_tol, orig_tolerance, -1, change_sym_count, no_scan);
  uint sgroup = (*this).SpaceGroup_ITC(use_tol, -1, 0, no_scan); //DX20180806 - added setting; 0 means default
  return GetSpaceGroupName(sgroup,(*this).directory) + " #" + aurostd::utype2string(sgroup); //DX20180526 - add directory
}

// ********************************************************************************************************************************
// AFLOW2SG:: Prints space group in string format (Tolerance Input and Manual Input)
// ********************************************************************************************************************************
string xstructure::aflow2sg(double& use_tol, const int& manual_it) {
  //DX20170905 - [OBSOLETE] double orig_tolerance = use_tol;
  //DX20170905 - [OBSOLETE] int change_sym_count = 1;
  bool no_scan = false;
  //DX20170905 - [OBSOLETE] uint sgroup = (*this).SpaceGroup_ITC(use_tol, orig_tolerance, manual_it, change_sym_count, no_scan);
  uint sgroup = (*this).SpaceGroup_ITC(use_tol, manual_it, 0, no_scan); //DX20180806 - added setting; 0 means default
  return GetSpaceGroupName(sgroup,(*this).directory) + " #" + aurostd::utype2string(sgroup); //DX20180526 - add directory
}

// ********************************************************************************************************************************
// FUNCTION OVERLOADING: Void input
// ********************************************************************************************************************************
uint xstructure::SpaceGroup_ITC(void) {
  double use_tol = SYM::defaultTolerance((*this));
  //DX20170905 - [OBSOLETE] double orig_tolerance = use_tol;
  //DX20170905 - [OBSOLETE] int change_sym_count = 1;
  bool no_scan = false;
  //DX20170905 - [OBSOLETE] return (*this).SpaceGroup_ITC(use_tol, orig_tolerance, -1, change_sym_count, no_scan);
  return (*this).SpaceGroup_ITC(use_tol, -1, 0, no_scan); //DX20180806 - added setting; 0 means default
}

// ********************************************************************************************************************************
// FUNCTION OVERLOADING: Tolerance Input
// ********************************************************************************************************************************
uint xstructure::SpaceGroup_ITC(double& use_tol) {
  //DX20170905 - [OBSOLETE] double orig_tolerance = use_tol;
  //DX20170905 - [OBSOLETE] int change_sym_count = 1;
  bool no_scan = false;
  //DX20170905 - [OBSOLETE] return (*this).SpaceGroup_ITC(use_tol, orig_tolerance, -1, change_sym_count, no_scan);
  return (*this).SpaceGroup_ITC(use_tol, -1, 0, no_scan); //DX20180806 - added setting; 0 means default
}

// ********************************************************************************************************************************
// FUNCTION OVERLOADING: Tolerance Input
// ********************************************************************************************************************************
uint xstructure::SpaceGroup_ITC(double& use_tol, bool& no_scan) {
  //DX20170905 - [OBSOLETE] double orig_tolerance = use_tol;
  //DX20170905 - [OBSOLETE] int change_sym_count = 1;
  //DX20170905 - [OBSOLETE] return (*this).SpaceGroup_ITC(use_tol, orig_tolerance, -1, change_sym_count, no_scan);
  return (*this).SpaceGroup_ITC(use_tol, -1, 0, no_scan); //DX20180806 - added setting; 0 means default
}

//// ********************************************************************************************************************************
//// FUNCTION OVERLOADING: Tolerance Input
//// ********************************************************************************************************************************
//uint xstructure::SpaceGroup_ITC(double& use_tol, const int& manual_it) {
//  //DX20170905 - [OBSOLETE] double orig_tolerance = use_tol;
//  //DX20170905 - [OBSOLETE] int change_sym_count = 1;
//  bool no_scan = false;
//  //DX20170905 - [OBSOLETE] return (*this).SpaceGroup_ITC(use_tol, orig_tolerance, manual_it, change_sym_count, no_scan);
//  return (*this).SpaceGroup_ITC(use_tol, manual_it, 0, no_scan); //DX20180806 - added setting; 0 means default
//}

// ********************************************************************************************************************************
// FUNCTION OVERLOADING: Tolerance Input
// ********************************************************************************************************************************
uint xstructure::SpaceGroup_ITC(double& use_tol, const int& setting) {
  bool no_scan = false;
  return (*this).SpaceGroup_ITC(use_tol, -1, setting, no_scan); //DX20180806 - added setting
}

//// ********************************************************************************************************************************
//// FUNCTION OVERLOADING: Tolerance Input
//// ********************************************************************************************************************************
//uint xstructure::SpaceGroup_ITC(double& use_tol, const int& manual_it, bool& no_scan) {
//  return (*this).SpaceGroup_ITC(use_tol, manual_it, 0, no_scan); //DX20180806 - added setting; 0 means default
//}

// ********************************************************************************************************************************
// MAIN FUNCTION:: Space group consistent with the International Tables of Crystallography
// ********************************************************************************************************************************
//DX20170905 - [OBSOLETE] uint xstructure::SpaceGroup_ITC(double& use_tol, double& orig_tolerance, const int& manual_it, int& change_sym_count, bool& no_scan)
uint xstructure::SpaceGroup_ITC(double& use_tol, const int& manual_it, const int& setting, bool& no_scan) //DX20180806 - added setting
{ //CO20200106 - patching for auto-indenting
  bool LDEBUG = (FALSE || XHOST.DEBUG);
  stringstream message;
  if(use_tol < _ZERO_TOL_) {
    cerr << __AFLOW_FUNC__ << " ERROR: Tolerance cannot be zero (i.e. less than 1e-10) [dir=" << (*this).directory << "]." << endl;
    return 0;
  }

  double min_dist = 100.0;
  if((*this).dist_nn_min == AUROSTD_NAN) {
    (*this).dist_nn_min = SYM::minimumDistance((*this));
  }
  min_dist = (*this).dist_nn_min;
  if(use_tol > min_dist) {
    cerr << __AFLOW_FUNC__ << " ERROR: The tolerance cannot be larger than the minimum interatomic distance (" << min_dist << " Angstroms). [dir= " << (*this).directory << "]" << endl;
    return 0;
  }

  // If iomode is WYCKCAR, convert to IOVASP temporarily, otherwise, printing of ITC lattice parameters/coordinates
  // gets messed up (this is just a printing problem when debugging)
  // DX20221128
  int iomode_tmp = (*this).iomode;
  if((*this).iomode==IOVASP_WYCKCAR) { iomode_tmp = (*this).iomode; (*this).iomode=IOVASP_POSCAR; }

  // ===== Initialize symmetry elements and set tolerance ===== //
  SymmetryInformationITC ITC_sym_info; //DX20190215
  ITC_sym_info.initsymmats(); //DX20190215
  ITC_sym_info.initglides(); //DX20190215
  ITC_sym_info.initsymops(); //DX20190215
  //DX20190215 [OBSOLETE] SYM::SetTolerance(use_tol);
  (*this).sym_eps = use_tol; //DX20190215 - replace SYM::SetTolerance() 

  //For running in manual mode (crystal cell axes mode is specified by hand)
  int MAN_IT = 0;
  if(manual_it != -1) {
    MAN_IT = manual_it;
  }

  xstructure xstr_orig;
  int iterate = -1;
  xstructure CrystOut_prev;

  xstructure xstr = (*this);

  // ===== Scale the lattice vectors at the beginning ===== //
  xstr.ReScale(1.0);

  // ===== If atoms not labeled, generate fake names ===== //
  bool generatefakenames = true;
  uint num_species_label = 0;
  for (uint i = 0; i < xstr.atoms.size(); i++) {
    if(xstr.atoms.at(i).name.size() != 0) {
      generatefakenames = false;
      num_species_label++;
    }
  }
  if(num_species_label != xstr.atoms.size() && num_species_label != 0) {
    cerr << __AFLOW_FUNC__ << " POSCAR ERROR (Some atoms labeled, others are not ... please label all atoms) [dir= " << (*this).directory << "]." << endl;
    return 0;
  }
  // If iatoms output, assign fake names (useful for wyccar output)
  if(xstr.atoms.at(0).name.find("*") != std::string::npos) {
    generatefakenames = true;
  }
  if(generatefakenames == true) {
    if(LDEBUG) { cerr << "SYM::SpaceGroup_ITC: Atoms in POSCAR not labeled ... assigning fake names [dir=" << (*this).directory << "]." << endl; }
    xstr.DecorateWithFakeElements();
  }
  xstr_orig = xstr;
  xstr.MinkowskiBasisReduction();

  // ===== Declare variables ===== //
  ostringstream woss, woss_min;  //stringstream for wyccar //DX+KC20221115 - need to store minimum wyccar stringstream for aflow --wyckoff --print=json
  int spacegroup = 0, spacegroup_global = 0; //DX20220920
  bool foundspacegroup = false, foundspacegroup_global = false; //DX20220827 - added foundspacegroup_global
  deque<_atom> atomicbasis;
  vector<int> insym;  //vector of index of symmetry operations found inside cell
  vector<xvector<double> > translations;
  vector<vector<int> > all_atom_maps;
  vector<vector<int> > all_type_maps;
  string pointgroup = "";
  //DX20190215 [OBSOLETE]using SYM::generators;
  //DX20190215 [OBSOLETE]using SYM::sgindex;
  //DX20190215 [OBSOLETE]extern vector<vector<SYM::symop> > generators;
  //DX20190215 [OBSOLETE]extern vector<int> sgindex;
  ITC_sym_info.initgenerators(""); //DX20190215
  vector<int> reduction_by_generator_set;
  vector<vector<vector<double> > > RHS_shifts;
  xvector<double> OriginShift, OriginShift_Wyckoff_min;
  //DX20190215 [OBSOLETE]using SYM::sym_mats;
  //DX20190215 [OBSOLETE]using SYM::symbol;
  //DX20190215 [OBSOLETE]using SYM::sym_mats_direction;
  //DX20190215 [OBSOLETE]extern vector<xmatrix<double> > sym_mats;
  //DX20190215 [OBSOLETE]extern vector<string> symbol;
  //DX20190215 [OBSOLETE]extern SYM::hash sym_mats_direction;
  xstructure CCell, CCell_Wyckoff_min; //DX20220827 - added Wyckoff_min
  char latticetypechar = '\0'; //DX20180514 - added initialization

  int cell_choice = setting; //DX20180806 - use setting instead of 0

  xmatrix<double> Ident;
  Ident(1, 1) = 1;
  Ident(2, 2) = 1;
  Ident(3, 3) = 1;

  // === Store operators from previous iteration (sometimes necessary in conventional cell reformation === //
  vector<xmatrix<double> > candidate_lattice_vectors_prev;
  vector<char> candidate_lattice_chars_prev;

  vector<int> lattice_pgroups;
  vector<xmatrix<double> > lattice_sym_mats;
  vector<xmatrix<double> > crystal_sym_mats;
  vector<int> reduced_insym;
  vector<int> sg_search;

  //WYCKOFF VARIABLE DECLARATIONS
  //  bool found_wyckoff = true;
  //  bool first_wyckoff = true;
  bool found_all_wyckoff = false;
  string spacegroupstring = "";
                  
  int min_wyckoff_sum = 1e9, minimum_scheme_sum = 1e9, global_min_wyckoff_sum = 1e9; //DX20220827 - track minimum Wyckoff sum
  vector<int> min_Wyckoff_numbers, global_min_Wyckoff_numbers; //DX20221103
  vector<int> wyckoffmult, wyckoffmult_min;                   // Wyckoff multiplicity //DX20220827 - added min variant
  vector<string> general_wyckoff_position, general_wyckoff_position_min;   // General wyckoff positions //DX20220827 - added min variant
  vector<string> wyckoffsymbols, wyckoffsymbols_min;             // Wyckoff symbol //DX20220827 - added min variant
  deque<_atom> wyckoffPositionsVector, wyckoffPositionsVector_min;       // Stored Wyckoff coordinate //DX20220827 - added min variant
  vector<string> wyckoffSymbols, wyckoffSymbols_min;             // Stored Wyckoff symbol //DX20220827 - added min variant
  vector<wyckoffsite_ITC> wyckoffVariables, wyckoffVariables_min;  //variable values + atom chemical label + wyckoff symbol //DX20220827 - added min variant
  vector<vector<vector<string> > > tmpvvvstring, tmpvvvstring_min; //DX20220827 - added min variant

  // === Tolerance boolean === //
  bool last_orientation = false;
  bool first_run_or_new_tol = true;

  bool lattice_reformed = false;
  bool symmetryfound = false;
  bool tested_all_lattices = false;

  string crystalsystem_prev = "";  //Correct conventional cell

  double sym_eps = xstr.sym_eps; //DX20190314 - added sym eps temp variable
  uint sym_eps_change_count = 0; //DX20180226 - added sym eps change count

  //********************************************************************************************************************************************
  // Space Group Loop
  //********************************************************************************************************************************************
  // This loop is used to determine the space group which is consistent with the ITC (generators/origin shift/Wyckoff positions)
  while (foundspacegroup == false) {
    uint gencentops = 0;
    vector<string> pointgroupops;
    uint centeringgroups = 1;
      

    while(!tested_all_lattices) {
      if(LDEBUG) { cerr << "SYM::SpaceGroup_ITC: [BEGINNING OF LOOP] iterate=" << iterate << " : min_wyckoff_sum=" << min_wyckoff_sum << " vs global=" << global_min_wyckoff_sum << endl; }
      ///// PUT IN STANDARD PRIMITIVE FORM /////
      // ========== Make the xstructure in primitive form if this is the first iteration, or the tolerance has been changed ========== //
      if(first_run_or_new_tol == true) {
        sym_eps = xstr.sym_eps; //DX20180226 - added sym eps temp variable
        sym_eps_change_count = xstr.sym_eps_change_count; //DX20180226 - added sym eps change count
        no_scan = xstr.sym_eps_no_scan; //DX20210430 - added no_scan
        xstr = xstr_orig;
        xstr.sym_eps = sym_eps; //DX20190314 - need to update sym_eps, since the structure was overwritten 
        xstr.sym_eps_change_count = sym_eps_change_count; //DX20180423 - need to update change count, since the structure was overwritten 
        xstr.sym_eps_no_scan = no_scan; //DX20210430 - added no_scan
        CCell.free();
        xstr.MinkowskiBasisReduction();
        //DX20210401 [OBSOLETE] xstr.GetPrimitiveCell();
        xstr.GetPrimitive(); //DX20210401

        // ===== Ensure Cartesian coordinates are consistent with fractional coordinates ===== //
        for (uint i = 0; i < xstr.atoms.size(); i++) {
          xstr.atoms[i].cpos = F2C(xstr.lattice, xstr.atoms[i].fpos);
        }
        xstr.BringInCell();
        first_run_or_new_tol = false;
      }

      // ===== If the last possible orientation of the conventional cell was attempted, and was incommensurate, change tolerance ===== //
      if(last_orientation == true) {
        // If we found the space group in an earlier iteration
        if(foundspacegroup_global) {
          if(LDEBUG) { cerr << "SYM::SpaceGroup_ITC: Explored all iterations, use the lattice and origin choice that found a valid space group." << endl; }
          foundspacegroup = true;
          spacegroup = spacegroup_global;
          CCell = CCell_Wyckoff_min;
          // update Wyckoff info
          wyckoffmult = wyckoffmult_min;
          wyckoffsymbols = wyckoffsymbols_min;
          wyckoffPositionsVector = wyckoffPositionsVector_min;
          wyckoffSymbols = wyckoffSymbols_min;
          wyckoffVariables = wyckoffVariables_min;
          tmpvvvstring = tmpvvvstring_min;
          woss.str(""); woss << woss_min.str(); //DX+KC20221103 - added woss_min for aflow --wyckoff --print=json
          // update shifts
          OriginShift = OriginShift_Wyckoff_min;
          break;
        }
        lattice_reformed = false;
        lattice_pgroups.clear();
        lattice_sym_mats.clear();
        crystal_sym_mats.clear();
        if(!no_scan){
          SYM::change_tolerance(xstr,xstr.sym_eps, min_dist, no_scan); //DX20190215 - _SYM_TOL_ to xstr.sym_eps
          sym_eps_change_count = xstr.sym_eps_change_count; //DX20180226 - added sym eps change count
          xstr.sym_eps_no_scan = no_scan; //DX20210430 - added no_scan
        }
        else {
          (*this).sym_eps_no_scan = no_scan; //DX20210430 - added no_scan
          return 0;
        }
        first_run_or_new_tol = true;
        iterate = -1;
        cell_choice = setting; //DX20180806 - use setting instead of 0
        last_orientation = false;
        symmetryfound = false;
        continue;
      }
      iterate++;
      insym.clear();
      all_atom_maps.clear();
      all_type_maps.clear();
      reduction_by_generator_set.clear();
      RHS_shifts.clear();
      OriginShift.clear();
      atomicbasis.clear();
      translations.clear();
      pointgroup = "";

      SYM::symfolder checkops;

      ///// FIND CONVENTIONAL CELL /////
      //manual or not manual mode
      if(MAN_IT != -1) {
        int IT = iterate;
        // ===== Conventional Cell Routine ===== //
        CCell = SYM::ConventionalCell(xstr, IT, cell_choice, last_orientation, crystalsystem_prev, CrystOut_prev, candidate_lattice_vectors_prev,
            candidate_lattice_chars_prev, checkops, ITC_sym_info, lattice_reformed, lattice_pgroups, lattice_sym_mats, crystal_sym_mats,
            symmetryfound, foundspacegroup_global); //DX20220927 - added foundspacegroup_global
        // ===== If the conventional cell was not found, change the tolerance ===== //
        if(symmetryfound == false) {
          if(foundspacegroup_global && last_orientation) {
          if(LDEBUG) { cerr << "SYM::SpaceGroup_ITC: Last lattice choice/orientation was explored. Consistent symmetry was found with a different choice/orientation, so we will use that conventional cell." << endl; }
            continue;
        }
          lattice_reformed = false;
          lattice_pgroups.clear();
          lattice_sym_mats.clear();
          crystal_sym_mats.clear();
          if(!no_scan){
            SYM::change_tolerance(xstr,xstr.sym_eps, min_dist, no_scan); //DX20190215 - _SYM_TOL_ to xstr.sym_eps
            sym_eps_change_count = xstr.sym_eps_change_count; //DX20180226 - added sym eps change count
            xstr.sym_eps_no_scan = no_scan; //DX20180226 - added sym eps change count
          }
          else {
            (*this).sym_eps_no_scan = no_scan; //DX20210430 - added no_scan
            return 0;
          }
          first_run_or_new_tol = true;
          iterate = -1;
          continue;
        }
        iterate = IT;
      } else {
        iterate = MAN_IT;
        // ===== Conventional Cell Routine ===== //
        CCell = SYM::ConventionalCell(xstr, iterate, cell_choice, last_orientation, crystalsystem_prev, CrystOut_prev, candidate_lattice_vectors_prev,
            candidate_lattice_chars_prev, checkops, ITC_sym_info, lattice_reformed, lattice_pgroups, lattice_sym_mats, crystal_sym_mats,
            symmetryfound, foundspacegroup_global); //DX20220927 - added foundspacegroup_global
        // ===== If the conventional cell was not found, change the tolerance ===== //
        if(symmetryfound == false) {
          lattice_reformed = false;
          lattice_pgroups.clear();
          lattice_sym_mats.clear();
          crystal_sym_mats.clear();
          if(!no_scan){
            SYM::change_tolerance(xstr,xstr.sym_eps, min_dist, no_scan); //DX20190215 - _SYM_TOL_ to xstr.sym_eps
            sym_eps_change_count = xstr.sym_eps_change_count; //DX20180226 - added sym eps change count
            xstr.sym_eps_no_scan = no_scan; //DX20210430 - added no_scan
          }
          else {
            (*this).sym_eps_no_scan = no_scan; //DX20210430 - added no_scan
            return 0;
          }
          first_run_or_new_tol = true;
          iterate = -1;
          continue;
        }
      }
      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " Conventional cell: " << endl;
        cerr << CCell << endl;
      }


      // === Find ratio of lattice vectors === //
      xmatrix<double> f2c = trasp(CCell.lattice);
      xmatrix<double> c2f = inverse(trasp(CCell.lattice));

      atomicbasis = CCell.atoms;

      // === Obtain symmetry elements from the symmetry check (checkops) from the Conventional Cell routine === //
      insym = checkops.insym;
      vector<xvector<double> > centeringops = checkops.centeringops;
      pointgroupops = checkops.pointgroupops;
      translations = checkops.translations;
      all_atom_maps = checkops.all_atom_maps;
      all_type_maps = checkops.all_type_maps;

      vector<xmatrix<double> > pointgroup_mats;
      vector<string> secondpgops;
      reduced_insym.clear();
      for (uint i = 0; i < insym.size(); i++) {
        pointgroup_mats.push_back(ITC_sym_info.sym_mats[insym[i]]);
        if(i % centeringops.size() == 0) {
          secondpgops.push_back(ITC_sym_info.symbol[insym[i]]);
          reduced_insym.push_back(insym[i]);
        }
      }

      // ===== Determine the pointgroup and crystal system based on the symmetry operations found ===== //
      pointgroup = SYM::point_group_library_search(secondpgops, CCell.crystal_system_ITC, centeringops.size());
      sg_search.clear();

      // ===== Use the point group and the lattice type/lattice centering to narrow down the space group ===== //
      sg_search = SYM::PointGroup_SpaceGroup(pointgroup, CCell.bravais_label_ITC);
      //std::reverse(sg_search.begin(), sg_search.end()); //DX20221230 - reverse to give preference to higher space groups
      // If we are looking at the rhombohedral setting, we need to use bravais label 'R' instead of 'P' //DX20180807 
      //DX20190130 - add anrl cell choice : if(cell_choice == 1 && CCell.lattice_label_ITC == 'r') //DX20180807
      if((cell_choice == SG_SETTING_1 || cell_choice == SG_SETTING_ANRL) && CCell.lattice_label_ITC == 'r') //DX20190131 - add anrl setting
      { //CO20200106 - patching for auto-indenting
        sg_search = SYM::PointGroup_SpaceGroup(pointgroup, 'R'); //DX20180807
      } //DX20180807
      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " Possible space groups based on point group (PG=" << pointgroup << "): ";
        for(uint s=0;s<sg_search.size();s++){
          cerr << sg_search[s] << " "; 
        }
        cerr << "[dir=" << xstr.directory << "]." << endl;
      }
      //DEBUG print(sg_search);
      //DEBUG xb();
      // === DIVIDE SYMMETRY BY CENTERING OPERATIONS === //
      gencentops = 0;
      centeringgroups = 1;
      if(CCell.bravais_label_ITC == 'I' || CCell.bravais_label_ITC == 'C') {
        gencentops = 1;
        centeringgroups = 2;
      }
      if(CCell.bravais_label_ITC == 'F') {
        gencentops = 2;
        centeringgroups = 4;
      }
      if(CCell.bravais_label_ITC == 'R') {
        gencentops = 1;
        centeringgroups = 3;
      }

      latticetypechar = CCell.lattice_label_ITC;

      // *********************************************************************************************************************************
      // Find the correct origin with respect to the International Tables of Crystallography (ITC)
      // *********************************************************************************************************************************
      vector<int> FCG;  //first_centered_group;
      FCG.clear();

      // DEBUG
      //for(uint i=0;i<pointgroupops.size();i++){
      //   cerr << pointgroupops[i] << endl;
      //}
      //for(uint i=0;i<insym.size();i++){
      //   cerr << insym[i] << endl;
      //}
      // DEBUG

      if(gencentops == 0) {  //No centeringops means that FCG should simply be an index for finalsymbolset
        for (uint i = 0; i < pointgroupops.size(); i++) {
          FCG.push_back(i);
          //cerr << "FCG.at(i): " << FCG.at(i) << endl;
        }
      } else {
        //CANT DO PRIMITIVE SET WITHOUT KNOWING ORIGIN (?) (SEE BRUTE FORCE METHOD AT END OF LOOP)

        // === Initialize containers for centered groups. === //
        vector<int> tmp;
        for (uint i = 0; i < pointgroupops.size(); i++) {
          //METHOD 1 (PRIMITIVE CELL):
          //atom tmp;
          //tmp.coord = mod_one_xvec(ITC_sym_info.sym_mats[insym[i]]*primitive_set_in_conventional[0].coord + translations[i]);
          //tmp.type = primitive_set_in_conventional[0].type;
          //
          //if(invec_atom(primitive_set_in_conventional,tmp,1e-4)){ //double tol
          //cerr << ITC_sym_info.symbol[insym[i]] << " " << dirparam[insym[i]] << "  " <<  translations[i] << endl;
          //FCG.push_back(i);
          //}

          //METHOD 2 (DIFFERENCE OF SHIFTS):
          //bool add = true;
          //for(int j=1;j<centeringops.size();j++){
          //	if((aurostd::modulus(translations[i])-aurostd::modulus(centeringops[j])) > -tol ){
          //	  add = false;
          //	}
          //}
          //if(add == true){
          //	cerr << ITC_sym_info.symbol[insym[i]] << " " << dirparam[insym[i]] << "  " <<  translations[i] << endl;
          //	FCG.push_back(i);
          //}

          //METHOD 3: BRUTE FORCE:
          FCG.clear();
          for (uint i = 0; i < pointgroupops.size(); i++) {
            FCG.push_back(i);
            //cerr << "FCG.at(i): " << FCG.at(i) << endl;
          }
        }
      }
      //cerr << "FCG: " << endl;
      //print(FCG);
      //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);

      stringstream axis_cell;
      axis_cell.str(std::string());
      //DX20180806 [OBSOLETE] axis_cell << CCell.lattice_label_ITC << cell_choice;
      axis_cell << cell_choice; //DX20180806 - use setting
      // ===== Construct the generator set based on the cell choice/unique axis ===== //
      ITC_sym_info.initgenerators(axis_cell.str());

      // ===== Check if space group = 1 ===== //
      if(sg_search.size() == 1 && sg_search[0] == 1) {
        spacegroup = 1;
        foundspacegroup = true;
        if(LDEBUG) { cerr << "Space Group " << sg_search[0] << endl; }
      }
      // ===== Check other space groups (2-230), and find origin ===== //
      for (uint j = 0; j < ITC_sym_info.sgindex.size(); j++) {
      //for (int j = ITC_sym_info.sgindex.size()-1; j>=0; j--) //DX20221230 - reverse order to give preference to higher space groups
        if(aurostd::WithinList(sg_search, ITC_sym_info.sgindex[j])) { //DX20210422 - SYM::invec() -> aurostd::WithinList()
          //DEBUGGER:
          if(LDEBUG) {
            cerr << __AFLOW_FUNC__ << " Checking generators for space group: " << ITC_sym_info.sgindex[j] << " (with cell choice = " << axis_cell.str() << ")." <<  endl;
          }
          //cerr << "////////////////////////////////////////" << endl;

          uint match_count = 0;
          // === Obtain origin shift set from the ITC === //
          vector<xvector<double> > ITCshiftset;

          // === For certian space groups, the cell choice and unique axis are important (i.e. MCL). === //
          ITCshiftset = SYM::ReturnITCGenShift(ITC_sym_info.sgindex[j], axis_cell.str());
          //cerr << "library shift set: "<< endl;
          //print(ITCshiftset);
          //xb();

          // === GET COMBINATIONS OF POSSIBLE GENERATOR OPERATIONS === //
          vector<vector<int> > op_ind;  // operation indices divided by type eg: 2+ : (1,3) -1: (2,4)
          for (uint i = 0; i < ITC_sym_info.generators[j].size(); i++) {
            //cerr << ITC_sym_info.generators[j][i].direction << endl;
            vector<int> sub;
            for (uint k = 0; k < FCG.size(); k++) {
              string symboltmp = "";
              //For the purpose of the reduction by generators, and due to the representation of symmetry
              //operations as operations about the origin plus a shift, must equate the glides
              //with a mirror operation
              if(ITC_sym_info.generators[j][i].symbol == "c" || ITC_sym_info.generators[j][i].symbol == "b" ||
                  ITC_sym_info.generators[j][i].symbol == "a" || ITC_sym_info.generators[j][i].symbol == "n" ||
                  ITC_sym_info.generators[j][i].symbol == "d") {
                symboltmp = "m";
              } else {
                symboltmp = ITC_sym_info.generators[j][i].symbol;
              }
              if(pointgroupops[FCG[k]] == symboltmp && ITC_sym_info.sym_mats_direction[insym[FCG[k]] + 1] == ITC_sym_info.generators[j][i].direction) {
                sub.push_back(FCG[k]);
                match_count++;
              }
            }
            //xb();
            //print(sub);
            op_ind.push_back(sub);
          }

          vector<vector<int> > CG;  //centered groups
          CG.clear();
          vector<int> CGtmp;

          // === GET COMBINATIONS OF GENERATORS SINCE THE OPERATIONS CANNOT BE SPLIT BY CENTERING GROUPS WITHOUT KNOWLEDGE OF ORIGIN === //
          int numofgens = ITC_sym_info.generators[j].size();
          if(numofgens == 5) {
            for (uint ix = 0; ix < op_ind[0].size(); ix++) {
              for (uint jx = 0; jx < op_ind[1].size(); jx++) {
                for (uint kx = 0; kx < op_ind[2].size(); kx++) {
                  for (uint lx = 0; lx < op_ind[3].size(); lx++) {
                    for (uint mx = 0; mx < op_ind[4].size(); mx++) {
                      CGtmp.push_back(op_ind[0][ix]);
                      CGtmp.push_back(op_ind[1][jx]);
                      CGtmp.push_back(op_ind[2][kx]);
                      CGtmp.push_back(op_ind[3][lx]);
                      CGtmp.push_back(op_ind[4][mx]);
                      CG.push_back(CGtmp);
                      CGtmp.clear();
                    }
                  }
                }
              }
            }
          } else if(numofgens == 4) {
            for (uint ix = 0; ix < op_ind[0].size(); ix++) {
              for (uint jx = 0; jx < op_ind[1].size(); jx++) {
                for (uint kx = 0; kx < op_ind[2].size(); kx++) {
                  for (uint lx = 0; lx < op_ind[3].size(); lx++) {
                    CGtmp.push_back(op_ind[0][ix]);
                    CGtmp.push_back(op_ind[1][jx]);
                    CGtmp.push_back(op_ind[2][kx]);
                    CGtmp.push_back(op_ind[3][lx]);
                    CG.push_back(CGtmp);
                    CGtmp.clear();
                  }
                }
              }
            }
          } else if(numofgens == 3) {
            for (uint ix = 0; ix < op_ind[0].size(); ix++) {
              for (uint jx = 0; jx < op_ind[1].size(); jx++) {
                for (uint kx = 0; kx < op_ind[2].size(); kx++) {
                  CGtmp.push_back(op_ind[0][ix]);
                  CGtmp.push_back(op_ind[1][jx]);
                  CGtmp.push_back(op_ind[2][kx]);
                  CG.push_back(CGtmp);
                  CGtmp.clear();
                }
              }
            }
          } else if(numofgens == 2) {
            for (uint ix = 0; ix < op_ind[0].size(); ix++) {
              for (uint jx = 0; jx < op_ind[1].size(); jx++) {
                CGtmp.push_back(op_ind[0][ix]);
                CGtmp.push_back(op_ind[1][jx]);
                CG.push_back(CGtmp);
                CGtmp.clear();
              }
            }
          } else if(numofgens == 1) {
            for (uint ix = 0; ix < op_ind[0].size(); ix++) {
              CGtmp.push_back(op_ind[0][ix]);
              CG.push_back(CGtmp);
              CGtmp.clear();
            }
          }
          // ===== Need to determine the reconcile the differences in the translations and ITC shift set ===== //
          //cerr << "match_count: " << match_count << endl;
          //cerr << "ITC_sym_info.generators[j].size()*(centeringgroups): " << ITC_sym_info.generators[j].size()*(centeringgroups) << endl;
          if(match_count == ITC_sym_info.generators[j].size() * (centeringgroups)) {
            if(LDEBUG) {
              cerr << __AFLOW_FUNC__ << " Generator reduction matches for space group: " << ITC_sym_info.sgindex[j] << endl;
            }
            if(LDEBUG) {
              cerr << __AFLOW_FUNC__ << "GENERATOR REDUCTION MATCH: SG" <<  ITC_sym_info.sgindex[j] << endl;
              for(uint ix=0;ix<CG.size();ix++){
                cerr << __AFLOW_FUNC__ << " CG: " << endl;
                print(CG[ix]);
              }
            }
            vector<xvector<double> > oneshiftset;
            vector<xmatrix<double> > onerotset;
            for (uint ix = 0; ix < CG.size(); ix++) {
              onerotset.clear();
              oneshiftset.clear();
              for (uint k = 0; k < CG[ix].size(); k++) {
                onerotset.push_back(Ident - ITC_sym_info.sym_mats[insym[CG[ix][k]]]);
                oneshiftset.push_back(ITCshiftset[k] - translations[CG[ix][k]]);
              }
              if(LDEBUG) {
                cerr << __AFLOW_FUNC__ << "Set of possible shifts: " << endl;
                print(oneshiftset);
                xb();
                cerr << __AFLOW_FUNC__ << "Set of possible rotations: " << endl;
                print(onerotset);
                xb();
              }

              // ===== It becomes a linear algebra problem to solve for a consistent origin shift ===== //
              vector<double> RHS;
              vector<xvector<double> > LHS;
              xvector<double> SOL;
              for (uint l = 0; l < oneshiftset.size(); l++) {
                for (int k = 1; k < 4; k++) {
                  RHS.push_back(oneshiftset[l][k]);
                  LHS.push_back(SYM::extract_row(onerotset[l], k));
                }
              }
              reduction_by_generator_set.push_back(j);
              if(SYM::solve_overdetermined_system(LHS, RHS, OriginShift, CCell.lattice, CCell.dist_nn_min, CCell.sym_eps)) { //DX20190215 - added sym_eps
                if(LDEBUG) {
                  cerr << __AFLOW_FUNC__ << " Generators and origin shift of space group #" << ITC_sym_info.sgindex[j] << " matched to ITC. Testing Wyckoff positions [dir=" << xstr.directory << "]." << endl;
                }

                //cerr << "Space Group " << ITC_sym_info.sgindex[j] << endl;
                //cerr << "ORIGIN SHIFT: " ;
                //cerr << OriginShift << endl;
                //xb();
                spacegroup = ITC_sym_info.sgindex[j];
                //for(int s=0;s<insym.size();s++){
                //cerr << "#" << s << " (" << insym[s]+1 << ") " <<  ITC_sym_info.symbol[insym[s]] << " " << dirparam[insym[s]] << "  " <<  translations[s] << endl;
                //}
                //cerr << "Indices of centered group operations (refer to above): ";
                //print(CG[ix]); //indices that define a centered group as in ITC.
                foundspacegroup = true;
                //DX20180613 - do not break loop in case we need to check other generators - break;
              }
              else {
                if(LDEBUG) {
                  cerr << __AFLOW_FUNC__ << " Generators and origin shift of space group #" << ITC_sym_info.sgindex[j] << " COULD NOT be matched to ITC [dir=" << xstr.directory << "]." << endl;
                }
              }
              //DX20180613 - do not break loop in case we need to check other generators - if(foundspacegroup == true) {
              //DX20180613 - do not break loop in case we need to check other generators -   break;
              //DX20180613 - do not break loop in case we need to check other generators - }
              // ***********************************************************************************************************************************
              // Determine the Wyckoff positions (and verify the space group)
              // ***********************************************************************************************************************************
              // This section maps the atoms in the cell to their corresponding Wyckoff positions.  
              // It is possible that multiple space groups may be consistent with the origin shift.  
              // Thus, the Wyckoff position seach is an additional check to verify the space group.  If the
              // Wyckoff positions are inconsistent, then we search the "sgsearch" list for another possible space group candidate.
              if(foundspacegroup == true) {
                bool obverse_force_transformed = false;
                // ===== Group atoms into equivalent atoms (i.e. mapped onto one another by a symmetry operator) ===== //
                deque<deque<_atom> > equivalent_atoms = SYM::groupSymmetryEquivalentAtoms(atomicbasis, CCell.lattice, pointgroup_mats, translations, CCell.dist_nn_min, CCell.sym_eps); //DX20190215 - added sym_eps
                // ===== Once grouped into equivalent atoms, apply the origin shift ===== //
                deque<deque<_atom> > equivalent_atoms_shifted = SYM::shiftSymmetryEquivalentAtoms(equivalent_atoms, CCell.lattice, OriginShift, CCell.dist_nn_min, CCell.sym_eps); //DX20190215 - added sym_eps

                // If Rhombohedral, a shift may reveal it is not in the obverse setting
                if(CCell.lattice_label_ITC == 'R' && !SYM::isObverseSetting(CCell.lattice, equivalent_atoms_shifted, CCell.dist_nn_min, CCell.sym_eps)) { //DX20190215 - added sym_eps
                  // To store the unshifted, simply shift back
                  xvector<double> ReverseOriginShift = -OriginShift;
                  equivalent_atoms = SYM::shiftSymmetryEquivalentAtoms(equivalent_atoms_shifted, CCell.lattice, ReverseOriginShift, CCell.dist_nn_min, CCell.sym_eps); //DX20190215 - added sym_eps
                  deque<_atom> atoms_tmp;
                  for (uint i = 0; i < equivalent_atoms.size(); i++) {
                    for (uint j = 0; j < equivalent_atoms[i].size(); j++) {
                      atoms_tmp.push_back(equivalent_atoms[i][j]);
                    }
                  }
                  atomicbasis = atoms_tmp;
                  CCell.atoms = atoms_tmp;
                }

                // =============== Obtain Wyckoff positions from the ITC =============== //
                stringstream axis_cell;
                axis_cell.str(std::string());
                //DX20180806 [OBSOLETE] axis_cell << CCell.lattice_label_ITC << cell_choice;
                axis_cell << cell_choice; //DX20180806 - use setting
                if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Lattice type (unique axis if specified) and cell choice  : " << axis_cell.str() << endl; }
                //DX20190215 [OBSOLETE] SYM::initsgs(axis_cell.str());
                ITC_sym_info.initsgs(axis_cell.str()); //DX20190215
                //DX20190215 [OBSOLETE]extern vector<string> gl_sgs;
                //DX20190215 [OBSOLETE]using SYM::gl_sgs;
                spacegroupstring = ITC_sym_info.gl_sgs[spacegroup - 1];
                wyckoffmult = SYM::get_multiplicities(spacegroupstring);
                //print(wyckoffmult);
                general_wyckoff_position = SYM::findGeneralWyckoffPosition(spacegroupstring, wyckoffmult[1]); //DX20170831
                wyckoffsymbols = SYM::get_symmetry_symbols(spacegroupstring);
                //cerr << "WYCKOFFSYMBOLS: " << endl;
                //print(wyckoffsymbols);
                //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);
                //xb();

                // ========== Match equivalent atoms to Wyckoff positions consistent with ITC ========== //
                vector<int> wyckoffsymbols_mult_index;  //wyckoff symbols with specified multiplicity index

                // ===== Determine minimum enumerated Wyckoff letter set ===== //
                // There may be variability in the choice of the Wyckoff letters if the multiplicity
                // and the site symmetry are the same.  We choose the Wyckoff Position scheme with the
                // smallest Wyckoff lettering combination possible (i.e, sum of enumerated Wyckoff letters
                // are minimized) --David Hicks.

                vector<xvector<double> > possible_shifts;
                vector<bool> shift_valid;
                bool orig_origin_shift = true;
                bool other_shifts_explored = false;
                bool final_shift = false;
                uint origin_shift_index = 0;
                vector<int> sum_wyckoff_letters;
                vector<vector<int> > all_enumerated_wyckoff_letters; //DX20221103

                while (orig_origin_shift == true || final_shift == false) {
                  // ========== Pick origin shift to test for Wyckoff Positions ========== //

                  // === Clear vector variables from previous loop === //
                  wyckoffsymbols_mult_index.clear();
                  wyckoffVariables.clear();
                  wyckoffPositionsVector.clear();
                  woss.str("");
                  wyckoffSymbols.clear();

                  // === Scan through all possible origin shifts (original origin shift + ones to minimize wyckoff letters) === //
                  if(other_shifts_explored == false) {
                    // === If not the original origin shift trial === //
                    if(possible_shifts.size() > 0) {
                      if(LDEBUG) {
                        cerr << __AFLOW_FUNC__ << " Exploring possible shifts to find minimum enumerated Wyckoff set. Testing shift "
                          << origin_shift_index << " (" << origin_shift_index+1 << " of " << possible_shifts.size() << "): "
                          << possible_shifts[origin_shift_index] << "." << endl;
                      }
                      //[DX20221001 - OBSOLETE]// if shift is zero, skip
                      //[DX20221001 - OBSOLETE]if(aurostd::iszero(possible_shifts[origin_shift_index])) {
                      //[DX20221001 - OBSOLETE]  origin_shift_index++;
                      //[DX20221001 - OBSOLETE]  if(origin_shift_index == possible_shifts.size()) {
                      //[DX20221001 - OBSOLETE]    other_shifts_explored = true;
                      //[DX20221001 - OBSOLETE]  }
                      //[DX20221001 - OBSOLETE]  if(LDEBUG) {
                      //[DX20221001 - OBSOLETE]    cerr << "SYM::SpaceGroup_ITC: Skipping shift of zero (already checked)." << endl;
                      //[DX20221001 - OBSOLETE]  }
                      //[DX20221001 - OBSOLETE]  continue;
                      //[DX20221001 - OBSOLETE]}

                      xvector<double> previous_shift;
                      if(origin_shift_index > 0) {
                        previous_shift = possible_shifts[origin_shift_index - 1];
                      }
                      SYM::shiftWyckoffPositions(equivalent_atoms_shifted, previous_shift, possible_shifts[origin_shift_index]);
                      //DX20180928 test - SYM::shiftWyckoffPositions(equivalent_atoms_shifted, possible_shifts, origin_shift_index);
                      origin_shift_index++;

                      // == Signifies when all possible origin shifts have been tested == //
                      if(origin_shift_index == possible_shifts.size()) {
                        other_shifts_explored = true;
                      }
                    }
                  } 
                  else {
                    if(LDEBUG) {
                      cerr << __AFLOW_FUNC__ << " Tested (" << possible_shifts.size() << ") possible origin choices. Finding minimum enumerated Wyckoff letters from valid shifts." << endl;
                    }
                    //If all possible origin shifts scanned, pick the one with the minimum sum of enumerated Wyckoff letters
                    int min_wyckoff_config = 0;
                    //DX20220828 - int min_wyckoff_sum = sum_wyckoff_letters[0];
                    min_wyckoff_sum = sum_wyckoff_letters[0]; //DX20220828
                    min_Wyckoff_numbers = all_enumerated_wyckoff_letters[0]; //DX20220828
                    for (uint s = 0; s < sum_wyckoff_letters.size(); s++) {
                      if(LDEBUG) { cerr << __AFLOW_FUNC__ << "shift_valid[s]=" << shift_valid[s] << " sum_wyckoff_letters[s]=" << aurostd::joinWDelimiter(all_enumerated_wyckoff_letters[s], ",") << endl; }
                      if((min_wyckoff_sum > sum_wyckoff_letters[s] && shift_valid[s]) || 
                        (min_wyckoff_sum==sum_wyckoff_letters[s] && all_enumerated_wyckoff_letters[s]<min_Wyckoff_numbers)) //DX20221110 - added tie-breaker
                      {
                        min_wyckoff_config = s;
                        min_wyckoff_sum = sum_wyckoff_letters[s];
                        min_Wyckoff_numbers = all_enumerated_wyckoff_letters[s];
                      }
                    }
                    xvector<double> shift_for_min_wyckoff_config;
                    if(min_wyckoff_config > 0) {
                      shift_for_min_wyckoff_config = possible_shifts[min_wyckoff_config-1];
                      OriginShift = OriginShift + shift_for_min_wyckoff_config; //DX NEW
                    }
                    SYM::shiftWyckoffPositions(equivalent_atoms_shifted, possible_shifts[origin_shift_index-1], shift_for_min_wyckoff_config);
                    final_shift = true;
                  }
                  // Find the Wyckoff positions in the ITC standard representation
                  //perhaps change findWyckoffPositions to mapEquivalentAtomsToWyckoffPositions
                  found_all_wyckoff = SYM::findWyckoffPositions(CCell, atomicbasis, tmpvvvstring, equivalent_atoms, equivalent_atoms_shifted, 
                      foundspacegroup, spacegroupstring, orig_origin_shift, OriginShift, 
                      wyckoffmult, wyckoffsymbols, wyckoffVariables, wyckoffPositionsVector, 
                      wyckoffSymbols, woss, obverse_force_transformed);
                  shift_valid.push_back(found_all_wyckoff);
                  if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Shift is valid ? : " << found_all_wyckoff << endl; }
                  if(found_all_wyckoff == false && orig_origin_shift == true) {
                    foundspacegroup=false;
                    break;
                  } else if(found_all_wyckoff == false && orig_origin_shift == false) {
                    sum_wyckoff_letters.push_back(AUROSTD_MAX_INT); //DX20220827 - changed from 100 to AUROSTD_MAX_DOUBLE
                    vector<int> vec_int_tmp;
                    all_enumerated_wyckoff_letters.push_back(vec_int_tmp); //DX20221103 
                    continue;
                  }
                  // ===== Determine multiplicity, letters, and site symmetry in WYCCAR ===== //
                  // This allows us to investigate if there are other Wyckoff positions corresponding to the
                  // ones in the POSCAR which will find the "lowest" Wyckoff letter enumeration scheme
                  // ===== Extract Wyckoff info =====
                  vector<GroupedWyckoffPosition> grouped_Wyckoff_positions_structure;
                  vector<wyckoffsite_ITC> sorted_Wyckoff_sites = wyckoffVariables;
                  std::sort(sorted_Wyckoff_sites.begin(), sorted_Wyckoff_sites.end(), sortWyckoffByLetter); //DX20221201 - need to sort for tie-breaker
                  std::sort(sorted_Wyckoff_sites.begin(), sorted_Wyckoff_sites.end(), sortWyckoffByType); //DX20221201 - need to sort for tie-breaker
                  compare::groupWyckoffPositions(sorted_Wyckoff_sites, grouped_Wyckoff_positions_structure);
                  if(LDEBUG) {
                    string Wyckoff_string_structure = anrl::groupedWyckoffPosition2ANRLString(grouped_Wyckoff_positions_structure, true);
                    cerr << __AFLOW_FUNC__ << "Wyckoff_string_structure=" << Wyckoff_string_structure << endl;
                  }
                  vector<int> wyckoff_mult;
                  vector<string> wyckoff_letter;
                  vector<string> wyckoff_site_sym;

                  // ---------------------------------------------------------------------------
                  // Store the wyckoff information into vectors
                  // DX20221201 - updated the loop below to use the sorted/grouped Wyckoff positions
                  // (sorting is required to find the minimum enumerated Wyckoff set)
                  for (uint w = 0; w < grouped_Wyckoff_positions_structure.size(); w++) {
                    for (uint x = 0; x < grouped_Wyckoff_positions_structure[w].multiplicities.size(); x++) {
                      wyckoff_mult.push_back(grouped_Wyckoff_positions_structure[w].multiplicities[x]);
                      wyckoff_letter.push_back(grouped_Wyckoff_positions_structure[w].letters[x]);
                      wyckoff_site_sym.push_back(grouped_Wyckoff_positions_structure[w].site_symmetries[x]);
                    }
                  }

                  // ===== Calculate the sum of the enumerated wyckoff letters ===== //
                  vector<int> enumerated_wyckoff_letters = SYM::enumerate_wyckoff_letters(wyckoff_letter); //DX20221103
                  sum_wyckoff_letters.push_back(aurostd::sum(enumerated_wyckoff_letters));
                  all_enumerated_wyckoff_letters.push_back(enumerated_wyckoff_letters); //DX20221103

                  // ===== If we have performed the original origin shift, now check the other possibilites ===== //
                  if(orig_origin_shift == true) {
                    orig_origin_shift = false;

                    // ===== Identify minimum enumerated Wyckoff letter occuring in structure ===== // 
                    int min_enumerated_letter = 1e9;
                    int min_multiplicity = 1e9;
                    string min_letter;
                    string min_site_sym;
                    for (uint w = 0; w < wyckoff_mult.size(); w++) {
                      int enumerated_letter = SYM::enumerate_wyckoff_letter(wyckoff_letter[w]);
                      if(enumerated_letter < min_enumerated_letter){
                        min_multiplicity = wyckoff_mult[w];
                        min_letter = wyckoff_letter[w];
                        min_site_sym = wyckoff_site_sym[w];
                        min_enumerated_letter = enumerated_letter;
                      }
                    }

                    // ===== Determine current Wyckoff lettering enumeration ===== //
                    vector<int> found_enumerated_Wyckoff_numbers = SYM::enumerate_wyckoff_letters(wyckoff_letter);
                    int sum_enumerated_found_letters = aurostd::sum(found_enumerated_Wyckoff_numbers);

                    // ===== Determine minimum Wyckoff lettering enumeration (only need to do this once) ===== //
                    if(minimum_scheme_sum==1e9) {
                      vector<string> minimum_enumerated_Wyckoff_letters = SYM::get_minimum_enumerated_Wyckoff_letters(spacegroupstring, wyckoff_mult, wyckoff_site_sym);
                      vector<int> minimum_enumerated_Wyckoff_numbers = SYM::enumerate_wyckoff_letters(minimum_enumerated_Wyckoff_letters);
                      //int minimum_scheme_sum = aurostd::sum(minimum_enumerated_Wyckoff_numbers);
                      minimum_scheme_sum = aurostd::sum(minimum_enumerated_Wyckoff_numbers); //DX20220828
                    }

                    // ===== If we do not have the minimum enumerated Wyckoff letter set, check possible origin shifts ===== //
                    if(LDEBUG){ cerr << __AFLOW_FUNC__ << "SYM::SpaceGroup_ITC: minimum_scheme_sum=" << minimum_scheme_sum << " vs current=" << sum_enumerated_found_letters << endl; }
                    if(sum_enumerated_found_letters >= minimum_scheme_sum) { //DX20221127 - added = condition, keep this condition in case we want to skip it later
                      if(LDEBUG) {
                        cerr << __AFLOW_FUNC__ << " Enumerated Wyckoff letters is not minimized; check origin shifts." << endl;
                      }
                      possible_shifts = SYM::get_possible_origin_shifts(spacegroupstring, min_multiplicity, min_site_sym);
                      min_wyckoff_sum = sum_enumerated_found_letters; //DX20230112
                      min_Wyckoff_numbers = found_enumerated_Wyckoff_numbers; //DX20230110
                    }
                    else {
                      final_shift = true;
                      min_wyckoff_sum = sum_enumerated_found_letters; //DX20220828
                      min_Wyckoff_numbers = found_enumerated_Wyckoff_numbers; //DX20230110
                    }
                  }
                }
              }
              if(foundspacegroup == true && found_all_wyckoff == true) { //DX20180613 - moved generator loop
                break; //DX20180613 - moved generator loop
              } //DX20180613 - moved generator loop
            } //DX20180613 - moved generator loop
          } //DX20180613 - moved generator loop
          if(foundspacegroup == true && found_all_wyckoff == true) {
            break;
          }
        }
        if(foundspacegroup == true && found_all_wyckoff == true) {
          break;
        }
      }
      //cerr << "SPACEGROUP:: " <<spacegroup << endl;
      if(foundspacegroup == true) {
        if(LDEBUG) {
          cerr << "SYM::SpaceGroup_ITC: spacegroup=" << spacegroup << " iterate=" << iterate << " : min_wyckoff_sum=" << min_wyckoff_sum << " vs global=" << global_min_wyckoff_sum << " (spacegroup_global=" << spacegroup_global << ")" <<  endl;
          cerr << "SYM::SpaceGroup_ITC: min_Wyckoff_numbers: " << min_Wyckoff_numbers.size() << " " << aurostd::joinWDelimiter(min_Wyckoff_numbers,",") << endl;
          cerr << "SYM::SpaceGroup_ITC: global_min_Wyckoff_numbers: size=" << global_min_Wyckoff_numbers.size() << " enumerations=" << aurostd::joinWDelimiter(global_min_Wyckoff_numbers,",") << endl;
        }
        foundspacegroup_global = true;
        //spacegroup_global = spacegroup;
        if(spacegroup_global == 0) { spacegroup_global = spacegroup; } //DX20221229
        foundspacegroup = false;
        //if(min_wyckoff_sum<global_min_wyckoff_sum || 
        if(spacegroup >= spacegroup_global &&
            (min_wyckoff_sum<global_min_wyckoff_sum || 
           (min_wyckoff_sum==global_min_wyckoff_sum && global_min_Wyckoff_numbers.size()>0 && min_Wyckoff_numbers<global_min_Wyckoff_numbers))
        ) { //DX20221103 - added tie-breaker
          // update minimum Wyckoff sum
          spacegroup_global = spacegroup; //DX20221229
          global_min_wyckoff_sum = min_wyckoff_sum;
          global_min_Wyckoff_numbers = min_Wyckoff_numbers; //DX20221103
          // update CCell
          CCell_Wyckoff_min = CCell;
          // update Wyckoff info
          wyckoffmult_min = wyckoffmult;
          wyckoffsymbols_min = wyckoffsymbols;
          wyckoffPositionsVector_min = wyckoffPositionsVector;
          wyckoffSymbols_min = wyckoffSymbols;
          wyckoffVariables_min = wyckoffVariables;
          tmpvvvstring_min = tmpvvvstring;
          woss_min.str(""); woss_min << woss.str(); //DX+KC20221103 - added woss_min for aflow --wyckoff --print=json
          // update shifts
          OriginShift_Wyckoff_min = OriginShift;
        }
        if(iterate == (int)candidate_lattice_vectors_prev.size()-1) {//DX20221201 [do not stop searching once minimum Wyckoff sequence is found - need to keep searching to get most alphabetical from left to right] 
          tested_all_lattices = true;
          foundspacegroup = true;
          spacegroup = spacegroup_global;
          CCell = CCell_Wyckoff_min;
          // update Wyckoff info
          wyckoffmult = wyckoffmult_min;
          wyckoffsymbols = wyckoffsymbols_min;
          wyckoffPositionsVector = wyckoffPositionsVector_min;
          wyckoffSymbols = wyckoffSymbols_min;
          wyckoffVariables = wyckoffVariables_min;
          tmpvvvstring = tmpvvvstring_min;
          woss.str(""); woss << woss_min.str(); //DX+KC20221103 - added woss_min for aflow --wyckoff --print=json
          // update shifts
          OriginShift = OriginShift_Wyckoff_min;
        }
        // loop through conventional lattices
        // else continue looking
      }
    } // end of tested_all_lattices
  }
  if(foundspacegroup_global == false) { //DX20220827 - changed from foundspacegroup to foundspacegroup_global
    message << "Failed to find WYCKOFF POSITIONS, ORIGIN SHIFT, or inconsistent number of GENERATORS [dir=" << (*this).directory << "].";
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_RUNTIME_ERROR_);
  }

  if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Tolerance used to obtain this space group: " << CCell.sym_eps << " [dir=" << xstr.directory << "]." << endl; } //DX20190215 - _SYM_TOL_ to CCell.sym_eps
  if(CCell.sym_eps != use_tol) { //DX20190215 - _SYM_TOL_ to CCell.sym_eps
    //cerr << "TOL: " << CCell.sym_eps << endl; //DX20190215 - _SYM_TOL to CCell.sym_eps
    use_tol = CCell.sym_eps; //DX20190215 - _SYM_TOL to CCell.sym_eps
  }
  (*this).sym_eps = CCell.sym_eps; //DX20190215 - _SYM_TOL to CCell.sym_eps, redundant, but safe
  (*this).sym_eps_change_count = sym_eps_change_count; //DX20180226 - added sym eps change count
  (*this).sym_eps_no_scan = no_scan; //DX20210430 - added no_scan

  // ***********************************************************************************************************************************
  // Prepare/Create Wyccar
  // ***********************************************************************************************************************************
  vector<int> wyckofftypescount = SYM::count_types(wyckoffPositionsVector);

  //Get angles alpha beta gamma and magnitudes of basis vectors
  xmatrix<double> latmat = CCell.lattice;
  xvector<double> latmat1 = SYM::extract_row(latmat, 1);
  xvector<double> latmat2 = SYM::extract_row(latmat, 2);
  xvector<double> latmat3 = SYM::extract_row(latmat, 3);

  double alpha = SYM::get_angle(latmat2, latmat3, "D");
  double beta = SYM::get_angle(latmat3, latmat1, "D");
  double gamma = SYM::get_angle(latmat1, latmat2, "D");

  double amag = aurostd::modulus(latmat1);
  double bmag = aurostd::modulus(latmat2);
  double cmag = aurostd::modulus(latmat3);

  // === ITC Cell Choice === //
  stringstream axis_cell;
  axis_cell.str(std::string());
  //DX20180806 [OBSOLETE] axis_cell << CCell.bravais_label_ITC << cell_choice;
  axis_cell << cell_choice; //DX20180806 - use setting

  // default setting choice
  int setting_choice = 1;
  if(CCell.bravais_label_ITC == 'R'){
    setting_choice = 2;  //DX: SC has HEX setting as option 2 in aflow_wyckoff.cpp
  }
  //check if preferred setting choice
  if(cell_choice!=0){
    setting_choice = cell_choice;
    //if the spacegroup only has option 1, change to option 1
    if(!spacegroup::SpaceGroupOptionRequired(spacegroup)){
      setting_choice = 1;
    }
    //DX20190131 - handle ANRL setting - START
    else if(cell_choice==SG_SETTING_ANRL){
      setting_choice = anrl::getANRLSettingChoice(spacegroup);
    }
    //DX20190131 - handle ANRL setting - END
  }

  // ===== WYCCAR OUTPUT: ===== //
  stringstream wss;

  // === Add space group information to the title line of the Wyccar === //
  wss << CCell.title << "| SG: " << GetSpaceGroupName(spacegroup,CCell.directory) << " " << spacegroup << " PG: " << pointgroup << " BL: " << CCell.bravais_label_ITC << " | sym_eps: " << (*this).sym_eps << endl;
  // === Add lattice parameter info === //
  wss << setprecision(4) << fixed << 1.00 << endl;  // MODIFIED (RHT)
  if(axis_cell.str().find("b") != std::string::npos) {
    wss << setprecision(4) << fixed << amag << " " << bmag << " " << cmag << " " << alpha << " " << beta << " " << gamma << " " << spacegroup << " "
      << "1" << endl;  //If b is unique axis ==> first choice (1)
  } else if(axis_cell.str().find("m") != std::string::npos) {
    wss << setprecision(4) << fixed << amag << " " << bmag << " " << cmag << " " << alpha << " " << beta << " " << gamma << " " << spacegroup << " "
      << "2" << endl;  //If c is unique axis ==> second choice (2)
  } else {
    wss << setprecision(4) << fixed << amag << " " << bmag << " " << cmag << " " << alpha << " " << beta << " " << gamma << " " << spacegroup << " " << setting_choice << endl;  // CODE PREFERS OPTION 1 (ie, first ITC origin choice)
  }
  // === Add number of types of Wyckoff positions === ..
  for (uint i = 0; i < wyckofftypescount.size(); i++) {
    wss << wyckofftypescount[i] << " ";
  }
  // === Add Wyckoff positions === //
  wss << endl
    << "Direct(WYCCAR)" << endl;
  wss << woss.str();  //<< endl;

  vector<string> wyccar_string_vec;
  aurostd::stream2vectorstring(wss, wyccar_string_vec);

  //DX20190212 START
  // lattice type/centering fixes
  // rhombohedral fixes
  if(latticetypechar == 'r' && CCell.bravais_label_ITC == 'P'){CCell.bravais_label_ITC = 'R'; } 
  if(latticetypechar == 'R' || latticetypechar == 'r'){latticetypechar = 'h'; }
  // monoclinic fixes
  if(latticetypechar == 'b'){latticetypechar = 'm'; }  //DX20190311 - turns b to m for standard monoclinic designation
  //DX20190212 END

  // ========== Update the characteristics for the xstructure ========== //
  (*this).crystal_system_ITC = CCell.crystal_system_ITC;
  (*this).point_group_ITC = pointgroup;
  (*this).bravais_label_ITC = CCell.bravais_label_ITC;
  (*this).lattice_label_ITC = latticetypechar;
  (*this).space_group_ITC = (uint)spacegroup;
  (*this).wyckoff_library_entry_ITC = spacegroupstring;
  (*this).wyccar_ITC = wyccar_string_vec;
  (*this).standard_lattice_ITC = CCell.lattice;
  (*this).wyckoff_sites_ITC = wyckoffVariables;
  (*this).wyckoff_symbols_ITC = wyckoffSymbols;
  (*this).setting_ITC = setting_choice; //DX20170830 - SGDATA
  (*this).origin_ITC = OriginShift; //DX20170830 - SGDATA
  (*this).general_position_ITC = general_wyckoff_position; //DX20170830 - SGDATA
  //cerr << "TOLERANCE USED: " << (*this).sym_eps << endl; //DX20190215 - _SYM_TOL to (*this).sym_eps
  
  // convert iomode to original format
  (*this).iomode = iomode_tmp;
  
  return (uint)spacegroup;
}

// ***************************************************************************
// Calculate orthogonality defect
// ***************************************************************************
namespace SYM{
  string OrthoDefect(istream& cin) {
    _aflags aflags;
    aflags.QUIET = TRUE;
    ofstream File("/dev/null");
    //bool verbose=FALSE;
    //bool WRITE=FALSE;
    stringstream oss;
    xstructure a(cin, IOVASP_AUTO);
    oss << a << endl;
    //  cerr << "";

    //  cerr << a.lattice(1) << endl;
    //  cerr << a.lattice(2) << endl;
    //  cerr << a.lattice(3) << endl;

    xmatrix<double> latticemat = xvec2xmat(a.lattice(1), a.lattice(2), a.lattice(3));

    oss << "ortho defect: " << orthogonality_defect(latticemat) << endl;
    oss << "volume of cell: " << aurostd::det(latticemat) << endl;
    oss << "-----------------------------------------" << endl;
    return oss.str();
  }
} //namespace SYM

// ***************************************************************************
// Latex plot lattice
// ***************************************************************************
namespace SYM {
  stringstream* latex_plot_lattice(xmatrix<double> L, string color, int a, int b, int c) {
    stringstream* ss = new stringstream;

    *ss << "\\coordinate (O) at (0,0,0);" << endl;
    *ss << "\\draw[dashed,gray,->] (O) -- (0.3,0,0) node[anchor=north east]{$x$};" << endl;
    *ss << "\\draw[dashed,gray,->] (O) -- (0,0.3,0) node[anchor=north west]{$y$};" << endl;
    *ss << "\\draw[dashed,gray,->] (O) -- (0,0,0.3) node[anchor=south]{$z$};" << endl;
    *ss << "\\draw[dashed," << color << ",->] (O) -- (" << L(1, 1) << "," << L(1, 2) << "," << L(1, 3) << ") node[anchor=north east]{$L1$};" << endl;
    *ss << "\\draw[dashed," << color << ",->] (O) -- (" << L(2, 1) << "," << L(2, 2) << "," << L(2, 3) << ") node[anchor=north east]{$L2$};" << endl;
    *ss << "\\draw[dashed," << color << ",->] (O) -- (" << L(3, 1) << "," << L(3, 2) << "," << L(3, 3) << ") node[anchor=north east]{$L3$};" << endl;

    *ss << "\\foreach \\x in {" << 0 << ",...," << a + 1 << "}{" << endl;
    *ss << "\\foreach \\y in {" << 0 << ",...," << b + 1 << "}{" << endl;
    *ss << "\\foreach \\z in {" << 0 << ",...," << c + 1 << "}{" << endl;
    *ss << "\\node[draw,circle,inner sep=5pt] at (";
    *ss << L(1, 1) << "*\\x+" << L(2, 1) << "*\\y+" << L(3, 1) << "*\\z," << L(1, 2) << "*\\x+" << L(2, 2) << "*\\y+" << L(3, 2) << "*\\z," << L(1, 3) << "*\\x+" << L(2, 3) << "*\\y+" << L(3, 3) << "*\\z) {};" << endl;
    *ss << "}" << endl
      << "}" << endl
      << "}" << endl;

    *ss << "\\foreach \\x in {" << 0 << ",...," << a << "}{" << endl;
    *ss << "\\foreach \\y in {" << 0 << ",...," << b << "}{" << endl;
    *ss << "\\foreach \\z in {" << 0 << ",...," << c << "}{" << endl;
    *ss << "\\coordinate (O) at (" << L(1, 1) << "*\\x+" << L(2, 1) << "*\\y+" << L(3, 1) << "*\\z," << L(1, 2) << "*\\x+" << L(2, 2) << "*\\y+" << L(3, 2) << "*\\z," << L(1, 3) << "*\\x+" << L(2, 3) << "*\\y+" << L(3, 3) << "*\\z);" << endl;
    *ss << "\\coordinate (P1) at (" << L(1, 1) << "+" << L(1, 1) << "*\\x+" << L(2, 1) << "*\\y+" << L(3, 1) << "*\\z," << L(1, 2) << "+" << L(1, 2) << "*\\x+" << L(2, 2) << "*\\y+" << L(3, 2) << "*\\z," << L(1, 3) << "+" << L(1, 3) << "*\\x+" << L(2, 3) << "*\\y+" << L(3, 3) << "*\\z);" << endl;
    *ss << "\\coordinate (P2) at (" << L(2, 1) << "+" << L(1, 1) << "*\\x+" << L(2, 1) << "*\\y+" << L(3, 1) << "*\\z," << L(2, 2) << "+" << L(1, 2) << "*\\x+" << L(2, 2) << "*\\y+" << L(3, 2) << "*\\z," << L(2, 3) << "+" << L(1, 3) << "*\\x+" << L(2, 3) << "*\\y+" << L(3, 3) << "*\\z);" << endl;
    *ss << "\\coordinate (P3) at (" << L(3, 1) << "+" << L(1, 1) << "*\\x+" << L(2, 1) << "*\\y+" << L(3, 1) << "*\\z," << L(3, 2) << "+" << L(1, 2) << "*\\x+" << L(2, 2) << "*\\y+" << L(3, 2) << "*\\z," << L(3, 3) << "+" << L(1, 3) << "*\\x+" << L(2, 3) << "*\\y+" << L(3, 3) << "*\\z);" << endl;
    *ss << "\\coordinate (P12) at (" << L(1, 1) + L(2, 1) << "+" << L(1, 1) << "*\\x+" << L(2, 1) << "*\\y+" << L(3, 1) << "*\\z," << L(1, 2) + L(2, 2) << "+" << L(1, 2) << "*\\x+" << L(2, 2) << "*\\y+" << L(3, 2) << "*\\z," << L(1, 3) + L(2, 3) << "+" << L(1, 3) << "*\\x+" << L(2, 3) << "*\\y+" << L(3, 3) << "*\\z);" << endl;
    *ss << "\\coordinate (P13) at (" << L(1, 1) + L(3, 1) << "+" << L(1, 1) << "*\\x+" << L(2, 1) << "*\\y+" << L(3, 1) << "*\\z," << L(1, 2) + L(3, 2) << "+" << L(1, 2) << "*\\x+" << L(2, 2) << "*\\y+" << L(3, 2) << "*\\z," << L(1, 3) + L(3, 3) << "+" << L(1, 3) << "*\\x+" << L(2, 3) << "*\\y+" << L(3, 3) << "*\\z);" << endl;
    *ss << "\\coordinate (P23) at (" << L(2, 1) + L(3, 1) << "+" << L(1, 1) << "*\\x+" << L(2, 1) << "*\\y+" << L(3, 1) << "*\\z," << L(2, 2) + L(3, 2) << "+" << L(1, 2) << "*\\x+" << L(2, 2) << "*\\y+" << L(3, 2) << "*\\z," << L(2, 3) + L(3, 3) << "+" << L(1, 3) << "*\\x+" << L(2, 3) << "*\\y+" << L(3, 3) << "*\\z);" << endl;
    *ss << "\\coordinate (P123) at (" << L(1, 1) + L(2, 1) + L(3, 1) << "+" << L(1, 1) << "*\\x+" << L(2, 1) << "*\\y+" << L(3, 1) << "*\\z," << L(1, 2) + L(2, 2) + L(3, 2) << "+" << L(1, 2) << "*\\x+" << L(2, 2) << "*\\y+" << L(3, 2) << "*\\z," << L(1, 3) + L(2, 3) + L(3, 3) << "+" << L(1, 3) << "*\\x+" << L(2, 3) << "*\\y+" << L(3, 3) << "*\\z);" << endl;

    *ss << "\\if{\\x = 0 \\AND  \\y = 0 \\AND \\z = 0}{" << endl;
    *ss << "\\draw[dashed,color=" << color << " ] (O) -- (P1);" << endl;
    *ss << "\\draw[dashed,color=" << color << " ] (O) -- (P2);" << endl;
    *ss << "\\draw[dashed,color=" << color << " ] (O) -- (P3);" << endl;
    *ss << "}" << endl;
    *ss << "\\else {" << endl;
    //  *ss << "\\draw[thick,color=" << color << " ] (O) -- (P1);" << endl;
    //  *ss << "\\draw[thick,color=" << color << " ] (O) -- (P2);" << endl;
    //  *ss << "\\draw[thick,color=" << color << " ] (O) -- (P3);" << endl;
    *ss << "\\draw[thick,color=" << color << " ] (P1) -- (P12);" << endl;
    *ss << "\\draw[thick,color=" << color << " ] (P1) -- (P13);" << endl;
    *ss << "\\draw[thick,color=" << color << " ] (P2) -- (P12);" << endl;
    *ss << "\\draw[thick,color=" << color << " ] (P2) -- (P23);" << endl;
    *ss << "\\draw[thick,color=" << color << " ] (P3) -- (P13);" << endl;
    *ss << "\\draw[thick,color=" << color << " ] (P3) -- (P23);" << endl;
    *ss << "\\draw[thick,color=" << color << " ] (P123) -- (P12);" << endl;
    *ss << "\\draw[thick,color=" << color << " ] (P123) -- (P13);" << endl;
    *ss << "\\draw[thick,color=" << color << " ] (P123) -- (P23);" << endl;
    *ss << "}\\fi";
    *ss << "}" << endl
      << "}" << endl
      << "}" << endl;
    return ss;
  }
} //namespace SYM

// ***************************************************************************
// Plot lattice
// ***************************************************************************
//DX20161222 - Is this function needed?  If so, need to switch over to xstructure
//string plot_lattice(vector<string> files){
//
//  int a = atoi(files[2].c_str());
//  int b = atoi(files[3].c_str());
//  int c = atoi(files[4].c_str());
//
//  CrystalStructure C;
//  xmatrix<double> L;
//  stringstream * ss = new stringstream;
//  *ss << "\\documentclass{article}" << endl;
//  *ss << "\\usepackage{verbatim}" << endl;
//  *ss << "\\usepackage{tikz}" << endl;
//  *ss << "\\usepackage{3dplot}" << endl;
//  *ss << "\\usepackage[active,tightpage]{preview}" << endl;
//  *ss << "\\PreviewEnvironment{tikzpicture}" << endl;
//  *ss << "\\setlength\\PreviewBorder{2mm}" << endl;
//  *ss << "\\begin{document}" << endl;
//  *ss << "\\tdplotsetmaincoords{60}{110}" << endl;
//  *ss << "\\begin{tikzpicture}[scale=5,tdplot_main_coords]" << endl;
//
//  vector<string> colors;
//  colors.push_back("red");
//  colors.push_back("green");
//  colors.push_back("blue");
//  colors.push_back("yellow");
//  colors.push_back("purple");
//  for(uint i=0;i<(files.size()-3)/5+1;i++){
//     colors.insert(colors.end(),colors.begin(),colors.end());
//  }
//
//  ifstream fin;
//  for(uint i=5;i<files.size();i++){
//    //    char* fn;
//    //    sprintf(fn, "fin.%n",i);
//    //    cerr << (*fn) << endl;
//
//    fin.open(files[i].c_str());
//    //GET LATTICE
//    L=C.get_latticexmat(fin);
//    //CrystalStructure C(cin);
//    //xmatrix<double> L = C.prim_lattice;
//    *ss << (*latex_plot_lattice(L,colors[i-5],a,b,c)).str() << endl;
//    fin.close();
//  }
//  *ss <<"\\end{tikzpicture}" << endl;
//  *ss <<"\\end{document}" << endl;
//  //cout << (*ss).str() << endl;
//  return ss->str();
//}

// ***************************************************************************
// //DX20210118 [OBSOLETE - now GetNumEachType in xatom] SYM::arrange_atoms
// ***************************************************************************

#endif

// Written by Richard H. Taylor
// Updated by David Hicks
// d.hicks@duke.edu
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
