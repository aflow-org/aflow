// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo
// Dane Morgan

#ifndef _AFLOW_PFLOW_PRINT_CPP_
#define _AFLOW_PFLOW_PRINT_CPP_

#include "aflow.h"
#include "aflow_pflow.h"
#include "aflow_symmetry_spacegroup.h" //DX20180823

// ***************************************************************************
string _CleanElementName(const string& name) {
  // Dane Morgan - Stefano Curtarolo
  string name1,name2,nameout;
  name1=name+" ";
  name2=name+" ";
  name1=name1.substr(0,1);
  if(!((name1[0]>=65 && name1[0]<=90)||(name1[0]>=97 && name1[0]<=122))) name1="";
  if((name1[0]>=97 && name1[0]<=122)) name1[0]-=-97+65;
  name2=name2.substr(1,1);
  if(!((name2[0]>=65 && name2[0]<=90)||(name2[0]>=97 && name2[0]<=122))) name2="";
  if((name2[0]>=65 && name2[0]<=90)) name2[0]+=-97+65;
  nameout=name1+name2;
  return nameout;
}

// ***************************************************************************
// pflow::PrintACE
// ***************************************************************************
// This funtion prints out structural data in a ace format (for CaRIne).
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintACE(const xstructure& str, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "pflow::PrintACE: BEGIN" << endl;
    // Print out data in ace format
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(10);
    oss << "CaRIne Crystallography 3.0" << endl;
    oss << "Cell Standard ASCII File=*.ACE" << endl;
    oss << "File Version Number=1" << endl;
    oss << endl;
    //    oss << "Colors Refered to: Mendele" << (char) 239 << "ev Table" << endl;
    oss << "Colors Refered to: RGB Components" << endl;
    oss << endl;
    oss << "----------------------" << endl;
    oss << "Cell Definition" << endl;
    oss << "----------------------" << endl;
    oss << endl;
    oss << "Cell Parameters(" << (char) 197 << " and " << (char) 176 << ")"  << endl;
    // Get cell data
    xvector<double> data(6);
    data=Getabc_angles(str.lattice,DEGREES);data(1)*=str.scale;data(2)*=str.scale;data(3)*=str.scale;
    oss << "a=" << data(1) << '\t'
      << "b=" << data(2) << '\t'
      << "c=" << data(3) << endl;
    oss << "alpha=" << data(4) << '\t'
      << "beta=" << data(5) << '\t'
      << "gamma=" << data(6) << endl;
    oss << endl;
    oss << "System=not used" << endl;
    oss << "Space Group Number=not used" << endl;
    oss << "Number of positions in Cell=" << str.atoms.size() << endl;
    oss << endl;
    // If cartesian convert to direct
    oss << "Atom" << '\t' << "Oxi." << '\t' << "X" << '\t' << "Y" << '\t' << "Z" << '\t'
      << "R(A)" << '\t' << "Occ." << endl;
    for(uint i=0;i<str.atoms.size();i++) {
      oss << (str.atoms.at(i)).cleanname << '\t' << "0" << '\t' << str.atoms.at(i).fpos(1) << '\t' << str.atoms.at(i).fpos(2) << '\t' << str.atoms.at(i).fpos(3) << '\t'
        << "1.00" << '\t' << "1.0000" << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::PrintAngles
// ***************************************************************************
// This funtion prints out 3 lattice vectors in cartesian coordinates.
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintAngles(xstructure str, const double& cutoff, ostream& oss) {  
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "pflow::PrintAngles: BEGIN" << endl;
    uint MAX_NUM_ANGLE=21;
    oss.setf(std::ios::left, std::ios::adjustfield);
    oss.setf(std::ios::fixed, std::ios::floatfield);
    //  oss.setf(std::ios::fixed,std::ios::floatfield);
    //  oss.precision(10);
    deque<deque<_atom> > neigh_mat;
    // [OBSOLETE]  GetNeighData(str,cutoff,neigh_mat);
    str.GetNeighData(cutoff,neigh_mat);   // once PrintAngles goes in xstructure I can remove the copy  //CO20220623 - using new GetNeighData()
    double tol=1e-15;

    // Output header
    oss << "All atoms identified by" << endl;
    oss << "[basis number]  [name]  [n1 n2 n3 (unit cell)]" << endl;

    // Print out angles
    oss << endl;
    oss << "********************" << " ANGLES " << "********************" << endl;
    oss << "Atoms of increasing indentation are A, B, and C." << endl;
    oss << "The given angles are those made by the lines A to B and B to C (atom B is at the vertex)." << endl;
    oss << "The given angles are always taken 0<=angle<=180." << endl;
    oss << "Angles are only calculated for the first " << MAX_NUM_ANGLE-1 << " neighbors within " << cutoff <<  endl;
    xmatrix<int> ijk(3,3);
    xvector<int> pflag(3);
    for(uint ia1=0;ia1<neigh_mat.size();ia1++) {
      pflag(1)=1;
      _atom a1=neigh_mat.at(ia1).at(0);
      xvector<double> a1pos=a1.cpos;
      for(int i=1;i<=3;i++) ijk(1,i)=a1.ijk(i);
      for(uint ia2=0;(ia2<neigh_mat.at(ia1).size() && ia2<MAX_NUM_ANGLE);ia2++) {
        pflag(2)=1;
        _atom a2=neigh_mat.at(ia1).at(ia2);
        xvector<double> a2pos=a2.cpos;
        for(int i=1;i<=3;i++) ijk(2,i)=a2.ijk(i);
        for(uint ia3=0;(ia3<neigh_mat.at(a2.basis).size() && ia3<MAX_NUM_ANGLE);ia3++) //[CO20200130 - number->basis]for(uint ia3=0;(ia3<neigh_mat.at(a2.number).size() && ia3<MAX_NUM_ANGLE);ia3++) 
          {
          pflag(3)=1;
          _atom a3=neigh_mat.at(a2.basis).at(ia3);       //[CO20200130 - number->basis]_atom a3=neigh_mat.at(a2.number).at(ia3);      
          xvector<double> a3pos=a3.cpos;
          for(int i=1;i<=3;i++) ijk(3,i)=a3.ijk(i);
          xvector<double> a1toa2(3),a2toa3(3),a1toa3(3);
          for(int ic=1;ic<=3;ic++) {
            a1toa2(ic)=a2pos(ic)-a1pos(ic);
            a2toa3(ic)=a3pos(ic)-a2pos(ic);
            a1toa3(ic)=a3pos(ic)-a1pos(ic);
          }
          if(modulus(a1toa2)>tol && modulus(a2toa3)>tol && modulus(a1toa3)>tol) {
            double angle = 180.0-rad2deg*acos(getcos(a1toa2,a2toa3));
            if(aurostd::abs(getcos(a1toa2,a2toa3)-1)<tol) angle=180;
            if(aurostd::abs(getcos(a1toa2,a2toa3)+1)<tol) angle=0;
            if(pflag(1)) {
              oss << setw(4) << a1.basis+1 //[CO20200130 - number->basis]oss << setw(4) << a1.number+1
                << " " << setw(4) << a1.name.c_str()
                << " " << setw(3) << ijk(1,1)
                << " " << setw(3) << ijk(1,2)
                << " " << setw(3) << ijk(1,3)
                << "   ";
              oss << endl;
            } // if plag(1)
            if(pflag(2)) {
              oss << "     " << setw(4) << a2.basis+1  //[CO20200130 - number->basis]oss << "     " << setw(4) << a2.number+1
                << " " << setw(4) << a2.name.c_str()
                << " " << setw(3) << ijk(2,1)
                << " " << setw(3) << ijk(2,2)
                << " " << setw(3) << ijk(2,3)
                << "   ";
              oss << endl;
            } // if plag(2)
            if(pflag(3)) {
              oss << "          " << setw(4) << a3.basis+1 //[CO20200130 - number->basis]oss << "          " << setw(4) << a3.number+1
                << " " << setw(4) << a3.name.c_str()
                << " " << setw(3) << ijk(3,1)
                << " " << setw(3) << ijk(3,2)
                << " " << setw(3) << ijk(3,3)
                << "   ";
              oss << "     "  << std::setprecision(4) << angle << endl;
            } // if plag(3)
            // Set print flags not to print until next looping.
            for(int i=1;i<=3;i++) pflag(i)=0;	  
          } // if norms>tol
        } // a3
      } // a2
    } // a1
  } // end routine
} // namespace pflow

// ***************************************************************************
// pflow::PrintBands
// ***************************************************************************
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void CalcPathLength(const pflow::projdata& pd, vector<double>& pathl,
      vector<double>& dpathl) {
    aurostd::matrix<double> rlat=RecipLat(pd.lat); //CO20200404 pflow::matrix()->aurostd::matrix()
    pathl[0]=0;
    dpathl[0]=0;
    for(int ik=1;ik<pd.nkpts;ik++) {
      vector<double> disp(3);
      disp=pflow::VVdiff(pd.kpts[ik],pd.kpts[ik-1]);
      dpathl[ik]=dpathl[ik-1]+pflow::norm(disp);
      disp=pflow::vecF2C(rlat,disp);
      pathl[ik]=pathl[ik-1]+pflow::norm(disp);
    }
  }// end
} // namespace pflow

namespace pflow {
  void PrintBands(const pflow::projdata& pd) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "pflow::PrintBands: BEGIN" << endl;
    ofstream outf_up("band.up.out");
    outf_up.precision(5);
    outf_up.setf(std::ios::fixed,std::ios::floatfield);
    outf_up.setf(std::ios::left,std::ios::adjustfield);
    vector<double> pathl(pd.nkpts);
    vector<double> dpathl(pd.nkpts);
    CalcPathLength(pd,pathl,dpathl);
    for(int ik=0;ik<pd.nkpts;ik++) {
      outf_up << pathl[ik] << " " << dpathl[ik];
      for(int ib=0;ib<pd.nbands;ib++) {
        outf_up << " " << pd.ener_k_b_u[ik][ib];
      }
      outf_up << "          " << ik;
      for(int ic=0;ic<3;ic++) {
        outf_up << " " << pd.kpts[ik][ic];
      }
      outf_up<<endl;
    }//kpt
    if(pd.sp) {
      ofstream outf_dn("band.dn.out");
      outf_dn.precision(5);
      outf_dn.setf(std::ios::fixed,std::ios::floatfield);
      outf_dn.setf(std::ios::left,std::ios::adjustfield);
      for(int ik=0;ik<pd.nkpts;ik++) {
        outf_dn << pathl[ik] << " " << dpathl[ik];
        for(int ib=0;ib<pd.nbands;ib++) {
          outf_dn << " " << pd.ener_k_b_d[ik][ib];
        }
        outf_dn << "          " << ik;
        for(int ic=0;ic<3;ic++) {
          outf_dn << " " << pd.kpts[ik][ic];
        }
        outf_dn<<endl;
      }//kpt
    }// if sp
  }// end
} // namespace pflow

// ***************************************************************************
// pflow::PrintCHGCAR
// ***************************************************************************
// Dane Morgan - Stefano Curtarolo
// edited by CO20140108
namespace pflow {
  bool PrintCHGCAR(const xstructure& str,
      const stringstream& chgcar_header,
      const vector<int>& ngrid,
      const vector<int>& format_dim,
      const vector<double>& chg_tot,
      const vector<double>& chg_diff,
      const string& output_name,
      ostream& oss) {  


    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) oss << __AFLOW_FUNC__ << "BEGIN" << endl;
    stringstream CHGCARout_ss;
    uint npts,natoms=pflow::GetNumAtoms(str),indx=0,numcolumns;

    // print POSCAR
    if(LDEBUG) oss << __AFLOW_FUNC__ << "PRINTING POSCAR" << endl;
    CHGCARout_ss << chgcar_header.str();
    // empty line
    CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
    // grid
    if(ngrid.size()!=3) {
      oss << endl;
      oss << __AFLOW_FUNC__ << "ERROR: ngrid must be length 3." << endl;
      oss << __AFLOW_FUNC__ << "Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    if(LDEBUG) oss << __AFLOW_FUNC__ << "PRINTING GRID" << endl;
    for(uint i=0;i<ngrid.size();i++) {
      CHGCARout_ss << "   " << ngrid.at(i);
    }
    CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
    npts=ngrid.at(0)*ngrid.at(1)*ngrid.at(2);
    // chg_tot
    if(LDEBUG) oss << __AFLOW_FUNC__ << "PRINTING CHG_TOT" << endl;
    for(uint i=0;i<npts;i+=format_dim.at(indx)) {
      // make sure we get exactly format_dim.at(indx) 0's and not more
      if(npts-i>=(uint)format_dim.at(indx)) {
        numcolumns=format_dim.at(indx);
      } else {
        numcolumns=npts-i;
      }
      for(uint j=0;j<numcolumns;j++) {    
        //CO START
        //chg_tot CAN be negative too 
        //http://cms.mpi.univie.ac.at/vasp-forum/viewtopic.php?t=4194&p=4194
        CHGCARout_ss << " " << std::setw(18) << std::setprecision(11) << std::scientific << std::uppercase << chg_tot.at(i+j);
        //CO END
      }
      CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
    }
    indx++;
    // augmentation occupancies, all set to 0
    if(format_dim.size()>1) {   // we have CHGCAR, not AECCAR
      if(LDEBUG) oss << __AFLOW_FUNC__ << "THIS IS A CHGCAR FILE, NOT AN AECCAR FILE" << endl;
      if(LDEBUG) oss << __AFLOW_FUNC__ << "PRINTING AUGMENTATION OCCUPANCIES (ALL SET TO 0)" << endl;
      for(uint i=1;i<(natoms+1);i++) {
        CHGCARout_ss << "augmentation occupancies " << i << " " << format_dim.at(indx) << endl;
        for(uint j=0;j<(uint)format_dim.at(indx);j+=(uint)format_dim.at(indx+1)) {
          // make sure we get exactly format_dim.at(indx) 0's and not more
          if((uint)format_dim.at(indx)-j>=(uint)format_dim.at(indx+1)) {
            numcolumns=format_dim.at(indx+1);
          } else {
            numcolumns=(uint)format_dim.at(indx)-j;
          }
          for(uint k=0;k<numcolumns;k++) {
            CHGCARout_ss << "  " << std::setprecision(7) << std::scientific << std::uppercase << 0.0;
          }
          CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
        }
        indx+=2;
      }
    }
    if(chg_diff.size()>0) {     // we have spin-polarized
      if(LDEBUG) oss << __AFLOW_FUNC__ << "THIS IS A SPIN-POLARIZED CALCULATION" << endl;
      // add line of 0's
      for(uint i=0;i<(uint)format_dim.at(indx);i++) { //mimic formatting of chg_diff for 0's line
        // standard shows +1 more 0
        CHGCARout_ss << "  " << std::setprecision(12) << std::scientific << std::uppercase << 0.0;
      }
      CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
      // grid
      if(LDEBUG) oss << __AFLOW_FUNC__ << "PRINTING GRID" << endl;
      for(uint i=0;i<ngrid.size();i++) {
        CHGCARout_ss << "   " << ngrid.at(i);
      }
      CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
      // chg_diff
      if(LDEBUG) oss << __AFLOW_FUNC__ << "PRINTING CHG_DIFF" << endl;
      for(uint i=0;i<npts;i+=format_dim.at(indx)) {
        // make sure we get exactly format_dim.at(indx) 0's and not more
        if(npts-i>=(uint)format_dim.at(indx)) {
          numcolumns=format_dim.at(indx);
        } else {
          numcolumns=npts-i;
        }
        for(uint j=0;j<numcolumns;j++) {    
          // chg_diff needs more space for minus sign
          CHGCARout_ss << " " << std::setw(18) << std::setprecision(11) << std::scientific << std::uppercase << chg_diff.at(i+j);
        }
        CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
      }
      indx++;
      // check for augmentation occupancies
      if(format_dim.size()>indx) { // we have more augmentation occupancies to write out, all set to 0
        if(LDEBUG) oss << __AFLOW_FUNC__ << "PRINTING AUGMENTATION OCCUPANCIES (ALL SET TO 0)" << endl;
        for(uint i=1;i<(natoms+1);i++) {
          CHGCARout_ss << "augmentation occupancies " << i << " " << format_dim.at(indx) << endl;
          for(uint j=0;j<(uint)format_dim.at(indx);j+=(uint)format_dim.at(indx+1)) {
            // make sure we get exactly format_dim.at(indx) 0's and not more
            if((uint)format_dim.at(indx)-j>=(uint)format_dim.at(indx+1)) {
              numcolumns=format_dim.at(indx+1);
            } else {
              numcolumns=(uint)format_dim.at(indx)-j;
            }
            for(uint k=0;k<numcolumns;k++) {
              CHGCARout_ss << "  " << std::setprecision(7) << std::scientific << std::uppercase << 0.0;
            }
            CHGCARout_ss << " " << endl;    //put space so AFLOW pick this up as a line
          }
          indx+=2;
        }

      }
    }

    // check if output file exists, delete it first
    if(aurostd::FileExist(output_name)) {
      if(LDEBUG) oss << __AFLOW_FUNC__ << "DELETING EXISTING OUTPUT FILE" << endl;
      aurostd::RemoveFile(output_name);
    }

    // write out to file
    if(LDEBUG) oss << __AFLOW_FUNC__ << "WRITING OUTPUT FILE" << endl;
    aurostd::stringstream2file(CHGCARout_ss,output_name);
    if(LDEBUG) oss << __AFLOW_FUNC__ << "DONE" << endl;
    return TRUE;
  }
} // namespace pflow

// ***************************************************************************
// pflow::PrintChgInt
// ***************************************************************************
// This funtion prints out 3 lattice vectors in cartesian coordinates.
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintChgInt(vector<aurostd::matrix<double> >& rad_chg_int,
      aurostd::matrix<double>& vor_chg_int, ostream& oss) {  //CO20200404 pflow::matrix()->aurostd::matrix()  
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "pflow::PrintChgInt: BEGIN" << endl;
    oss.setf(std::ios::left, std::ios::adjustfield);
    oss.setf(std::ios::fixed, std::ios::floatfield);

    // Output Voronoi data
    oss << "******************** Charge in Voronoi volume: Tot Diff Up Dn ********************" << endl;
    for(uint iat=0;iat<vor_chg_int.size();iat++) {
      oss << iat+1;
      for(int i=0;i<4;i++) {
        oss << " " << vor_chg_int[iat][i];
      }
      oss << endl;
    }

    // Output Radial data
    oss << "******************** Charge in sphere: Tot Diff Up Dn ********************" << endl;
    for(uint iat=0;iat<rad_chg_int.size();iat++) {
      oss << "Atom " << iat+1 << endl;
      for(uint ib=0;ib<rad_chg_int[iat].size();ib++) {
        oss << "  " << rad_chg_int[iat][ib][0];
        for(int i=1;i<5;i++) {
          oss << " " << rad_chg_int[iat][ib][i];
        }
        oss << endl;
      }
    }
  } // end routine
} // namespace pflow

// ***************************************************************************
// pflow::PrintCIF
// ***************************************************************************
// This funtion prints out structural data in a cif format.
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintCIF(ostream& oss,const xstructure& str,int _spacegroupnumber, int setting) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "pflow::PrintCIF: BEGIN" << endl;
    // Print out data in cif format
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(10);
    int spacegroupnumber=_spacegroupnumber;
    if(spacegroupnumber<1 || spacegroupnumber>230) spacegroupnumber=1;
    //title
    //CO MODS START 20170613 - fixes for BH cif
    string compound="",pd_name="";
    string tmp;
    if(!str.title.empty()){
      vector<string> tokens;
      aurostd::string2tokens(str.title,tokens," ");
      if(tokens.size()&&!tokens[0].empty()){
        if(aurostd::substring2bool(tokens[0],"/")){ //something out of LIBX: Ru_pvTa_pvY_sv/TFCC012.ABC, ICSD looks like this: Ag1F6Sb1
          tmp=tokens[0];
          aurostd::string2tokens(tmp,tokens,"/");
          compound=tokens[0];
          if(tokens.size()>1){pd_name=tokens[1];}  //just safety
        } else {
          compound=pd_name=tokens[0];
        }
      }
    }
    if(compound.empty()){
      if(!pd_name.empty()){compound=pd_name;}
      else {compound="1";}
    }
    //oss << "data_" << str.title << endl;
    oss << "# AFLOW.org Repositories" << endl;  //CO fix 20170606
    oss << "# " << str.title << endl;  //CO fix 20170606
    oss << "data_" << compound << endl;  //CO fix 20170606
    //CO try to get phase name, usually written as first part of title in AFLOW
    if(!pd_name.empty()){oss << "_pd_phase_name " << pd_name << endl;}  //just safety
    //CO MODS END 20170613 - fixes for BH cif
    xvector<double> data(6);
    data=Getabc_angles(str.lattice,DEGREES);data(1)*=str.scale;data(2)*=str.scale;data(3)*=str.scale;
    oss << "_cell_length_a  " << setw(12) << data(1) <<endl;
    oss << "_cell_length_b  " << setw(12) << data(2) <<endl;
    oss << "_cell_length_c  " << setw(12) << data(3) <<endl;
    oss << "_cell_angle_alpha  " << setw(9) << data(4) <<endl;
    oss << "_cell_angle_beta  " << setw(9) << data(5) <<endl;
    oss << "_cell_angle_gamma  " << setw(9) << data(6) <<endl;
    if(str.wyccar_ITC.size()!=0){
      oss << "_symmetry_space_group_name_H-M  '" << GetSpaceGroupName(spacegroupnumber,str.directory) << "'" << endl; //DX20180526 - add directory
      oss << "_symmetry_Int_Tables_Number  " << spacegroupnumber << endl;
    }
    else {
      // Output space group - Use P1 for now
      oss << "_symmetry_space_group_name_H-M  '" << GetSpaceGroupName(1,str.directory) << "'" << endl; //DX20180526 - add directory
      oss << "_symmetry_Int_Tables_Number  " << 1 << endl;
    }
    oss << "loop_" << endl;
    oss << " _symmetry_equiv_pos_site_id" << endl;
    //oss << " _symmetry_equiv_pos_as_xyz_" << endl;
    oss << " _symmetry_equiv_pos_as_xyz" << endl;  //CO fix 20170606

    // if the Wyckoff positions are stored/calculated
    if(str.wyccar_ITC.size()!=0){
      //DX20180801 - grab general Wyckoff equations
      string setting_string = aurostd::utype2string<uint>(setting);
      int general_wyckoff_multiplicity=0; // general Wyckoff position multiplicity
      vector<string> general_wyckoff_position; // general Wyckoff position equations
      SYM::getGeneralWyckoffMultiplicityAndPosition(spacegroupnumber, setting_string, general_wyckoff_multiplicity, general_wyckoff_position);

      for(uint i=0;i<general_wyckoff_position.size();i++){
        oss << "  " << i+1 << "  " << general_wyckoff_position[i] << endl;
      }

      //atom posns header
      oss << "loop_" << endl;
      oss << " _atom_site_label" << endl;
      oss << " _atom_site_occupancy" << endl;
      oss << " _atom_site_fract_x" << endl;
      oss << " _atom_site_fract_y" << endl;
      oss << " _atom_site_fract_z" << endl;
      //[DX+CO20220715 - OBSOLETE]oss << " _atom_site_thermal_displace_type" << endl; //https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_thermal_displace_type.html
      oss << " _atom_site_adp_type" << endl;  //DX+CO20220715
      oss << " _atom_site_B_iso_or_equiv" << endl;
      oss << " _atom_site_type_symbol" << endl;
      oss << " _atom_site_symmetry_multiplicity" << endl;
      oss << " _atom_site_Wyckoff_label" << endl;

      oss.precision(10);

      for(uint i=0;i<str.wyckoff_sites_ITC.size();i++){
        vector<string> tokens;
        aurostd::string2tokens(str.wyckoff_sites_ITC[i].wyckoffSymbol,tokens," ");
        string wyckoff_multiplicity = tokens[0];
        string wyckoff_letter = tokens[1];
        //DX20191010 [OBOSLETE] oss << str.wyckoff_sites_ITC[i].type << i+1 << " 1.0 "
        oss << str.wyckoff_sites_ITC[i].type << i+1 << " " << str.wyckoff_sites_ITC[i].site_occupation << " " //DX20191010
          << setw(12) << str.wyckoff_sites_ITC[i].coord(1) <<" "
          << setw(12) << str.wyckoff_sites_ITC[i].coord(2) <<" "
          << setw(12) << str.wyckoff_sites_ITC[i].coord(3) <<" "
          << " Biso 1.0 " << str.wyckoff_sites_ITC[i].type <<" "
          << wyckoff_multiplicity <<" "
          << wyckoff_letter << endl;
      }

      //      for(uint i=0;i<str.wyccar_ITC.size();i++){
      //        cerr << "str.wyccar_ITC[i]: " << str.wyccar_ITC[i] << endl;
      //        vector<string> tokens;
      //        aurostd::string2tokens(str.wyccar_ITC[i],tokens," ");
      //        cerr << "tokens: " << tokens.size() << endl;
      //        if(tokens.size()==7){
      //          oss << tokens[3] << i+1 << " 1.0 " 
      //              << setw(12) << tokens[0] <<" "
      //              << setw(12) << tokens[1] <<" "
      //              << setw(12) << tokens[2] <<" "
      //              << " Biso 1.0 " << tokens[3] <<" "
      //              << tokens[4] <<" "
      //              << tokens[5] << endl;
      //        }
      //      }
    }
    else {
      oss << "  1  x,y,z" << endl;

      //atom posns header
      oss << "loop_" << endl;
      oss << " _atom_site_label" << endl;
      oss << " _atom_site_occupancy" << endl;
      oss << " _atom_site_fract_x" << endl;
      oss << " _atom_site_fract_y" << endl;
      oss << " _atom_site_fract_z" << endl;
      //[DX+CO20220715 - OBSOLETE]oss << " _atom_site_thermal_displace_type" << endl; //https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_thermal_displace_type.html
      oss << " _atom_site_adp_type" << endl;  //DX+CO20220715
      oss << " _atom_site_B_iso_or_equiv" << endl;
      oss << " _atom_site_type_symbol" << endl;

      // Coordinates in direct

      oss.precision(10);
      for(uint i=0;i<str.atoms.size();i++) {
        //DX20191010 [OBSOLETE] oss << str.atoms.at(i).cleanname.c_str() << i+1 << " 1.0 "
        oss << str.atoms.at(i).cleanname.c_str() << i+1 << " " << str.atoms.at(i).partial_occupation_value << " " //DX20191010
          << setw(12) << str.atoms.at(i).fpos(1) <<" "
          << setw(12) << str.atoms.at(i).fpos(2) <<" "
          << setw(12) << str.atoms.at(i).fpos(3) <<" "
          << " Biso 1.0 " << str.atoms.at(i).cleanname.c_str() <<endl;
      }
    }
  }
} // namespace pflow

// ***************************************************************************
// pflow::PrintClat
// ***************************************************************************
// This funtion prints out 3 lattice vectors in cartesian coordinates.
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintClat(const xvector<double>& data, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    if(data.rows!=6) {
      init::ErrorOption("",__AFLOW_FUNC__,aurostd::liststring2string("data.size()=",aurostd::utype2string(data.rows)));
    }
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(10);
    xvector<double> abc_angles(data);
    xmatrix<double> lattice = GetClat(abc_angles); // angles in radiants
    roundoff(lattice,1.0e-4);
    oss << lattice(1,1) << " " << lattice(1,2) << " " << lattice(1,3) << endl;
    oss << lattice(2,1) << " " << lattice(2,2) << " " << lattice(2,3) << endl;
    oss << lattice(3,1) << " " << lattice(3,2) << " " << lattice(3,3) << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::PrintCmpStr
// ***************************************************************************
// Prints out the comparison between two structures.
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintCmpStr(const xstructure& str1, const xstructure& str2, const double& rcut, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "pflow::PrintCmpStr: BEGIN" << endl;
    oss.setf(std::ios::left, std::ios::adjustfield);
    oss.setf(std::ios::fixed, std::ios::floatfield);
    oss.precision(4);
    double cutoff=rcut;
    //  oss.width(10);
    // Rescale to avoid scale problems
    xstructure sstr1=str1;
    xstructure sstr2=str2;
    sstr1=ReScale(sstr1,1.0);
    sstr2=ReScale(sstr2,1.0);
    int nat1=sstr1.atoms.size();
    int nat2=sstr2.atoms.size();
    std::deque<int> num_each_type1=sstr1.num_each_type;
    std::deque<int> num_each_type2=sstr2.num_each_type;
    int ntyp_min;
    int ntyp_max;
    int nprs_min;
    int nprs_max;

    oss << "*************************  Structure Comparison *************************" << endl;

    // Nun atoms
    oss << "Number of atoms" << endl;
    oss << "  Num Atom Str1: " << nat1 << endl;
    oss << "  Num Atom Str2: " << nat2 << endl;
    oss << "  Num Atom %Err: " << 100*2*(nat2-nat1)/(nat2+nat1) << endl;

    // Num types
    int ntyp1=sstr1.num_each_type.size();
    int ntyp2=sstr2.num_each_type.size();
    oss << "Number of types" << endl;
    oss << "  Num Type Str1: " << ntyp1 << endl;
    oss << "  Num Type Str2: " << ntyp2 << endl;
    oss << "  Num Type Err : " << ntyp2-ntyp1 << endl;

    // Stoichiometry
    int wid=8;
    ntyp_max=max(ntyp1,ntyp2);
    vector<double> tmp_stoich1(ntyp_max,0);
    vector<double> tmp_stoich2(ntyp_max,0);
    oss << "Stoichiometry" << endl;
    oss << "  Type:        ";
    for(int it=0;it<ntyp_max;it++) {
      oss << setw(wid) << it << " ";
    }
    oss << endl;
    oss << "  Stoich Str1: ";
    for(int it=0;it<ntyp1;it++) {
      tmp_stoich1[it]=(double) num_each_type1[it]/ (double) nat1;
      oss << setw(wid) << tmp_stoich1[it] << " ";
    }
    oss << endl;
    oss << "  Stoich Str2: ";
    for(int it=0;it<ntyp2;it++) {
      tmp_stoich2[it]=(double) num_each_type2[it]/ (double) nat2;
      oss << setw(wid) << tmp_stoich2[it] << " ";
    }
    oss << endl;
    oss << "  Stoich Diff: ";
    for(int it=0;it<ntyp_max;it++) {
      oss << setw(wid) << tmp_stoich2[it]-tmp_stoich1[it] << " ";
    }
    oss << endl;

    // Volume
    double v1=GetVol(sstr1.lattice)/nat1;
    double v2=GetVol(sstr2.lattice)/nat2;
    oss << "Volume per atom" << endl;
    oss << "  Vol Str1: " << v1 << endl;
    oss << "  Vol Str2: " << v2 << endl;
    oss << "  Vol %Err: " << 100*2*(v2-v1)/(v2+v1) << endl;

    // lat par and angles
    // vector<double> ldat1=pflow::Getabc_angles(aurostd::xmatrix2matrix(sstr1.lattice));  //CO20200404 pflow::matrix()->aurostd::matrix()
    vector<double> ldat1=xvector2vector(Getabc_angles(sstr1.lattice,DEGREES));
    // vector<double> ldat2=pflow::Getabc_angles(aurostd::xmatrix2matrix(sstr2.lattice));  //CO20200404 pflow::matrix()->aurostd::matrix()
    vector<double> ldat2=xvector2vector(Getabc_angles(sstr2.lattice,DEGREES));

    ldat1=pflow::Sort_abc_angles(ldat1);  
    ldat2=pflow::Sort_abc_angles(ldat2);  
    vector<double> ldat_err(7,0.0);
    ldat_err[0]=100*2*(v2*nat2-v1*nat1)/(v2*nat2+v1*nat1);
    for(int i=0;i<6;i++) {
      double tmp=(abs(ldat2[i])+abs(ldat1[i]));
      if(tmp>0) {
        ldat_err[i+1]=100*2*(ldat2[i]-ldat1[i])/(ldat2[i]+ldat1[i]);
      }
      else {
        ldat_err[i+1]=0.0;
      }
    }
    oss << "Unit cell: tot_vol a b c alpha beta gamma" << endl;
    oss << "  Cell Str1: ";
    oss << v1*nat1 << " ";
    pflow::Vout(ldat1,oss);
    oss << "  Cell Str2: ";
    oss << v2*nat2 << " ";
    pflow::Vout(ldat2,oss);
    oss << "  Cell %Err: ";
    pflow::Vout(ldat_err,oss);

    // Volume adjusted distances
    double vavg=(v1+v2)/2.0;
    double rescale1=std::pow((double) vavg/v1,(double) 1.0/3.0);
    double rescale2=std::pow((double) vavg/v2,(double) 1.0/3.0);
    xstructure str1_newvol=sstr1;
    xstructure str2_newvol=sstr2;
    str1_newvol=SetScale(str1_newvol,rescale1);
    str2_newvol=SetScale(str2_newvol,rescale2);
    aurostd::matrix<double> dist1; //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> dist2; //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> dist_diff; //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> dist_diff_n; //CO20200404 pflow::matrix()->aurostd::matrix()
    // If cutoff < 0 use default of 4*(avg rad of an atom).  Note
    // that in a rescaled str where each atom has vavg=1, rad=.6204.
    // In the true str, where each atom has vol vavg, rad=0.6204*vavg^(1/3).
    if(rcut<0) cutoff=4*0.6204*std::pow((double) vavg,(double) 1.0/3.0);
    pflow::CmpStrDist(str1_newvol,str2_newvol,cutoff,dist1,dist2,dist_diff,dist_diff_n);
    ntyp_min=min(ntyp1,ntyp2);
    nprs_min=dist_diff.size();
    oss << "Avg of Distance Magnitude Differences" << endl;
    oss << "  Rescaling vol/atom of Str1 and Str2 to the average vol/atom: " << vavg << endl;
    oss << "  Working with the all the types up to num types: " << ntyp_min << endl;
    oss << "  This corresponds to the following number of distinct pairs: " << nprs_min << endl;
    double tadiff=0;
    double tadiff_n=0;
    int tnn=0;
    for(int it1=0;it1<ntyp_min;it1++) {
      for(int it2=it1;it2<ntyp_min;it2++) {
        int id=it2-it1+ntyp_min*it1-max(0,it1*(it1-1)/2);
        int nn=dist_diff[id].size();      
        oss << "  P" << id << " Pair between types: "
          << it1 << " " << it2 <<endl;
        oss << "    P" << id << " Comparing neighbors within cutoff (number): "
          << cutoff << " (" << nn << ")" <<  endl;
        double adiff=0;
        double adiff_n=0;
        // double perc_err=0;  //DM not used
        for(int in=0;in<nn;in++) {
          adiff=adiff+abs(dist_diff[id][in]);
          adiff_n=adiff_n+abs(dist_diff_n[id][in]);
        }
        tadiff=tadiff+adiff;
        tadiff_n=tadiff_n+adiff_n;
        if(nn>0) {
          adiff=adiff/nn;
          adiff_n=adiff_n/nn;
        }
        else {
          adiff=0.0;
          adiff_n=0.0;
        }
        tnn=tnn+nn;
        oss << "    P" << id << " Avg Difference: " << adiff << endl;
        oss << "    P" << id << " Avg Diff/vol_avg^1/3 %: "
          << 100*adiff/std::pow((double) vavg,(double) 1.0/3.0) << endl;
        oss << "    P" << id << " Avg Difference/bond length %: " << 100*adiff_n << endl;
      }// it2
    }// it1
    if(tnn>0) {
      tadiff=tadiff/tnn;
      tadiff_n=tadiff_n/tnn;
    }
    else {
      tadiff=0;
      tadiff_n=0;
    }
    oss << "  Totals for all pairs, compared by type, within cutoff (number): "
      << cutoff << " (" << tnn << ")" <<  endl;
    oss << "    Tot Avg Difference: " << tadiff << endl;
    oss << "    Tot Avg Diff/vol_avg^1/3 %: "
      << 100*tadiff/std::pow((double) vavg,(double) 1.0/3.0) << endl;
    oss << "    Tot Avg Diff/bond_length %: "
      << 100*tadiff_n << endl;

    /* THIS NEEDS WORK - need to sum the rdfsh over all atoms of one type,
       then output average numbers.  This is still not perfect, since you could
       get cancelling errors.  Could try to do site to site comparison.
    // Compare radial distribution functions
    ntyp_max=max(ntyp1,ntyp2);
    nprs_max=(ntyp_max+1)*(ntyp_max+2)/2;
    vector<double> tmpsh1(nprs_max,0);
    vector<double> tmpsh2(nprs_max,0);
    double bin_width=0.01;
    int nbins=Nint(cutoff/bin_width);
    int smooth_width=5;
    int nsh_max=2;
    aurostd::matrix<double> rdf_all_1; //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> rdf_all_2; //CO20200404 pflow::matrix()->aurostd::matrix()
    GetRDF_20220101(sstr1,cutoff,nbins,rdf_all_1);
    GetRDF_20220101(sstr2,cutoff,nbins,rdf_all_2);
    aurostd::matrix<double> rdf_all_1_sm=GetSmoothRDF(rdf_all_1,smooth_width); //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> rdf_all_2_sm=GetSmoothRDF(rdf_all_2,smooth_width); //CO20200404 pflow::matrix()->aurostd::matrix()
    // Get shells
    aurostd::matrix<double> rdfsh_all_1; //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> rdfsh_loc_1; // Radial location of rdf shells. //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> rdfsh_all_2; //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> rdfsh_loc_2; // Radial location of rdf shells. //CO20200404 pflow::matrix()->aurostd::matrix()
    GetRDFShells(sstr1,cutoff,nbins,smooth_width,
    rdf_all_1_sm,rdfsh_all_1,rdfsh_loc_1);
    GetRDFShells(sstr2,cutoff,nbins,smooth_width,
    rdf_all_2_sm,rdfsh_all_2,rdfsh_loc_2);
    // print out position, number for each radial shell up to nsh_max.
    // Print out for each type of pair
    oss << "Radial shell comparison" << endl;
    for(int ish=0;ish<nsh_max;ish++) { // Radial shells loop.
    oss << "  Shell " << ish << endl;
    oss << "                    ";
    for(int it1=0;it1<ntyp_max;it1++) {
    for(int it2=it1;it2<ntyp_max;it2++) {
    oss << it1 << "-" << it2 << "  ";
    }
    }
    oss << endl;
    oss << "    Shell Occ Str1: ";
    for(int it1=0;it1<ntyp1;it1++) {
    for(int it2=0;it2<ntyp1;it2++) {
    int id=it2+(ntyp1+1)*it1;
    int idb=ntyp1*nat1+it2+(ntyp1+1)*it1;
    cout << id<<" "<<idb<< " "<<rdfsh_all_1.size()<<" "<<ntyp1<<" "<< nat1<< endl;
    tmpsh1[id]=rdfsh_all_1[ntyp1*nat1+id][ish];
    oss << tmpsh1[id] << " ";
    }
    }
    oss << endl;
    }    

*/

    // Number of each type up to rcut.
    nprs_max=ntyp_max*(ntyp_max+1)/2;
    vector<int> n1(nprs_max,0);
    vector<int> n2(nprs_max,0);
    oss << "Number neighbors of each type within cutoff for rescaled volume" << endl;
    oss << "  Cutoff: " << cutoff << endl;
    oss << "  Rescaling vol/atom of Str1 and Str2 to the average vol/atom: " << vavg << endl;
    oss << "           ";
    for(int it1=0;it1<ntyp_max;it1++) {
      for(int it2=it1;it2<ntyp_max;it2++) {
        oss << it1 << "-" << it2 << "  ";
      }
    }
    oss << endl;
    oss << "  NN Str1  ";
    for(int it1=0;it1<ntyp_max;it1++) {
      for(int it2=it1;it2<ntyp_max;it2++) {
        if(it1<ntyp1 && it2<ntyp1) {
          int id1=it2-it1+ntyp1*it1-max(0,it1*(it1-1)/2);
          n1[id1]=dist1[id1].size();
          oss << setw(4) << n1[id1] << " ";
        }
      }
    }
    oss << endl;
    oss << "  NN Str2  ";
    for(int it1=0;it1<ntyp_max;it1++) {
      for(int it2=it1;it2<ntyp_max;it2++) {
        if(it1<ntyp2 && it2<ntyp2) {
          int id2=it2-it1+ntyp2*it1-max(0,it1*(it1-1)/2);
          n2[id2]=dist2[id2].size();
          oss << setw(4) << n2[id2] << " ";
        }
      }
    }
    oss << endl;  
    oss << "  NN Err   ";
    for(int it1=0;it1<ntyp_max;it1++) {
      for(int it2=it1;it2<ntyp_max;it2++) {
        int id=it2-it1+ntyp_max*it1-max(0,it1*(it1-1)/2);
        oss << setw(4) << n2[id]-n1[id] << " ";
      }
    }
    oss << endl;  

    // Space group
    bool Platon_EQUAL=FALSE,Platon_EXACT=FALSE;
    double Platon_ang=1.0,Platon_d1=0.25,Platon_d2=0.25,Platon_d3=0.25;
    system("touch pflow_tmp_in"); // Create input file to connect to stream.
    ofstream outfile("pflow_tmp_out");
    ifstream infile("pflow_tmp_in");
    string spcgrp1a,spcgrp2a,spcgrp1b,spcgrp2b;
    platon2print(str1,Platon_EQUAL,Platon_EXACT,Platon_ang,Platon_d1,Platon_d2,Platon_d3,outfile);
    system("cat pflow_tmp_out | platonSG > pflow_tmp_in");
    infile >> spcgrp1a >> spcgrp1b;
    // Puts file position pointers back to beginning.
    outfile.seekp(0);
    infile.seekg(0);
    platon2print(sstr2,Platon_EQUAL,Platon_EXACT,Platon_ang,Platon_d1,Platon_d2,Platon_d3,outfile);
    system("cat pflow_tmp_out | platonSG > pflow_tmp_in");
    infile >> spcgrp2a >> spcgrp2b;
    system("/bin/rm pflow_tmp_out");
    system("/bin/rm pflow_tmp_in");
    oss << "Space groups " << endl;
    oss << "  SG Str1: " << spcgrp1a << " " << spcgrp1b << endl;
    oss << "  SG Str2: " << spcgrp2a << " " << spcgrp2b << endl;

    // Center of mass  
  } // end routine
} // namespace pflow

// ***************************************************************************
// pflow::PrintData()
// ***************************************************************************
// This funtion prints out structural data.
// Stefano Curtarolo
// Modified by David Hicks (DX) - DX20210301
//  - The symmetry descriptors are now calculated in
//    GetExtendedCrystallographicData() (in aflow_xatom).
//  - Return type changed from void to string
//  - Furthermore, each symmetry descriptor (real, reciprocal, superlattice)
//    has their own printing function, enabling them to be called separately.
//    The xoptions (EDATA::<PROPERTY>) are also updated in these separate
//    functions (necessary for lib2raw runs).
namespace pflow {
  string PrintData(const xstructure& xstr, const string& smode, filetype ftype, bool already_calculated, double sym_eps, bool no_scan, int setting){
    aurostd::xoption vpflow;
    return PrintData(xstr, vpflow, smode, ftype, already_calculated, sym_eps, no_scan, setting);
  }
}

namespace pflow {
  string PrintData(const xstructure& xstr, aurostd::xoption& vpflow, const string& smode, filetype ftype, bool already_calculated, double sym_eps, bool no_scan, int setting){
    xstructure str_sym,str_sp,str_sc;
    return PrintData(xstr, str_sym, str_sp, str_sc, vpflow, smode, ftype, already_calculated, sym_eps, no_scan, setting);
  }
}
namespace pflow {
  string PrintData(const xstructure& xstr,
      xstructure& str_sp,
      xstructure& str_sc,
      aurostd::xoption& vpflow,
      const string& smode,
      filetype ftype,
      bool already_calculated,
      double sym_eps,
      bool no_scan,
      int setting){
    xstructure str_sym;
    return PrintData(xstr, str_sym, str_sp, str_sc, vpflow, smode, ftype, already_calculated, sym_eps, no_scan, setting);
  }
}

namespace pflow {
  string PrintData(const xstructure& xstr,
      xstructure& str_sym,
      xstructure& str_sp,
      xstructure& str_sc,
      const string& smode,
      filetype ftype,
      bool already_calculated,
      double sym_eps,
      bool no_scan,
      int setting){
    aurostd::xoption vpflow;
    return PrintData(xstr, str_sym, str_sp, str_sc, vpflow, smode, ftype, already_calculated, sym_eps, no_scan, setting);
  }
}

namespace pflow {
  string PrintData(const xstructure& xstr,
      xstructure& str_sym,
      xstructure& str_sp,
      xstructure& str_sc,
      aurostd::xoption& vpflow,
      const string& smode,
      filetype ftype,
      bool already_calculated,
      double sym_eps,
      bool no_scan,
      int setting){

    bool LDEBUG=(FALSE || XHOST.DEBUG);

    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;

    stringstream oss;
    oss.setf(std::ios::fixed,std::ios::floatfield);

    bool standalone=false; // edata combines all symmetries

    // ---------------------------------------------------------------------------
    // if just lattice data (DATA, i.e., no symmetry), perform and return
    if(smode=="DATA"){
      standalone = true; // if smode==DATA, then it is standalone
      oss << PrintRealLatticeData(xstr, smode, ftype, standalone, already_calculated, sym_eps);
      return oss.str();
    }

    xstructure str_aus(xstr);
    str_aus.ReScale(1.0);

    // ---------------------------------------------------------------------------
    // calculate the extended crystallographic data information
    // (real, space group, reciprocal, and superlattice symmetries)
    if(!already_calculated){
      str_aus.GetExtendedCrystallographicData(str_sp, str_sc, sym_eps, no_scan, setting);
      already_calculated=true;
      if(str_aus.sym_eps != sym_eps){ sym_eps = str_aus.sym_eps; }
    }

    // ---------------------------------------------------------------------------
    // now print the symmetry information
    // save to vector of strings so we can combine based on txt or json format
    // later (e.g., jsons need to be joined with commas)
    vector<string> content;

    // real lattice data
    content.push_back(PrintRealLatticeData(str_aus, vpflow, smode, ftype, standalone, already_calculated, sym_eps));
    // point group data
    content.push_back(PrintCrystalPointGroupData(str_aus, vpflow, ftype, standalone, already_calculated, sym_eps));
    // space group data
    content.push_back(PrintSGData(str_aus, vpflow, ftype, standalone, already_calculated, sym_eps, no_scan, setting));
    // lattice (no atoms) data
    content.push_back(PrintLatticeLatticeData(str_aus, vpflow, ftype, standalone, already_calculated, sym_eps));
    // superlattice data
    content.push_back(PrintSuperlatticeData(str_aus, vpflow, ftype, standalone, already_calculated, sym_eps));
    // reciprocal lattice data
    content.push_back(PrintReciprocalLatticeData(str_aus, vpflow, ftype, standalone, already_calculated, sym_eps));

    // control printing
    if(ftype == txt_ft){
      oss << aurostd::joinWDelimiter(content,"");
      oss << "SPRIM" << endl;
      oss << str_sp;
      oss << "SCONV" << endl;
      oss << str_sc;
    }
    else if(ftype == json_ft){
      aurostd::JSONwriter json;
      json.mergeRawJSON(aurostd::joinWDelimiter(content,","));
      json.addRawJSON("standard_primitive_structure", xstructure2json(str_sp)); // add sprim
      json.addRawJSON("standard_conventional_structure", xstructure2json(str_sc)); // add sconv
      oss << json.toString(true) << endl;
    }

    // ---------------------------------------------------------------------------
    // update vpflow flag (note: other edata attributes are updated in the
    // individual lattice-data functions)
    if(vpflow.flag("EDATA::CALCULATED")){ //DX20180823 - if calculated remove
      vpflow.pop_attached("EDATA::SPRIM");
      vpflow.pop_attached("EDATA::SCONV");
    }
    vpflow.flag("EDATA::CALCULATED",TRUE);
    stringstream ss_str_sp; ss_str_sp << str_sp << endl;
    vpflow.push_attached("EDATA::SPRIM",ss_str_sp.str());
    stringstream ss_str_sc; ss_str_sc << str_sc << endl;
    vpflow.push_attached("EDATA::SCONV",ss_str_sc.str());

    // ---------------------------------------------------------------------------
    // save xstructure with symmetry info to str_sym
    str_sym=str_aus;

    return oss.str();
  }
}

//DX20210301 [OBSOLETE] namespace pflow {
//DX20210301 [OBSOLETE]   void PrintData(const xstructure& str,xstructure& str_sym,xstructure& str_sp,xstructure& str_sc, ostream& oss,const string& smode, const string& format,bool already_calculated) { //CO20171027 //DX20180806 - added space group setting
//DX20210301 [OBSOLETE]     //DX20181215 [OBSOLETE] double tolerance = AUROSTD_NAN;
//DX20210301 [OBSOLETE]     double tolerance = SYM::defaultTolerance(str_sym); //DX20180215 START with default tolerance //DX20180226 - str_sym not str_sp
//DX20210301 [OBSOLETE]     if(already_calculated){tolerance=str_sp.sym_eps;} //CO20171025
//DX20210301 [OBSOLETE]     bool no_scan = false;
//DX20210301 [OBSOLETE]     int sg_setting = 1;
//DX20210301 [OBSOLETE]     PrintData(str,str_sym,str_sp,str_sc,oss,tolerance,smode,no_scan,sg_setting,format,already_calculated); //CO20171027
//DX20210301 [OBSOLETE]   }
//DX20210301 [OBSOLETE] }
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE] //DX20180822 - add xoption - START
//DX20210301 [OBSOLETE] namespace pflow {
//DX20210301 [OBSOLETE]   void PrintData(const xstructure& str,xstructure& str_sym,xstructure& str_sp,xstructure& str_sc, ostream& oss,aurostd::xoption& vpflow,const string& smode,const string& format,bool already_calculated) {
//DX20210301 [OBSOLETE]     //DX20181215 [OBSOLETE] double tolerance = AUROSTD_NAN;
//DX20210301 [OBSOLETE]     double tolerance = SYM::defaultTolerance(str_sym); //DX20180215 START with default tolerance //DX20180226 - str_sym not str_sp
//DX20210301 [OBSOLETE]     if(already_calculated){tolerance=str_sp.sym_eps;} //CO20171025
//DX20210301 [OBSOLETE]     bool no_scan = false;
//DX20210301 [OBSOLETE]     int sg_setting = 1;
//DX20210301 [OBSOLETE]     PrintData(str,str_sym,str_sp,str_sc,oss,vpflow,tolerance,smode,no_scan,sg_setting,format,already_calculated); //CO20171027
//DX20210301 [OBSOLETE]   }
//DX20210301 [OBSOLETE] }
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE] namespace pflow {
//DX20210301 [OBSOLETE]   void PrintData(const xstructure& str,xstructure& str_sym,xstructure& str_sp,xstructure& str_sc, ostream& oss,double tolerance,const string& smode, bool no_scan, int sg_setting, const string& format,bool already_calculated) {
//DX20210301 [OBSOLETE]     aurostd::xoption vpflow;
//DX20210301 [OBSOLETE]     PrintData(str,str_sym,str_sp,str_sc,oss,vpflow,tolerance,smode,no_scan,sg_setting,format,already_calculated);
//DX20210301 [OBSOLETE]   }
//DX20210301 [OBSOLETE] }
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE] //DX20180822 - add xoption - END
//DX20210301 [OBSOLETE] namespace pflow {
//DX20210301 [OBSOLETE]   void PrintData(const xstructure& str,xstructure& str_sym,xstructure& str_sp,xstructure& str_sc, ostream& oss_final,aurostd::xoption& vpflow,double tolerance,const string& smode, bool no_scan, int sg_setting, const string& format,bool already_calculated) { //CO20171027 //DX20180822 -added xoption
//DX20210301 [OBSOLETE]     bool LDEBUG=(FALSE || XHOST.DEBUG);
//DX20210301 [OBSOLETE]     if(LDEBUG) cerr << XPID << "pflow::PrintData: BEGIN" << endl;
//DX20210301 [OBSOLETE]     // smode=="DATA" or "EDATA"
//DX20210301 [OBSOLETE]     // Print out structural data
//DX20210301 [OBSOLETE]     stringstream oss;
//DX20210301 [OBSOLETE]     oss.setf(std::ios::fixed,std::ios::floatfield);
//DX20210301 [OBSOLETE]     double vol;
//DX20210301 [OBSOLETE]     vol=GetVol(str.lattice)*str.scale*str.scale*str.scale;
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]     bool symmetry_commensurate = false;
//DX20210301 [OBSOLETE]     bool force_perform = true;  // Force a result to be found at the end, even if incommensurate
//DX20210301 [OBSOLETE]     double orig_tolerance = tolerance;
//DX20210301 [OBSOLETE]     uint sym_eps_change_count = 0; //DX20180226 - added sym eps change count
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]     xstructure str_aus(str);//,str_sp,str_sc; //CO20171027
//DX20210301 [OBSOLETE]     str_aus.ReScale(1.0); //DX20210211
//DX20210301 [OBSOLETE]     while(symmetry_commensurate == false){
//DX20210301 [OBSOLETE]       if(smode!="EDATA" || no_scan){
//DX20210301 [OBSOLETE]         symmetry_commensurate=true;
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       // FORMAT = TXT
//DX20210301 [OBSOLETE]       if(format=="txt"){
//DX20210301 [OBSOLETE]         oss << "REAL LATTICE" << endl;
//DX20210301 [OBSOLETE]         oss << " Real space a b c alpha beta gamma: ";
//DX20210301 [OBSOLETE]         oss.precision(10);
//DX20210301 [OBSOLETE]         xvector<double> data(6);
//DX20210301 [OBSOLETE]         data=Getabc_angles(str.lattice,DEGREES);data(1)*=str.scale;data(2)*=str.scale;data(3)*=str.scale;
//DX20210301 [OBSOLETE]         oss.precision(10);oss << data(1) << " " << data(2) << " " << data(3) << " ";
//DX20210301 [OBSOLETE]         oss.precision(3); oss << data(4) << " " << data(5) << " " << data(6) << endl;
//DX20210301 [OBSOLETE]         oss.precision(4);
//DX20210301 [OBSOLETE]         if(smode=="EDATA") {
//DX20210301 [OBSOLETE]           oss << " Real space a b c alpha beta gamma: ";
//DX20210301 [OBSOLETE]           oss.precision(10);oss << data(1)*angstrom2bohr << " " << data(2)*angstrom2bohr << " " << data(3)*angstrom2bohr << " ";
//DX20210301 [OBSOLETE]           oss.precision(3); oss << data(4) << " " << data(5) << " " << data(6) << "   Bohrs/Degs " << endl;
//DX20210301 [OBSOLETE]           oss.precision(4);
//DX20210301 [OBSOLETE]         }
//DX20210301 [OBSOLETE]         oss << " Real space Volume: " << vol << endl;
//DX20210301 [OBSOLETE]         oss << " Real space c/a = " << data(3)/data(1) << endl;
//DX20210301 [OBSOLETE]         str_aus=str; //,str_sp,str_sc; //CO20171027
//DX20210301 [OBSOLETE]         //DX20170901 -Add tolerance - START
//DX20210301 [OBSOLETE]         str_aus.sym_eps = str_sp.sym_eps = str_sc.sym_eps = tolerance;
//DX20210301 [OBSOLETE]         str_aus.sym_eps_change_count = str_sp.sym_eps_change_count = str_sc.sym_eps_change_count = sym_eps_change_count; //DX20180226 - added sym eps change count
//DX20210301 [OBSOLETE]         //DX20170901 -Add tolerance - END
//DX20210301 [OBSOLETE]         if(smode=="EDATA") {
//DX20210301 [OBSOLETE]           // str_aus.CalculateSymmetryPointGroup(FALSE);
//DX20210301 [OBSOLETE]           // str_aus.CalculateSymmetryFactorGroup(FALSE);
//DX20210301 [OBSOLETE]           // str_aus.CalculateSymmetryPointGroupCrystal(FALSE);
//DX20210301 [OBSOLETE]           if(!already_calculated){  //CO20171025
//DX20210301 [OBSOLETE]             str_aus.GetLatticeType(str_sp,str_sc);
//DX20210301 [OBSOLETE]             if(orig_tolerance == AUROSTD_NAN){
//DX20210301 [OBSOLETE]               orig_tolerance = str_sp.sym_eps;  // Use new tolerance calculated here
//DX20210301 [OBSOLETE]             }
//DX20210301 [OBSOLETE]             tolerance = str_sp.sym_eps;
//DX20210301 [OBSOLETE]             sym_eps_change_count = str_aus.sym_eps_change_count = str_sp.sym_eps_change_count; //DX20180226 - added sym eps change count
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           oss << "BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal)" << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Bravais Lattice Primitive        = " << str_aus.bravais_lattice_type << endl;// " " << str.title << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Lattice Variation                = " << str_aus.bravais_lattice_variation_type << endl;//WSETYAWAN mod
//DX20210301 [OBSOLETE]           oss << " Real space: Lattice System                   = " << str_aus.bravais_lattice_system << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Pearson Symbol                   = " << str_aus.pearson_symbol << endl;
//DX20210301 [OBSOLETE]           oss << "POINT GROUP CRYSTAL" << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Crystal Family                   = " << str_aus.crystal_family << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Crystal System                   = " << str_aus.crystal_system << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Crystal Class                    = " << str_aus.point_group_crystal_class << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Point Group (Hermann Mauguin)    = " << str_aus.point_group_Hermann_Mauguin << endl; // "      PGXTAL" << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Point Group (Schoenflies)        = " << str_aus.point_group_Shoenflies << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Point Group Orbifold             = " << str_aus.point_group_orbifold << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Point Group Type                 = " << str_aus.point_group_type << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Point Group Order                = " << str_aus.point_group_order << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Point Group Structure            = " << str_aus.point_group_structure << endl;
//DX20210301 [OBSOLETE]           PrintSGData(str_aus,tolerance,oss,vpflow,no_scan,sg_setting,false,format,already_calculated);  //DX20170831 - SGDATA //CO20171025 //CO20171027 //DX20180823 - added xoption
//DX20210301 [OBSOLETE]           //DX20180608 - recalculate GetLatticeType() before changing tolerance - START
//DX20210301 [OBSOLETE]           if(sym_eps_change_count != str_aus.sym_eps_change_count){
//DX20210301 [OBSOLETE]             if(LDEBUG) {
//DX20210301 [OBSOLETE]               cerr << XPID << "pflow::PrintData: WARNING: PrintSGData() changed the tolerance, need to recalculate GetLatticeType() at the new tolerance (i.e., sym_eps=" << str_aus.sym_eps << ") [dir=" << str.directory << "]" << endl; //DX20180526 - add directory
//DX20210301 [OBSOLETE]             }
//DX20210301 [OBSOLETE]             sym_eps_change_count = str_aus.sym_eps_change_count;
//DX20210301 [OBSOLETE]             str_aus.ClearSymmetry();
//DX20210301 [OBSOLETE]             oss.str("");
//DX20210301 [OBSOLETE]             continue;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           //DX20180608 - recalculate GetLatticeType() before changing tolerance - END
//DX20210301 [OBSOLETE]           sym_eps_change_count = str_aus.sym_eps_change_count; //DX20180226 - added sym eps change count
//DX20210301 [OBSOLETE]           //DX20170901 - Add consistency check for input symmetry method and ITC method - START
//DX20210301 [OBSOLETE]           int multiplicity_of_primitive=str_sp.fgroup.size()/str_sp.pgroup_xtal.size();
//DX20210301 [OBSOLETE]           bool derivative_structure=false;
//DX20210301 [OBSOLETE]           if(SYM::ComparePointGroupAndSpaceGroupString(str_aus,multiplicity_of_primitive,derivative_structure) || no_scan){
//DX20210301 [OBSOLETE]             symmetry_commensurate=true;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           else {
//DX20210301 [OBSOLETE]             if(LDEBUG) {
//DX20210301 [OBSOLETE]               cerr << XPID << "pflow::PrintData: WARNING: Space group symbol and point group symbol do not match. (sg=" << GetSpaceGroupName(str_aus.space_group_ITC,str_aus.directory) << ", pg=" << str_aus.point_group_Hermann_Mauguin << ") [dir=" << str.directory << "]" << endl; //DX20180526 - add directory
//DX20210301 [OBSOLETE]             }
//DX20210301 [OBSOLETE]             str_aus.ClearSymmetry();
//DX20210301 [OBSOLETE]             oss.str("");
//DX20210301 [OBSOLETE]             //if(!SYM::change_tolerance(str_aus,tolerance,orig_tolerance,change_sym_count,str_aus.dist_nn_min,no_scan))
//DX20210301 [OBSOLETE]             if(!SYM::change_tolerance(str_aus,tolerance,str_aus.dist_nn_min,no_scan)){
//DX20210301 [OBSOLETE]               if(force_perform){
//DX20210301 [OBSOLETE]                 if(LDEBUG) {
//DX20210301 [OBSOLETE]                   cerr << XPID << "pflow::PrintData: WARNING: Scan failed. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies). [dir=" << str.directory << "]" << endl;
//DX20210301 [OBSOLETE]                 }
//DX20210301 [OBSOLETE]                 no_scan=true; //Force it to continue
//DX20210301 [OBSOLETE]                 tolerance = orig_tolerance;
//DX20210301 [OBSOLETE]                 str_aus.ClearSymmetry();
//DX20210301 [OBSOLETE]                 oss.str("");
//DX20210301 [OBSOLETE]                 continue;
//DX20210301 [OBSOLETE]               }
//DX20210301 [OBSOLETE]             }
//DX20210301 [OBSOLETE]             else {
//DX20210301 [OBSOLETE]               str_aus.ClearSymmetry();
//DX20210301 [OBSOLETE]               oss.str("");
//DX20210301 [OBSOLETE]               sym_eps_change_count = str_aus.sym_eps_change_count; //DX20180226 - added sym eps change count
//DX20210301 [OBSOLETE]               continue;
//DX20210301 [OBSOLETE]             }
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           //DX20170901 - Add consistency check for input symmetry method and ITC method - END
//DX20210301 [OBSOLETE]           oss << "BRAVAIS LATTICE OF THE LATTICE (pgroup)" << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Bravais Lattice Primitive        = " << str_aus.bravais_lattice_lattice_type << endl;// " " << str.title << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Lattice Variation                = " << str_aus.bravais_lattice_lattice_variation_type << endl; //WSETYAWAN mod
//DX20210301 [OBSOLETE]           oss << " Real space: Lattice System                   = " << str_aus.bravais_lattice_lattice_system << endl;
//DX20210301 [OBSOLETE]           oss << "SUPERLATTICE (equally decorated)" << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Bravais Superlattice Primitive   = " << str_aus.bravais_superlattice_type << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Superlattice Variation           = " << str_aus.bravais_superlattice_variation_type << endl;
//DX20210301 [OBSOLETE]           //DX oss << " Real space: Superlattice System              = " << str_aus.bravais_lattice_system << endl;
//DX20210301 [OBSOLETE]           oss << " Real space: Superlattice System              = " << str_aus.bravais_superlattice_system << endl; //DX - fixed mistake in line above
//DX20210301 [OBSOLETE]           oss << " Real space: Pearson Symbol Superlattice      = " << str_aus.pearson_symbol_superlattice << endl;
//DX20210301 [OBSOLETE]         }
//DX20210301 [OBSOLETE]         oss << "RECIPROCAL LATTICE" << endl;
//DX20210301 [OBSOLETE]         oss << " Reciprocal space lattice:" << endl;
//DX20210301 [OBSOLETE]         oss.precision(10);
//DX20210301 [OBSOLETE]         oss <<"   "<<str.klattice(1,1)<<"  "<<str.klattice(1,2) <<" "<<str.klattice(1,3)<<endl;
//DX20210301 [OBSOLETE]         oss <<"   "<<str.klattice(2,1)<<"  "<<str.klattice(2,2) <<" "<<str.klattice(2,3)<<endl;
//DX20210301 [OBSOLETE]         oss <<"   "<<str.klattice(3,1)<<"  "<<str.klattice(3,2) <<" "<<str.klattice(3,3)<<endl;
//DX20210301 [OBSOLETE]         oss << " Reciprocal space a b c alpha beta gamma: ";
//DX20210301 [OBSOLETE]         data=Getabc_angles(str.klattice,DEGREES);
//DX20210301 [OBSOLETE]         oss.precision(10);oss << data(1) << " " << data(2) << " " << data(3) << " ";
//DX20210301 [OBSOLETE]         oss.precision(3); oss << data(4) << " " << data(5) << " " << data(6) << endl;
//DX20210301 [OBSOLETE]         oss.precision(4);
//DX20210301 [OBSOLETE]         double kvol=GetVol(str.klattice);
//DX20210301 [OBSOLETE]         oss << " Reciprocal space Volume: " << kvol << endl;
//DX20210301 [OBSOLETE]         if(smode=="EDATA") {
//DX20210301 [OBSOLETE]           oss << " Reciprocal lattice primitive            = " << str_aus.reciprocal_lattice_type << endl;
//DX20210301 [OBSOLETE]           oss << " Reciprocal lattice variation            = " << str_aus.reciprocal_lattice_variation_type << endl; //WSETYAWAN mod
//DX20210301 [OBSOLETE]           //oss << " Reciprocal conventional lattice         = " << str_aus.reciprocal_conventional_lattice_type << endl;
//DX20210301 [OBSOLETE]           oss << "SPRIM" << endl;
//DX20210301 [OBSOLETE]           oss << str_sp << endl;
//DX20210301 [OBSOLETE]           oss << "SCONV" << endl;
//DX20210301 [OBSOLETE]           oss << str_sc << endl;
//DX20210301 [OBSOLETE]         }
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       // FORMAT = JSON
//DX20210301 [OBSOLETE]       if(format=="json"){
//DX20210301 [OBSOLETE]         string eendl="";
//DX20210301 [OBSOLETE]         bool roff=true; //round off
//DX20210301 [OBSOLETE]         stringstream sscontent_json;
//DX20210301 [OBSOLETE]         vector<string> vcontent_json;
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]         // Real lattice parameters
//DX20210301 [OBSOLETE]         xvector<double> data(6);
//DX20210301 [OBSOLETE]         data=Getabc_angles(str.lattice,DEGREES);data(1)*=str.scale;data(2)*=str.scale;data(3)*=str.scale;
//DX20210301 [OBSOLETE]         if(data.rows){
//DX20210301 [OBSOLETE]           sscontent_json << "\"lattice_parameters\":[" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(data,_AFLOWLIB_DATA_GEOMETRY_PREC_,roff),",") << "]" << eendl;  //CO20200731 - precision
//DX20210301 [OBSOLETE]         } else {
//DX20210301 [OBSOLETE]           if(PRINT_NULL_JSON){ sscontent_json << "\"lattice_parameters\":null" << eendl;}
//DX20210301 [OBSOLETE]         }
//DX20210301 [OBSOLETE]         vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]         // Real lattice parameters (Bohr/Deg)
//DX20210301 [OBSOLETE]         if(smode=="EDATA") {
//DX20210301 [OBSOLETE]           xvector<double> Bohr_Degs_data(6);
//DX20210301 [OBSOLETE]           Bohr_Degs_data(1) =  data(1)*angstrom2bohr;
//DX20210301 [OBSOLETE]           Bohr_Degs_data(2) =  data(2)*angstrom2bohr;
//DX20210301 [OBSOLETE]           Bohr_Degs_data(3) =  data(3)*angstrom2bohr;
//DX20210301 [OBSOLETE]           Bohr_Degs_data(4) =  data(4);
//DX20210301 [OBSOLETE]           Bohr_Degs_data(5) =  data(5);
//DX20210301 [OBSOLETE]           Bohr_Degs_data(6) =  data(6);
//DX20210301 [OBSOLETE]           if(data.rows){
//DX20210301 [OBSOLETE]             sscontent_json << "\"lattice_parameters_Bohr_deg\":[" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(Bohr_Degs_data,_AFLOWLIB_DATA_GEOMETRY_PREC_,roff),",") << "]" << eendl; //CO20200731 - precision
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             if(PRINT_NULL_JSON){ sscontent_json << "\"lattice_parameters_Bohr_deg\":null" << eendl;}
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]         }
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]         // Real space volume
//DX20210301 [OBSOLETE]         sscontent_json << "\"volume\":" << vol << eendl;
//DX20210301 [OBSOLETE]         vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]         // Real space c/a
//DX20210301 [OBSOLETE]         sscontent_json << "\"c_over_a\":" << data(3)/data(1) << eendl;
//DX20210301 [OBSOLETE]         vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]         str_aus=str; //,str_sp,str_sc; //CO20171027
//DX20210301 [OBSOLETE]         str_aus.sym_eps = str_sp.sym_eps = str_sc.sym_eps = tolerance;
//DX20210301 [OBSOLETE]         if(smode=="EDATA") {
//DX20210301 [OBSOLETE]           if(!already_calculated){  //CO20171025
//DX20210301 [OBSOLETE]             str_aus.GetLatticeType(str_sp,str_sc);
//DX20210301 [OBSOLETE]             if(orig_tolerance == AUROSTD_NAN){
//DX20210301 [OBSOLETE]               orig_tolerance = str_sp.sym_eps;  // Use new tolerance calculated here
//DX20210301 [OBSOLETE]             }
//DX20210301 [OBSOLETE]             tolerance = str_sp.sym_eps;
//DX20210301 [OBSOLETE]             sym_eps_change_count = str_aus.sym_eps_change_count = str_sp.sym_eps_change_count; //DX20180226 - added sym eps change count
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: bravais lattice primitive
//DX20210301 [OBSOLETE]           if(str_aus.bravais_lattice_type.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_lattice_type\":\"" << str_aus.bravais_lattice_type << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_lattice_type\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: bravais lattice variation
//DX20210301 [OBSOLETE]           if(str_aus.bravais_lattice_variation_type.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_lattice_variation_type\":\"" << str_aus.bravais_lattice_variation_type << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_lattice_variation_type\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: lattice system
//DX20210301 [OBSOLETE]           if(str_aus.bravais_lattice_system.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_lattice_system\":\"" << str_aus.bravais_lattice_system << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_lattice_system\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: Pearson symbol
//DX20210301 [OBSOLETE]           if(str_aus.pearson_symbol.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"Pearson_symbol\":\"" << str_aus.pearson_symbol << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"Pearson_symbol\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: crystal family
//DX20210301 [OBSOLETE]           if(str_aus.crystal_family.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"crystal_family\":\"" << str_aus.crystal_family << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"crystal_family\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: crystal system
//DX20210301 [OBSOLETE]           if(str_aus.crystal_system.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"crystal_system\":\"" << str_aus.crystal_system << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"crystal_system\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: point group crystal class
//DX20210301 [OBSOLETE]           if(str_aus.point_group_crystal_class.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"point_group_crystal_class\":\"" << str_aus.point_group_crystal_class << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"point_group_crystal_class\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: point group Hermann Mauguin
//DX20210301 [OBSOLETE]           if(str_aus.point_group_Hermann_Mauguin.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"point_group_Hermann_Mauguin\":\"" << str_aus.point_group_Hermann_Mauguin << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"point_group_Hermann_Mauguin\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: point group Schoenflies
//DX20210301 [OBSOLETE]           if(str_aus.point_group_Shoenflies.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"point_group_Schoenflies\":\"" << str_aus.point_group_Shoenflies << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"point_group_Schoenflies\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: point group orbifold
//DX20210301 [OBSOLETE]           if(str_aus.point_group_orbifold.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"point_group_orbifold\":\"" << str_aus.point_group_orbifold << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"point_group_orbifold\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: point group type
//DX20210301 [OBSOLETE]           if(str_aus.point_group_type.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"point_group_type\":\"" << str_aus.point_group_type << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"point_group_type\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: point group order
//DX20210301 [OBSOLETE]           if(str_aus.point_group_order.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"point_group_order\":\"" << str_aus.point_group_order << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"point_group_order\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: point group structure
//DX20210301 [OBSOLETE]           if(str_aus.point_group_structure.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"point_group_structure\":\"" << str_aus.point_group_structure << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"point_group_structure\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // SGDATA
//DX20210301 [OBSOLETE]           PrintSGData(str_aus,tolerance,sscontent_json,vpflow,no_scan,sg_setting,false,format,already_calculated);  //DX20170831 - SGDATA //CO20171027 //DX20180823 - added xoption
//DX20210301 [OBSOLETE]           //DX20180608 - recalculate GetLatticeType() before changing tolerance - START
//DX20210301 [OBSOLETE]           if(sym_eps_change_count != str_aus.sym_eps_change_count){
//DX20210301 [OBSOLETE]             if(LDEBUG) {
//DX20210301 [OBSOLETE]               cerr << XPID << "pflow::PrintData: WARNING: PrintSGData() changed the tolerance, need to recalculate GetLatticeType() at the new tolerance (i.e., sym_eps=" << str_aus.sym_eps << ") [dir=" << str.directory << "]" << endl; //DX20180526 - add directory
//DX20210301 [OBSOLETE]             }
//DX20210301 [OBSOLETE]             sym_eps_change_count = str_aus.sym_eps_change_count;
//DX20210301 [OBSOLETE]             str_aus.ClearSymmetry();
//DX20210301 [OBSOLETE]             sscontent_json.str("");
//DX20210301 [OBSOLETE]             vcontent_json.clear();
//DX20210301 [OBSOLETE]             continue;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           //DX20180608 - recalculate GetLatticeType() before changing tolerance - END
//DX20210301 [OBSOLETE]           sym_eps_change_count = str_aus.sym_eps_change_count; //DX20180226 - added sym eps change count
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]           //DX20170901 - Add consistency check for input symmetry method and ITC method - START
//DX20210301 [OBSOLETE]           int multiplicity_of_primitive =str_sp.fgroup.size()/str_sp.pgroup_xtal.size();
//DX20210301 [OBSOLETE]           bool derivative_structure=false;
//DX20210301 [OBSOLETE]           if(SYM::ComparePointGroupAndSpaceGroupString(str_aus,multiplicity_of_primitive,derivative_structure) || no_scan == true){
//DX20210301 [OBSOLETE]             symmetry_commensurate=true;
//DX20210301 [OBSOLETE]           } else { //CO20171027
//DX20210301 [OBSOLETE]             if(LDEBUG) {
//DX20210301 [OBSOLETE]               cerr << XPID << "pflow::PrintData: WARNING: Space group symbol and point group symbol do not match." << endl;
//DX20210301 [OBSOLETE]             }
//DX20210301 [OBSOLETE]             str_aus.ClearSymmetry();
//DX20210301 [OBSOLETE]             sscontent_json.str("");
//DX20210301 [OBSOLETE]             vcontent_json.clear();
//DX20210301 [OBSOLETE]             //if(!SYM::change_tolerance(str_aus,tolerance,orig_tolerance,change_sym_count,str_aus.dist_nn_min,no_scan))
//DX20210301 [OBSOLETE]             if(!SYM::change_tolerance(str_aus,tolerance,str_aus.dist_nn_min,no_scan))
//DX20210301 [OBSOLETE]             {   //CO20200106 - patching for auto-indenting
//DX20210301 [OBSOLETE]               if(force_perform){
//DX20210301 [OBSOLETE]                 if(LDEBUG) {
//DX20210301 [OBSOLETE]                   cerr << XPID << "pflow::PrintData: WARNING: Scan failed. Reverting back to original tolerance and recalculating as is (with aforementioned inconsistencies)." << str.directory << endl;
//DX20210301 [OBSOLETE]                 }
//DX20210301 [OBSOLETE]                 no_scan=true; //Force it to continue
//DX20210301 [OBSOLETE]                 tolerance = orig_tolerance;
//DX20210301 [OBSOLETE]                 str_aus.ClearSymmetry();
//DX20210301 [OBSOLETE]                 sscontent_json.str("");
//DX20210301 [OBSOLETE]                 vcontent_json.clear();
//DX20210301 [OBSOLETE]                 continue;
//DX20210301 [OBSOLETE]               }
//DX20210301 [OBSOLETE]             }
//DX20210301 [OBSOLETE]             else {
//DX20210301 [OBSOLETE]               str_aus.ClearSymmetry();
//DX20210301 [OBSOLETE]               sscontent_json.str("");
//DX20210301 [OBSOLETE]               vcontent_json.clear();
//DX20210301 [OBSOLETE]               sym_eps_change_count = str_aus.sym_eps_change_count; //DX20180226 - added sym eps change count
//DX20210301 [OBSOLETE]               continue;
//DX20210301 [OBSOLETE]             }
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           //DX20170901 - Add consistency check for input symmetry method and ITC method - END
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: bravais lattice lattice type
//DX20210301 [OBSOLETE]           if(str_aus.bravais_lattice_lattice_type.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_lattice_lattice_type\":\"" << str_aus.bravais_lattice_lattice_type << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_lattice_lattice_type\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: bravais lattice lattice variation type
//DX20210301 [OBSOLETE]           if(str_aus.bravais_lattice_lattice_variation_type.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_lattice_lattice_variation_type\":\"" << str_aus.bravais_lattice_lattice_variation_type << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_lattice_lattice_variation_type\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: bravais lattice lattice system
//DX20210301 [OBSOLETE]           if(str_aus.bravais_lattice_lattice_system.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_lattice_lattice_system\":\"" << str_aus.bravais_lattice_lattice_system << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_lattice_lattice_system\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: bravais superlattice type
//DX20210301 [OBSOLETE]           if(str_aus.bravais_superlattice_type.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_superlattice_type\":\"" << str_aus.bravais_superlattice_type << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_superlattice_type\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: bravais superlattice variation type
//DX20210301 [OBSOLETE]           if(str_aus.bravais_superlattice_variation_type.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_superlattice_variation_type\":\"" << str_aus.bravais_superlattice_variation_type << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_superlattice_variation_type\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: bravais superlattice system
//DX20210301 [OBSOLETE]           if(str_aus.bravais_superlattice_system.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_superlattice_system\":\"" << str_aus.bravais_superlattice_system << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"Bravais_superlattice_system\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Real space: pearson symbol superlattice
//DX20210301 [OBSOLETE]           if(str_aus.pearson_symbol_superlattice.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"Pearson_symbol_superlattice\":\"" << str_aus.pearson_symbol_superlattice << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"Pearson_symbol_superlattice\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]         }
//DX20210301 [OBSOLETE]         //RECIP
//DX20210301 [OBSOLETE]         // Reciprocal space lattice
//DX20210301 [OBSOLETE]         if(str_aus.klattice.rows){
//DX20210301 [OBSOLETE]           sscontent_json << "\"reciprocal_lattice_vectors\":[" << aurostd::xmatDouble2String(str_aus.klattice,_AFLOWLIB_DATA_GEOMETRY_PREC_,roff) << "]" << eendl;  //CO20200731 - precision
//DX20210301 [OBSOLETE]         } else {
//DX20210301 [OBSOLETE]           sscontent_json << "\"reciprocal_lattice_vectors\":null" << eendl;
//DX20210301 [OBSOLETE]         }
//DX20210301 [OBSOLETE]         vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]         // Reciprocal lattice parameters
//DX20210301 [OBSOLETE]         data=Getabc_angles(str.klattice,DEGREES);
//DX20210301 [OBSOLETE]         if(data.rows){
//DX20210301 [OBSOLETE]           sscontent_json << "\"reciprocal_lattice_parameters\":[" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(data,_AFLOWLIB_DATA_GEOMETRY_PREC_,roff),",") << "]" << eendl; //CO20200731 - precision
//DX20210301 [OBSOLETE]         } else {
//DX20210301 [OBSOLETE]           sscontent_json << "\"reciprocal_lattice_parameters\":null" << eendl;
//DX20210301 [OBSOLETE]         }
//DX20210301 [OBSOLETE]         vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]         // Reciprocal space volume
//DX20210301 [OBSOLETE]         double kvol=GetVol(str.klattice);
//DX20210301 [OBSOLETE]         sscontent_json << "\"reciprocal_volume\":" << kvol << eendl;
//DX20210301 [OBSOLETE]         vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]         if(smode=="EDATA") {
//DX20210301 [OBSOLETE]           // Reciprocal space: reciprocal lattice type
//DX20210301 [OBSOLETE]           if(str_aus.reciprocal_lattice_type.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"reciprocal_lattice_type\":\"" << str_aus.reciprocal_lattice_type << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"reciprocal_lattice_type\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // Reciprocal space: reciprocal lattice variation type
//DX20210301 [OBSOLETE]           if(str_aus.reciprocal_lattice_variation_type.size()){
//DX20210301 [OBSOLETE]             sscontent_json << "\"reciprocal_lattice_variation_type\":\"" << str_aus.reciprocal_lattice_variation_type << "\"" << eendl;
//DX20210301 [OBSOLETE]           } else {
//DX20210301 [OBSOLETE]             sscontent_json << "\"reciprocal_lattice_variation_type\":null" << eendl;
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // SPRIM
//DX20210301 [OBSOLETE]           sscontent_json << "\"standard_primitive_structure\":" << xstructure2json(str_sp) << eendl;
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]           // SCONV
//DX20210301 [OBSOLETE]           sscontent_json << "\"standard_conventional_structure\":" << xstructure2json(str_sc) << eendl;
//DX20210301 [OBSOLETE]           vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]         }
//DX20210301 [OBSOLETE]         oss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}" << endl;
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]     }
//DX20210301 [OBSOLETE]     oss_final << oss.str();
//DX20210301 [OBSOLETE]     str_sym=str_aus; //CO20171027
//DX20210301 [OBSOLETE]     //DX20180822 - put attributes into vpflow - START
//DX20210301 [OBSOLETE]     if(vpflow.flag("EDATA::CALCULATED")){ //DX20180823 - if calculated remove
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::LATTICE_PARAMETERS");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::COVERA");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::VOLUME");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::BRAVAIS_LATTICE_TYPE");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::BRAVAIS_LATTICE_VARIATION_TYPE");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::BRAVAIS_LATTICE_SYSTEM");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::PEARSON_SYMBOL");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::CRYSTAL_FAMILY");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::CRYSTAL_SYSTEM");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::POINT_GROUP_CRYSTAL_CLASS");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::POINT_GROUP_HERMANN_MAUGUIN");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::POINT_GROUP_SHOENFLIES");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::POINT_GROUP_ORBIFOLD");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::POINT_GROUP_TYPE");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::POINT_GROUP_ORDER");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::POINT_GROUP_STRUCTURE");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::BRAVAIS_LATTICE_LATTICE_TYPE");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::BRAVAIS_LATTICE_LATTICE_SYSTEM");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::BRAVAIS_SUPERLATTICE_TYPE");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::BRAVAIS_SUPERLATTICE_VARIATION_TYPE");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::BRAVAIS_SUPERLATTICE_SYSTEM");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::PEARSON_SYMBOL_SUPERLATTICE");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::RECIPROCAL_LATTICE_PARAMETERS");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::REAL_SPACE_VOLUME");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::RECIPROCAL_LATTICE_TYPE");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::RECIPROCAL_LATTICE_VARIATION_TYPE");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::SPRIM");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("EDATA::SCONV");
//DX20210301 [OBSOLETE]     }
//DX20210301 [OBSOLETE]     //DX20180823 - now populate xoption
//DX20210301 [OBSOLETE]     vpflow.flag("EDATA::CALCULATED",TRUE);
//DX20210301 [OBSOLETE]     xvector<double> data(6);
//DX20210301 [OBSOLETE]     data=Getabc_angles(str.lattice,DEGREES);data(1)*=str.scale;data(2)*=str.scale;data(3)*=str.scale;
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::LATTICE_PARAMETERS",aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(data,_AFLOWLIB_DATA_GEOMETRY_PREC_,true),","));  //CO20200731 - precision
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::COVERA",aurostd::utype2string<double>(data(3)/data(1)));
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::VOLUME",aurostd::utype2string<double>(vol));
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::BRAVAIS_LATTICE_TYPE",str_aus.bravais_lattice_type);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::BRAVAIS_LATTICE_VARIATION_TYPE",str_aus.bravais_lattice_variation_type);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::BRAVAIS_LATTICE_SYSTEM",str_aus.bravais_lattice_system);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::PEARSON_SYMBOL",str_aus.pearson_symbol);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::CRYSTAL_FAMILY",str_aus.crystal_family);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::CRYSTAL_SYSTEM",str_aus.crystal_system);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::POINT_GROUP_CRYSTAL_CLASS",str_aus.point_group_crystal_class);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::POINT_GROUP_HERMANN_MAUGUIN",str_aus.point_group_Hermann_Mauguin);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::POINT_GROUP_SCHOENFLIES",str_aus.point_group_Shoenflies);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::POINT_GROUP_ORBIFOLD",str_aus.point_group_orbifold);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::POINT_GROUP_TYPE",aurostd::StringSubst(str_aus.point_group_type, " ", "_"));
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::POINT_GROUP_ORDER",str_aus.point_group_order);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::POINT_GROUP_STRUCTURE",aurostd::StringSubst(str_aus.point_group_structure, " ", "_"));
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::BRAVAIS_LATTICE_LATTICE_TYPE",str_aus.bravais_lattice_lattice_type);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE",str_aus.bravais_lattice_lattice_variation_type);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::BRAVAIS_LATTICE_LATTICE_SYSTEM",str_aus.bravais_lattice_lattice_system);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::BRAVAIS_SUPERLATTICE_TYPE",str_aus.bravais_superlattice_type);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::BRAVAIS_SUPERLATTICE_VARIATION_TYPE",str_aus.bravais_superlattice_variation_type);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::BRAVAIS_SUPERLATTICE_SYSTEM",str_aus.bravais_superlattice_system);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::PEARSON_SYMBOL_SUPERLATTICE",str_aus.pearson_symbol_superlattice);
//DX20210301 [OBSOLETE]     data=Getabc_angles(str.klattice,DEGREES);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::RECIPROCAL_LATTICE_PARAMETERS",aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(data,_AFLOWLIB_DATA_GEOMETRY_PREC_,true),","));   //CO20200731 - precision
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::RECIPROCAL_SPACE_VOLUME",aurostd::utype2string<double>(GetVol(str.klattice)));
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::RECIPROCAL_LATTICE_TYPE",str_aus.reciprocal_lattice_type);
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::RECIPROCAL_LATTICE_VARIATION_TYPE",str_aus.reciprocal_lattice_variation_type);
//DX20210301 [OBSOLETE]     stringstream ss_str_sp; ss_str_sp << str_sp << endl;
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::SPRIM",ss_str_sp.str());
//DX20210301 [OBSOLETE]     stringstream ss_str_sc; ss_str_sc << str_sc << endl;
//DX20210301 [OBSOLETE]     vpflow.push_attached("EDATA::SCONV",ss_str_sc.str());
//DX20210301 [OBSOLETE]     //DX20180822 - put attributes into vpflow - END
//DX20210301 [OBSOLETE]   }
//DX20210301 [OBSOLETE] } // namespace pflow
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE] namespace pflow {
//DX20210301 [OBSOLETE]   void PrintData(const xstructure& str,ostream& oss,const string& smode,const string& format) {
//DX20210301 [OBSOLETE]     xstructure str_sym,str_sp,str_sc; //CO20171027
//DX20210301 [OBSOLETE]     double tolerance = AUROSTD_NAN;
//DX20210301 [OBSOLETE]     bool no_scan = false;
//DX20210301 [OBSOLETE]     int sg_setting = 1;
//DX20210301 [OBSOLETE]     return PrintData(str,str_sym,str_sp,str_sc,oss,tolerance,smode,no_scan,sg_setting,format); //CO20171027
//DX20210301 [OBSOLETE]   }
//DX20210301 [OBSOLETE]   void PrintData(const xstructure& str,ostream& oss,aurostd::xoption& vpflow,const string& smode,const string& format) {  //CO20200731
//DX20210301 [OBSOLETE]     xstructure str_sym,str_sp,str_sc; //CO20171027
//DX20210301 [OBSOLETE]     double tolerance = AUROSTD_NAN;
//DX20210301 [OBSOLETE]     bool no_scan = false;
//DX20210301 [OBSOLETE]     int sg_setting = 1;
//DX20210301 [OBSOLETE]     return PrintData(str,str_sym,str_sp,str_sc,oss,vpflow,tolerance,smode,no_scan,sg_setting,format); //CO20171027
//DX20210301 [OBSOLETE]   }
//DX20210301 [OBSOLETE]   void PrintData(const xstructure& str,xstructure& str_sp,xstructure& str_sc,ostream& oss,aurostd::xoption& vpflow,const string& smode,const string& format) {  //CO20200731
//DX20210301 [OBSOLETE]     xstructure str_sym; //CO20171027
//DX20210301 [OBSOLETE]     double tolerance = AUROSTD_NAN;
//DX20210301 [OBSOLETE]     bool no_scan = false;
//DX20210301 [OBSOLETE]     int sg_setting = 1;
//DX20210301 [OBSOLETE]     return PrintData(str,str_sym,str_sp,str_sc,oss,vpflow,tolerance,smode,no_scan,sg_setting,format); //CO20171027
//DX20210301 [OBSOLETE]   }
//DX20210301 [OBSOLETE] } // namespace pflow
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE] namespace pflow {
//DX20210301 [OBSOLETE]   void PrintData(const xstructure& str,ostream& oss,double tolerance,const string& smode, bool no_scan, int sg_setting, const string& format) {
//DX20210301 [OBSOLETE]     xstructure str_sym,str_sp,str_sc; //CO20171027
//DX20210301 [OBSOLETE]     return PrintData(str,str_sym,str_sp,str_sc,oss,tolerance,smode,no_scan,sg_setting,format); //CO20171027
//DX20210301 [OBSOLETE]   }
//DX20210301 [OBSOLETE] } // namespace pflow

// ***************************************************************************
// pflow::PrintRealLatticeData() //DX20210209
// ***************************************************************************
// Print the real lattice data (lattice, type/variation, Pearson, etc.)
// Also updates the xoption (EDATA::<PROPERTY>) for lib2raw runs
namespace pflow {
  string PrintRealLatticeData(const xstructure& xstr, const string& smode, filetype ftype, bool standalone, bool already_calculated, double sym_eps){
    aurostd::xoption vpflow;
    return PrintRealLatticeData(xstr, vpflow, smode, ftype, standalone, already_calculated, sym_eps);
  }
}
namespace pflow {
  string PrintRealLatticeData(const xstructure& xstr, aurostd::xoption& vpflow, const string& smode, filetype ftype, bool standalone, bool already_calculated, double sym_eps){


    // ---------------------------------------------------------------------------
    // calculate the real lattice symmetry if not already calculated
    xstructure str_aus(xstr);
    str_aus.ReScale(1.0);
    if(smode == "EDATA" && !already_calculated){
      str_aus.GetRealLatticeType(sym_eps);
    }

    // get lattice and angles
    xvector<double> data(6);
    data=Getabc_angles(xstr.lattice,DEGREES);data(1)*=xstr.scale;data(2)*=xstr.scale;data(3)*=xstr.scale;
    double c_over_a = data(3)/data(1);

    // get lattice and angles (Bohr)
    xvector<double> data_Bohr(6);
    data_Bohr(1) =  data(1)*angstrom2bohr;
    data_Bohr(2) =  data(2)*angstrom2bohr;
    data_Bohr(3) =  data(3)*angstrom2bohr;
    data_Bohr(4) =  data(4);
    data_Bohr(5) =  data(5);
    data_Bohr(6) =  data(6);

    // get volume
    double vol = GetVol(xstr.lattice)*xstr.scale*xstr.scale*xstr.scale;

    stringstream ss_output;
    // ---------------------------------------------------------------------------
    // text output
    if(ftype == txt_ft){

      string info_prefix = "";
      if(!standalone){ info_prefix = "Real space "; }

      ss_output << "REAL LATTICE" << endl;
      ss_output << " Real space lattice:" << endl;
      ss_output << "  " << aurostd::roundoff(str_aus.lattice(1),AUROSTD_ROUNDOFF_TOL) << endl;
      ss_output << "  " << aurostd::roundoff(str_aus.lattice(2),AUROSTD_ROUNDOFF_TOL) << endl;
      ss_output << "  " << aurostd::roundoff(str_aus.lattice(3),AUROSTD_ROUNDOFF_TOL) << endl;
      ss_output << " " << info_prefix << "a b c alpha beta gamma: ";
      ss_output.precision(10);ss_output << data(1) << " " << data(2) << " " << data(3) << " ";
      ss_output.precision(3); ss_output << data(4) << " " << data(5) << " " << data(6) << endl;
      ss_output.precision(4);
      if(smode == "EDATA"){
        ss_output << " " << info_prefix << "a b c alpha beta gamma: ";
        ss_output.precision(10);ss_output << data_Bohr(1) << " " << data_Bohr(2) << " " << data_Bohr(3) << " ";
        ss_output.precision(3); ss_output << data(4) << " " << data(5) << " " << data(6) << "   Bohrs/Degs " << endl;
        ss_output.precision(4);
      }
      ss_output << " " << info_prefix << "Volume: " << vol << endl;
      ss_output << " " << info_prefix << "c/a = " << c_over_a << endl;
      if(smode == "EDATA"){
        ss_output << "BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal)" << endl;
        ss_output << " " << info_prefix << "Bravais Lattice Primitive        = " << str_aus.bravais_lattice_type << endl;// " " << str.title << endl;
        ss_output << " " << info_prefix << "Lattice Variation                = " << str_aus.bravais_lattice_variation_type << endl;//WSETYAWAN mod
        ss_output << " " << info_prefix << "Lattice System                   = " << str_aus.bravais_lattice_system << endl;
        ss_output << " " << info_prefix << "Pearson Symbol                   = " << str_aus.pearson_symbol << endl;
      }
    }
    // ---------------------------------------------------------------------------
    // json output
    else if(ftype == json_ft){

      aurostd::JSONwriter json;
      bool roff = true;

      // Real space lattice
      if(str_aus.lattice.rows != 0){
        json.addMatrix("lattice_vectors",str_aus.lattice,_AFLOWLIB_DATA_GEOMETRY_PREC_,roff);
      } else if(PRINT_NULL_JSON) {
        json.addNull("lattice_vectors");
      }

      // Real lattice parameters
      if(data.rows != 0){
        json.addVector("lattice_parameters",data,_AFLOWLIB_DATA_GEOMETRY_PREC_,roff);
      } else if (PRINT_NULL_JSON){
        json.addNull("lattice_parameters");
      }

      if(smode == "EDATA"){
        // Real lattice parameters (Bohr/Deg)
        if(data_Bohr.rows != 0){
          json.addVector("lattice_parameters_Bohr_deg",data_Bohr,_AFLOWLIB_DATA_GEOMETRY_PREC_,roff);
        } else if (PRINT_NULL_JSON){
          json.addNull("lattice_parameters_Bohr_deg");
        }
      }

      // Real space volume
      json.addNumber("volume", vol);

      // Real space c/a
      json.addNumber("c_over_a", c_over_a);

      if(smode == "EDATA"){
        // Real space: bravais lattice primitive
        if(!str_aus.bravais_lattice_type.empty()){
          json.addString("Bravais_lattice_type", str_aus.bravais_lattice_type);
        } else if (PRINT_NULL_JSON){
          json.addNull("Bravais_lattice_type");
        }

        // Real space: bravais lattice variation
        if(!str_aus.bravais_lattice_variation_type.empty()){
          json.addString("Bravais_lattice_variation_type", str_aus.bravais_lattice_variation_type);
        } else if (PRINT_NULL_JSON){
          json.addNull("Bravais_lattice_variation_type");
        }

        // Real space: lattice system
        if(!str_aus.bravais_lattice_system.empty()){
          json.addString("Bravais_lattice_system", str_aus.bravais_lattice_system);
        } else if (PRINT_NULL_JSON){
          json.addNull("Bravais_lattice_system");
        }

        // Real space: Pearson symbol
        if(!str_aus.pearson_symbol.empty()){
          json.addString("Pearson_symbol", str_aus.pearson_symbol);
        } else if (PRINT_NULL_JSON){
          json.addNull("Pearson_symbol");
        }
      }

      ss_output << json.toString(standalone); //standalone: determines if we enclose in brackets
      if(standalone) { ss_output << endl; }
    }
    else{
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Format type is not supported.",_INPUT_ILLEGAL_);
    }

    // ---------------------------------------------------------------------------
    // update vpflow, needed for library runs
    // clear first, if needed
    if(smode == "EDATA"){
      if(vpflow.flag("EDATA::REAL_LATTICE::CALCULATED")){
        vpflow.pop_attached("EDATA::LATTICE_PARAMETERS");
        vpflow.pop_attached("EDATA::COVERA");
        vpflow.pop_attached("EDATA::VOLUME");
        vpflow.pop_attached("EDATA::BRAVAIS_LATTICE_TYPE");
        vpflow.pop_attached("EDATA::BRAVAIS_LATTICE_VARIATION_TYPE");
        vpflow.pop_attached("EDATA::BRAVAIS_LATTICE_SYSTEM");
        vpflow.pop_attached("EDATA::PEARSON_SYMBOL");
      }
      // update vpflow
      vpflow.flag("EDATA::REAL_LATTICE::CALCULATED",TRUE);
      vpflow.push_attached("EDATA::LATTICE_PARAMETERS",aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(data,_AFLOWLIB_DATA_GEOMETRY_PREC_,true),","));  //CO20200731 - precision
      vpflow.push_attached("EDATA::COVERA",aurostd::utype2string<double>(c_over_a));
      vpflow.push_attached("EDATA::VOLUME",aurostd::utype2string<double>(vol));
      vpflow.push_attached("EDATA::BRAVAIS_LATTICE_TYPE",str_aus.bravais_lattice_type);
      vpflow.push_attached("EDATA::BRAVAIS_LATTICE_VARIATION_TYPE",str_aus.bravais_lattice_variation_type);
      vpflow.push_attached("EDATA::BRAVAIS_LATTICE_SYSTEM",str_aus.bravais_lattice_system);
      vpflow.push_attached("EDATA::PEARSON_SYMBOL",str_aus.pearson_symbol);
    }

    return ss_output.str();
  }
}

// ***************************************************************************
// pflow::PrintLatticeLatticeData() //DX20210209
// ***************************************************************************
// Print the lattice lattice data (lattice type/variation/system)
// Also updates the xoption (EDATA::<PROPERTY>) for lib2raw runs
namespace pflow {
  string PrintLatticeLatticeData(const xstructure& xstr, filetype ftype, bool standalone, bool already_calculated, double sym_eps){
    aurostd::xoption vpflow;
    return PrintLatticeLatticeData(xstr, vpflow, ftype, standalone, already_calculated, sym_eps);
  }
}
namespace pflow {
  string PrintLatticeLatticeData(const xstructure& xstr, aurostd::xoption& vpflow, filetype ftype, bool standalone, bool already_calculated, double sym_eps){


    // ---------------------------------------------------------------------------
    // calculate the real lattice symmetry if not already calculated
    xstructure str_aus(xstr);
    str_aus.ReScale(1.0);
    if(!already_calculated){
      str_aus.GetRealLatticeType(sym_eps);
    }

    stringstream ss_output;
    // ---------------------------------------------------------------------------
    // text output
    if(ftype == txt_ft){

      string info_prefix = "";
      if(!standalone){ info_prefix = "Real space "; }

      ss_output << "BRAVAIS LATTICE OF THE LATTICE (pgroup)" << endl;
      ss_output << " " << info_prefix << "Bravais Lattice Primitive        = " << str_aus.bravais_lattice_lattice_type << endl;// " " << str.title << endl;
      ss_output << " " << info_prefix << "Lattice Variation                = " << str_aus.bravais_lattice_lattice_variation_type << endl; //WSETYAWAN mod
      ss_output << " " << info_prefix << "Lattice System                   = " << str_aus.bravais_lattice_lattice_system << endl;
    }
    // ---------------------------------------------------------------------------
    // json output
    else if(ftype == json_ft){

      aurostd::JSONwriter json;

      // Real space: bravais lattice lattice type
      if(!str_aus.bravais_lattice_lattice_type.empty()){
        json.addString("Bravais_lattice_lattice_type", str_aus.bravais_lattice_lattice_type);
      } else if (PRINT_NULL_JSON){
        json.addNull("Bravais_lattice_lattice_type");
      }

      // Real space: bravais lattice lattice variation type
      if(!str_aus.bravais_lattice_lattice_variation_type.empty()){
        json.addString("Bravais_lattice_lattice_variation_type", str_aus.bravais_lattice_lattice_variation_type);
      } else if (PRINT_NULL_JSON){
        json.addNull("Bravais_lattice_lattice_variation_type");
      }

      // Real space: bravais lattice lattice system
      if(!str_aus.bravais_lattice_lattice_system.empty()){
        json.addString("Bravais_lattice_lattice_system", str_aus.bravais_lattice_lattice_system);
      } else if (PRINT_NULL_JSON){
        json.addNull("Bravais_lattice_lattice_system");
      }

      ss_output << json.toString(standalone); //standalone: determines if we enclose in brackets
      if(standalone) { ss_output << endl; }

    }
    else{
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Format type is not supported.",_INPUT_ILLEGAL_);
    }

    // ---------------------------------------------------------------------------
    // update vpflow, needed for library runs
    // clear first, if needed
    if(vpflow.flag("EDATA::LATTICE_LATTICE::CALCULATED")){
      vpflow.pop_attached("EDATA::BRAVAIS_LATTICE_LATTICE_TYPE");
      vpflow.pop_attached("EDATA::BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE");
      vpflow.pop_attached("EDATA::BRAVAIS_LATTICE_LATTICE_SYSTEM");
    }
    // update vpflow
    vpflow.flag("EDATA::LATTICE_LATTICE::CALCULATED",TRUE);
    vpflow.push_attached("EDATA::BRAVAIS_LATTICE_LATTICE_TYPE",str_aus.bravais_lattice_lattice_type);
    vpflow.push_attached("EDATA::BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE",str_aus.bravais_lattice_lattice_variation_type);
    vpflow.push_attached("EDATA::BRAVAIS_LATTICE_LATTICE_SYSTEM",str_aus.bravais_lattice_lattice_system);

    return ss_output.str();
  }
}

// ***************************************************************************
// pflow::PrintCrystalPointGroupData() //DX20210209
// ***************************************************************************
// Print the crystallographic point group data (family, class, type, etc.)
// Also updates the xoption (EDATA::<PROPERTY>) for lib2raw runs
namespace pflow {
  string PrintCrystalPointGroupData(const xstructure& xstr, filetype ftype, bool standalone, bool already_calculated, double sym_eps){
    aurostd::xoption vpflow;
    return PrintCrystalPointGroupData(xstr, vpflow, ftype, standalone, already_calculated, sym_eps);
  }
}
namespace pflow {
  string PrintCrystalPointGroupData(const xstructure& xstr, aurostd::xoption& vpflow, filetype ftype, bool standalone, bool already_calculated, double sym_eps){


    // ---------------------------------------------------------------------------
    // calculate the real lattice symmetry (contains point group information)
    // if not already calculated
    xstructure str_aus(xstr);
    str_aus.ReScale(1.0);
    if(!already_calculated){
      str_aus.GetRealLatticeType(sym_eps);
    }

    stringstream ss_output;
    // ---------------------------------------------------------------------------
    // text output
    if(ftype == txt_ft){

      string info_prefix = "";
      if(!standalone){ info_prefix = "Real space "; }

      ss_output << "POINT GROUP CRYSTAL" << endl;
      ss_output << " " << info_prefix << "Crystal Family                   = " << str_aus.crystal_family << endl;
      ss_output << " " << info_prefix << "Crystal System                   = " << str_aus.crystal_system << endl;
      ss_output << " " << info_prefix << "Crystal Class                    = " << str_aus.point_group_crystal_class << endl;
      ss_output << " " << info_prefix << "Point Group (Hermann Mauguin)    = " << str_aus.point_group_Hermann_Mauguin << endl; // "      PGXTAL" << endl;
      ss_output << " " << info_prefix << "Point Group (Schoenflies)        = " << str_aus.point_group_Shoenflies << endl;
      ss_output << " " << info_prefix << "Point Group Orbifold             = " << str_aus.point_group_orbifold << endl;
      ss_output << " " << info_prefix << "Point Group Type                 = " << str_aus.point_group_type << endl;
      ss_output << " " << info_prefix << "Point Group Order                = " << str_aus.point_group_order << endl;
      ss_output << " " << info_prefix << "Point Group Structure            = " << str_aus.point_group_structure << endl;
    }
    // ---------------------------------------------------------------------------
    // json output
    else if(ftype == json_ft){

      aurostd::JSONwriter json;

      // Real space: crystal family
      if(!str_aus.crystal_family.empty()){
        json.addString("crystal_family", str_aus.crystal_family);
      } else if (PRINT_NULL_JSON){
        json.addNull("crystal_family");
      }

      // Real space: crystal system
      if(!str_aus.crystal_system.empty()){
        json.addString("crystal_system", str_aus.crystal_system);
      } else if (PRINT_NULL_JSON){
        json.addNull("crystal_system");
      }

      // Real space: point group crystal class
      if(!str_aus.point_group_crystal_class.empty()){
        json.addString("point_group_crystal_class", str_aus.point_group_crystal_class);
      } else if (PRINT_NULL_JSON){
        json.addNull("point_group_crystal_class");
      }

      // Real space: point group Hermann Mauguin
      if(!str_aus.point_group_Hermann_Mauguin.empty()){
        json.addString("point_group_Hermann_Mauguin", str_aus.point_group_Hermann_Mauguin);
      } else if (PRINT_NULL_JSON){
        json.addNull("point_group_Hermann_Mauguin");
      }

      // Real space: point group Schoenflies
      if(!str_aus.point_group_Shoenflies.empty()){
        json.addString("point_group_Schoenflies", str_aus.point_group_Shoenflies);
      } else if (PRINT_NULL_JSON){
        json.addNull("point_group_Schoenflies");
      }

      // Real space: point group orbifold
      if(!str_aus.point_group_orbifold.empty()){
        json.addString("point_group_orbifold", str_aus.point_group_orbifold);
      } else if (PRINT_NULL_JSON){
        json.addNull("point_group_orbifold");
      }

      // Real space: point group type
      if(!str_aus.point_group_type.empty()){
        json.addString("point_group_type", str_aus.point_group_type);
      } else if (PRINT_NULL_JSON){
        json.addNull("point_group_type");
      }

      // Real space: point group order
      if(!str_aus.point_group_order.empty()){
        json.addNumber("point_group_order", str_aus.point_group_order);
      } else if (PRINT_NULL_JSON){
        json.addNull("point_group_order");
      }

      // Real space: point group structure
      if(!str_aus.point_group_structure.empty()){
        json.addString("point_group_structure", str_aus.point_group_structure);
      } else if (PRINT_NULL_JSON){
        json.addNull("point_group_structure");
      }

      ss_output << json.toString(standalone); //standalone: determines if we enclose in brackets
      if(standalone) { ss_output << endl; }

    }
    else{
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Format type is not supported.",_INPUT_ILLEGAL_);
    }

    // ---------------------------------------------------------------------------
    // update vpflow, needed for library runs
    // clear first, if needed
    if(vpflow.flag("EDATA::POINT_GROUP::CALCULATED")){
      vpflow.pop_attached("EDATA::CRYSTAL_FAMILY");
      vpflow.pop_attached("EDATA::CRYSTAL_SYSTEM");
      vpflow.pop_attached("EDATA::POINT_GROUP_CRYSTAL_CLASS");
      vpflow.pop_attached("EDATA::POINT_GROUP_HERMANN_MAUGUIN");
      vpflow.pop_attached("EDATA::POINT_GROUP_SHOENFLIES");
      vpflow.pop_attached("EDATA::POINT_GROUP_ORBIFOLD");
      vpflow.pop_attached("EDATA::POINT_GROUP_TYPE");
      vpflow.pop_attached("EDATA::POINT_GROUP_ORDER");
      vpflow.pop_attached("EDATA::POINT_GROUP_STRUCTURE");
    }
    // update vpflow
    vpflow.flag("EDATA::POINT_GROUP::CALCULATED",TRUE);
    vpflow.push_attached("EDATA::CRYSTAL_FAMILY",str_aus.crystal_family);
    vpflow.push_attached("EDATA::CRYSTAL_SYSTEM",str_aus.crystal_system);
    vpflow.push_attached("EDATA::POINT_GROUP_CRYSTAL_CLASS",str_aus.point_group_crystal_class);
    vpflow.push_attached("EDATA::POINT_GROUP_HERMANN_MAUGUIN",str_aus.point_group_Hermann_Mauguin);
    vpflow.push_attached("EDATA::POINT_GROUP_SCHOENFLIES",str_aus.point_group_Shoenflies);
    vpflow.push_attached("EDATA::POINT_GROUP_ORBIFOLD",str_aus.point_group_orbifold);
    vpflow.push_attached("EDATA::POINT_GROUP_TYPE",aurostd::StringSubst(str_aus.point_group_type, " ", "_"));
    vpflow.push_attached("EDATA::POINT_GROUP_ORDER",str_aus.point_group_order);
    vpflow.push_attached("EDATA::POINT_GROUP_STRUCTURE",aurostd::StringSubst(str_aus.point_group_structure, " ", "_"));

    return ss_output.str();
  }
}

// ***************************************************************************
// pflow::PrintReciprocalLatticeData() //DX20210209
// ***************************************************************************
// Print the reciprocal lattice data (lattice, type/variation, Pearson, etc.)
// Also updates the xoption (EDATA::<PROPERTY>) for lib2raw runs
namespace pflow {
  string PrintReciprocalLatticeData(const xstructure& xstr, filetype ftype, bool standalone, bool already_calculated, double sym_eps){
    aurostd::xoption vpflow;
    return PrintReciprocalLatticeData(xstr, vpflow, ftype, standalone, already_calculated, sym_eps);
  }
}
namespace pflow {
  string PrintReciprocalLatticeData(const xstructure& xstr, aurostd::xoption& vpflow, filetype ftype, bool standalone, bool already_calculated, double sym_eps){


    // ---------------------------------------------------------------------------
    // calculate the reciprocal lattice symmetry if not already calculated
    xstructure str_aus(xstr);
    str_aus.ReScale(1.0);
    if(!already_calculated){
      str_aus.GetReciprocalLatticeType(sym_eps);
    }

    xvector<double> data=Getabc_angles(str_aus.klattice,DEGREES);
    double kvol=GetVol(str_aus.klattice);

    stringstream ss_output;
    // ---------------------------------------------------------------------------
    // text output
    if(ftype == txt_ft){

      string info_prefix = "";
      if(!standalone){ info_prefix = "Reciprocal space "; }

      ss_output << "RECIPROCAL LATTICE" << endl;
      ss_output << " Reciprocal space lattice:" << endl;
      ss_output << "  " << aurostd::roundoff(str_aus.klattice(1),AUROSTD_ROUNDOFF_TOL) << endl;
      ss_output << "  " << aurostd::roundoff(str_aus.klattice(2),AUROSTD_ROUNDOFF_TOL) << endl;
      ss_output << "  " << aurostd::roundoff(str_aus.klattice(3),AUROSTD_ROUNDOFF_TOL) << endl;
      ss_output << " " << info_prefix << "a b c alpha beta gamma: ";
      ss_output.precision(10);ss_output << data(1) << " " << data(2) << " " << data(3) << " ";
      ss_output.precision(3); ss_output << data(4) << " " << data(5) << " " << data(6) << endl;
      ss_output.precision(4);
      ss_output << " " << info_prefix << "Volume: " << kvol << endl;
      ss_output << " Reciprocal lattice primitive            = " << str_aus.reciprocal_lattice_type << endl;
      ss_output << " Reciprocal lattice variation            = " << str_aus.reciprocal_lattice_variation_type << endl; //WSETYAWAN mod
    }
    // ---------------------------------------------------------------------------
    // json output
    else if(ftype == json_ft){

      aurostd::JSONwriter json;
      bool roff = true;

      // Reciprocal space lattice
      if(str_aus.klattice.rows != 0){
        json.addMatrix("reciprocal_lattice_vectors", str_aus.klattice, _AFLOWLIB_DATA_GEOMETRY_PREC_, roff);
      } else if(PRINT_NULL_JSON) {
        json.addNull("reciprocal_lattice_vectors");
      }

      // Reciprocal lattice parameters
      if(data.rows != 0){
        json.addVector("reciprocal_lattice_parameters", data, _AFLOWLIB_DATA_GEOMETRY_PREC_, roff);
      } else if(PRINT_NULL_JSON) {
        json.addNull("reciprocal_lattice_parameters");
      }

      // Reciprocal space volume
      json.addNumber("reciprocal_volume", kvol);

      // Reciprocal space: reciprocal lattice type
      if(!str_aus.reciprocal_lattice_type.empty()){
        json.addString("reciprocal_lattice_type", str_aus.reciprocal_lattice_type);
      } else if(PRINT_NULL_JSON) {
        json.addNull("reciprocal_lattice_type");
      }

      // Reciprocal space: reciprocal lattice variation type
      if(!str_aus.reciprocal_lattice_variation_type.empty()){
        json.addString("reciprocal_lattice_variation_type", str_aus.reciprocal_lattice_variation_type);
      } else if(PRINT_NULL_JSON) {
        json.addNull("reciprocal_lattice_variation_type");
      }
      ss_output << json.toString(standalone); //standalone: determines if we enclose in brackets
      if(standalone) { ss_output << endl; }
    }
    else{
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Format type is not supported.",_INPUT_ILLEGAL_);
    }

    // ---------------------------------------------------------------------------
    // update vpflow, needed for library runs
    // clear first, if needed
    if(vpflow.flag("EDATA::RECIPROCAL_LATTICE::CALCULATED")){
      vpflow.pop_attached("EDATA::RECIPROCAL_LATTICE_PARAMETERS");
      vpflow.pop_attached("EDATA::REAL_SPACE_VOLUME");
      vpflow.pop_attached("EDATA::RECIPROCAL_LATTICE_TYPE");
      vpflow.pop_attached("EDATA::RECIPROCAL_LATTICE_VARIATION_TYPE");
    }
    // update vpflow
    vpflow.flag("EDATA::RECIPROCAL_LATTICE::CALCULATED",TRUE);
    vpflow.push_attached("EDATA::RECIPROCAL_LATTICE_PARAMETERS",aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(data,_AFLOWLIB_DATA_GEOMETRY_PREC_,true),","));   //CO20200731 - precision
    vpflow.push_attached("EDATA::RECIPROCAL_SPACE_VOLUME",aurostd::utype2string<double>(kvol));
    vpflow.push_attached("EDATA::RECIPROCAL_LATTICE_TYPE",str_aus.reciprocal_lattice_type);
    vpflow.push_attached("EDATA::RECIPROCAL_LATTICE_VARIATION_TYPE",str_aus.reciprocal_lattice_variation_type);

    return ss_output.str();
  }
}

// ***************************************************************************
// pflow::PrintSuperlatticeData() //DX20210209
// ***************************************************************************
// Print the superlattice data (lattice, type/variation, Pearson, etc.)
// Also updates the xoption (EDATA::<PROPERTY>) for lib2raw runs
namespace pflow {
  string PrintSuperlatticeData(const xstructure& xstr, filetype ftype, bool standalone, bool already_calculated, double sym_eps){
    aurostd::xoption vpflow;
    return PrintSuperlatticeData(xstr, vpflow, ftype, standalone, already_calculated, sym_eps);
  }
}
namespace pflow {
  string PrintSuperlatticeData(const xstructure& xstr, aurostd::xoption& vpflow, filetype ftype, bool standalone, bool already_calculated, double sym_eps){


    // ---------------------------------------------------------------------------
    // calculate the superlattice symmetry if not already calculated
    xstructure str_aus(xstr);
    str_aus.ReScale(1.0);
    if(!already_calculated){ 
      str_aus.GetSuperlatticeType(sym_eps);
    }

    // get lattice and angles
    xvector<double> data(6);
    data=Getabc_angles(str_aus.bravais_superlattice_lattice,DEGREES);//DX20210201 - scaling factor from data(1)*=str_aus.scale;data(2)*=str_aus.scale;data(3)*=str_aus.scale;
    double c_over_a = data(3)/data(1);

    // get lattice and angles (Bohr)
    xvector<double> data_Bohr(6);
    data_Bohr(1) =  data(1)*angstrom2bohr;
    data_Bohr(2) =  data(2)*angstrom2bohr;
    data_Bohr(3) =  data(3)*angstrom2bohr;
    data_Bohr(4) =  data(4);
    data_Bohr(5) =  data(5);
    data_Bohr(6) =  data(6);

    // get volume
    double vol = GetVol(str_aus.bravais_superlattice_lattice)*str_aus.scale*str_aus.scale*str_aus.scale;

    stringstream ss_output;
    // ---------------------------------------------------------------------------
    // text output
    if(ftype == txt_ft){

      string info_prefix = "";
      if(!standalone){ info_prefix = "Real space "; }

      ss_output << "SUPERLATTICE (equally decorated)" << endl;
      ss_output << " Superlattice lattice:" << endl;
      ss_output << "  " << aurostd::roundoff(str_aus.bravais_superlattice_lattice(1),AUROSTD_ROUNDOFF_TOL) << endl;
      ss_output << "  " << aurostd::roundoff(str_aus.bravais_superlattice_lattice(2),AUROSTD_ROUNDOFF_TOL) << endl;
      ss_output << "  " << aurostd::roundoff(str_aus.bravais_superlattice_lattice(3),AUROSTD_ROUNDOFF_TOL) << endl;
      ss_output << " " << info_prefix << "a b c alpha beta gamma: ";
      ss_output.precision(10);ss_output << data(1) << " " << data(2) << " " << data(3) << " ";
      ss_output.precision(3); ss_output << data(4) << " " << data(5) << " " << data(6) << endl;
      ss_output.precision(4);
      ss_output << " " << info_prefix << "a b c alpha beta gamma: ";
      ss_output.precision(10);ss_output << data_Bohr(1) << " " << data_Bohr(2) << " " << data_Bohr(3) << " ";
      ss_output.precision(3); ss_output << data(4) << " " << data(5) << " " << data(6) << "   Bohrs/Degs " << endl;
      ss_output.precision(4);
      ss_output << " " << info_prefix << "Volume: " << vol << endl;
      ss_output << " " << info_prefix << "c/a = " << c_over_a << endl;
      ss_output << " " << info_prefix << "Bravais Superlattice Primitive   = " << str_aus.bravais_superlattice_type << endl;
      ss_output << " " << info_prefix << "Superlattice Variation           = " << str_aus.bravais_superlattice_variation_type << endl;
      ss_output << " " << info_prefix << "Superlattice System              = " << str_aus.bravais_superlattice_system << endl; //DX - fixed mistake in line above
      ss_output << " " << info_prefix << "Pearson Symbol Superlattice      = " << str_aus.pearson_symbol_superlattice << endl;
    }
    // ---------------------------------------------------------------------------
    // json output
    else if(ftype == json_ft){

      aurostd::JSONwriter json;
      bool roff = true;

      // Real space: bravais superlattice lattice
      if(!str_aus.bravais_superlattice_type.empty()){
        json.addMatrix("Bravais_superlattice_lattice", str_aus.bravais_superlattice_lattice, _AFLOWLIB_DATA_GEOMETRY_PREC_, roff);
      } else if(PRINT_NULL_JSON){
        json.addNull("Bravais_superlattice_lattice");
      }

      // Real lattice parameters
      if(data.rows != 0){
        json.addVector("lattice_parameters_superlattice", data, _AFLOWLIB_DATA_GEOMETRY_PREC_, roff);
      } else if(PRINT_NULL_JSON){
        json.addNull("lattice_parameters_superlattice");
      }

      // Real lattice parameters (Bohr/Deg)
      if(data_Bohr.rows != 0){
        json.addVector("lattice_parameters_superlattice_Bohr_deg", data_Bohr, _AFLOWLIB_DATA_GEOMETRY_PREC_, roff);
      } else if(PRINT_NULL_JSON){
        json.addNull("lattice_parameters_superlattice_Bohr_deg");
      }

      // Real space volume
      json.addNumber("volume_superlattice", vol);

      // Real space c/a
      json.addNumber("c_over_a_superlattice", c_over_a);

      // Real space: bravais superlattice type
      if(!str_aus.bravais_superlattice_type.empty()){
        json.addString("Bravais_superlattice_type", str_aus.bravais_superlattice_type);
      } else if(PRINT_NULL_JSON){
        json.addNull("Bravais_superlattice_type");
      }

      // Real space: bravais superlattice variation type
      if(!str_aus.bravais_superlattice_type.empty()){
        json.addString("Bravais_superlattice_variation_type", str_aus.bravais_superlattice_variation_type);
      } else if(PRINT_NULL_JSON){
        json.addNull("Bravais_superlattice_variation_type");
      }

      // Real space: bravais superlattice system
      if(!str_aus.bravais_superlattice_type.empty()){
        json.addString("Bravais_superlattice_system", str_aus.bravais_superlattice_system);
      } else if(PRINT_NULL_JSON){
        json.addNull("Bravais_superlattice_system");
      }

      // Real space: bravais superlattice type
      if(!str_aus.bravais_superlattice_type.empty()){
        json.addString("Pearson_symbol_superlattice", str_aus.pearson_symbol_superlattice);
      } else if(PRINT_NULL_JSON){
        json.addNull("Pearson_symbol_superlattice");
      }
      ss_output << json.toString(standalone); //standalone: determines if we enclose in brackets
      if(standalone) { ss_output << endl; }
    }
    else{
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Format type is not supported.",_INPUT_ILLEGAL_);
    }

    // ---------------------------------------------------------------------------
    // update vpflow, needed for library runs
    // clear first, if needed
    if(vpflow.flag("EDATA::SUPERLATTICE::CALCULATED")){
      vpflow.pop_attached("EDATA::BRAVAIS_SUPERLATTICE_TYPE");
      vpflow.pop_attached("EDATA::BRAVAIS_SUPERLATTICE_VARIATION_TYPE");
      vpflow.pop_attached("EDATA::BRAVAIS_SUPERLATTICE_SYSTEM");
      vpflow.pop_attached("EDATA::PEARSON_SYMBOL_SUPERLATTICE");
    }
    // update vpflow
    vpflow.flag("EDATA::SUPERLATTICE::CALCULATED",TRUE);
    vpflow.push_attached("EDATA::BRAVAIS_SUPERLATTICE_TYPE",str_aus.bravais_superlattice_type);
    vpflow.push_attached("EDATA::BRAVAIS_SUPERLATTICE_VARIATION_TYPE",str_aus.bravais_superlattice_variation_type);
    vpflow.push_attached("EDATA::BRAVAIS_SUPERLATTICE_SYSTEM",str_aus.bravais_superlattice_system);
    vpflow.push_attached("EDATA::PEARSON_SYMBOL_SUPERLATTICE",str_aus.pearson_symbol_superlattice);

    return ss_output.str();
  }
}

// ***************************************************************************
// pflow::PrintData1
// ***************************************************************************
// Prints out the comparison between two xstructures.
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintData1(const xstructure& str1, const double& rcut, ostream& oss) {
    oss << PrintData1(str1,rcut);
  }
} // namespace pflow

namespace pflow {
  string PrintData1(const xstructure& str1, const double& rcut) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "pflow::PrintData1: BEGIN" << endl;
    // LDEBUG=TRUE;
    stringstream oss;
    oss.setf(std::ios::left, std::ios::adjustfield);
    oss.setf(std::ios::fixed, std::ios::floatfield);
    oss.precision(4);
    double cutoff=rcut;
    //  oss.width(10);
    // Rescale to avoid scale problems
    xstructure sstr1=str1;
    sstr1=ReScale(sstr1,1.0);
    int nat1=sstr1.atoms.size();
    int ntyp1=sstr1.num_each_type.size();
    std::deque<int> num_each_type1=sstr1.num_each_type;

    if(LDEBUG) cerr << "PrintData1 [2]" << endl;
    oss << "*************************  Structure Comparison *************************" << endl;

    // Nun atoms
    oss << "Number of atoms" << endl;
    oss << "  Num Atom Str1: " << nat1 << endl;

    // Num types
    oss << "Number of types" << endl;
    oss << "  Num Type Str1: " << ntyp1 << endl;

    // Stoichiometry
    int wid=8;
    vector<double> tmp_stoich1(ntyp1,0);
    oss << "Stoichiometry" << endl;
    oss << "  Type:        ";
    for(int it=0;it<ntyp1;it++) {
      oss << setw(wid) << it << " ";
    }
    oss << endl;
    oss << "  Stoich Str1: ";
    for(int it=0;it<ntyp1;it++) {
      tmp_stoich1[it]=(double) num_each_type1[it]/ (double) nat1;
      oss << setw(wid) << tmp_stoich1[it] << " ";
    }
    oss << endl;

    if(LDEBUG) cerr << "PrintData1 [3]" << endl;
    // Volume
    double v1=sstr1.scale*sstr1.scale*sstr1.scale*GetVol(sstr1.lattice)/nat1;
    oss << "Volume per atom" << endl;
    oss << "  Vol Str1: " << v1 << endl;

    // lat par and angles
    //  vector<double> ldat1=pflow::Getabc_angles(aurostd::xmatrix2matrix(sstr1.lattice)); //CO20200404 pflow::matrix()->aurostd::matrix()
    vector<double> ldat1=xvector2vector(Getabc_angles(sstr1.lattice,DEGREES));
    ldat1=pflow::Sort_abc_angles(ldat1);  
    oss << "Unit cell: tot_vol a b c alpha beta gamma" << endl;
    oss << "  Cell Str1: ";
    oss << v1*nat1 << " ";
    pflow::Vout(ldat1,oss);

    // Volume adjusted distances

    double rescale1=std::pow((double) 1.0/v1,(double) 1.0/3.0);
    xstructure str1_newvol=sstr1;
    // str1_newvol=SetScale(str1_newvol,rescale1); // Sets vol/atom=1
    str1_newvol.scale=rescale1;
    aurostd::matrix<double> dist1; //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> dist2; //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> dist_diff; //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> dist_diff_n; //CO20200404 pflow::matrix()->aurostd::matrix()
    // If cutoff < 0 use default of 4*(avg rad of an atom).  Note
    // that in the rescaled str, where each atom has vavg=1, rad=.6204
    // In the true str, where each atom has vol = v1, rad=0.6204*v1^(1/3).
    if(LDEBUG) cerr << "PrintData1 [3]" << endl;
    if(rcut<0) cutoff=4*0.6204*std::pow((double) v1,(double) 1.0/3.0);  // THAT`S WRONG
    //  cerr << cutoff << endl;
    if(rcut<0) cutoff=1.01*RadiusSphereLattice(str1_newvol.lattice); // better way to nail it for sure I`ll hit something. 1% to be sure of roundoff
    //  cerr << cutoff << endl;
    //  if(rcut<0) cutoff=max(modulus(str1_newvol.lattice(1)),modulus(str1_newvol.lattice(2)),modulus(str1_newvol.lattice(3))); // better way to nail it for sure I`ll hit something.
    //  cerr << cutoff << endl;
    if(LDEBUG) {cerr << "PrintData1 [4]  " << 4*0.6204*std::pow((double) v1,(double) 1.0/3.0)<< " " << RadiusSphereLattice(str1_newvol.lattice) << " " << max(modulus(str1_newvol.lattice(1)),modulus(str1_newvol.lattice(2)),modulus(str1_newvol.lattice(3))) << endl;}
    pflow::CmpStrDist(str1_newvol,str1_newvol,cutoff*rescale1,dist1,dist2,dist_diff,dist_diff_n);
    int nprs1=dist1.size();
    oss << "Avg of Distance Magnitude Differences" << endl;
    oss << "  Rescaling vol/atom of Str1 to the average vol/atom (return to orig scale factor): "
      << str1_newvol.scale*str1_newvol.scale*str1_newvol.scale*GetVol(str1_newvol.lattice)/nat1 << " (" << 1/rescale1 << ")" << endl;
    oss << "  Working with the all the types up to num types (num pair types): " << ntyp1 << " (" << nprs1 << ")" << endl;
    oss << "  Including all pairs within cutoff: " << cutoff << endl;

    // double tadiff=0;   //DM not used
    // double tadiff_n=0; //DM not used
    // int tnn=0;  //DM not used
    for(int it1=0;it1<ntyp1;it1++) {
      for(int it2=it1;it2<ntyp1;it2++) {
        int id=it2-it1+ntyp1*it1-std::max(0,it1*(it1-1)/2);
        int nn=dist1[id].size();      
        oss << "  P" << id << "  " << nn << " Pairs between types: "
          << it1 << " " << it2 << "   ";
        for(int in=0;in<nn;in++) {
          oss << " " << dist1[id][in];
        }
        oss << endl;
      }// it2
    }// it1

    //    // Number of each type up to rcut.
    //    nprs_max=ntyp_max*(ntyp_max+1)/2;
    //    vector<int> n1(nprs_max,0);
    //    vector<int> n2(nprs_max,0);
    //    oss << "Number neighbors of each type within cutoff for rescaled volume" << endl;
    //    oss << "  Cutoff: " << cutoff << endl;
    //    oss << "  Rescaling vol/atom of Str1 and Str2 to the average vol/atom: " << vavg << endl;
    //    oss << "           ";
    //    for(int it1=0;it1<ntyp_max;it1++) {
    //      for(int it2=it1;it2<ntyp_max;it2++) {
    //        oss << it1 << "-" << it2 << "  ";
    //      }
    //    }
    //    oss << endl;
    //    oss << "  NN Str1  ";
    //    for(int it1=0;it1<ntyp_max;it1++) {
    //      for(int it2=it1;it2<ntyp_max;it2++) {
    //        if(it1<ntyp1 && it2<ntyp1) {
    //  	int id1=it2-it1+ntyp1*it1-std::max(0,it1*(it1-1)/2);
    //  	n1[id1]=dist1[id1].size();
    //  	oss << setw(4) << n1[id1] << " ";
    //        }
    //      }
    //    }
    //    oss << endl;
    //    oss << "  NN Str2  ";
    //    for(int it1=0;it1<ntyp_max;it1++) {
    //      for(int it2=it1;it2<ntyp_max;it2++) {
    //        if(it1<ntyp2 && it2<ntyp2) {
    //  	int id2=it2-it1+ntyp2*it1-std::max(0,it1*(it1-1)/2);
    //  	n2[id2]=dist2[id2].size();
    //  	oss << setw(4) << n2[id2] << " ";
    //        }
    //      }
    //    }
    //    oss << endl;  
    //    oss << "  NN Err   ";
    //    for(int it1=0;it1<ntyp_max;it1++) {
    //      for(int it2=it1;it2<ntyp_max;it2++) {
    //        int id=it2-it1+ntyp_max*it1-std::max(0,it1*(it1-1)/2);
    //        oss << setw(4) << n2[id]-n1[id] << " ";
    //      }
    //    }
    //    oss << endl;  

    //    // Space group
    //    bool Platon_EQUAL=FALSE,Platon_EXACT=FALSE;
    //    double Platon_ang=1.0,Platon_d1=0.25,Platon_d2=0.25,Platon_d3=0.25;
    //    system("touch pflow_tmp_in"); // Create input file to connect to stream.
    //    ofstream outfile("pflow_tmp_out");
    //    ifstream infile("pflow_tmp_in");
    //    string spcgrp1a,spcgrp2a,spcgrp1b,spcgrp2b;
    //    platon2print(str1,Platon_EQUAL,Platon_EXACT,Platon_ang,Platon_d1,Platon_d2,Platon_d3,outfile);
    //    system("cat pflow_tmp_out | platonSG > pflow_tmp_in");
    //    infile >> spcgrp1a >> spcgrp1b;
    //    // Puts file position pointers back to beginning.
    //    outfile.seekp(0);
    //    infile.seekg(0);
    //    platon2print(sstr2,Platon_EQUAL,Platon_EXACT,Platon_ang,Platon_d1,Platon_d2,Platon_d3,outfile);
    //    system("cat pflow_tmp_out | platonSG > pflow_tmp_in");
    //    infile >> spcgrp2a >> spcgrp2b;
    //    system("/bin/rm pflow_tmp_out");
    //    system("/bin/rm pflow_tmp_in");
    //    oss << "Space groups " << endl;
    //    oss << "  SG Str1: " << spcgrp1a << " " << spcgrp1b << endl;
    //    oss << "  SG Str2: " << spcgrp2a << " " << spcgrp2b << endl;

    // Center of mass
    return oss.str();
  } // end routine
} // namespace pflow

// ***************************************************************************
// pflow::PrintData2
// ***************************************************************************
// This funtion prints out structural data.
// Stefano Curtarolo
namespace pflow {
  void PrintData2(const xstructure& str, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "pflow::PrintData2: BEGIN" << endl;
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(6);
    xvector<double> data(6);
    oss << endl;
    oss << "DIRECT LATTICE" << endl;
    data=Getabc_angles(str.lattice,DEGREES);
    oss << "volume = " << str.scale*str.scale*str.scale*GetVol(str.lattice) << endl;
    oss << "a,b,c,alpha,beta,gamma = "
      << str.scale*data(1) << " " << str.scale*data(2) << " " << str.scale*data(3) << " "
      << data(4) << " " << data(5) << " " << data(6) << endl;
    oss << "b/a = " << data(2)/data(1) << "  c/a = " << data(3)/data(1) << endl;
    oss << "Lattice =" << endl;
    oss << str.scale*str.lattice << endl;
    oss << endl;      
    oss << "RECIPROCAL LATTICE" << endl;
    data=Getabc_angles(str.klattice,DEGREES);
    oss << "volume = " << GetVol(str.klattice) << endl;
    oss << "a,b,c,alpha,beta,gamma = "
      << str.scale*data(1) << " " << str.scale*data(2) << " " << str.scale*data(3) << " "
      << data(4) << " " << data(5) << " " << data(6) << endl;
    oss << "K-Lattice =" << endl;
    oss << str.klattice << endl;
    oss << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::PrintDisplacements
// ***************************************************************************
//  This funtion prints out displacements
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintDisplacements(xstructure str, const double cutoff, ostream& oss) {  
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "pflow::PrintDisplacements: BEGIN" << endl;
    oss.setf(std::ios::left, std::ios::adjustfield);
    oss.setf(std::ios::fixed, std::ios::floatfield);
    //  oss.setf(std::ios::fixed,std::ios::floatfield);
    //  oss.precision(10);
    deque<deque<_atom> > neigh_mat;
    // [OBSOLETE]  pflow::GetNeighData(str,cutoff,neigh_mat);
    str.GetNeighData(cutoff,neigh_mat); //CO20220623 - using new GetNeighData()
    // double tol=1e-15;   //DM not used
    xmatrix<double> lattice(3);
    lattice=str.lattice;
    int coord_flag=str.coord_flag;

    // Output header
    oss << "All neighbor data given as" << endl;
    oss << "[basis number]  [name]  [n1 n2 n3 (unit cell)] [d1 d2 d3 (displacement)] [distance]" << endl;
    oss << "Distances are only calculated for the neighbors within " << cutoff <<  endl;

    // Print out dislpacements and distances
    oss << endl;
    oss << "********************" << " DISPLACEMENTS " << "********************" << endl;
    for(uint ia=0;ia<neigh_mat.size();ia++) {
      _atom a = neigh_mat.at(ia).at(0);

      // Output reference atom info.
      oss << setw(4) << a.basis+1 << " " << setw(4) << a.name.c_str(); //[CO20200130 - number->basis]oss << setw(4) << a.number+1 << " " << setw(4) << a.name.c_str();
      xvector<double> pos(3);
      if(str.coord_flag==FALSE) { // direct
        pos=a.fpos;
      } else {
        pos=a.cpos;
      }
      oss << "   " << setw(9) << setprecision(4) << pos(1)
        << setw(9) << setprecision(4) << pos(2)
        << setw(9) << setprecision(4) << pos(3)
        << endl;

      for(uint in=1;in<neigh_mat.at(ia).size();in++) {
        _atom an = neigh_mat.at(ia).at(in);
        xvector<int> ijk(3);
        ijk=an.ijk;
        oss << "        ";
        oss << setw(4) << an.basis+1 << " "; //[CO20200130 - number->basis]oss << setw(4) << an.number+1 << " ";
        oss << setw(4) << an.name.c_str() << "   ";
        oss << setw(3) << ijk(1) << " " << setw(3) << ijk(2) << " " << setw(3) << ijk(3) << "   ";
        xvector<double> disp(3);
        disp=AtomCDisp(a,an);
        if(coord_flag==0) { // direct
          disp=C2F(lattice,disp);
        }
        oss << setw(9) << setprecision(4) << disp(1) << " " << setw(9) << setprecision(4) << disp(2) << " " << setw(9) << setprecision(4) << disp(3);
        oss << "   " << setw(9) << setprecision(4) << AtomDist(a,an);
        oss << endl;
      } // in
    } // ia

  } // end routine
} // namespace pflow

// ***************************************************************************
// pflow::PrintDistances
// ***************************************************************************
// This funtion prints out distances
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintDistances(xstructure str, const double cutoff, ostream& oss) {  
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "pflow::PrintDistances: BEGIN" << endl;
    oss.setf(std::ios::left, std::ios::adjustfield);
    oss.setf(std::ios::fixed, std::ios::floatfield);
    //  oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(10);
    deque<deque<_atom> > neigh_mat;
    // [OBSOLETE] pflow::GetNeighData(str,cutoff,neigh_mat);  //CO20220623 - using new GetNeighData()
    str.GetNeighData(cutoff,neigh_mat);                       //CO20220623 - using new GetNeighData()
    //  cerr << "DEBUG: neigh_mat.size()=" << neigh_mat.size() << endl;
    //  for(uint i=0;i<neigh_mat.size();i++)
    //   cerr << "DEBUG: neigh_mat.at(i).size()=" << neigh_mat.at(i).size() << endl;

    // double tol=1e-15;  //DM not used
    uint xtra=0;//3;

    // Output header
    oss << "All atoms identified by" << endl;
    oss << "[basis number]  [name]  [n1 n2 n3 (unit cell)]" << endl;
    oss << "Distances are only calculated for the neighbors within " << cutoff <<  endl;

    // Print out distances
    oss << endl;
    oss << "********************" << " DISTANCES " << "********************" << endl;
    for(uint ia=0;ia<neigh_mat.size();ia++) {
      _atom a = neigh_mat.at(ia).at(0);

      // Output reference atom info.
      oss << setw(4) << a.basis+1 << " " << setw(4) << a.name.c_str(); //[CO20200130 - number->basis]oss << setw(4) << a.number+1 << " " << setw(4) << a.name.c_str();
      xvector<double> pos(3);
      if(str.coord_flag==FALSE) { // direct
        pos=a.fpos;
      } else {
        pos=a.cpos;
      }
      oss << "   "
        << setw(9+xtra) << setprecision(4+xtra) << pos(1)
        << setw(9+xtra) << setprecision(4+xtra) << pos(2)
        << setw(9+xtra) << setprecision(4+xtra) << pos(3)
        << endl;

      for(uint in=1;in<neigh_mat.at(ia).size();in++) {
        xvector<int> ijk(3);ijk=neigh_mat.at(ia).at(in).ijk;
        xvector<double> fpos(3);fpos=neigh_mat.at(ia).at(in).fpos;
        xvector<double> cpos(3);cpos=neigh_mat.at(ia).at(in).cpos;
        oss << "      ";
        oss << setw(4) << neigh_mat.at(ia).at(in).basis+1 << " ";  //[CO20200130 - number->basis]oss << setw(4) << neigh_mat.at(ia).at(in).number+1 << " ";
        oss << setw(4) << neigh_mat.at(ia).at(in).name.c_str() << "   ";
        oss << setw(3+xtra) <<  ijk(1) << " " << setw(3+xtra) <<  ijk(2) << " " << setw(3+xtra) <<  ijk(3) << "   ";
        //     oss << setw(3+xtra) << "F " << fpos(1) << " " << setw(3+xtra) << fpos(2) << " " << setw(3+xtra) << fpos(3) << "   ";
        //   oss << setw(3+xtra) << "C " << cpos(1) << " " << setw(3+xtra) << cpos(2) << " " << setw(3+xtra) << cpos(3) << "   ";
        oss << setprecision(4+xtra) << AtomDist(a,neigh_mat.at(ia).at(in));
        oss << endl;
      } // in
    } // ia

  } // end routine
} // namespace pflow

// ***************************************************************************
// pflow::PrintEwald
// ***************************************************************************
// Dane Morgan - Stefano Curtarolo
namespace pflow {
  void PrintEwald(const xstructure& in_str, double& epoint,
      double& ereal, double& erecip, double& eewald,
      double& eta, const double& SUMTOL, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "pflow::PrintEwald: BEGIN" << endl;
    oss.setf(std::ios::left, std::ios::adjustfield);
    oss.setf(std::ios::fixed, std::ios::floatfield);
    aurostd::matrix<double> lat=pflow::GetScaledLat(in_str); //CO20200404 pflow::matrix()->aurostd::matrix()
    int nat = pflow::GetNumAtoms(in_str);
    double vol=pflow::GetVol(lat);
    double d=std::pow((double) vol,(double) 1.0/3.0)*1E-10; // In meters
    double eV2NoUnit=2*TWOPI*EPS_VACUUM*d/(E_ELECTRON*nat);
    int wid;
    int prec;
    wid=15;
    prec=6;
    oss << "*************************** Ewald Results ***********************" << endl;
    oss << "E(NoUnit)=PerAtom(eV)*4*pi*eps_0*vol^1/3/e_electron=Charge/Distance unit independent."<< endl;
    oss << "Using eta= " << setprecision(2*prec) << eta << endl;
    oss << "All real and k-space sums stop when maximum term in shell is less than SUMTOL= " << std::setiosflags(std::ios::scientific) << SUMTOL << endl;
    oss << "SUMTOL is the following percentage of your total Ewald energy: " << std::setiosflags(std::ios::scientific) << 100*SUMTOL/eewald << endl;
    oss.setf(std::ios::fixed, std::ios::floatfield);
    oss << "                         "
      << setw(wid) << "Total(eV)"
      << setw(wid) << "Per atom(eV)"
      << setw(wid) << "Total(NoUnit)" << endl;
    oss << "Point Energy:            "
      << " " << setw(wid) << setprecision(prec) << epoint
      << " " << setw(wid) << setprecision(prec) << epoint/nat
      << " " << setw(wid) << setprecision(prec) << eV2NoUnit*epoint
      << endl;
    oss << "Real Space Energy:       "
      << " " << setw(wid) << setprecision(prec) << ereal
      << " " << setw(wid) << setprecision(prec) << ereal/nat
      << " " << setw(wid) << setprecision(prec) << eV2NoUnit*ereal
      << endl;
    oss << "Recipricol Space Energy: "
      << " " << setw(wid) << setprecision(prec) << erecip
      << " " << setw(wid) << setprecision(prec) << erecip/nat
      << " " << setw(wid) << setprecision(prec) << eV2NoUnit*erecip
      << endl;
    oss << "Total Ewald Energy:      "
      << " " << setw(wid) << setprecision(prec) << eewald
      << " " << setw(wid) << setprecision(prec) << eewald/nat
      << " " << setw(wid) << setprecision(prec) << eV2NoUnit*eewald
      << endl;
  }
} // namespace pflow

// ***************************************************************************
// pflow::PrintGulp
// ***************************************************************************
// This funtion prints out structural data in a format for
// gulp (file will cause gulp to calculate distances between
// atoms when used as input).
// Stefano Curtarolo
namespace pflow {
  void PrintGulp(const xstructure& str, ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "pflow::PrintGulp: BEGIN" << endl;
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(10);
    xstructure sstr=str;
    // Set scale to 1 so you don't need to rescale coordinates.
    sstr=ReScale(sstr,1.0);

    // Print out data in Gulp format
    oss << "distance" << endl;
    oss << "cutd 5" << endl;
    oss << "vectors" << endl;
    oss << sstr.lattice(1,1) << " " << sstr.lattice(1,2) << " " << sstr.lattice(1,3) << "    Bravais lattice" << endl;
    oss << sstr.lattice(2,1) << " " << sstr.lattice(2,2) << " " << sstr.lattice(2,3) << "    Bravais lattice" << endl;
    oss << sstr.lattice(3,1) << " " << sstr.lattice(3,2) << " " << sstr.lattice(3,3) << "    Bravais lattice" << endl;

    for(uint i=0;i<sstr.num_each_type.size();i++)
      oss <<  sstr.num_each_type.at(i) << " ";

    oss << endl;
    if(sstr.coord_flag==0) oss << "Fractional" << endl;
    if(sstr.coord_flag==1) oss << "Cartesian" << endl;
    for(uint i=0;i<sstr.atoms.size();i++) {
      if(sstr.coord_flag==0)
        oss << " " << sstr.atoms.at(i).fpos(1) << " " << sstr.atoms.at(i).fpos(2) << " " << sstr.atoms.at(i).fpos(3) << endl;
      if(sstr.coord_flag==1) oss << " " << sstr.atoms.at(i).cpos(1) << " " << sstr.atoms.at(i).cpos(2) << " " << sstr.atoms.at(i).cpos(3) << endl;
    }
  }
} // namespace pflow

// ***************************************************************************
// PrintKmesh
// ***************************************************************************
//  This funtion prints a set of Monkhurst-Pack kpoint mesh values
// (divisions along each reciprical lattice parameter)
// based on the desired kpt density and the lattice parameters.  The
// values assure the most even distribution of kpts along the lattice
// params consistent with the kpt density.  Returns 4 kmesh vectors of length
// 3
//    1. Rounded to Nint
//    2. Exact
//    3. Rounded down
//    4. Rounded up
// Dane Morgan - Stefano Curtarolo
void PrintKmesh(const xmatrix<double>& kmesh, ostream& oss) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(3);
  vector<string> s(4);
  s[0]="RoundNint:";
  s[1]="Exact:    ";
  s[2]="RoundDown:";
  s[3]="RoundUp:  ";
  for(int i=1;i<=4;i++)
    oss << s[i-1] << "  " << kmesh(i,1) << "  " << kmesh(i,2) << "  " << kmesh(i,3) << endl;
}

// ***************************************************************************
// PrintImages
// ***************************************************************************
//  This function prints out all the image files.
// Dane Morgan - Stefano Curtarolo
void PrintImages(xstructure strA, xstructure strB, const int& ni, const string& path_flag) {  
  // Get cm of strA. Get image at i/(n-1) between
  // the A and B.  Set mom1 of image to be that of strA.  
  // Output that image as image i.
  xvector<double> mom1A=GetMom1(strA);
  for(int im=0;im<ni+2;im++) {
    double f;
    if((ni+2)>1) {
      f=(double)im/double(ni+1);
    } else {
      f=0;
    }
    xstructure imstr=GetIntpolStr(strA,strB,f,path_flag);
    imstr=SetMom1(imstr,mom1A);
    // Set title to include image name.
    ostringstream titlestream;
    if(im<10) {
      titlestream << "POSCAR_0" << im << ":" << imstr.title << ends;
    } else {
      titlestream << "POSCAR_" << im << ":" << imstr.title << ends;
    }
    imstr.title=titlestream.str();
    // Set output file name.
    ostringstream namestream;
    if(im<10) {
      namestream << "POSCAR_0" << im << ends;
    } else {
      namestream << "POSCAR_" << im << ends;
    }
    //  char* namechar = namestream.str().c_str();
    //  ofstream outfile(namechar);
    ofstream outfile(namestream.str().c_str());
    outfile << imstr;
    // Create directory
    ostringstream aus;
    if(im<10) {
      aus << "mkdir 0" << im << ends;
      aurostd::execute(aus);
      aus << "cp " << namestream.str().c_str() << " " << "./" << "0" << im << "/POSCAR" << endl;
      aurostd::execute(aus);
    } else {
      aus << "mkdir " << im << ends;
      aurostd::execute(aus);
      aus << "cp " << namestream.str().c_str() << " " << "./" << im << "/POSCAR" << endl;
      aurostd::execute(aus);
    }
  }
} // end routine

// ***************************************************************************
// PrintMSI
// ***************************************************************************
//  This funtion prints out structural data in a msi format (for cerius).
// Dane Morgan - Stefano Curtarolo
void PrintMSI(const xstructure& str, ostream& oss) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  xstructure sstr=str;

  // Set scale to 1 so you don't need to rescale coordinates.
  sstr=ReScale(sstr,1.0);

  // Print out data in msi format

  // Tolerance
  double Tol=1e-10;

  // First we must determine the atomic numbers.
  for(uint i=0;i<sstr.atoms.size();i++) {
    sstr.atoms.at(i).CleanName();
    if(sstr.atoms.at(i).atomic_number<1) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"atomic_number not found: sstr.atoms.at("+aurostd::utype2string(i)+").cleanname="+sstr.atoms.at(i).cleanname,_INPUT_ILLEGAL_);  //CO20200624
    }
  }

  //   Now we have to orient the lattice vectors for msi.
  //   The desired orientation is obtained by
  //   rotating c around Z into the XZ plane (R1a), then rotating c around
  //   Y to align it with Z (R1b), then rotating b around Z (leaving c fixed)
  //   until it is in the YZ plane (R2).  The total rotation needed is then
  //   Rtot=R2*R1b*R1a.

  xvector<double> v(3);
  xmatrix<double> R1a(3,3);  
  xmatrix<double> R1b(3,3);  
  xmatrix<double> R2(3,3);  
  xmatrix<double> Rtot(3,3);  
  xvector<double> a(3);a=sstr.lattice(1);//double an=modulus(a); //DM not used
  xvector<double> b(3);b=sstr.lattice(2);//double bn=modulus(b); //DM not used
  xvector<double> c(3);c=sstr.lattice(3);double cn=modulus(c);

  double cxyn=sqrt(cn*cn-c(3)*c(3));
  // double bxyn=sqrt(bn*bn-b(3)*b(3)); //DM not used
  // double cosa=cos(a,b); //DM not used
  // double cosb=cos(a,c); //DM not used
  // double cosg=cos(b,c); //DM not used

  // Build R1a.
  double cosp,sinp;
  if(cxyn>Tol && (abs(c(1))>Tol || abs(c(2))>Tol)) {
    cosp=c(1)/cxyn;
    sinp=c(2)/cxyn;
  } else {
    cosp=1.0;
    sinp=0.0;
  }
  clear(R1a);
  R1a(1,1)=cosp;
  R1a(1,2)=sinp;
  R1a(2,1)=-sinp;
  R1a(2,2)=cosp;
  R1a(3,3)=1.0;

  // tpx
  //  Mout(R1a,oss);
  //  oss << "cosp, sinp" << cosp << " " << sinp << endl;
  // Build R1b.
  double costh=c(3)/cn;
  double sinth=sqrt(1-costh*costh); // This is OK since sinth>0.
  R1b(1,1)=costh;
  R1b(1,3)=-sinth;
  R1b(3,1)=sinth;
  R1b(3,3)=costh;
  R1b(2,2)=1.0;
  xmatrix<double> RR(3,3);
  RR=R1b*R1a;
  // Build R2.  
  // Using rotated b vector, called bp, (bp=R1b*R1a*b) to define the angles.
  xvector<double> bp(3);
  bp=RR*b;
  double bpn=modulus(bp);
  double bpxyn=sqrt(bpn*bpn-bp(3)*bp(3));
  double cospp,sinpp;
  if(bpxyn>Tol && (b(1)>Tol || b(2)>Tol)) {
    cospp=bp(2)/bpxyn;
    sinpp=bp(1)/bpxyn;
  } else {
    cospp=1.0;
    sinpp=0.0;
  }
  clear(R2);
  R2(1,1)=cospp;
  R2(1,2)=-sinpp;
  R2(2,1)=sinpp;
  R2(2,2)=cospp;
  R2(3,3)=1.0;
  // Calculate final rotation xmatrix.
  // tpx
  //  for(int i=1;i<=3;i++)
  //  for(int j=1;j<=3;j++)
  //  oss << "R2 " << R2(i,j) << endl;
  Rtot=R2*RR;

  //tpx
  //Mout(Rtot,oss);
  //oss << "HERE " << cxyn << " " << cosa << " " << cosb << " " << cosg << " " << cosp << " " << sinp << endl;
  // Now rotate all the lattice vectors and the positions.
  xmatrix<double> c2lat(3,3);
  std::vector<xvector<double> > c2pos(sstr.atoms.size());

  v=Rtot*a;
  for(int i=1;i<=3;i++)
    c2lat(1,i)=v(i);
  v=Rtot*b;
  for(int i=1;i<=3;i++)
    c2lat(2,i)=v(i);
  v=Rtot*c;
  for(int i=1;i<=3;i++)
    c2lat(3,i)=v(i);

  for(uint i=0;i<sstr.atoms.size();i++) {
    v=Rtot*sstr.atoms.at(i).cpos;
    c2pos.at(i)(1)=v(1);
    c2pos.at(i)(2)=v(2);
    c2pos.at(i)(3)=v(3);
  }

  // Write output.
  // Title.
  oss << "# MSI CERIUS2 DataModel File Version 3 5" << endl;

  // Stuff.
  oss << "(1 Model" << endl;

  // Lattice vectors.
  oss << " (A D A3 (" << c2lat(1,1) << " " << c2lat(1,2) << " " << c2lat(1,3) << " " << "))" << endl;
  oss << " (A D B3 (" << c2lat(2,1) << " " << c2lat(2,2) << " " << c2lat(2,3) << " " << "))" << endl;
  oss << " (A D C3 (" << c2lat(3,1) << " " << c2lat(3,2) << " " << c2lat(3,3) << " " << "))" << endl;

  // Stuff.
  oss << " (A I Id 1)" << endl;
  oss << " (A C Label \"Model1\")" << endl;

  // The atoms.
  for(uint i=0;i<sstr.atoms.size();i++) {
    oss << " (" << i+2 << " Atom" << endl;
    oss << "  (A C ACL \"" << sstr.atoms.at(i).atomic_number << " " << sstr.atoms.at(i).cleanname << "\")" << endl;
    oss << "  (A D XYZ (" << c2pos.at(i)(1) << " " << c2pos.at(i)(2) << " " << c2pos.at(i)(3) << " " << "))" << endl;
    oss << "  (A I Id " << i+1 << ")" << endl;
    oss << " )" << endl;
  }
  // Final parenthesis
  oss << ")" << endl;

}

// ***************************************************************************
// PrintNData
// ***************************************************************************
// This funtion prints out structural data.
// STefano Curtarolo
void PrintNdata(const xstructure& str, ostream& oss) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  xstructure sstr=str;

  // Set scale to 1 so you don't need to rescale coordinates.
  sstr=ReScale(sstr,1.0);
  // Print out data in Platon format
  double volume=GetVol(sstr.lattice);
  double normalization=std::pow((double) volume,(double) 1.0/3.0);

  printf(" %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f \n",
      modulus(sstr.lattice(1))/normalization,
      modulus(sstr.lattice(2))/normalization,
      modulus(sstr.lattice(3))/normalization,
      aurostd::angle(sstr.lattice(2,1),sstr.lattice(2,2),sstr.lattice(2,3),sstr.lattice(3,1),sstr.lattice(3,2),sstr.lattice(3,3)),
      aurostd::angle(sstr.lattice(1,1),sstr.lattice(1,2),sstr.lattice(1,3),sstr.lattice(3,1),sstr.lattice(3,2),sstr.lattice(3,3)),
      aurostd::angle(sstr.lattice(1,1),sstr.lattice(1,2),sstr.lattice(1,3),sstr.lattice(2,1),sstr.lattice(2,2),sstr.lattice(2,3)));//,volume);
}

// ***************************************************************************
// PrintNeatProj
// ***************************************************************************
// Dane Morgan - Stefano Curtarolo
void PrintNeatProj(pflow::projdata& pd, ostream& oss) {
  //  namespace pflow
  pd.Print(oss);
}

// ***************************************************************************
// PrintPDB
// ***************************************************************************
// This funtion prints out structural data in a PDB format.
// Note that I don't really know anything about PDB but it seems
// to be an annoyingly column formatted.  Therefore, all the widths
// below must tbe kept as is.  This means that if we have width W
// for variable X, and prec P, then X takes up P+1 for for the decimal
// part and decimal point, and we have only R=W-(P+1)-1=W-P-2 remaining digits
// (the -1 is because if X takes up all of W then you run into
// the previous field).  So we have the constraints
// Cell vectors: W=9,P=3,R=4 => <=10^5
// Cell angles: W=7,P=2,R=3 => <=10^3 (which always works since angles are
// given as >=0 and <=360)
// Cartesian positions: W=12,8,P=3,R=7,3 => <=10^8,10^4 (>0)
// and <=10^7,10^3 (<0 since you need a
// space for the - sign)
// This is all probably fine unless an atom makes it to more negative than
// -999.999.
// Dane Morgan - Stefano Curtarolo
void PrintPDB(const xstructure& str, ostream& oss) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  xstructure sstr=str;

  // Set scale to 1 so you don't need to rescale coordinates.
  sstr=ReScale(sstr,1.0);

  // Print out data in pdb format
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(3);
  oss << "HEADER    " << sstr.title << endl;
  oss << "REMARK   1" << endl;
  oss << "REMARK   PDB Format: Created by convasp from POSCAR file format." << endl;
  oss << "REMARK   1" << endl;
  xvector<double> data(6);
  data=Getabc_angles(sstr.lattice,DEGREES);//  no need to rescale
  oss << "CRYST1" << setw(9) << data(1)
    << setw(9) << data(2)
    << setw(9) << data(3);
  oss.precision(2);
  oss << setw(7) << data(4)
    << setw(7) << data(5)
    << setw(7) << data(6);
  // Output space group - Use P1 for now
  oss << setw(11) << "P1" << endl;
  // oss << "              P1" << endl;
  // Coordinates in cartesian
  oss.precision(3);
  for(uint i=0;i<sstr.atoms.size();i++) {
    oss << "ATOM  " << setw(5) << i+1 << setw(4) << sstr.atoms.at(i).cleanname.c_str()
      << " " << "MOL" << " H" << setw(4) << 0 << " " << "   "
      << setw(8) << sstr.atoms.at(i).cpos(1)
      << setw(8) << sstr.atoms.at(i).cpos(2)
      << setw(8) << sstr.atoms.at(i).cpos(3)
      << "  1.00" << "  0.00"
      << setw(12) << sstr.atoms.at(i).cleanname.c_str() << "0" << endl;
  }
  oss << "TER" << endl;
}

// ***************************************************************************
// platon2print
// ***************************************************************************
// This funtion prints out structural data in a format for PLATON
// Stefano Curtarolo
void platon2print(xstructure str,bool P_EQUAL,bool P_EXACT,double P_ang,double P_d1,double P_d2,double P_d3, ostream& oss) {
  oss << str.platon2print(P_EQUAL,P_EXACT,P_ang,P_d1,P_d2,P_d3);
}

// **************************************************************************
//  Function PrintRBAnal
// **************************************************************************
// This function prints out analysis of the rubber band run.
// Dane Morgan - Stefano Curtarolo
void PrintRBAnal(const int& nim, const string& path_flag, ostream& oss) {
  int DEFAULT_NPT_FOR_SPLINE=50;
  oss.setf(std::ios::left, std::ios::adjustfield);
  oss.setf(std::ios::fixed, std::ios::floatfield);
  int w=10;
  int p=4;

  // Get data
  vector<string> dirvec=pflow::GetRBDir(nim);
  vector<double> enervec=pflow::GetRBEner(nim);
  vector<xstructure> strvec=pflow::GetRBStruct(nim);
  vector<double> dcum(nim+2);
  vector<double> d_from_init(nim+2);
  vector<double> d_from_final(nim+2);

  // Print data
  oss << "********* Analysis of rubber band run *********" << endl;

  // Exact flag
  string path_flag_temp="E";
  dcum=pflow::GetRBDistCum(strvec,path_flag_temp);
  d_from_init=pflow::GetRBDistFromStrI(strvec,strvec[0],path_flag_temp);
  d_from_final=pflow::GetRBDistFromStrI(strvec,strvec[nim+1],path_flag_temp);
  oss << "** Here we calculate with distance method (eE=exact, nN=nearest): "
    << path_flag_temp << " **" << endl;
  oss << "Dir    Dist(Cum)   Dist(00)    Dist(END)   Ener" <<endl;  
  for(int im=0;im<nim+2;im++) {
    oss << setprecision(p) << setw(5) << dirvec[im]
      << "  " << setw(w) << dcum[im]
      << "  " << setw(w) << d_from_init[im]
      << "  " << setw(w) << d_from_final[im]
      << "  " << setw(w) << enervec[im] << endl;
  }
  // Nearest flag
  path_flag_temp="N";
  dcum=pflow::GetRBDistCum(strvec,path_flag_temp);
  d_from_init=pflow::GetRBDistFromStrI(strvec,strvec[0],path_flag_temp);
  d_from_final=pflow::GetRBDistFromStrI(strvec,strvec[nim+1],path_flag_temp);
  oss << "** Here we calculate with distance method (eE=exact, nN=nearest): "
    << path_flag_temp << " **" << endl;
  oss << "Dir    Dist(Cum)   Dist(00)    Dist(END)   Ener" <<endl;  
  for(int im=0;im<nim+2;im++) {
    oss << setprecision(p) << setw(5) << dirvec[im]
      << "  " << setw(w) << dcum[im]
      << "  " << setw(w) << d_from_init[im]
      << "  " << setw(w) << d_from_final[im]
      << "  " << setw(w) << enervec[im] << endl;
  }

  oss << endl;
  oss << "*******  SPLINE FIT *******" << endl;
  oss << "Distance used in spline (eE=exact, nN=nearest): "
    << path_flag << endl;  
  oss << endl;
  pflow::PrintSpline(dcum,enervec,DEFAULT_NPT_FOR_SPLINE,oss);
}  

// ***************************************************************************
// Function PrintRBPoscarDisp
// ***************************************************************************
// This function prints out the displacement info for 2 POSCAR.
// Dane Morgan - Stefano Curtarolo
void PrintRBPoscarDisp(const xstructure& diffstr, double& totdist,
    aurostd::matrix<double>& cm, const string& path_flag,  //CO20200404 pflow::matrix()->aurostd::matrix()
    ostream& oss) {
  oss << "*******  BEGIN:  Structure 2 - Structure 1  *******" << endl;
  oss << diffstr;
  oss << "*******  END:  Structure 2 - Structure 1  *******" << endl;
  oss << endl;
  oss << "Center of mass of POSCARs: "<< endl;
  oss << "  POSCAR1  ";
  pflow::Vout(cm[0],oss);
  oss << "  POSCAR2  ";
  pflow::Vout(cm[1],oss);
  oss << "Total Distance: " << totdist << endl;
  oss << "Method (eE=exact, nN=nearest): " << path_flag << endl;
}

// **************************************************************************
// PrintRDF
// **************************************************************************
string RDF2string(const xstructure& xstr,const double rmax,const int nbins,const xmatrix<double>& rdf_all){ //CO20220627 - new procedure
  int padding=8,padding_extra=5;  //4 because of 10,100,1000,neg below
  string eset="";
  vector<vector<int> > vvitypes=pflow::getvvitypesRDF(xstr);
  double rad=0.0;
  double drad=rmax/(double)nbins; //sphere shell volume=4*pi*r^2*dr
  uint k=0;
  int ibin=0,itype=0;
  stringstream oss;
  oss << aurostd::PaddedPOST("#bin",padding) << " ";
  oss << aurostd::PaddedPRE("rad",padding+padding_extra) << " ";
  for(k=0;k<vvitypes.size();k++){
    eset=xstr.species[vvitypes[k][0]]+"-"+xstr.species[vvitypes[k][1]];
    oss << aurostd::PaddedPRE(eset,padding+padding_extra) << " ";
  }
  oss << aurostd::PaddedPRE("total",padding+padding_extra) << " ";
  oss << endl;
  for(ibin=rdf_all.lrows;ibin<=rdf_all.urows;ibin++){
    oss << aurostd::PaddedPOST(aurostd::utype2string(ibin,padding-3,FIXED_STREAM),padding) << " ";
    rad=drad*(ibin-rdf_all.lrows);
    oss << (rad<10?" ":"");
    oss << (rad<100?" ":"");
    oss << (rad<1000?" ":"");
    oss << (rad<10000?" ":"");
    oss << (std::signbit(rad)?"":" ");
    oss << aurostd::PaddedPOST(aurostd::utype2string(rad,padding-2,FIXED_STREAM),padding) << " ";
    for(itype=rdf_all.lcols;itype<=rdf_all.ucols;itype++){
      oss << (rdf_all[ibin][itype]<10?" ":"");
      oss << (rdf_all[ibin][itype]<100?" ":"");
      oss << (rdf_all[ibin][itype]<1000?" ":"");
      oss << (rdf_all[ibin][itype]<10000?" ":"");
      oss << (std::signbit(rdf_all[ibin][itype])?"":" ");
      oss << aurostd::PaddedPOST(aurostd::utype2string(rdf_all[ibin][itype],padding-2,FIXED_STREAM),padding) << " ";
    }
    oss << endl;
  }
  return oss.str();
}
void PrintRDF(const xstructure& xstr,const double rmax,const int nbins,const xmatrix<double>& rdf_all,ostream& oss){ //CO20220627 - new procedure
  oss << RDF2string(xstr,rmax,nbins,rdf_all);
}
// This function prints out the RDF information.
// Dane Morgan - Stefano Curtarolo
void PrintRDF_20220101(const xstructure& str, const double& rmax,
    const int& nbins,
    const int& smooth_width,
    const aurostd::matrix<double>& rdf_all,  //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double>& rdfsh_all,  //CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double>& rdfsh_loc,  //CO20200404 pflow::matrix()->aurostd::matrix()
    ostream& oss) {
  int natoms=str.atoms.size();
  //  vector<int> types=str.GetTypes();
  //  vector<int> num_each_type=str.GetNumEachType();
  int ntypes=str.num_each_type.size();
  oss.setf(std::ios::left, std::ios::adjustfield);
  oss.setf(std::ios::fixed, std::ios::floatfield);
  int w1=5;
  int _rdi=rdfsh_loc.size();
  int _rdj=nbins;
  aurostd::matrix<double> rdfsh_bin(_rdi,_rdj);pflow::VVset(rdfsh_bin,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()

  // Put rdfsh_all into bins so you can write out in same
  // way as rdf_all
  for(uint i=0;i<rdfsh_loc.size();i++) {
    for(uint j=0;j<rdfsh_loc[i].size();j++) {
      int ib=int((rdfsh_loc[i][j]/rmax)*nbins); // first bin is [0,rmax/nbins).
      rdfsh_bin[i][ib]+=rdfsh_all[i][j];
    }
  }

  // RDF
  oss <<"# Radial Distribution Function (RDF)"<< endl;
  oss << "# rmax= " << rmax << "     nbins = " << nbins << "     sigma= " << smooth_width << endl;
  oss << "# Format: bin  rad  type1:rdf shells  type2:rdf shells ... " << endl;
  double bw=rmax/(double)nbins;
  for(int ia=0;ia<natoms;ia++) {
    oss << "# atom / type: " << ia+1 << " / " << str.atoms.at(ia).type+1 << endl;  // CONVASP_MODE
    oss << "# bin    rad ";
    oss << "       type" << 1;
    for(int it=1;it<ntypes;it++) {
      oss << "          type" << it+1;
    }
    oss <<   "          typeAll";
    oss << endl;
    for(int ib=0;ib<nbins;ib++) {
      oss << "  " << setw(w1) << ib+1;
      oss << "  " << ib*bw;
      oss << "  ";
      for(int it=0;it<ntypes;it++) {
        oss <<" " << setw(w1) << rdf_all[(ntypes+1)*ia+it][ib];
        oss <<" " << setw(w1) << Nint(rdfsh_bin[(ntypes+1)*ia+it][ib]);
      }
      oss <<" " << setw(w1) << rdf_all[(ntypes+1)*ia+ntypes][ib];
      oss <<" " << setw(w1) << Nint(rdfsh_bin[(ntypes+1)*ia+ntypes][ib]);
      oss << endl;
    }

    // RDF SHELLS
    // Header line
    oss << "# atom / type: " << ia+1 << " / " << str.atoms.at(ia).type+1 << endl;  // CONVASP_MODE
    oss << "# NN   ";
    for(int it=0;it<ntypes;it++) {
      oss << "  type" << it+1;
      oss << ":rad/size";
    }
    oss <<   "  typeAll";
    oss << ":rad/size";
    oss << endl;
    // Get nshmx for this atom
    int nshmx=0;
    for(int i=0;i<ntypes+1;i++) {
      int temp=rdfsh_all[(ntypes+1)*ia+i].size();
      if(temp>nshmx) nshmx=temp;
    }
    // Output shell data
    for(int is=0;is<nshmx;is++) {
      oss << "  " << setw(w1) << is+1;
      for(int it=0;it<ntypes+1;it++) {
        int id=(ntypes+1)*ia+it;
        if(is<(int) rdfsh_all[id].size()) {
          oss << "  " << setw(w1+2) << setprecision(4) << rdfsh_loc[id][is];
          oss << "  " << setw(w1) << Nint(rdfsh_all[id][is]);
        }
        else {
          oss << "  " << setw(w1) << "NA";
          oss << "  " << setw(w1) << "NA";
        }
      } // for it
      oss << endl;
    } // for is
  } // for ia

}

// **************************************************************************
// PrintRDFCmp
// **************************************************************************
// This function prints out the RDF shells comparison information.
// Dane Morgan - Stefano Curtarolo
void PrintRDFCmp(const xstructure& str_A, const xstructure& str_B,
    const double& rmax, const int nbins,
    const double& smooth_width,const int nsh,
    const aurostd::matrix<double>& rdfsh_all_A,
    const aurostd::matrix<double>& rdfsh_all_B,
    const vector<int>& best_match,
    const aurostd::matrix<double>& rms_mat, ostream& oss) {  //CO20200404 pflow::matrix()->aurostd::matrix()
  int nat=str_A.atoms.size();   if(str_B.atoms.size()) {;} // phony to keep str_B busy
  std::deque<int> num_each_type=str_A.num_each_type;
  int nt=num_each_type.size();
  int w1=4;
  oss.setf(std::ios::left, std::ios::adjustfield);
  oss.setf(std::ios::fixed, std::ios::floatfield);

  // Info on RDF calculation parameters.
  oss <<"# Comparison of structures based on radial distribution functions (RDFs)"<< endl;
  oss << "# rmax= " << rmax << "     nbins = " << nbins << endl;
  oss << "# sigma= " << smooth_width << "     nsh = " << nsh << endl;
  oss << endl;
  oss << "# The matrix of rms error between atoms in str A and str B" << endl;

  // The RMS matrix.

  // header line
  oss << setw(w1+1) << " ";
  for(int i=0;i<nat;i++) {
    oss << " " << setw(w1) << i+1;
  }
  oss << endl;

  // The RMS matrix as a row/column matrix.
  for(int i=0;i<nat;i++) {
    oss << setw(w1) << i+1 << " ";
    for(int j=0;j<nat;j++) {
      oss << " " << setprecision(2) << rms_mat[i][j];
    }
    oss << endl;
  }

  // The fit data for each atom's best fit.
  double rms_tot=0;
  oss << endl;

  // Basic atoms and rms.
  for(int iA=0;iA<nat;iA++) {
    int iB=best_match[iA];
    rms_tot+=rms_mat[iA][iB];
    oss << "A atom / Best matching B atom / RMS: "
      << iA+1 << " " << iB+1 << " " << rms_mat[iA][iB] << endl;

    // Detailed data on the A and B shells
    // Header line
    oss << "# CN ";
    oss << "   type" << 1;
    for(int it=1;it<nt;it++) {
      oss << "             type" << it+1;
    }
    oss <<   "             typeAll";
    oss << endl;

    // Get nshmx for this atom
    int nshmx=0;
    int temp;
    for(int i=0;i<nt+1;i++) {
      temp=rdfsh_all_A[(nt+1)*iA+i].size();
      if(temp>nshmx) nshmx=temp;
      temp=rdfsh_all_B[(nt+1)*iB+i].size();
      if(temp>nshmx) nshmx=temp;
    }

    // Output shell data
    for(int is=0;is<nshmx;is++) {
      oss << "  " << setw(w1) << is+1;
      for(int it=0;it<nt+1;it++) {
        int idA=(nt+1)*iA+it;
        int idB=(nt+1)*iB+it;
        if(is<(int) rdfsh_all_A[idA].size()) {
          oss << " " << setw(w1) << Nint(rdfsh_all_A[idA][is]);
        }
        else {
          oss << " " << setw(w1) << "NA";
        }
        if(is<(int) rdfsh_all_B[idB].size()) {
          oss << " " << setw(w1) << Nint(rdfsh_all_B[idB][is]);
        }
        else {
          oss << " " << setw(w1) << "NA";
        }
        if(is<(int) rdfsh_all_A[idA].size() && is<(int) rdfsh_all_B[idB].size()) {
          oss << " " << setw(w1) << Nint(rdfsh_all_A[idA][is]-rdfsh_all_B[idB][is]);
        }
        else {
          oss << " " << setw(w1) << "NA";
        }
        oss << "   ";
      } // for it
      oss << endl;
    } // for is

  } // for iA

  // Total RMS
  oss << "RMS TOTAL: " << rms_tot/nat;
}

// ***************************************************************************
// PrintRSM
// ***************************************************************************
// This funtion prints out structural data in a Rasmol (RSM) format with unit cell wireframe
// Wahyu Setyawan
void PrintRSM(const xstructure& str, ostream& oss) {
  string atName[]={"H    ","He   ","Li   ","Be   ","B    ","C    ","N    ","O    ","F    ","Ne   ",
    "Na   ","Mg   ","Al   ","Si   ","P    ","S    ","Cl   ","Ar   ","K    ","Ca   ",
    "Sc   ","Ti   ","V    ","Cr   ","Mn   ","Fe   ","Co   ","Ni   ","Cu   ","Zn   ",
    "Ga   ","Ge   ","As   ","Se   ","Br   ","Kr   ","Rb   ","Sr   ","Y    ","Zr   ",
    "Nb   ","Mo   ","Tc   ","Ru   ","Rh   ","Pd   ","Ag   ","Cd   ","In   ","Sn   ",
    "Sb   ","Te   ","I    ","Xe   ","Cs   ","Ba   ","La   ","Ce   ","Pr   ","Nd   ",
    "Pm   ","Sm   ","Eu   ","Gd   ","Tb   ","Dy   ","Ho   ","Er   ","Tm   ","Yb   ",
    "Lu   ","Hf   ","Ta   ","W    ","Re   ","Os   ","Ir   ","Pt   ","Au   ","Hg   ",
    "Tl   ","Pb   ","Bi   ","Po   ","At   ","Rn   ","Fr   ","Ra   ","Ac   ","Th   ",
    "Pa   ","U    ","Np   ","Pu   ","Am   ","Cm   ","Bk   ","Cf   ","Es   ","Fm   "};
  oss<<"#Title"<<endl
    <<"#23456#8901##4567#9012#456#890#2345678#0123456#8901234#67890#23456789"<<endl;
  int num_atoms;
  num_atoms=0;
  for(uint i=0;i<str.num_each_type.size();i++)
    num_atoms+=str.num_each_type.at(i);

  oss.setf(std::ios::fixed);
  for(int i=0;i<num_atoms;i++) {
    oss<<"ATOM  "
      <<setw(5)<<i+9
      <<" "
      <<atName[-1+str.atoms.at(i).type]
      <<"UNK  "
      <<setw(4)<<str.atoms.at(i).type
      <<"    "
      <<setprecision(3)<<setw(8)<<str.scale*str.atoms.at(i).cpos(1)
      <<setprecision(3)<<setw(8)<<str.scale*str.atoms.at(i).cpos(2)
      <<setprecision(3)<<setw(8)<<str.scale*str.atoms.at(i).cpos(3)<<endl;
    //        fprintf(fid,'%6s%5d%1s%5s%5s %3d%4s%8.4f%8.4f%8.4f\n', ...^M
    //            'ATOM  ',ptr+8,' ',atName((label(i)-1)*5+1:label(i)*5),'UNK  ',label(i),'    ',x,y,z);^M
  }

  double latt1[3]={str.scale*str.lattice(1,1),str.scale*str.lattice(1,2),str.scale*str.lattice(1,3)};
  double latt2[3]={str.scale*str.lattice(2,1),str.scale*str.lattice(2,2),str.scale*str.lattice(2,3)};
  double latt3[3]={str.scale*str.lattice(3,1),str.scale*str.lattice(3,2),str.scale*str.lattice(3,3)};

  double cellframe[8][3]={{0,0,0},
    {latt1[0],latt1[1],latt1[2]}, //latt(1,:)
    {latt2[0],latt2[1],latt2[2]}, //latt(2,:)
    {latt3[0],latt3[1],latt3[2]}, //latt(3,:)
    {latt1[0]+latt2[0],latt1[1]+latt2[1],latt1[2]+latt2[2]}, //latt(1,:)+latt(2,:);
    {latt1[0]+latt3[0],latt1[1]+latt3[1],latt1[2]+latt3[2]}, //latt(1,:)+latt(3,:);
    {latt3[0]+latt2[0],latt3[1]+latt2[1],latt3[2]+latt2[2]}, //latt(3,:)+latt(2,:);
    {latt1[0]+latt2[0]+latt3[0],latt1[1]+latt2[1]+latt3[1],latt1[2]+latt2[2]+latt3[2]}};

  for(int i=0;i<8;i++) {
    oss<<"ATOM  "
      <<setw(5)<<i+1
      <<" "
      <<"X    "
      <<"UNK  "
      <<setw(4)<<0
      <<"    "
      <<setprecision(3)<<setw(8)<<cellframe[i][0]
      <<setprecision(3)<<setw(8)<<cellframe[i][1]
      <<setprecision(3)<<setw(8)<<cellframe[i][2]<<endl;
    //        fprintf(fid,'%6s%5d%1s%5s%5s %3d%4s%8.4f%8.4f%8.4f\n', ...^M
    //            'ATOM  ',i,' ','X    ','UNK  ',0,'    ',x,y,z);^M
  }

  int lineframe[12][2]={{1,2},  {1,3},  {1,4},  {8,7},  {8,6},  {8,5},  {2,5},  {2,6},  {3,5},  {3,7},  {4,6},  {4,7}};
  for(int i=0;i<12;i++) {
    oss<<"CONECT "
      <<setw(4)<<lineframe[i][0]
      <<" "
      <<setw(4)<<lineframe[i][1]<<endl;
    //    fprintf(fid,'%6s %4d %4d\n','CONECT',lineframe(i,1),lineframe(i,2));^M
  }

  oss<<"END"<<endl
    <<"#!rasmol -rsm file"<<endl
    <<"load pdb inline"<<endl
    <<"spacefill 150"<<endl
    <<"select 0"<<endl
    <<"spacefill 10"<<endl
    <<"exit"<<endl;

  //%#23456#8901##4567#9012#456#890#2345678#0123456#8901234#67890#23456789^M
  //%COLO      1 Al   UNK  0001       0.500   0.500   0.000   0.2^M
  //%ATOM      1 Al   UNK  0001       0.000   0.000   0.000^M
  //%ATOM      2 O    UNK  0002       2.000   0.000   0.000^M
  //%CONECTxxxxxyyyyy^M
}

// ***************************************************************************
// pflow::PrintSGData()
// ***************************************************************************
// This funtion prints out space group data.
// David Hicks (DX) 
// DX20210301 - cleaned, use filetype, use JSON writer, and other minor mods
// Note: xstructure is modified in this function (so do not use const)
namespace pflow {
  string PrintSGData(xstructure& xstr,
      filetype ftype,
      bool standalone,
      bool already_calculated,
      double sym_eps,
      bool no_scan,
      int setting,
      bool suppress_Wyckoff) {
    aurostd::xoption vpflow;
    return PrintSGData(xstr, vpflow, ftype, standalone, already_calculated, sym_eps, no_scan, setting, suppress_Wyckoff);
  }
}

namespace pflow {
  string PrintSGData(xstructure& xstr,
      aurostd::xoption& vpflow,
      filetype ftype,
      bool standalone,
      bool already_calculated,
      double sym_eps,
      bool no_scan,
      int setting,
      bool suppress_Wyckoff) {

    bool LDEBUG=(FALSE || XHOST.DEBUG);

    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " BEGIN" << endl; }

    // ---------------------------------------------------------------------------
    // calculate the space group symmetry if not already calculated
    xstructure str_sg(xstr);
    if(!already_calculated){
      str_sg.SpaceGroup_ITC(sym_eps, -1, setting, no_scan);
    }

    // ---------------------------------------------------------------------------
    // get symmetry labels ("try" allows the function to continue if no_scan is true)
    string space_group_HM="", space_group_Hall="", space_group_Schoenflies="", class_Laue="";
    try {
      space_group_HM = GetSpaceGroupName(str_sg.space_group_ITC, str_sg.directory);
      space_group_Hall = GetSpaceGroupHall(str_sg.space_group_ITC, setting, str_sg.directory);
      space_group_Schoenflies = GetSpaceGroupSchoenflies(str_sg.space_group_ITC, str_sg.directory);
      class_Laue = GetLaueLabel(str_sg.point_group_ITC);
    }
    catch(aurostd::xerror& re){
      if(!xstr.sym_eps_no_scan){ throw; } //DX20210430 - only throw if symmetry scan is implemented
    }

    stringstream ss_output;
    ss_output.setf(std::ios::fixed,std::ios::floatfield);
    // ---------------------------------------------------------------------------
    // text output
    // ME20210402 - print both for web
    if((ftype == txt_ft) || (XHOST.vflag_control.flag("WWW") && standalone)) { //DX20210521 - add standalone
      ss_output << "SPACE GROUP OF THE CRYSTAL" << endl;
      ss_output << " Space group number                           = " << str_sg.space_group_ITC << endl;
      ss_output << " Space group label (Hermann Mauguin)          = " << space_group_HM << endl;
      ss_output << " Space group label (Hall)                     = " << space_group_Hall << endl;
      ss_output << " Space group label (Schoenflies)              = " << space_group_Schoenflies << endl;
      ss_output << " Laue class                                   = " << class_Laue << endl;
      ss_output << " Crystal class                                = " << str_sg.point_group_ITC << endl;
      ss_output << "ITC REPRESENTATION OF THE CRYSTAL" << endl;
      ss_output << " Setting                                      = " << str_sg.setting_ITC << endl;
      ss_output << " Origin                                       = " << roundoff(str_sg.origin_ITC,1e-8) << endl;

      // if printing all Wyckoff info
      if(!suppress_Wyckoff){
        ss_output << " General Wyckoff position" << endl;
        for(uint i=0;i<str_sg.general_position_ITC.size();i++){
          ss_output << "  " << i+1 << " " << str_sg.general_position_ITC[i] << endl;
        }
        ss_output << " Representative Wyckoff positions" << endl;
        for(uint i=0;i<str_sg.wyccar_ITC.size();i++){
          if(i>4 && i!=str_sg.wyccar_ITC.size()-1){ // skip title, scale, lattice parameters, number of atoms, and coordinate type, and last newline
            ss_output << str_sg.wyccar_ITC[i] << endl;
          }
        }
        ss_output << "WYCCAR" << endl;
        // ---------------------------------------------------------------------------
        // expanded form of WYCCAR (i.e., poscar, not just the wyccar) //DX20210708
        xstructure str_expanded;
        str_expanded.title = str_sg.title; // transfer title over
        str_expanded.sym_eps = str_sg.sym_eps; // transfer symmetry eps over
        // ---------------------------------------------------------------------------
        // add representative atoms only //DX20220921
        _atom atom_tmp;
        for(uint i=0;i<str_sg.wyckoff_sites_ITC.size();i++){
          atom_tmp.fpos = str_sg.wyckoff_sites_ITC[i].coord;
          atom_tmp.type = str_sg.wyckoff_sites_ITC[i].index;
          atom_tmp.name = str_sg.wyckoff_sites_ITC[i].type;
          atom_tmp.partial_occupation_flag = aurostd::isdifferent(str_sg.wyckoff_sites_ITC[i].site_occupation,1.0);
          atom_tmp.partial_occupation_value = str_sg.wyckoff_sites_ITC[i].site_occupation;
          str_expanded.AddAtom(atom_tmp,false);  //CO20230319 - add by species
        }
        xvector<double> data = Getabc_angles(str_sg.standard_lattice_ITC,DEGREES);
        str_expanded = WyckoffPOSITIONS(str_sg.space_group_ITC, str_sg.setting_ITC, str_expanded); //DX20220921 - changed str_sg to str_expanded in last argument
        str_expanded.lattice=GetClat(data(1), data(2), data(3), data(4), data(5), data(6)); //DX20221128 - need to use GetClat to get standard representation of the ITC unit cell (rather than using std_lattice_ITC(), which could be rotated)
        str_expanded.ReScale(1.0);
        str_expanded.neg_scale=FALSE;
        str_expanded.iomode = IOVASP_POSCAR;
        ss_output << str_expanded;
        //DX20210621 [OBSOLETE] convert vector<string> of WYCCAR to xstructure
        //DX20210621 [OBSOLETE] stringstream wss;
        //DX20210621 [OBSOLETE] for(uint i=0;i<str_sg.wyccar_ITC.size();i++){
        //DX20210621 [OBSOLETE]   wss << str_sg.wyccar_ITC[i] << endl;
        //DX20210621 [OBSOLETE] }
        //DX20210621 [OBSOLETE] if(str_sg.wyccar_ITC.size()!=0){ //DX20180526 - skip if fails
        //DX20210621 [OBSOLETE]   xstructure xstr_wyccar(wss,IOVASP_WYCKCAR);
        //DX20210621 [OBSOLETE]   ss_output << xstr_wyccar;
        //DX20210621 [OBSOLETE] }
      }
      // ME20210402 - Convert to array of strings for web
      if (XHOST.vflag_control.flag("WWW") && standalone) { //DX20210521 - add standalone
        vector<string> voutput;
        aurostd::stream2vectorstring(ss_output, voutput);
        ss_output.clear();
        ss_output.str("");
        ss_output << "{\"txt\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(voutput, "\"", "\""), ",") << "],";
      }
    }
    // ---------------------------------------------------------------------------
    // json output
    // ME20210402 - print both for web
    if((ftype == json_ft) || (XHOST.vflag_control.flag("WWW") && standalone)){ //DX20210521 - add standalone

      aurostd::JSONwriter json;
      bool roff = true;

      // space group number
      if(str_sg.space_group_ITC > 0 && str_sg.space_group_ITC < 231){
        json.addNumber("space_group_number", str_sg.space_group_ITC);
      } else if (PRINT_NULL_JSON){
        json.addNull("space_group_number");
      }

      // space group label (HM)
      if(!space_group_HM.empty()){
        json.addString("space_group_Hermann_Mauguin", space_group_HM);
      } else if (PRINT_NULL_JSON){
        json.addNull("space_group_Hermann_Mauguin");
      }

      // space group label (Hall)
      if(!space_group_Hall.empty()){
        json.addString("space_group_Hall", space_group_Hall);
      } else if (PRINT_NULL_JSON){
        json.addNull("space_group_Hall");
      }

      // space group label (Schoenflies)
      if(!space_group_Schoenflies.empty()){
        json.addString("space_group_Schoenflies", space_group_Schoenflies);
      } else if (PRINT_NULL_JSON){
        json.addNull("space_group_Schoenflies");
      }

      // Laue
      if(!class_Laue.empty()){
        json.addString("Laue", class_Laue);
      } else if (PRINT_NULL_JSON){
        json.addNull("Laue");
      }

      // crystal class
      if(!str_sg.point_group_ITC.empty()){
        json.addString("crystal_class", str_sg.point_group_ITC);
      } else if (PRINT_NULL_JSON){
        json.addNull("crystal_class");
      }

      // ITC setting
      if(str_sg.setting_ITC == 1 || str_sg.setting_ITC == 2){
        json.addNumber("setting_ITC", str_sg.setting_ITC);
      } else if (PRINT_NULL_JSON){
        json.addNull("setting_ITC");
      }

      // ITC origin
      if(str_sg.origin_ITC.rows != 0){
        json.addVector("origin_ITC", str_sg.origin_ITC, _AFLOWLIB_DATA_GEOMETRY_PREC_, roff);
      } else if (PRINT_NULL_JSON){
        json.addNull("origin_ITC");
      }

      if(!suppress_Wyckoff){
        // general Wyckoff position
        if(!str_sg.general_position_ITC.empty()){
          vector<vector<string> > positions;
          vector<string> tokens_positions;
          for(uint i=0;i<str_sg.general_position_ITC.size();i++){
            tokens_positions.clear();
            aurostd::string2tokens(str_sg.general_position_ITC[i],tokens_positions,",");
            positions.push_back(tokens_positions);
          }
          json.addMatrix("general_position_ITC",positions);
        } else if (PRINT_NULL_JSON){
          json.addNull("general_position_ITC");
        }

        // representative Wyckoff positions
        if(str_sg.wyccar_ITC.size()){

          aurostd::JSONwriter Wyckoff_json;
          vector<aurostd::JSONwriter> vset_json;

          xvector<double> position;
          for(uint i=0;i<str_sg.wyccar_ITC.size();i++){
            if(i>4 && i!=str_sg.wyccar_ITC.size()-1){ //Skip title, scale, lattice parameters, number of atoms, and coordinate type, and last newline
              Wyckoff_json.clear();
              position.clear();

              vector<string> tokens;
              aurostd::string2tokens(str_sg.wyccar_ITC[i],tokens," ");
              if(tokens.size() == 7){
                position(1) = aurostd::string2utype<double>(tokens[0]);
                position(2) = aurostd::string2utype<double>(tokens[1]);
                position(3) = aurostd::string2utype<double>(tokens[2]);
                Wyckoff_json.addVector("position", position, _AFLOWLIB_DATA_GEOMETRY_PREC_, roff);
                Wyckoff_json.addString("name", tokens[3]);
                Wyckoff_json.addNumber("mulitiplicity", tokens[4]);
                Wyckoff_json.addString("Wyckoff_letter", tokens[5]);
                Wyckoff_json.addString("site_symmetry", tokens[6]);
                vset_json.push_back(Wyckoff_json);
              }
            }
          }

          json.addVector("Wyckoff_positions", vset_json);
        } else if (PRINT_NULL_JSON){
          json.addNull("Wyckoff_positions");
        }

        // WYCCAR
        if(!str_sg.wyccar_ITC.empty()){
          // convert vector<string> of WYCCAR to xstructure
          stringstream wss;
          for(uint i=0;i<str_sg.wyccar_ITC.size();i++){
            wss << str_sg.wyccar_ITC[i] << endl;
          }
          xstructure xstr_wyccar(wss,IOVASP_WYCKCAR);
          json.addRawJSON("wyccar", xstructure2json(xstr_wyccar));
        } else if (PRINT_NULL_JSON){
          json.addNull("wyccar");
        }
      }
      // ME20210402 - added web output
      if (XHOST.vflag_control.flag("WWW") && standalone) ss_output << "\"json\":" << json.toString(standalone) << "}" << std::endl; //DX20210521 - add standalone
      else ss_output << json.toString(standalone);
      if(standalone) { ss_output << endl; }
    }

    //DX20180822 - put attributes into vpflow - START
    if(vpflow.flag("SGDATA::CALCULATED")){ //DX20180823 - if calculated, remove
      vpflow.pop_attached("SGDATA::SPACE_GROUP_NUMBER");
      vpflow.pop_attached("SGDATA::SPACE_GROUP_HERMANN_MAUGUIN");
      vpflow.pop_attached("SGDATA::SPACE_GROUP_HALL");
      vpflow.pop_attached("SGDATA::SPACE_GROUP_SCHOENFLIES");
      vpflow.pop_attached("SGDATA::LAUE");
      vpflow.pop_attached("SGDATA::CRYSTAL_CLASS");
      vpflow.pop_attached("SGDATA::SETTING_ITC");
      vpflow.pop_attached("SGDATA::ORIGIN_ITC");
      vpflow.pop_attached("SGDATA::GENERAL_POSITION_ITC");
      vpflow.pop_attached("SGDATA::WYCKOFF_LETTERS");
      vpflow.pop_attached("SGDATA::WYCKOFF_MULTIPLICITIES");
      vpflow.pop_attached("SGDATA::WYCKOFF_SITE_SYMMETRIES");
      vpflow.pop_attached("SGDATA::WYCKOFF_POSITIONS");
      vpflow.pop_attached("SGDATA::WYCCAR");
    }
    //DX20180823 - now populate xoption
    vpflow.flag("SGDATA::CALCULATED");
    vpflow.push_attached("SGDATA::SPACE_GROUP_NUMBER",aurostd::utype2string<uint>(str_sg.space_group_ITC));
    vpflow.push_attached("SGDATA::SPACE_GROUP_HERMANN_MAUGUIN",space_group_HM);
    vpflow.push_attached("SGDATA::SPACE_GROUP_HALL",space_group_Hall);
    vpflow.push_attached("SGDATA::SPACE_GROUP_SCHOENFLIES",space_group_Schoenflies);
    vpflow.push_attached("SGDATA::LAUE",class_Laue);
    vpflow.push_attached("SGDATA::CRYSTAL_CLASS",str_sg.point_group_ITC);
    vpflow.push_attached("SGDATA::SETTING_ITC",aurostd::utype2string<uint>(str_sg.setting_ITC));
    vpflow.push_attached("SGDATA::ORIGIN_ITC",aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(str_sg.origin_ITC,5,true),","));
    // turn vector of strings into single string
    stringstream general_position;
    general_position << "[";
    for(uint i=0;i<str_sg.general_position_ITC.size();i++){
      vector<string> tokens;
      aurostd::string2tokens(str_sg.general_position_ITC[i],tokens,",");
      general_position << "[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(tokens,"\""),",") << "]";
      if(i != str_sg.general_position_ITC.size()-1){
        general_position << ",";
      }
    }
    general_position << "]";
    vpflow.push_attached("SGDATA::GENERAL_POSITION_ITC",general_position.str());
    stringstream wyccar;
    stringstream Wyckoff_positions;
    for(uint i=0;i<str_sg.wyccar_ITC.size();i++){
      if(i>4 && i!=str_sg.wyccar_ITC.size()-1){ //Skip title, scale, lattice parameters, number of atoms, and coordinate type, and last newline
        Wyckoff_positions << str_sg.wyccar_ITC[i] << endl;
      }
      wyccar << str_sg.wyccar_ITC[i] << endl;
    }
    vpflow.push_attached("SGDATA::WYCKOFF_LETTERS",SYM::ExtractWyckoffLettersString(str_sg.wyccar_ITC));
    vpflow.push_attached("SGDATA::WYCKOFF_MULTIPLICITIES",SYM::ExtractWyckoffMultiplicitiesString(str_sg.wyccar_ITC));
    vpflow.push_attached("SGDATA::WYCKOFF_SITE_SYMMETRIES",SYM::ExtractWyckoffSiteSymmetriesString(str_sg.wyccar_ITC));
    vpflow.push_attached("SGDATA::WYCKOFF_POSITIONS",Wyckoff_positions.str());
    vpflow.push_attached("SGDATA::WYCCAR",wyccar.str());
    //DX20180822 - put attributes into vpflow - END
    return ss_output.str();
  }
}

// ***************************************************************************
// pflow::PrintWyckoffData()
// ***************************************************************************
// This funtion prints out Wyckoff position data.
// David Hicks (DX) 
// DX20210611 - cleaned, use filetype, use JSON writer, and other minor mods
// Note: xstructure is modified in this function (so do not use const)
namespace pflow {
  string PrintWyckoffData(xstructure& xstr,
      filetype ftype,
      bool standalone,
      bool already_calculated,
      double sym_eps,
      bool no_scan,
      int setting) {

    bool LDEBUG=(FALSE || XHOST.DEBUG);

    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " BEGIN" << endl; }

    // ---------------------------------------------------------------------------
    // calculate the space group symmetry if not already calculated
    xstructure str_sg(xstr);
    if(!already_calculated){
      str_sg.SpaceGroup_ITC(sym_eps, -1, setting, no_scan);
    }

    xvector<double> data = Getabc_angles(str_sg.standard_lattice_ITC,DEGREES);
    stringstream ss_output;
    if(ftype == txt_ft || (XHOST.vflag_control.flag("WWW") && standalone)){
      ss_output << " Space group number                           = " << str_sg.space_group_ITC << endl; 
      ss_output << " Setting                                      = " << str_sg.setting_ITC << endl;
      ss_output << " ITC cell parameters - a b c alpha beta gamma = ";
      ss_output.precision(10); ss_output << data(1) << " " << data(2) << " " << data(3) << " ";
      ss_output.precision(3); ss_output << data(4) << " " << data(5) << " " << data(6) << endl;
      ss_output << " Wyckoff positions" << endl;
      uint nWyckoff_sites = str_sg.wyckoff_sites_ITC.size();
      for(uint i=0;i<nWyckoff_sites;i++){
        ss_output << std::setprecision(14) << std::left << std::fixed
          << "   " << setw(5) << str_sg.wyckoff_sites_ITC[i].type         // set(5): species have max of two characters + pseudopoential (up to three characters)
          << " " << setw(3) << str_sg.wyckoff_sites_ITC[i].multiplicity   // set(3): Wyckoff positions have a max of three digits
          << " " << setw(5) << str_sg.wyckoff_sites_ITC[i].letter         // set(5): space group #47 has Wyckoff letter "alpha" (i.e, five letters)
          << " " << setw(8) << str_sg.wyckoff_sites_ITC[i].site_symmetry; // set(8): Wyckoff site symmetry generally lists three/four directions, each with one to two characters (e.g., inversion)
        double _coord = AUROSTD_MAX_DOUBLE;
        for(uint j=1;j<=3;j++) {
          _coord=aurostd::roundoff(str_sg.wyckoff_sites_ITC[i].coord(j),pow(10.0,-(double)14));
          if(abs(_coord)<10.0) ss_output << " ";
          if(!std::signbit(_coord)) ss_output << " ";
          ss_output << _coord << " ";
        }
        ss_output << endl;
      }
      // ME20210402 - Convert to array of strings for web
      if (XHOST.vflag_control.flag("WWW") && standalone) { //DX20210521 - add standalone
        vector<string> voutput;
        aurostd::stream2vectorstring(ss_output, voutput);
        ss_output.clear();
        ss_output.str("");
        ss_output << "{\"txt\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(voutput, "\"", "\""), ",") << "],";
      }
    }
    if(ftype == json_ft || (XHOST.vflag_control.flag("WWW") && standalone)){
      aurostd::JSONwriter json;
      bool roff = true;

      // space group number
      if(str_sg.space_group_ITC > 0 && str_sg.space_group_ITC < 231){
        json.addNumber("space_group_number", str_sg.space_group_ITC);
      } else if (PRINT_NULL_JSON){
        json.addNull("space_group_number");
      }

      // ITC setting
      if(str_sg.setting_ITC == 1 || str_sg.setting_ITC == 2){
        json.addNumber("setting_ITC", str_sg.setting_ITC);
      } else if (PRINT_NULL_JSON){
        json.addNull("setting_ITC");
      }
      // Real lattice parameters
      if(data.rows != 0){
        json.addVector("lattice_parameters_ITC",data,_AFLOWLIB_DATA_GEOMETRY_PREC_,roff);
      } else if (PRINT_NULL_JSON){
        json.addNull("lattice_parameters_ITC");
      }
      // representative Wyckoff positions
      if(str_sg.wyccar_ITC.size()){

        aurostd::JSONwriter Wyckoff_json;
        vector<aurostd::JSONwriter> vset_json;

        xvector<double> position;
        for(uint i=0;i<str_sg.wyccar_ITC.size();i++){
          if(i>4 && i!=str_sg.wyccar_ITC.size()-1){ //Skip title, scale, lattice parameters, number of atoms, and coordinate type, and last newline
            Wyckoff_json.clear();
            position.clear();

            vector<string> tokens;
            aurostd::string2tokens(str_sg.wyccar_ITC[i],tokens," ");
            if(tokens.size() == 7){
              position(1) = aurostd::string2utype<double>(tokens[0]);
              position(2) = aurostd::string2utype<double>(tokens[1]);
              position(3) = aurostd::string2utype<double>(tokens[2]);
              Wyckoff_json.addVector("position", position, _AFLOWLIB_DATA_GEOMETRY_PREC_, roff);
              Wyckoff_json.addString("name", tokens[3]);
              Wyckoff_json.addNumber("mulitiplicity", tokens[4]);
              Wyckoff_json.addString("Wyckoff_letter", tokens[5]);
              Wyckoff_json.addString("site_symmetry", tokens[6]);
              vset_json.push_back(Wyckoff_json);
            }
          }
        }

        json.addVector("Wyckoff_positions", vset_json);
      } else if (PRINT_NULL_JSON){
        json.addNull("Wyckoff_positions");
      }
      if (XHOST.vflag_control.flag("WWW") && standalone) ss_output << "\"json\":" << json.toString(standalone) << "}"; //DX20210521 - add standalone
      else ss_output << json.toString(standalone);
      if(standalone) { ss_output << endl; }
    }
    return ss_output.str(); 
  }
}

//DX20210301 [OBSOLETE] namespace pflow {
//DX20210301 [OBSOLETE]   bool PrintSGData(xstructure& str_sg, ostream& oss, bool standalone, const string& format,bool already_calculated) {
//DX20210301 [OBSOLETE]     double tolerance=SYM::defaultTolerance(str_sg);
//DX20210301 [OBSOLETE]     if(already_calculated){tolerance=str_sg.sym_eps;} //CO20171025
//DX20210301 [OBSOLETE]     bool no_scan=false;
//DX20210301 [OBSOLETE]     int setting=1;
//DX20210301 [OBSOLETE]     return PrintSGData(str_sg,tolerance,oss,no_scan,setting,standalone,format,already_calculated);
//DX20210301 [OBSOLETE]   }
//DX20210301 [OBSOLETE] }
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE] //DX20180822 - add xoption - START
//DX20210301 [OBSOLETE] namespace pflow {
//DX20210301 [OBSOLETE]   bool PrintSGData(xstructure& str_sg, double& tolerance, ostream& oss, bool no_scan, int sg_setting, bool standalone, const string& format,bool already_calculated) {
//DX20210301 [OBSOLETE]     aurostd::xoption vpflow;
//DX20210301 [OBSOLETE]     return PrintSGData(str_sg,tolerance,oss,vpflow,no_scan,sg_setting,standalone,format,already_calculated);
//DX20210301 [OBSOLETE]   }
//DX20210301 [OBSOLETE] }
//DX20210301 [OBSOLETE] //DX20180822 - add xoption - END
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE] namespace pflow {
//DX20210301 [OBSOLETE]   bool PrintSGData(xstructure& str_sg, double& tolerance, ostream& oss_final, aurostd::xoption& vpflow, bool no_scan, int setting, bool standalone, const string& format,bool already_calculated) { //DX20180226 - added & to tolerance
//DX20210301 [OBSOLETE]     bool LDEBUG=(FALSE || XHOST.DEBUG);
//DX20210301 [OBSOLETE]     if(LDEBUG) cerr << XPID << "pflow::PrintSGData: BEGIN" << endl;
//DX20210301 [OBSOLETE]     stringstream oss;
//DX20210301 [OBSOLETE]     // smode=="DATA" or "EDATA"
//DX20210301 [OBSOLETE]     // Print out structural data
//DX20210301 [OBSOLETE]     oss.setf(std::ios::fixed,std::ios::floatfield);
//DX20210301 [OBSOLETE]     //DX TEST str_sg.SpaceGroup_ITC(tolerance);
//DX20210301 [OBSOLETE]     if(!already_calculated){str_sg.SpaceGroup_ITC(tolerance, -1, setting, no_scan);} //CO20171025
//DX20210301 [OBSOLETE]     // FORMAT = TXT
//DX20210301 [OBSOLETE]     if(format=="txt"){
//DX20210301 [OBSOLETE]       oss << "SPACE GROUP OF THE CRYSTAL" << endl;
//DX20210301 [OBSOLETE]       oss << " Space group number                           = " << str_sg.space_group_ITC << endl;
//DX20210301 [OBSOLETE]       oss << " Space group label (Hermann Mauguin)          = " << GetSpaceGroupName(str_sg.space_group_ITC, str_sg.directory) << endl; //DX20180526 - add directory
//DX20210301 [OBSOLETE]       oss << " Space group label (Hall)                     = " << GetSpaceGroupHall(str_sg.space_group_ITC, setting, str_sg.directory) << endl; //DX20180526 - add directory
//DX20210301 [OBSOLETE]       oss << " Space group label (Schoenflies)              = " << GetSpaceGroupSchoenflies(str_sg.space_group_ITC, str_sg.directory) << endl; //DX20180526 - add directory
//DX20210301 [OBSOLETE]       oss << " Laue class                                   = " << GetLaueLabel(str_sg.point_group_ITC) << endl;
//DX20210301 [OBSOLETE]       oss << " Crystal class                                = " << str_sg.point_group_ITC << endl;
//DX20210301 [OBSOLETE]       oss << "ITC REPRESENTATION OF THE CRYSTAL" << endl;
//DX20210301 [OBSOLETE]       oss << " Setting                                      = " << str_sg.setting_ITC << endl;
//DX20210301 [OBSOLETE]       oss << " Origin                                       = " << roundoff(str_sg.origin_ITC,1e-8) << endl;
//DX20210301 [OBSOLETE]       oss << " General Wyckoff position" << endl;
//DX20210301 [OBSOLETE]       for(uint i=0;i<str_sg.general_position_ITC.size();i++){
//DX20210301 [OBSOLETE]         oss << "  " << i+1 << " " << str_sg.general_position_ITC[i] << endl;
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       oss << " Representative Wyckoff positions" << endl;
//DX20210301 [OBSOLETE]       for(uint i=0;i<str_sg.wyccar_ITC.size();i++){
//DX20210301 [OBSOLETE]         if(i>4 && i!=str_sg.wyccar_ITC.size()-1){ //Skip title, scale, lattice parameters, number of atoms, and coordinate type, and last newline
//DX20210301 [OBSOLETE]           oss << str_sg.wyccar_ITC[i] << endl;
//DX20210301 [OBSOLETE]         }
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       oss << "WYCCAR" << endl;
//DX20210301 [OBSOLETE]       // convert vector<string> of WYCCAR to xstructure
//DX20210301 [OBSOLETE]       stringstream wss;
//DX20210301 [OBSOLETE]       for(uint i=0;i<str_sg.wyccar_ITC.size();i++){
//DX20210301 [OBSOLETE]         wss << str_sg.wyccar_ITC[i] << endl;
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       if(str_sg.wyccar_ITC.size()!=0){ //DX20180526 - skip if fails
//DX20210301 [OBSOLETE]         xstructure xstr_wyccar(wss,IOVASP_WYCKCAR);
//DX20210301 [OBSOLETE]         oss << xstr_wyccar << endl;
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]     }
//DX20210301 [OBSOLETE]     // FORMAT = JSON
//DX20210301 [OBSOLETE]     else if(format=="json"){
//DX20210301 [OBSOLETE]       string eendl="";
//DX20210301 [OBSOLETE]       bool roff=true; //round off
//DX20210301 [OBSOLETE]       stringstream sscontent_json;
//DX20210301 [OBSOLETE]       vector<string> vcontent_json;
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]       // space group number
//DX20210301 [OBSOLETE]       if(str_sg.space_group_ITC){
//DX20210301 [OBSOLETE]         sscontent_json << "\"space_group_number\":\"" << str_sg.space_group_ITC << "\"" << eendl;
//DX20210301 [OBSOLETE]       } else {
//DX20210301 [OBSOLETE]         if(PRINT_NULL_JSON){ sscontent_json << "\"space_group_number\":null" << eendl;}
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]       // space group label (HM)
//DX20210301 [OBSOLETE]       string sg_HM = GetSpaceGroupName(str_sg.space_group_ITC, str_sg.directory); //DX20180526 - add directory
//DX20210301 [OBSOLETE]       if(sg_HM.size()){
//DX20210301 [OBSOLETE]         sscontent_json << "\"space_group_Hermann_Mauguin\":\"" << sg_HM << "\"" << eendl;
//DX20210301 [OBSOLETE]       } else {
//DX20210301 [OBSOLETE]         if(PRINT_NULL_JSON){ sscontent_json << "\"space_group_Hermann_Mauguin\":null" << eendl;}
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]       // space group label (Hall)
//DX20210301 [OBSOLETE]       string sg_Hall = GetSpaceGroupHall(str_sg.space_group_ITC, setting, str_sg.directory); //DX20180526 - add directory
//DX20210301 [OBSOLETE]       if(sg_Hall.size()){
//DX20210301 [OBSOLETE]         sscontent_json << "\"space_group_Hall\":\"" << sg_Hall << "\"" << eendl;
//DX20210301 [OBSOLETE]       } else {
//DX20210301 [OBSOLETE]         if(PRINT_NULL_JSON){ sscontent_json << "\"space_group_Hall\":null" << eendl;}
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]       // space group label (Schoenflies)
//DX20210301 [OBSOLETE]       string sg_Schoenflies = GetSpaceGroupSchoenflies(str_sg.space_group_ITC, str_sg.directory); //DX20180526 - add directory
//DX20210301 [OBSOLETE]       if(sg_Schoenflies.size()){
//DX20210301 [OBSOLETE]         sscontent_json << "\"space_group_Schoenflies\":\"" << sg_Schoenflies << "\"" << eendl;
//DX20210301 [OBSOLETE]       } else {
//DX20210301 [OBSOLETE]         if(PRINT_NULL_JSON){ sscontent_json << "\"space_group_Schoenflies\":null" << eendl;}
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]       // Laue
//DX20210301 [OBSOLETE]       string laue = GetLaueLabel(str_sg.point_group_ITC);
//DX20210301 [OBSOLETE]       if(laue.size()){
//DX20210301 [OBSOLETE]         sscontent_json << "\"Laue\":\"" << laue << "\"" << eendl;
//DX20210301 [OBSOLETE]       } else {
//DX20210301 [OBSOLETE]         if(PRINT_NULL_JSON){ sscontent_json << "\"Laue\":null" << eendl;}
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]       // crystal class
//DX20210301 [OBSOLETE]       if(str_sg.point_group_ITC.size()){
//DX20210301 [OBSOLETE]         sscontent_json << "\"crystal_class\":\"" << str_sg.point_group_ITC << "\"" << eendl;
//DX20210301 [OBSOLETE]       } else {
//DX20210301 [OBSOLETE]         if(PRINT_NULL_JSON){ sscontent_json << "\"crystal_class\":null" << eendl;}
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]       // ITC setting
//DX20210301 [OBSOLETE]       if(str_sg.setting_ITC){
//DX20210301 [OBSOLETE]         sscontent_json << "\"setting_ITC\":\"" << str_sg.setting_ITC << "\"" << eendl;
//DX20210301 [OBSOLETE]       } else {
//DX20210301 [OBSOLETE]         if(PRINT_NULL_JSON){ sscontent_json << "\"setting_ITC\":null" << eendl;}
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]       // ITC origin
//DX20210301 [OBSOLETE]       if(str_sg.origin_ITC.rows){
//DX20210301 [OBSOLETE]         sscontent_json << "\"origin_ITC\":[" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(str_sg.origin_ITC,5,roff),",") << "]" << eendl;
//DX20210301 [OBSOLETE]       } else {
//DX20210301 [OBSOLETE]         if(PRINT_NULL_JSON){ sscontent_json << "\"origin_ITC\":null" << eendl;}
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]       // general Wyckoff position
//DX20210301 [OBSOLETE]       if(str_sg.general_position_ITC.size()){
//DX20210301 [OBSOLETE]         sscontent_json << "\"general_position_ITC\":[";
//DX20210301 [OBSOLETE]         for(uint i=0;i<str_sg.general_position_ITC.size();i++){
//DX20210301 [OBSOLETE]           vector<string> tokens;
//DX20210301 [OBSOLETE]           aurostd::string2tokens(str_sg.general_position_ITC[i],tokens,",");
//DX20210301 [OBSOLETE]           sscontent_json << "[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(tokens,"\""),",") << "]" << eendl;
//DX20210301 [OBSOLETE]           if(i != str_sg.general_position_ITC.size()-1){
//DX20210301 [OBSOLETE]             sscontent_json << ",";
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]         }
//DX20210301 [OBSOLETE]         sscontent_json << "]";
//DX20210301 [OBSOLETE]       } else {
//DX20210301 [OBSOLETE]         if(PRINT_NULL_JSON){ sscontent_json << "\"setting_ITC\":null" << eendl;}
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]       // representative Wyckoff positions
//DX20210301 [OBSOLETE]       if(str_sg.wyccar_ITC.size()){
//DX20210301 [OBSOLETE]         sscontent_json << "\"Wyckoff_positions\":["<< eendl;
//DX20210301 [OBSOLETE]         vector<string> wyckoff_set;
//DX20210301 [OBSOLETE]         for(uint i=0;i<str_sg.wyccar_ITC.size();i++){
//DX20210301 [OBSOLETE]           if(i>4 && i!=str_sg.wyccar_ITC.size()-1){ //Skip title, scale, lattice parameters, number of atoms, and coordinate type, and last newline
//DX20210301 [OBSOLETE]             vector<string> tokens;
//DX20210301 [OBSOLETE]             aurostd::string2tokens(str_sg.wyccar_ITC[i],tokens," ");
//DX20210301 [OBSOLETE]             string name = "\"name\":";
//DX20210301 [OBSOLETE]             xvector<double> position;
//DX20210301 [OBSOLETE]             string multiplicity = "\"multiplicity\":";
//DX20210301 [OBSOLETE]             string wyckoff_letter = "\"Wyckoff_letter\":";
//DX20210301 [OBSOLETE]             string site_symmetry = "\"site_symmetry\":";
//DX20210301 [OBSOLETE]             stringstream sswyckoff;
//DX20210301 [OBSOLETE]             vector<string> vwyckoff_json;
//DX20210301 [OBSOLETE]             for(uint t=0;t<tokens.size();t++){
//DX20210301 [OBSOLETE]               if(t==0){
//DX20210301 [OBSOLETE]                 position(1) = aurostd::string2utype<double>(tokens[t]);
//DX20210301 [OBSOLETE]               }
//DX20210301 [OBSOLETE]               if(t==1){
//DX20210301 [OBSOLETE]                 position(2) = aurostd::string2utype<double>(tokens[t]);
//DX20210301 [OBSOLETE]               }
//DX20210301 [OBSOLETE]               if(t==2){
//DX20210301 [OBSOLETE]                 position(3) = aurostd::string2utype<double>(tokens[t]);
//DX20210301 [OBSOLETE]               }
//DX20210301 [OBSOLETE]               if(t==3){
//DX20210301 [OBSOLETE]                 name += "\""+tokens[t]+"\"";
//DX20210301 [OBSOLETE]               }
//DX20210301 [OBSOLETE]               if(t==4){
//DX20210301 [OBSOLETE]                 multiplicity += tokens[t];
//DX20210301 [OBSOLETE]               }
//DX20210301 [OBSOLETE]               if(t==5){
//DX20210301 [OBSOLETE]                 wyckoff_letter += "\""+tokens[t]+"\"";
//DX20210301 [OBSOLETE]               }
//DX20210301 [OBSOLETE]               if(t==6){
//DX20210301 [OBSOLETE]                 site_symmetry += "\""+tokens[t]+"\"";
//DX20210301 [OBSOLETE]               }
//DX20210301 [OBSOLETE]             }
//DX20210301 [OBSOLETE]             sswyckoff << "\"position\":[" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(position,5,roff),",") << "]" << eendl;
//DX20210301 [OBSOLETE]             vwyckoff_json.push_back(sswyckoff.str()); sswyckoff.str("");
//DX20210301 [OBSOLETE]             sswyckoff << name << eendl;
//DX20210301 [OBSOLETE]             vwyckoff_json.push_back(sswyckoff.str()); sswyckoff.str("");
//DX20210301 [OBSOLETE]             sswyckoff << multiplicity << eendl;
//DX20210301 [OBSOLETE]             vwyckoff_json.push_back(sswyckoff.str()); sswyckoff.str("");
//DX20210301 [OBSOLETE]             sswyckoff << wyckoff_letter << eendl;
//DX20210301 [OBSOLETE]             vwyckoff_json.push_back(sswyckoff.str()); sswyckoff.str("");
//DX20210301 [OBSOLETE]             sswyckoff << site_symmetry << eendl;
//DX20210301 [OBSOLETE]             vwyckoff_json.push_back(sswyckoff.str()); sswyckoff.str("");
//DX20210301 [OBSOLETE]             sswyckoff << "{" << aurostd::joinWDelimiter(vwyckoff_json,",") << "}" << eendl;
//DX20210301 [OBSOLETE]             wyckoff_set.push_back(sswyckoff.str()); sswyckoff.str("");
//DX20210301 [OBSOLETE]           }
//DX20210301 [OBSOLETE]         }
//DX20210301 [OBSOLETE]         sscontent_json << aurostd::joinWDelimiter(wyckoff_set,",") << eendl;  	
//DX20210301 [OBSOLETE]         sscontent_json << "]" << eendl;  	
//DX20210301 [OBSOLETE]         vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]       // WYCCAR
//DX20210301 [OBSOLETE]       if(str_sg.wyccar_ITC.size()){
//DX20210301 [OBSOLETE]         // convert vector<string> of WYCCAR to xstructure
//DX20210301 [OBSOLETE]         stringstream wss;
//DX20210301 [OBSOLETE]         for(uint i=0;i<str_sg.wyccar_ITC.size();i++){
//DX20210301 [OBSOLETE]           wss << str_sg.wyccar_ITC[i] << endl;
//DX20210301 [OBSOLETE]         }
//DX20210301 [OBSOLETE]         xstructure xstr_wyccar(wss,IOVASP_WYCKCAR);
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]         sscontent_json << "\"wyccar\":" << xstructure2json(xstr_wyccar) << eendl;
//DX20210301 [OBSOLETE]       } else {
//DX20210301 [OBSOLETE]         if(PRINT_NULL_JSON){ sscontent_json << "\"wyccar\":null" << eendl;}
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
//DX20210301 [OBSOLETE]
//DX20210301 [OBSOLETE]       if(standalone){
//DX20210301 [OBSOLETE]         oss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}" << endl;
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       else {
//DX20210301 [OBSOLETE]         oss << aurostd::joinWDelimiter(vcontent_json,",");
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]     }
//DX20210301 [OBSOLETE]     oss_final << oss.str();
//DX20210301 [OBSOLETE]     //DX20180822 - put attributes into vpflow - START
//DX20210301 [OBSOLETE]     if(vpflow.flag("SGDATA::CALCULATED")){ //DX20180823 - if calculated, remove
//DX20210301 [OBSOLETE]       vpflow.pop_attached("SGDATA::SPACE_GROUP_NUMBER");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("SGDATA::SPACE_GROUP_HERMANN_MAUGUIN");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("SGDATA::SPACE_GROUP_HALL");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("SGDATA::SPACE_GROUP_SCHOENFLIES");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("SGDATA::LAUE");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("SGDATA::CRYSTAL_CLASS");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("SGDATA::SETTING_ITC");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("SGDATA::ORIGIN_ITC");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("SGDATA::GENERAL_POSITION_ITC");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("SGDATA::WYCKOFF_LETTERS");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("SGDATA::WYCKOFF_MULTIPLICITIES");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("SGDATA::WYCKOFF_SITE_SYMMETRIES");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("SGDATA::WYCKOFF_POSITIONS");
//DX20210301 [OBSOLETE]       vpflow.pop_attached("SGDATA::WYCCAR");
//DX20210301 [OBSOLETE]     }
//DX20210301 [OBSOLETE]     //DX20180823 - now populate xoption
//DX20210301 [OBSOLETE]     vpflow.flag("SGDATA::CALCULATED");
//DX20210301 [OBSOLETE]     vpflow.push_attached("SGDATA::SPACE_GROUP_NUMBER",aurostd::utype2string<uint>(str_sg.space_group_ITC));
//DX20210301 [OBSOLETE]     vpflow.push_attached("SGDATA::SPACE_GROUP_HERMANN_MAUGUIN",GetSpaceGroupName(str_sg.space_group_ITC, str_sg.directory));
//DX20210301 [OBSOLETE]     vpflow.push_attached("SGDATA::SPACE_GROUP_HALL",GetSpaceGroupHall(str_sg.space_group_ITC, setting, str_sg.directory));
//DX20210301 [OBSOLETE]     vpflow.push_attached("SGDATA::SPACE_GROUP_SCHOENFLIES",GetSpaceGroupSchoenflies(str_sg.space_group_ITC, str_sg.directory));
//DX20210301 [OBSOLETE]     vpflow.push_attached("SGDATA::LAUE",GetLaueLabel(str_sg.point_group_ITC));
//DX20210301 [OBSOLETE]     vpflow.push_attached("SGDATA::CRYSTAL_CLASS",str_sg.point_group_ITC);
//DX20210301 [OBSOLETE]     vpflow.push_attached("SGDATA::SETTING_ITC",aurostd::utype2string<uint>(str_sg.setting_ITC));
//DX20210301 [OBSOLETE]     vpflow.push_attached("SGDATA::ORIGIN_ITC",aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(str_sg.origin_ITC,5,true),","));
//DX20210301 [OBSOLETE]     // turn vector of strings into single string
//DX20210301 [OBSOLETE]     stringstream general_position;
//DX20210301 [OBSOLETE]     general_position << "[";
//DX20210301 [OBSOLETE]     for(uint i=0;i<str_sg.general_position_ITC.size();i++){
//DX20210301 [OBSOLETE]       vector<string> tokens;
//DX20210301 [OBSOLETE]       aurostd::string2tokens(str_sg.general_position_ITC[i],tokens,",");
//DX20210301 [OBSOLETE]       general_position << "[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(tokens,"\""),",") << "]";
//DX20210301 [OBSOLETE]       if(i != str_sg.general_position_ITC.size()-1){
//DX20210301 [OBSOLETE]         general_position << ",";
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]     }
//DX20210301 [OBSOLETE]     general_position << "]";
//DX20210301 [OBSOLETE]     vpflow.push_attached("SGDATA::GENERAL_POSITION_ITC",general_position.str());
//DX20210301 [OBSOLETE]     stringstream wyccar;
//DX20210301 [OBSOLETE]     stringstream Wyckoff_positions;
//DX20210301 [OBSOLETE]     for(uint i=0;i<str_sg.wyccar_ITC.size();i++){
//DX20210301 [OBSOLETE]       if(i>4 && i!=str_sg.wyccar_ITC.size()-1){ //Skip title, scale, lattice parameters, number of atoms, and coordinate type, and last newline
//DX20210301 [OBSOLETE]         Wyckoff_positions << str_sg.wyccar_ITC[i] << endl;
//DX20210301 [OBSOLETE]       }
//DX20210301 [OBSOLETE]       wyccar << str_sg.wyccar_ITC[i] << endl;
//DX20210301 [OBSOLETE]     }
//DX20210301 [OBSOLETE]     vpflow.push_attached("SGDATA::WYCKOFF_LETTERS",SYM::ExtractWyckoffLettersString(str_sg.wyccar_ITC));
//DX20210301 [OBSOLETE]     vpflow.push_attached("SGDATA::WYCKOFF_MULTIPLICITIES",SYM::ExtractWyckoffMultiplicitiesString(str_sg.wyccar_ITC));
//DX20210301 [OBSOLETE]     vpflow.push_attached("SGDATA::WYCKOFF_SITE_SYMMETRIES",SYM::ExtractWyckoffSiteSymmetriesString(str_sg.wyccar_ITC));
//DX20210301 [OBSOLETE]     vpflow.push_attached("SGDATA::WYCKOFF_POSITIONS",Wyckoff_positions.str());
//DX20210301 [OBSOLETE]     vpflow.push_attached("SGDATA::WYCCAR",wyccar.str());
//DX20210301 [OBSOLETE]     //DX20180822 - put attributes into vpflow - END
//DX20210301 [OBSOLETE]     return true;
//DX20210301 [OBSOLETE]   }
//DX20210301 [OBSOLETE] }

// **************************************************************************
// PrintShell PrintShell
// **************************************************************************
// This function prints out all points in unit cell with
//  ns atoms between rmin and rmax.
// --------------------------------------------------
// GetAngleMat
// Given a point q and list of points, {p_i} returns the xmatrix
// M_ij = Angle from pi-q-pj.
// Dane Morgan - Stefano Curtarolo
aurostd::matrix<double> GetAngleMat(const xvector<double>& q, const vector<xvector<double> >& p) { //CO20200404 pflow::matrix()->aurostd::matrix()
  double TOL=1e-6;
  int np=p.size();
  aurostd::matrix<double> M(np,np);pflow::VVset(M,0.0);  //CO20200404 pflow::matrix()->aurostd::matrix()
  for(int i=0;i<np;i++) {
    for(int j=0;j<np;j++) {
      xvector<double> a(3);a=p.at(i)-q;
      xvector<double> b(3);b=p.at(j)-q;
      //  double na=aurostd::modulus(a); //DM not used
      //  double nb=aurostd::modulus(b); //DM not used
      double angle;
      double xcos=getcos(a,b);
      if(xcos>(1-TOL)) {
        angle=0;
      } else {
        if(xcos<(-1+TOL)) {
          angle=180;
        } else {
          angle=rad2deg*acos(xcos);
        }
      }
      M[i][j]=angle;
    }
  }
  return M;
}

// ***************************************************************************
// PrintShell GetGoodShellPoints
// ***************************************************************************
// This function collects the information about atoms
// within a shell around a set of points.  The points are
// input in fractional coordinates.  The points with the
// correct shells (ns atoms between rmin and rmax) are
// stored, with their neighbors, as atoms in neigh_mat.
// Dane Morgan - Stefano Curtarolo
void GetGoodShellPoints(vector<xvector<double> >& points, const xstructure& str,
    const int& ns, const double& rmin,const double& rmax,
    const string& sname,deque<deque<_atom> >& neigh_mat) {
  xstructure sstr=str;
  sstr=ReScale(sstr,1.0);
  xmatrix<double> lat(3,3);
  lat=sstr.lattice;

  // define fake atoms for each point
  deque<_atom> atvec;
  for(uint ipt=0;ipt<points.size();ipt++) {
    _atom a;
    // Set positions, other defaults are fine for a.
    a.fpos=points.at(ipt);
    a.cpos=F2C(lat,a.fpos);
    atvec.push_back(a);
  }

  // Send atoms to GetNeighData to get their neigh info.
  //[OBSOLETE]  pflow::GetNeighData(atvec,str,rmin,rmax,neigh_mat);
  sstr.GetNeighData_20220101(atvec,rmin,rmax,neigh_mat);  //CO20220623 - use orig GetNeighData() since input is NOT the atoms of the structure, this routine should be rewritten

  // Remove every point that does not have ns atoms in shell.
  deque<deque<_atom> > nnm;
  for(uint ip=0;ip<neigh_mat.size();ip++) {
    if((int) neigh_mat.at(ip).size()==(ns+1)) { // size is number of neighbors + the point.
      nnm.push_back(neigh_mat.at(ip));
    }
  }
  // Remove every point that does not have neighbors with the required name.
  // If sname==NONE then keep all the points.
  deque<deque<_atom> > nnnm;
  string none="NONE";
  if(!(sname==none)) {
    for(uint ip=0;ip<nnm.size();ip++) {
      int names_match=1;
      for(uint in=1;in<nnm.at(ip).size();in++) {
        if(!(sname==nnm.at(ip).at(in).name)) {
          names_match=0;
        }
      } // for in
      if(names_match) nnnm.push_back(nnm.at(ip));
    } // for ip
    neigh_mat=nnnm;
  } // if sname
  else {
    neigh_mat=nnm;
  } // if sname
}

// ***************************************************************************
// PrintShell new_point
// ***************************************************************************
// return true for a point with a different average
// neighbor position.
// Dane Morgan - Stefano Curtarolo
bool new_point(const vector<xvector<double> >& avnfpos, int na) {
  for(int ia=0;ia<na;ia++)
    if(identical(avnfpos.at(ia),avnfpos.at(na),1e-15)) return FALSE;
  return TRUE;
}

// --------------------------------------------------
// PrintShell
// Dane Morgan - Stefano Curtarolo
void PrintShell(const xstructure& str, const int& ns,const double& rmin, const double& rmax,const string& sname, const int lin_dens, ostream& oss) {  
  oss.setf(std::ios::left, std::ios::adjustfield);
  oss.setf(std::ios::fixed, std::ios::floatfield);
  xstructure sstr=str;
  // double scale=str.scale; //DM not used
  sstr=ReScale(sstr,1.0);
  xmatrix<double> lattice=sstr.lattice;
  //  oss.setf(std::ios::fixed,std::ios::floatfield);
  //  oss.precision(10);
  deque<deque<_atom> > neigh_mat;
  vector<xvector<double> > points;
  for(int ix=0;ix<lin_dens;ix++) {
    for(int iy=0;iy<lin_dens;iy++) {
      for(int iz=0;iz<lin_dens;iz++) {
        xvector<double> pt(3);
        pt(1)=ix/(double)lin_dens;
        pt(2)=iy/(double)lin_dens;
        pt(3)=iz/(double)lin_dens;
        points.push_back(pt);
      }
    }
  }
  GetGoodShellPoints(points,str,ns,rmin,rmax,sname,neigh_mat);
  // double tol=1e-15; //DM not used

  // Get average fpos of neighbors of each point.
  vector<xvector<double> > avnfpos;
  for(uint ia=0;ia<neigh_mat.size();ia++) {
    xvector<double> av(3);
    for(uint in=1;in<neigh_mat.at(ia).size();in++) {
      _atom a=neigh_mat.at(ia).at(in);
      //   xvector<double> pos(3);
      //  pos=a.fpos;
      av=av+a.fpos;
    } // for in
    for(int ic=1;ic<=3;ic++) {
      av(ic)=av(ic)/(neigh_mat.at(ia).size()-1);
    } // for ic
    avnfpos.push_back(av);
  }
  // Output header
  oss << "All atoms identified by:  ";
  oss << "[basis number]  [name]  [n1 n2 n3 (unit cell)]" << endl;
  oss << "Points must have " << ns << " atoms around them within a shell bounded by " << rmin << " and " << rmax << endl;
  oss << "Atoms must be of type: " << sname << ", ";
  oss << "Linear density of searched points is: " << lin_dens << endl;
  oss << "All coordinates are direct (cartesian) if the input POSCAR is direct (cartesian)." << endl;

  // Print out all points with unique environments.
  oss << endl;
  oss << "********** Unique environment points meeting shell constraints **********" << endl;
  int newcnt=0;
  vector<xvector<double> > all_new_pts;
  for(uint ia=0;ia<neigh_mat.size();ia++) {
    vector<xvector<double> > p; // For neighbors of averge point q.
    if(new_point(avnfpos,ia)) {
      newcnt++;
      _atom a;
      xvector<double> pt(3);pt=avnfpos.at(ia);
      all_new_pts.push_back(pt);
      a.fpos=pt;
      a.cpos=F2C(lattice,pt);
      if(sstr.coord_flag==FALSE) { // direct
        oss << setprecision(10) << setw(7) << ia+1 << "(" << newcnt << ")"
          << "  UAvg "
          << " " << setw(10) << pt(1)
          << " " << setw(10) << pt(2)
          << " " << setw(10) << pt(3) << endl;
      }
      if(sstr.coord_flag==TRUE) { // cart
        xvector<double> cv(3);cv=F2C(lattice,pt);
        oss << setprecision(10) << setw(7) << ia+1 << "(" << newcnt << ")"
          << "  UAvg "
          << " " << setw(10) << cv(1)
          << " " << setw(10) << cv(2)
          << " " << setw(10) << cv(3) << endl;
      }

      // Print out avg position of neighbors.
      //      oss << "        ";
      //      oss << setw(4) << "Avg ";
      //      oss << setw(4) << " ";
      //      if(sstr.coord_flag==0) { // direct
      //      oss << setprecision(10)<< setw(10) << avnfpos.at(ia)(1) << " " << setw(10) << avnfpos.at(ia)(2) << " " << setw(10) << avnfpos.at(ia)(3) << endl;
      //      }
      //      if(sstr.coord_flag==1) { // cart
      //      vector<double> cv=F2C(lattice,avnfpos.at(ia));
      //      oss << setprecision(10)<< setw(10) << cv(1) << " " << setw(10) << cv(2) << " " << setw(10) << cv(3) << endl;
      //      }
      // Print out neighbors.
      for(uint in=1;in<neigh_mat.at(ia).size();in++) {
        _atom an=neigh_mat.at(ia)[in];
        xvector<int> ijk(3);
        ijk=an.ijk;
        oss << "        ";
        oss << setw(4) << an.basis+1 << " "; //[CO20200130 - number->basis]oss << setw(4) << an.number+1 << " ";
        oss << setw(4) << an.name.c_str() << "   ";
        oss << setw(3) << ijk(1) << " " << setw(3) << ijk(2) << " " << setw(3) << ijk(3) << "   ";
        oss << setprecision(4) << AtomDist(a,an);
        oss << endl;
        p.push_back(an.cpos);
      } // in

      // Print out angles with neighbors.
      xvector<double> q(3);q=F2C(lattice,avnfpos.at(ia));
      aurostd::matrix<double> ang_mat=GetAngleMat(q,p);  //CO20200404 pflow::matrix()->aurostd::matrix()
      int np=p.size();
      oss << "        ";
      oss << "Angles"<<endl;
      oss << "           ";
      for(int i=0;i<np;i++) {
        oss << " " << setw(6) << i+1;
      }
      oss << endl;
      for(int i=0;i<np;i++) {
        oss << "        ";
        oss << setw(2) << i+1 ;
        for(int j=0;j<np;j++) {
          oss << " " << setprecision(2) << setw(6) << ang_mat[i][j];
        }
        oss << endl;
      }
    } // if newpoint
  } //ia

  // Print out all points.
  oss << endl;
  oss << "********** All points meeting shell constraints **********" << endl;
  for(uint ia=0;ia<neigh_mat.size();ia++) {
    _atom a=neigh_mat.at(ia).at(0);
    xvector<double> pt(3);pt=a.fpos;
    if(sstr.coord_flag==FALSE) { // direct
      oss << setprecision(10) << setw(7) << ia+1
        << " " << setw(10) << pt(1)
        << " " << setw(10) << pt(2)
        << " " << setw(10) << pt(3) << endl;
    }
    if(sstr.coord_flag==TRUE) { // cart
      xvector<double> cv(3);cv=F2C(lattice,pt);
      oss << setprecision(10) << setw(7) << ia+1
        << " " << setw(10) << cv(1)
        << " " << setw(10) << cv(2)
        << " " << setw(10) << cv(3) << endl;
    }
    for(uint in=1;in<neigh_mat.at(ia).size();in++) {
      _atom an=neigh_mat.at(ia).at(in);
      xvector<int> ijk(3);
      ijk=an.ijk;
      oss << "        ";
      oss << setw(4) << an.basis+1 << " "; //[CO20200130 - number->basis]oss << setw(4) << an.number+1 << " ";
      oss << setw(4) << an.name.c_str() << "   ";
      oss << setw(3) << ijk(1) << " " << setw(3) << ijk(2) << " " << setw(3) << ijk(3) << "   ";
      oss << setprecision(4) << AtomDist(a,an);
      oss << endl;
    } // in
    // Print out avg position of neighbors.
    oss << "        ";
    oss << setw(4) << "Avg ";
    oss << setw(4) << " ";
    if(sstr.coord_flag==FALSE) { // direct
      oss << setprecision(10)<< setw(10) << avnfpos.at(ia)(1) << " " << setw(10) << avnfpos.at(ia)(2) << " " << setw(10) << avnfpos.at(ia)(3) << endl;
    }
    if(sstr.coord_flag==TRUE) { // cart
      xvector<double> cv(3);cv=F2C(lattice,avnfpos.at(ia));
      oss << setprecision(10)<< setw(10) << cv(1) << " " << setw(10) << cv(2) << " " << setw(10) << cv(3) << endl;
    }
  } // ia

  // WORKING HERE ???  Finish this and check consistency of
  //  setting names, pos, nat, and particularly types - xstructure should
  //    keep all these consistent when any one is changed.
  //  Add new points to the POSCAR file and output.
  oss << endl;
  oss << "********** POSCAR file with unique environment points meeting shell constraints **********" << endl;
  // Number of atoms

  uint nnat=all_new_pts.size();
  uint nat=str.atoms.size();
  // uint tnat=nat+nnat;  //DM not used
  // Number of types
  sstr.num_each_type.push_back(nnat);
  sstr.comp_each_type.push_back((double) nnat);
  for(uint iat=0;iat<nnat;iat++) {
    _atom newatom;
    newatom.fpos=all_new_pts.at(iat);
    newatom.cpos=F2C(sstr.lattice,newatom.fpos);
    newatom.corigin=sstr.origin;
    newatom.type=sstr.num_each_type.size()-1;   // CONVASP_MODE
    ostringstream ostmp;ostmp<<iat+1;
    newatom.name="SHELL"+ostmp.str();
    newatom.name_is_given=TRUE;
    newatom.cleanname="XX";
    newatom.sd="";
    newatom.atomic_number=-1;
    //[CO20200130 - number->basis]newatom.number=nat+iat;                 // reference position for convasp
    newatom.basis=nat+iat;                  // position in the basis	
    clear(newatom.ijk);                     // position in the unit cell
    sstr.atoms.push_back(newatom);
  }
  //  sstr.write_DEBUG_flag=TRUE;
  oss << sstr;
} // end routine

// void PrintShell(const xstructure& str, const int& ns,const double& rmin,
// const double& rmax,const string& sname, const int lin_dens)

// ***************************************************************************
// PrintSumDOS
// ***************************************************************************
// Prints out the summed pdos in user friendly format.
// Dane Morgan - Stefano Curtarolo
void PrintSumPDOS(pflow::pdosdata& pdd, ostream& oss) {
  int spin=pdd.spin-1;
  pdd.PrintPDOS(oss,spin);
  if(pdd.print_params==1) {
    pflow::projdata prd;
    std::vector<string> LMnames=prd.LMnames;
    pdd.PrintParams(oss,LMnames);
  }
}

// ***************************************************************************
// PrintXYZ  PrintXYZws PrintXYZInSphere PrintXYZNanoparticle
// ***************************************************************************
// This funtion prints out structural data in an XYZ format
// Stefano Curtarolo
void PrintXYZ(const xstructure& a, const xvector<int>& n, ostream& oss) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  xstructure sstr=a;
  sstr=ReScale(sstr,1.0);
  // Print out data in XYZ format
  oss << sstr.atoms.size()*n(1)*n(2)*n(3) << endl;
  oss << sstr.title << endl;
  for(int i=0;i<n(1);i++) {
    for(int j=0;j<n(2);j++) {
      for(int k=0;k<n(3);k++) {
        for(uint iat=0;iat<sstr.atoms.size();iat++) {
          oss << sstr.atoms.at(iat).cleanname << " "
            << (sstr.atoms.at(iat).cpos(1)+i*sstr.lattice(1,1)+j*sstr.lattice(2,1)+k*sstr.lattice(3,1)) << " "
            << (sstr.atoms.at(iat).cpos(2)+i*sstr.lattice(1,2)+j*sstr.lattice(2,2)+k*sstr.lattice(3,2)) << " "
            << (sstr.atoms.at(iat).cpos(3)+i*sstr.lattice(1,3)+j*sstr.lattice(2,3)+k*sstr.lattice(3,3)) << " " << endl;
        } // iat
      } // k
    } // j
  } // i
}

// ***************************************************************************
// PrintXYZws
// ***************************************************************************
// This funtion moves atoms inside WS cell and prints out structural
// data in an XYZ format
// Stefano Curtarolo
void PrintXYZws(const xstructure& a, ostream& oss) {
  //SC20030110
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  xstructure sstr=a;
  // Set scale to 1 so you don't need to rescale coordinates.
  sstr=ReScale(sstr,1.0);
  // zoology of vectors
  xvector<double> rat(3),xoo(3),yoo(3),zoo(3),xyo(3),xzo(3),yzo(3),xyz(3);
  double axoo,ayoo,azoo,axyo,axzo,ayzo,axyz;
  xoo=sstr.lattice(1);axoo=aurostd::modulus(xoo);
  yoo=sstr.lattice(2);ayoo=aurostd::modulus(yoo);
  zoo=sstr.lattice(3);azoo=aurostd::modulus(zoo);
  xyo=sstr.lattice(1)+sstr.lattice(2);axyo=aurostd::modulus(xyo);
  xzo=sstr.lattice(1)+sstr.lattice(3);axzo=aurostd::modulus(xzo);
  yzo=sstr.lattice(2)+sstr.lattice(3);ayzo=aurostd::modulus(yzo);
  xyz=sstr.lattice(1)+sstr.lattice(2)+sstr.lattice(3);axyz=aurostd::modulus(xyz);
  double projxoo,projyoo,projzoo,projxyo,projxzo,projyzo,projxyz;

  // Print out data in XYZ format
  oss << sstr.atoms.size() << endl;
  oss << sstr.title << endl;
  for(uint iat=0;iat<sstr.atoms.size();iat++)
    for(int i=-1;i<=1;i++)
      for(int j=-1;j<=1;j++)
        for(int k=-1;k<=1;k++)
        {
          rat=sstr.atoms.at(iat).cpos+((double)i)*sstr.lattice(1)+((double)j)*sstr.lattice(2)+((double)k)*sstr.lattice(3);
          projxoo=scalar_product(rat,xoo)/axoo/axoo;
          projyoo=scalar_product(rat,yoo)/ayoo/ayoo;
          projzoo=scalar_product(rat,zoo)/azoo/azoo;
          projxyo=scalar_product(rat,xyo)/axyo/axyo;
          projxzo=scalar_product(rat,xzo)/axzo/axzo;
          projyzo=scalar_product(rat,yzo)/ayzo/ayzo;
          projxyz=scalar_product(rat,xyz)/axyz/axyz;
          //  if((projxoo>-0.5 && projxoo<=0.5))
          //  oss << projxoo << " " << projyoo << " " << projzoo << endl;
          if((projxoo>-0.5 && projxoo<=0.5) &&
              (projyoo>-0.5 && projyoo<=0.5) &&
              (projzoo>-0.5 && projzoo<=0.5) &&
              (projxyo>-0.5 && projxyo<=0.5) &&
              (projxzo>-0.5 && projxzo<=0.5) &&
              (projyzo>-0.5 && projyzo<=0.5) &&
              (projxyz>-0.5 && projxyz<=0.5))
            oss << sstr.atoms.at(iat).cleanname << " " << rat(1) << " " << rat(2) << " " << rat(3) << " " << endl;
        }
}

// ***************************************************************************
// PrintXYZInSphere
// ***************************************************************************
// This funtion prints out structural data inside a sphere with radius radius
// and 0,0,0 as origin, in an XYZ format
// Stefano Curtarolo
void PrintXYZInSphere(const xstructure& a, const double& radius, ostream& oss) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  xstructure sstr=a;
#define _SPHERE_CUT_OFF_ 5

  // Set scale to 1 so you don't need to rescale coordinates.
  sstr=ReScale(sstr,1.0);

  int atoms_inside=0;
  xvector<double> x(3);
  // Print out data in InSphere format
  for(int i=-_SPHERE_CUT_OFF_;i<_SPHERE_CUT_OFF_;i++) {
    for(int j=-_SPHERE_CUT_OFF_;j<_SPHERE_CUT_OFF_;j++) {
      for(int k=-_SPHERE_CUT_OFF_;k<_SPHERE_CUT_OFF_;k++) {
        for(uint iat=0;iat<sstr.atoms.size();iat++) {
          x(1)=(sstr.atoms.at(iat).cpos(1)+i*sstr.lattice(1,1)+j*sstr.lattice(2,1)+k*sstr.lattice(3,1));
          x(2)=(sstr.atoms.at(iat).cpos(2)+i*sstr.lattice(1,2)+j*sstr.lattice(2,2)+k*sstr.lattice(3,2));
          x(3)=(sstr.atoms.at(iat).cpos(3)+i*sstr.lattice(1,3)+j*sstr.lattice(2,3)+k*sstr.lattice(3,3));
          if(aurostd::modulus(x)<=radius) atoms_inside++;
        } // iat
      } // k
    } // j
  } // i

  oss << atoms_inside << endl;
  oss << sstr.title << endl;
  // Print out data in InSphere format
  for(int i=-_SPHERE_CUT_OFF_;i<_SPHERE_CUT_OFF_;i++) {
    for(int j=-_SPHERE_CUT_OFF_;j<_SPHERE_CUT_OFF_;j++) {
      for(int k=-_SPHERE_CUT_OFF_;k<_SPHERE_CUT_OFF_;k++) {
        for(uint iat=0;iat<sstr.atoms.size();iat++) {
          x(1)=(sstr.atoms.at(iat).cpos(1)+i*sstr.lattice(1,1)+j*sstr.lattice(2,1)+k*sstr.lattice(3,1));
          x(2)=(sstr.atoms.at(iat).cpos(2)+i*sstr.lattice(1,2)+j*sstr.lattice(2,2)+k*sstr.lattice(3,2));
          x(3)=(sstr.atoms.at(iat).cpos(3)+i*sstr.lattice(1,3)+j*sstr.lattice(2,3)+k*sstr.lattice(3,3));
          if(aurostd::modulus(x)<=radius)
            oss << sstr.atoms.at(iat).cleanname << " " << x(1) << " " << x(2) << " " << x(3) << " " << endl;
        } // iat
      } // k
    } // j
  } // i
}

// ***************************************************************************
// PrintXRAY CorrectionFactor
// ***************************************************************************
// This function contains any correction factors that need
// to be multiply the structure factors.  Theta should be in rads.
// Includes - Lorentz-Polarization
// Dane Morgan - Stefano Curtarolo
double CorrectionFactor(const double& th) {
  //  double th=theta*TWOPI/360.0; // rad->degrees
  double correction_factor=1;
  // Lorentz-polarization
  double lp;
  double c2t=cos(2.0*th);
  double st=sin(th);
  if(th>0 || th<0) {
    lp=(1.0+c2t*c2t)/(st*st);
  }
  else {
    lp=1;
  }
  correction_factor=correction_factor*lp;
  return lp;
}  

// ***************************************************************************
// PrintXRAY PrintXray
// ***************************************************************************
// Dane Morgan - Stefano Curtarolo
// Updated by Corey Oses 20190520
void PrintXray(const xstructure& str, double lambda, ostream& oss) {
  oss.setf(std::ios::left, std::ios::adjustfield);
  oss.setf(std::ios::fixed, std::ios::floatfield);
  //  oss.setf(std::ios::fixed,std::ios::floatfield);
  int PREC_DEFAULT=4;
  oss.precision(PREC_DEFAULT);

  int num_atoms=str.atoms.size();
  //  vector<string> names=str.GetNames();
  vector<double> dist,sf;
  vector<double> scatt_fact(num_atoms,0.0);
  vector<double> mass(num_atoms,0.0);
  vector<double> twoB_vec(num_atoms,0.0);
  //pflow::GetXray(str,dist,sf,lambda,scatt_fact,mass,twoB_vec);
  vector<vector<double> > ids;
  aurostd::matrix<double> data;  //CO20200404 pflow::matrix()->aurostd::matrix()
  pflow::GetXrayData(str,dist,sf,scatt_fact,mass,twoB_vec,ids,data,lambda); //CO20190620 - intmax can be grabbed later

  //get intmax
  double intmax=1e-8;
  for(uint i=0;i<data.size();i++){if(data[i][4]>intmax){intmax=data[i][4];}}

  double tol=XRAY_THETA_TOL; //1.0E-5;
  int w1=4; // int
  int w2=12; // some doubles
  int w3=20; // Integrated intensities
  int tlen=sf.size();
  int len=Nint(std::pow((double) tlen,(double) 1.0/3.0));
  int kmx=(len-1)/2; // len should be odd.

  // Print out all data.
  oss << "Wavelength(Ang) = " << lambda << endl;
  oss << "Atom_Name  ScattFact   Mass(amu)   B(Ang)(DW=exp(-B*sin(theta)^2/lambda^2))" << endl; //CO20190329
  for(uint iat=0;iat<(uint) num_atoms;iat++) {
    oss <<setw(4)<<iat+1<<" "<<setw(4)<<str.atoms.at(iat).cleanname << " ";
    if(str.atoms.at(iat).cleanname.length()>1) oss << " ";
    oss <<setw(w2)<<scatt_fact[iat]<<setw(w2)<<KILOGRAM2AMU*mass[iat]<<setw(w2)<<1E+20*twoB_vec[iat]/2.0<<endl;
  }
  oss << "******************** All data ********************" << endl;
  oss << "2*theta      Intensity            h    k    l    dist         keyword " << endl;
  for(int i0=-kmx;i0<=kmx;i0++) {
    for(int i1=-kmx;i1<=kmx;i1++) {
      for(int i2=-kmx;i2<=kmx;i2++) {
        int ii0=i0+kmx;
        int ii1=i1+kmx;
        int ii2=i2+kmx;
        int id1=ii2+ii1*len+ii0*len*len;
        int id=(int) ids[id1][1];
        ii0=(int) ids[id1][2];
        ii1=(int) ids[id1][3];
        ii2=(int) ids[id1][4];
        double theta=0;
        if(dist[id]>0) {
          double term=lambda/(2.0*dist[id]);
          if(term<=1) {
            theta=asin(term);
            theta=theta*360.0/TWOPI; // rad->degrees
            if(theta>tol) oss
              <<setw(w2)<<2.0*theta<<" " // angle
                <<setw(w3)<<sf[id]<<" " // sf
                <<setw(w1)<<ii0<<" " // h
                <<setw(w1)<<ii1<<" " // k
                <<setw(w1)<<ii2<<" " // l
                <<setw(w2)<<dist[id]<<" " // dist
                << "SINGLE"
                << endl;
          } // if term<=1
        } // if dist>0
      } // i2
    } // i1
  } // i0

  // Output grouped data
  oss << "******************** Grouped data ********************" << endl;
  oss << "2*theta      IntIntensity         %ofMaxInt    h    k    l    dist         mult. correction    keyword " << endl;
  for(uint i=0;i<data.size();i++) {
    double theta=0;
    if(data[i][3]>0) {
      double term=lambda/(2.0*data[i][3]);
      if(term<=1) {
        theta=asin(term);
        theta=theta*360.0/TWOPI; // rad->degrees
        if(theta>tol) oss
          <<setw(w2)<<2.0*theta<<" " // angle
            <<setw(w3)<<setprecision(2)<<data[i][4]<<setprecision(PREC_DEFAULT)<<" " // sf
            <<setw(w2)<<setprecision(2)<<100*data[i][4]/intmax<<setprecision(PREC_DEFAULT)<<" " // % max sf
            <<setw(w1)<<(int)data[i][0]<<" " // h
            <<setw(w1)<<(int)data[i][1]<<" " // k
            <<setw(w1)<<(int)data[i][2]<<" " // l
            <<setw(w2)<<data[i][3]<<" " // dist
            <<setw(5)<<(int)data[i][5]<<" " // mult.
            <<setw(w2)<<CorrectionFactor(theta*TWOPI/360.0)<<" " // correction.
            << " GROUP"
            << endl;
      } // if term<=1
    } // if dist>0
  } // for

  // Output data to plot
  oss << "******************** To Plot data ********************" << endl;
  oss << "2*theta      Amplitude    keyword " << endl;
  for(uint i=0;i<data.size();i++) {
    double theta=0;
    if(data[i][3]>0) {
      double term=lambda/(2.0*data[i][3]);
      if(term<=1) {
        theta=asin(term);
        theta=theta*360.0/TWOPI; // rad->degrees
        if(theta>tol) {
          // initial 0.
          oss <<setw(w2)<<2.0*theta<<" " // angle
            <<setw(w2)<<"0"<<" " // sf
            << "TOPLOT"
            << endl;
          // true value of sf/intmax.
          oss <<setw(w2)<<2.0*theta<<" " // angle
            <<setw(w2)<<setprecision(2)<<100*data[i][4]/intmax<<setprecision(PREC_DEFAULT)<<" " // sf
            << "TOPLOT"
            << endl;
          // final 0.
          oss <<setw(w2)<<2.0*theta<<" " // angle
            <<setw(w2)<<"0"<<" " // sf
            << "TOPLOT"
            << endl;
          // tpx
        } // if theta<tol
      } // if term<=1
    } // if dist>0
  } // for
}

#endif

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************
