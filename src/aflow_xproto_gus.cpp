// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2008 - 2009
// fixed for xz - SC 2018

#ifndef _AFLOW_XPROTO_GUS_CPP_
#define _AFLOW_XPROTO_GUS_CPP_
#include "aflow.h"
#include "aflow_xproto_gus_lib.cpp"

// ***************************************************************************
// GUS GUS GUS GUS GUS GUS GUS GUS GUS GUS GUS GUS GUS GUS GUS GUS GUS GUS GUS
// ***************************************************************************
namespace aflowlib {
  xstructure PrototypeBinaryGUS(ostream &oss,string label) {
    return aflowlib::PrototypeBinaryGUS(oss,label,"A",1.0,"B",1.0,0.0);
  }
} // namespace aflowlib

namespace aflowlib {
  xstructure PrototypeBinaryGUS(ostream &oss,string label,string atomA,string atomB) {
    double atomvolumeA,atomvolumeB;
    //[CO20181106]atomvolumeA=GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(atomA)); //CO20181128
    //[CO20181106]atomvolumeB=GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(atomB)); //CO20181128
    atomvolumeA=GetAtomVolume(atomA); //CO20181128
    atomvolumeB=GetAtomVolume(atomB); //CO20181128
    return aflowlib::PrototypeBinaryGUS(oss,label,atomA,atomvolumeA,atomB,atomvolumeB,0.0);
  }
} // namespace aflowlib

// ***************************************************************************
namespace aflowlib {
  string PrototypeBinaryGUS_Cache_LibraryS_Extract(ostream &oss,const string& auslat,const string& labelclean) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    ostringstream oaus;
    string FileLibrary="",structure_line="";
    for(uint i=0;aflowlib::PrototypeBinaryGUS_Cache_LibraryS[i]!="XXX"&&structure_line.empty();i++)
      if(aurostd::substring2bool(aflowlib::PrototypeBinaryGUS_Cache_LibraryS[i],"gus"+auslat))
        if(aurostd::substring2bool(aflowlib::PrototypeBinaryGUS_Cache_LibraryS[i],auslat+labelclean+" "))
          return aflowlib::PrototypeBinaryGUS_Cache_LibraryS[i];
    if(structure_line.empty()) {
      // NEED GUS LIBRARY
      // Search LIBRARY
      string FileLibrary="";
      for(uint j=0;j<vAFLOW_LIBRARY_DIRECTORIES.size()&&FileLibrary.empty();j++) {   // cycle through possible directories
        if(FileLibrary.empty() && aurostd::FileExist(vAFLOW_LIBRARY_DIRECTORIES.at(j)+"/aflow_library_gus.dat")) FileLibrary=vAFLOW_LIBRARY_DIRECTORIES.at(j)+"/aflow_library_gus.dat";
        if(FileLibrary.empty() && aurostd::FileExist(vAFLOW_LIBRARY_DIRECTORIES.at(j)+"/aflow_library_gus.dat.gz")) FileLibrary=vAFLOW_LIBRARY_DIRECTORIES.at(j)+"/aflow_library_gus.dat.gz";
        if(FileLibrary.empty() && aurostd::FileExist(vAFLOW_LIBRARY_DIRECTORIES.at(j)+"/aflow_library_gus.dat.bz2")) FileLibrary=vAFLOW_LIBRARY_DIRECTORIES.at(j)+"/aflow_library_gus.dat.bz2";
        if(FileLibrary.empty() && aurostd::FileExist(vAFLOW_LIBRARY_DIRECTORIES.at(j)+"/aflow_library_gus.dat.xz")) FileLibrary=vAFLOW_LIBRARY_DIRECTORIES.at(j)+"/aflow_library_gus.dat.xz";
      } // cycle through AFLOW LIBRARY as postfix
      if(FileLibrary!="") {
        if(LDEBUG) { oaus << "00000  MESSAGE AFLOW LIBRARY  Found library file = [" << FileLibrary << "]" << endl; }
        if(LDEBUG) { aurostd::PrintMessageStream(oaus,XHOST.QUIET,oss); } //CO20200624
      } else {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"AFLOW_LIBRARY not found.",_RUNTIME_ERROR_);
      }
      // FOUND
      if(aurostd::substring2bool(FileLibrary,".gz")) {
        structure_line=aurostd::execute2string("zcat "+FileLibrary+" | grep gus"+auslat+" | grep \""+auslat+labelclean+" \"");
      } else {
        if(aurostd::substring2bool(FileLibrary,".bz2")) {
          structure_line=aurostd::execute2string("bzcat "+FileLibrary+" | grep gus"+auslat+" | grep \""+auslat+labelclean+" \"");
        } else {
          if(aurostd::substring2bool(FileLibrary,".xz")) {
            structure_line=aurostd::execute2string("xzcat "+FileLibrary+" | grep gus"+auslat+" | grep \""+auslat+labelclean+" \"");
          } else {
            structure_line=aurostd::execute2string("cat "+FileLibrary+" | grep gus"+auslat+" | grep \""+auslat+labelclean+" \"");
          }
        }
      }
    }
    return structure_line;  // if found, otherwise return NIHIL;
  }
} // namespace aflowlib

namespace aflowlib {
  xstructure PrototypeBinaryGUS(ostream &oss,string label,
      string atomA,double volumeA,
      string atomB,double volumeB,
      double volume_in) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    // good for all
    xstructure str("");str.lattice.clear();
    double volumeC=0.0;string atomC="C";
    str.species.clear();str.species_pp.clear();;str.species_pp_type.clear();str.species_pp_version.clear();str.species_pp_ZVAL.clear();str.species_volume.clear();
    deque<string> species_tmp; species_tmp.push_back("");species_tmp.push_back("");species_tmp.push_back("");

    std::string title_database,bulksurf,title,prototype;
    vector<string> tokens;
    stringstream aus;
    ostringstream oaus;
    uint i,j,iat,natoms=0,natomsA=0,natomsB=0,natomsC,nspecies=0;
    // uint nlabels;
    uint nD=1,labelmax,labelnum;
    double volume=0.0;
    title=atomA+atomB+"/"+label;
    // *********************************************************************
    string labelclean=aflowlib::PrototypeCleanLatticeString(label);
    labelnum=aurostd::string2utype<uint>(labelclean);
    //  vector<string> vcache;
    //  aurostd::string2vectorstring(aflowlib::PrototypeBinaryGUS_Cache_Library(),vcache);
    // *********************************************************************
    xmatrix<double> p(3,3), sLV(3,3), sLVinv(3,3), Sinv(3,3), HNF(3,3);
    //  xmatrix<double> dvec(3,nD),sLVinvXdvec(3,nD);
    xmatrix<double> dvec(3,2),sLVinvXdvec(3,2); // up to 5 basis in the primitive
    dvec.clear();sLVinvXdvec.clear();
    string _strstream;
    string structure_line="";
    string auslat="";
    // --------- FCC found --------------
    if(label[0]=='F' || label[0]=='f') {
      if(LDEBUG) { cerr << "DEBUG  fcc found" << endl; }
      p(1,1)=0.50000000;p(2,1)=0.50000000;p(3,1)=0.00000000;  //   # a1 parent lattice vector
      p(1,2)=0.50000000;p(2,2)=0.00000000;p(3,2)=0.50000000;  //   # a2 parent lattice vector
      p(1,3)=0.00000000;p(2,3)=0.50000000;p(3,3)=0.50000000;  //   # a3 parent lattice vector
      nD=1;
      dvec(1,1)=0.0000000;dvec(2,1)=0.0000000;dvec(3,1)=0.0000000;           // FCC
      //    nlabels=2;
      labelmax=163373;
      title_database="fcc";
      bulksurf="bulk";
      if(labelnum==0 || labelnum>labelmax) {
        oaus << "WWWWW  AFLOW_LIBRARY GUS FCC, label (" << labelnum << ") out of boundary (1," << labelmax << ")" << endl;
        aurostd::PrintWarningStream(oaus,XHOST.QUIET);  //CO20200624
        return str;
      } else {
        structure_line=aflowlib::PrototypeBinaryGUS_Cache_LibraryS_Extract(oss,title_database,labelclean);
      }
    }
    // --------- BCC found --------------
    if(label[0]=='B' || label[0]=='b') {
      if(LDEBUG) { cerr << "DEBUG  bcc found" << endl; }
      p(1,1)=0.50000000;p(2,1)=0.50000000;p(3,1)=-0.5000000;  //   # a1 parent lattice vector
      p(1,2)=0.50000000;p(2,2)=-0.5000000;p(3,2)=0.50000000;  //   # a2 parent lattice vector
      p(1,3)=-0.5000000;p(2,3)=0.50000000;p(3,3)=0.50000000;  //   # a3 parent lattice vector
      nD=1;
      dvec(1,1)=0.0000000;dvec(2,1)=0.0000000;dvec(3,1)=0.0000000;           // BCC
      //   nlabels=2;
      labelmax=163373;
      title_database="bcc";
      bulksurf="bulk";
      if(labelnum==0 || labelnum>labelmax) {
        oaus << "WWWWW  AFLOW_LIBRARY GUS BCC, label (" << labelnum << ") out of boundary (1," << labelmax << ")" << endl;
        aurostd::PrintWarningStream(oaus,XHOST.QUIET);  //CO20200624
        return str;
      } else {
        structure_line=aflowlib::PrototypeBinaryGUS_Cache_LibraryS_Extract(oss,title_database,labelclean);
      }
    }
    // --------- HCP found --------------
    if(label[0]=='H' || label[0]=='h') {
      if(LDEBUG) { cerr << "DEBUG  hcp found" << endl; }
      p(1,1)=1.00000000;p(2,1)=0.00000000;p(3,1)=0.00000000;  //   # a1 parent lattice vector
      p(1,2)=0.50000000;p(2,2)=sqrt(3.0)/2.0;p(3,2)=0.00000000;  //   # a2 parent lattice vector
      p(1,3)=0.00000000;p(2,3)=0.00000000;p(3,3)=4.0/sqrt(6.0);  //   # a3 parent lattice vector
      nD=2;
      dvec(1,1)=0.0000000;dvec(2,1)=0.00000000;dvec(3,1)=0.0000000;          // HCP
      dvec(1,2)=0.5000000;dvec(2,2)=(1.0/sqrt(3))/2.0;dvec(3,2)=2/sqrt(6.0); // HCP
      //   nlabels=2;
      labelmax=1643380;
      title_database="hcp";
      bulksurf="bulk";
      if(labelnum==0 || labelnum>labelmax) {
        oaus << "WWWWW  AFLOW_LIBRARY GUS HCP, label (" << labelnum << ") out of boundary (1," << labelmax << ")" << endl;
        aurostd::PrintWarningStream(oaus,XHOST.QUIET);  //CO20200624
        return str;
      } else {
        structure_line=aflowlib::PrototypeBinaryGUS_Cache_LibraryS_Extract(oss,title_database,labelclean);
      }
    }
    // --------- sc found --------------
    if(label[0]=='S' || label[0]=='s') {
      if(LDEBUG) { cerr << "DEBUG  SC found" << endl; }
      p(1,1)=1.0000000;p(2,1)= 0.0000000;p(3,1)= 0.0000000;   //   # a1 parent lattice vector
      p(1,2)=0.0000000;p(2,2)= 1.0000000;p(3,2)= 0.0000000;   //   # a2 parent lattice vector
      p(1,3)=0.0000000;p(2,3)= 0.0000000;p(3,3)= 1.0000000;   //   # a3 parent lattice vector
      nD=1;
      dvec(1,1)=0.0000000;dvec(2,1)=0.0000000;dvec(3,1)=0.0000000;           //SC
      //    nlabels=2;
      labelmax=188729;
      title_database="sc";
      bulksurf="bulk";
      if(labelnum==0 || labelnum>labelmax) {
        oaus << "WWWWW  AFLOW_LIBRARY GUS HCP, label (" << labelnum << ") out of boundary (1," << labelmax << ")" << endl;
        aurostd::PrintWarningStream(oaus,XHOST.QUIET);  //CO20200624
        return str;
      } else {
        structure_line=aflowlib::PrototypeBinaryGUS_Cache_LibraryS_Extract(oss,title_database,labelclean);
      }
    }
    // *********************************************************************
    if(structure_line.empty()) {
      oaus << "EEEEE  aflowlib::PrototypeBinaryGUS: lattice not found, label=" << label << endl;
      aurostd::PrintErrorStream(oaus,XHOST.QUIET);  //CO20200624
    }
    // clear up the structure from the beginning
    aurostd::string2tokens(structure_line,tokens," ");
    structure_line.clear();
    for(i=0;i<tokens.size();i++)
      if(i>=2) structure_line+=tokens[i]+" "; // remove the 1st two entries

    int sizeN,nAt,pgOps,a,b,c,d,e,f;
    xvector<int> diag(3);
    xmatrix<int> L(3,3);
    string labeling="";
    aus.clear();aus.str(structure_line);
    aus >> sizeN >> nAt >> pgOps;
    aus >> diag(1) >> diag(2) >> diag(3);
    aus >> a >> b >> c >> d >> e >> f;
    for(i=1;i<=3;i++)
      for(j=1;j<=3;j++)
        aus >> L(i,j);
    aus >> labeling;
    if(LDEBUG) { cerr << "DEBUG sizeN=" << sizeN << " nAt=" << nAt << " pgOps=" << pgOps << endl; }
    if(LDEBUG) { cerr << "DEBUG diag=" << diag << endl; }
    if(LDEBUG) { cerr << "DEBUG a=" << a << " b=" << b << " c=" << c << " d=" << d << " e=" << e << " f=" << f << endl; }
    if(LDEBUG) { cerr << "DEBUG L=" << L << endl; }
    if(LDEBUG) { cerr << "DEBUG labeling=" << labeling << endl; }
    // Define the full HNF matrix
    HNF.clear();
    HNF(1,1) = a; HNF(2,1) = b; HNF(2,2) = c;
    HNF(3,1) = d; HNF(3,2) = e; HNF(3,3) = f;
    if(LDEBUG) { cerr << HNF << endl; }
    // Compute the superlattice vectors
    sLV = p*HNF;
    // Find the coordinates of the basis atoms
    xmatrix<double> aBas(3,nAt*nD);
    xvector<double> ausv(3),aussLVinvXdvec(3);
    int z1, z2, z3, ic,iD;
    //[CO20191201 - OBSOLETE]aurostd::inverse(HNF,Sinv); //**
    //[CO20191201 - OBSOLETE]aurostd::inverse(sLV,sLVinv);
    Sinv=aurostd::inverse(HNF); //**  //CO20191201
    sLVinv=aurostd::inverse(sLV);  //CO20191201
    sLVinvXdvec=sLVinv*dvec; //**
    ic = 0;
    for(iD=1;iD<=(int) nD;iD++) {//**
      for( z1=0; z1<=a-1; z1++) {
        for(j=1;j<=3;j++) aussLVinvXdvec(j)=sLVinvXdvec(j,iD);
        for( z2=(b*z1)/a; z2<=c+(b*z1)/a - 1; z2++) {
          for( z3 = z1*(d-(e*b)/c)/a+(e*z2)/c; z3<= f+z1*(d-(e*b)/c)/a+(e*z2)/c - 1; z3++) {
            ic++;
            if(ic>(int) (nAt*nD)) { throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Problem in basis atoms.",_RUNTIME_ERROR_); }
            // call inverse(real(HNF,dp),Sinv)
            //** Move this to outside the loop: aurostd::inverse(HNF,Sinv);
            ausv=Sinv*aurostd::xvector<double>({(double) z1,(double) z2,(double) z3})+aussLVinvXdvec; //**
            // aurostd::inverse(HNF,Sinv);
            // ausv=Sinv*aurostd::reshape((double) z1,(double) z2,(double) z3);
            for(j=1;j<=3;j++) aBas(j,ic)=ausv(j);
          } // z3
        } // z2
      } // z1
    } // iD
    if(LDEBUG) { cerr << "DEBUG ic=" << ic << endl; }
    if(LDEBUG) { cerr << "DEBUG nD=" << nD << endl; }
    if(ic!=(int) (nAt*nD)) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Not enough basis atoms.",_RUNTIME_ERROR_);}
    if(LDEBUG) { cerr << "DEBUG sLV=" << sLV << endl; }
    if(LDEBUG) { cerr << "DEBUG aBas=" << aBas << endl; }

    // *********************************************************************
    // FIX title
    str.scale=1.0;
    str.title=title+" - "+label+" - ("+title_database+","+bulksurf+") "+_HTQC_PROJECT_STRING_;
    //  str.prototype=prototype;
    // FIX lattice/basis
    str.lattice=trasp(sLV);
    if(det(str.lattice)<0) {
      //    cerr << "WWWWW  determinant(lattice)<0, reversing a2 with a3 (and fixing the basis)" << endl;
      double temp;
      if(LDEBUG) { cerr << "DEBUG str.lattice=" << str.lattice << endl; }
      for(i=1;i<=3;i++) SWAP(str.lattice(2,i),str.lattice(3,i));
      for(i=1;i<=(uint) nAt;i++) SWAP(aBas(2,i),aBas(3,i));
      if(LDEBUG) { cerr << "DEBUG str.lattice=" << str.lattice << endl; }
    }
    // FIX lattices and operators
    str.FixLattices();

    // FIX ATOMS
    natoms=0;natomsA=0;natomsB=0;natomsC=0;nspecies=0;
    for(i=0;i<labeling.length();i++) {
      natoms++;
      if(labeling[i]=='0') {natomsA++;species_tmp.at(0)=atomA;}
      if(labeling[i]=='1') {natomsB++;species_tmp.at(1)=atomB;}
      if(labeling[i]=='2') {natomsC++;species_tmp.at(2)=atomC;}
    }
    if(natoms!=(nAt*nD)) {
      stringstream message; message << "natoms!=nAt  natoms=" << natoms << "  nAt=" << nAt;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ERROR_);
    };
    if(species_tmp.at(0)!="") {str.num_each_type.push_back(natomsA);str.comp_each_type.push_back((double) natomsA);nspecies++;}
    if(species_tmp.at(1)!="") {str.num_each_type.push_back(natomsB);str.comp_each_type.push_back((double) natomsB);nspecies++;}
    if(species_tmp.at(2)!="") {str.num_each_type.push_back(natomsC);str.comp_each_type.push_back((double) natomsC);nspecies++;}
    if(LDEBUG) { cerr << "DEBUG natoms=" << natoms << endl; }
    if(LDEBUG) { cerr << "DEBUG species_tmp.at(0)=" << species_tmp.at(0) << " natomsA=" << natomsA << endl; }
    if(LDEBUG) { cerr << "DEBUG species_tmp.at(1)=" << species_tmp.at(1) << " natomsB=" << natomsB << endl; }
    if(LDEBUG) { cerr << "DEBUG species_tmp.at(2)=" << species_tmp.at(2) << " natomsC=" << natomsC << endl; }
    if(LDEBUG) { cerr << "DEBUG nspecies=" << nspecies << endl; }
    if(LDEBUG) { cerr << "DEBUG labeling.length()=" << labeling.length() << endl; }

    // Plug ATOMS
    volume=0.0;
    xvector<double> fpos(3),cpos(3);
    // perform atoms A
    for(iat=0;iat<natoms;iat++) {
      cpos.clear();fpos.clear(); // clear positions
      if(labeling[iat]=='0') { // atoms A
        _atom atom;
        atom.CleanName();
        atom.CleanSpin();
        fpos(1)=aBas(1,iat+1);
        fpos(2)=aBas(2,iat+1);
        fpos(3)=aBas(3,iat+1);
        cpos=F2C(str.lattice,fpos);
        atom.fpos=fpos;
        atom.cpos=cpos;
        atom.name=atomA;
        atom.type=0; // A is always type 0 !
        atom.name_is_given=TRUE;
        volume+=volumeA;
        str.atoms.push_back(atom);
      }
    }
    // perform atom B
    for(iat=0;iat<natoms;iat++) {
      cpos.clear();fpos.clear(); // clear positions
      if(labeling[iat]=='1') { // atoms B
        _atom atom;
        atom.CleanName();
        atom.CleanSpin();
        fpos(1)=aBas(1,iat+1);
        fpos(2)=aBas(2,iat+1);
        fpos(3)=aBas(3,iat+1);
        cpos=F2C(str.lattice,fpos);
        atom.fpos=fpos;
        atom.cpos=cpos;
        atom.name=atomB;
        atom.type=0;
        if(species_tmp.at(0)!="") { atom.type++; } // B is 0 or 1
        atom.name_is_given=TRUE;
        volume+=volumeB;
        str.atoms.push_back(atom);
      }
    }
    // perform atom C
    for(iat=0;iat<natoms;iat++) {
      cpos.clear();fpos.clear(); // clear positions
      if(labeling[iat]=='2') { // atoms C
        _atom atom;
        atom.CleanName();
        atom.CleanSpin();
        fpos(1)=aBas(1,iat+1);
        fpos(2)=aBas(2,iat+1);
        fpos(3)=aBas(3,iat+1);
        cpos=F2C(str.lattice,fpos);
        atom.fpos=fpos;
        atom.cpos=cpos;
        atom.name=atomC;
        atom.type=0;
        if(species_tmp.at(0)!="") { atom.type++; } // C is 0 or 1 or 2
        if(species_tmp.at(1)!="") { atom.type++; } // C is 0 or 1 or 2
        atom.name_is_given=TRUE;
        volume+=volumeB;
        str.atoms.push_back(atom);
      }
    }
    if(LDEBUG) { cerr << "DEBUG volumeA=" << volumeA << endl; }
    if(LDEBUG) { cerr << "DEBUG volumeB=" << volumeB << endl; }
    if(LDEBUG) { cerr << "DEBUG volumeC=" << volumeC << endl; }
    if(LDEBUG) { cerr << "DEBUG volume=" << volume << endl; }
    // FIX scale
    if(volume_in>0.0) { volume=natoms*volume_in; }
    //[CO20190205 - OBSOLETE]str.scale=std::pow((double) (abs(volume)/det(str.lattice)),(double) 1.0/3.0);
    str.SetVolume(volume);  //CO20190205 - more robust
    str.neg_scale=TRUE;
    // str.SetCoordinates(_COORDS_CARTESIAN_);
    str.MinkowskiBasisReduction();   // BY DEFINITION MAKE THEM MINKOSKWIAN  GUS DISCUSSION = [Mon Apr  8 11:04:47 EDT 2013]
    str.iomode=IOVASP_ABCCAR; // put in ABCCAR GUS DISCUSSION = [Mon Apr  8 11:04:47 EDT 2013]
    str.BringInCell();
    str.FixLattices();
    // make species
    for(uint i=0;i<species_tmp.size();i++) {
      if(species_tmp.at(i)!="") {
        str.species.push_back(species_tmp.at(i));
        str.species_pp.push_back(species_tmp.at(i));
        str.species_pp_type.push_back("");
        str.species_pp_version.push_back("");
        str.species_pp_ZVAL.push_back(0.0);
        str.species_volume.push_back(0.0);
      }
    }
    return str;
  }
} // namespace aflowlib

#endif  // _AFLOW_XPROTO_GUS_CPP_
// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************
