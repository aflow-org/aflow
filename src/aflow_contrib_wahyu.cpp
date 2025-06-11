// ***************************************************************************
// *                                                                         *
// *              AFlow WAHYU SETYAWAN - Duke University 2007-2011           *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_wahyu.cpp
// functions written by
// 2007-2011: wahyu@alumni.duke.edu
// Fixed for xz by SC - 2018


#ifndef _AFLOW_CONTRIB_WAHYU_CPP_
#define _AFLOW_CONTRIB_WAHYU_CPP_
#define tiny6 1e-6
#define DEFAULT_METAL_GAP_TOLERANCE 0.001  // < 1 meV for metal gap  used to be 0.026

// ***************************************************************************
#include "aflow.h"
#include "aflow_pflow.h"
// [OBSOLETE] #include "aflow_contrib_wahyu.h"

// ***************************************************************************
// AConvaspBandgap
// ***************************************************************************
void AConvaspBandgap(vector<string>& argv) {
  // tmp safe, checked SC+JX 0512
  //Calculate band gap from bands run.
  //example:
  //   aflow --bandgap=/common/DATA/ZnO/
  //in this case, the bands run directory is "/common/DATA/ZnO/"
  //  grep E-fermi from OUTCAR.bands.EXT
  //  determine the occupied/unocc bands in the EIGENVAL.bands.EXT
  //  Egap = CBM-VBM
  // in the bandsdir, there must be OUTCAR.bands.EXT and EIGENVAL.bands.EXT
  char gaptype;
  float Egap,Efermi;
  string directory,stmp,tag="";
  stringstream straus;
  bool found=FALSE;

  directory=argv.at(2); //directory
  //cerr << "directory: " << directory << endl;
  //OUTCAR.bands

  string file_tmp=aurostd::TmpFileCreate("AConvaspBandgap");
  aurostd::RemoveFile(file_tmp);

  if(XHOST.vext.size()!=XHOST.vcat.size()) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"XHOST.vext.size()!=XHOST.vcat.size(), aborting",_RUNTIME_ERROR_); //CO20200624
  }

  // OUTCAR.bands
  found=FALSE;
  for(uint iext=0;iext<XHOST.vext.size();iext++) {
    if(!found && aurostd::FileExist(directory+"/OUTCAR.bands"+XHOST.vext.at(iext))) {
      found=TRUE;
      aurostd::execute(XHOST.vcat.at(iext)+" "+directory+"/OUTCAR.bands"+XHOST.vext.at(iext)+" | grep E-fermi > "+file_tmp);
    }
  }
  if(!found) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"OUTCAR.bands[.EXT] not found in the directory, aborting",_FILE_CORRUPT_); //CO20200624
  }
  straus.clear();straus.str(std::string());
  aurostd::file2stringstream(file_tmp,straus);
  straus >> stmp >> stmp >> Efermi;
  aurostd::RemoveFile(file_tmp);

  // EIGENVAL.bands
  found=FALSE;
  for(uint iext=0;iext<XHOST.vext.size();iext++) {
    if(!found && aurostd::FileExist(directory+"/EIGENVAL.bands"+XHOST.vext.at(iext))) {
      found=TRUE;
      aurostd::execute(XHOST.vcat.at(iext)+" "+directory+"/EIGENVAL.bands"+XHOST.vext.at(iext)+" > "+file_tmp);
    }
  }
  if(!found) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"EIGENVAL.bands[.EXT] not found in the directory, aborting",_FILE_CORRUPT_); //CO20200624
  }
  straus.clear();straus.str(std::string());
  aurostd::file2stringstream(file_tmp,straus);
  Egap=-1.0;
  Egap=GetBandGap_WAHYU(straus,Efermi,gaptype); //straus is stringstream of EIGENVAL.bands
  aurostd::RemoveFile(file_tmp);
  //output to stdout
  cout << tag << " " << directory << " ";
  if(Egap<tiny6) cout << "0.0 (M)" << endl;
  else cout << Egap << " (" << gaptype << ")" << endl;
}
// ***************************************************************************
// AConvaspBandgaps
// ***************************************************************************
void AConvaspBandgaps(istream& bandsdir, ostringstream& oss);

void AConvaspBandgaps(istream& bandsdir, ostream& oss) {
  ostringstream tmposs;
  AConvaspBandgaps(bandsdir,tmposs);
  oss << tmposs.str();
}

void AConvaspBandgaps(istream& bandsdir, ostringstream& oss) {
  // tmp safe, checked SC+JX 0512
  //  void AConvaspBandgaps(istream& bandsdir)  //JX

  //Calculate band gap from bands run.
  //bandsdir contains directories of bands run, one directory per line
  //example:
  //   aflow --bandgaps < dirlist
  //in this case, bandsdir is istream of file dirlist containing e.g.
  //  "/common/DATA/ZnO/"
  //  "/common/DATA/ZnSe/"
  //   so on...
  //-grep E-fermi from OUTCAR.bands.EXT
  //-determine the occupied/unocc bands in the EIGENVAL.bands.EXT
  //-Eg=CBM-VBM
  //in the bandsdir, there must be OUTCAR.bands.EXT and EIGENVAL.bands.EXT
  char gaptype;
  float Egap,Efermi;
  string directory,stmp,tag="";
  stringstream straus;
  bool found=FALSE;

  if(XHOST.vext.size()!=XHOST.vcat.size()) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"XHOST.vext.size()!=XHOST.vcat.size(), aborting",_RUNTIME_ERROR_); //CO20200624
  }

  while(bandsdir.good()) {
    bandsdir >> directory; //directory
    //cerr << "directory: " << directory << endl;
    //OUTCAR.bands
    //Getting the fermi energy from OUTCAR.bands

    string file_tmp=aurostd::TmpFileCreate("AConvaspBandgap");

    // OUTCAR.bands
    found=FALSE;
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found && aurostd::FileExist(directory+"/OUTCAR.bands"+XHOST.vext.at(iext))) {
        found=TRUE;
        aurostd::execute(XHOST.vcat.at(iext)+" "+directory+"/OUTCAR.bands"+XHOST.vext.at(iext)+" | grep E-fermi > "+file_tmp);
      }
    }
    if(!found) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"OUTCAR.bands[.EXT] not found in the directory, aborting",_FILE_CORRUPT_); //CO20200624
    }
    straus.clear();straus.str(std::string());
    aurostd::file2stringstream(file_tmp,straus);
    straus >> stmp >> stmp >> Efermi;
    aurostd::RemoveFile(file_tmp);

    // EIGENVAL.bands
    found=FALSE;
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found && aurostd::FileExist(directory+"/EIGENVAL.bands"+XHOST.vext.at(iext))) {
        found=TRUE;
        aurostd::execute(XHOST.vcat.at(iext)+" "+directory+"/EIGENVAL.bands"+XHOST.vext.at(iext)+" > "+file_tmp);
      }
    }
    if(!found) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"EIGENVAL.bands[.EXT] not found in the directory, aborting",_FILE_CORRUPT_); //CO20200624
    }
    straus.clear();straus.str(std::string());
    aurostd::file2stringstream(file_tmp,straus);
    Egap=-1.0;
    Egap=GetBandGap_WAHYU(straus,Efermi,gaptype);
    aurostd::RemoveFile(file_tmp);
    //output to stdout
    oss << tag << " " << directory << " ";
    //cout << tag << " " << s << " ";   //JX
    if(Egap<tiny6) oss << "0.0 (M)" << endl;
    //if(Egap<tiny6) cout << "METAL" << endl; //JX
    else oss << Egap << " (" << gaptype << ")" << endl;
    //else cout << Egap << " (" << gaptype << ")" << endl; //JX
  }
}
// ***************************************************************************
// AConvaspBandgapFromDOS
// ***************************************************************************
void AConvaspBandgapFromDOS(istream& doscar) {
  //take a single DOSCAR and print out to stdout the band gap
  float f;
  f=GetBandgapFromDOS(doscar);
  if(f<tiny6) cout << "0.0 (M)" << endl;
  else cout << f << endl;
}
// ***************************************************************************
// AConvaspBandgapListFromDOS
// ***************************************************************************
void AConvaspBandgapListFromDOS(istream& doscar) {
  //take a list of DOSCAR files, one on each line
  //and print out to stdout a list of band gap
  float f;
  bool fgz=false,fbz=false,fbz2=false,fxz=false;
  size_t ibz,ibz2,igz,ixz,idoscar;
  string s,tag,stmp,scom,fext;
  s=tag=stmp=scom=fext="";
  ifstream dos;
  while(doscar.good()) {
    doscar >> s;
    fxz=fgz=fbz=fbz2=false;
    ixz=s.find(".xz");
    if(ixz!=s.npos)	{fxz=true; fbz2=false; fbz=false; fgz=false; fext=".xz";}
    igz=s.find(".gz");
    if(igz!=s.npos)	{fxz=false; fbz2=false; fbz=false; fgz=true; fext=".gz";}
    ibz=s.find(".bz");
    if(ibz!=s.npos)	{fxz=false; fbz2=false; fbz=true; fgz=false; fext=".bz";}
    ibz2=s.find(".bz2");
    if(ibz2!=s.npos)    {fxz=false; fbz2=true; fbz=false; fgz=false; fext=".bz2";}
    if(fgz || fbz || fbz2) {
      stmp="";
      idoscar=s.find("DOSCAR");
      stmp=s.substr(idoscar);
      if(stmp.size()>5) stmp="DOSCAR.tmp";
      aurostd::execute("cp "+s+" "+stmp+fext);
      if(fxz) aurostd::execute("xz -dfq "+stmp+fext);
      if(fgz) aurostd::execute("gzip -dfq "+stmp+fext);
      if(fbz) aurostd::execute("bzip2 -dfq "+stmp+fext);
      if(fbz2) aurostd::execute("bzip2 -dfq "+stmp+fext);
    }
    else stmp=s;
    dos.open(stmp.data());
    if(dos.good()) {
      f=GetBandgapFromDOS(dos);
      cout << tag << " " << s << " ";
      if(f<tiny6) cout << "0.0 (M)" << endl;
      else cout << f << endl;
      dos.close();
    }
    aurostd::execute("rm -f DOSCAR.tmp");
  }
}

// ***************************************************************************
// pflow::ICSD_
// ***************************************************************************
namespace pflow {
  void ICSD(vector<string> argv, istream& input) {
    //   Output to stdout in "ICSD-format" (the NIST's Inorganic Crystal Structure Database)
    //   all compounds containing particular elements according to the options:

    //   --icsd Pb              : extract compounds containing Pb
    //   --icsd Pb Sn Se        : extract compounds containing Pb, Sn, and Se
    //   we can also substitute elements with their atomic number
    //   --icsd_alllessthan Ra  : all elements in the compound MUST have Z<Z_Ra
    //   --icsd_allmorethan Ar
    //   --icsd_id 26675        : extract compound with icsd entry ID number 26675
    //   --icsd_morethan Co     : all compounds containing element with atomic number Z>Z_Co
    //   --icsd_morethan 27     : all compounds containing element with atomic number Z>27
    //   --icsd_lessthan z1     : ...
    //   --icsd_lessthan z1 --icsd_morethan z2 : ..

    //   --icsd_densmorethan 4.5   : all compounds with density > 4.5
    //   --icsd_sg 62              : all compounds with SG number 62
    //   --icsd_sglessthan 62      : all compounds with SG number < 62
    //   --icsd_cubic              : all compounds that belong to cubic system

    //   aflow --icsd_tri < input.icsd \n		
    //   aflow --icsd_mcl < input.icsd \n		
    //   aflow --icsd_mclc < input.icsd \n	
    //   aflow --icsd_orc < input.icsd \n		
    //   aflow --icsd_orcc < input.icsd \n	
    //   aflow --icsd_orcf < input.icsd \n	
    //   aflow --icsd_orci < input.icsd \n	
    //   aflow --icsd_tet < input.icsd \n		
    //   aflow --icsd_bct < input.icsd \n		
    //   aflow --icsd_rhl < input.icsd \n		
    //   aflow --icsd_hex < input.icsd \n		
    //   aflow --icsd_cub < input.icsd \n		
    //   aflow --icsd_fcc < input.icsd \n		
    //   aflow --icsd_bcc < input.icsd \n		

    //   --icsd_basislessthan #basis : all compounds with number of basis atoms in PRIMITIVE cell < #basis
    //   --icsd_basismorethan #basis : all compounds with number of basis atoms in PRIMITIVE cell > #basis
    //   --icsd_n_ary #species       : (#species=1: unary, #species=2: binary, so on)
    //   --icsd_nopartialocc   : discard compounds with partial occupancies
    //   --icsd_nobrokenbasis  : select only compounds with complete basis is chosen (some times there is a broken basis, meaning
    //   that some elements are missing from the basis in the database)
    //   --icsd_unique         : discard redundant compounds having the same chemical formula AND SG#, the first
    //   compound with complete basis is chosen (some times there is a broken basis, meaning
    //   that some elements are missing from the basis in the database)

    //   --icsd_chem  MgB4     : extract all MgB4 structures
    //   --icsd_makelabel	  : make label from icsd format
    //   --icsd_remove_and A B : remove compounds containing A AND B elements
    //   --icsd_remove_or A B  : remove compounds containing A OR B elements
    //   --icsd_removemetals   : remove compounds composed ENTIRELY by metallic elements only.
    //   --icsd_proto 1 2 7    : all compounds with prototype AB2C7

    //   The input is from standard istream input in "ICSD-format". Note that the heaviest
    //   element supported in this version is Thorium (Z=90).

    //  bool  ifounbd,ifound1, ifound2, ifound3,ifound4,
    bool iwrite,iread,ifoundstart,ifoundend,AtomLineFound,
         ALLLESSTHAN,ALLMORETHAN,
         ICSD,ICSD_ID,ICSD_LESSTHAN,ICSD_MORETHAN,DENSLESSTHAN,DENSMORETHAN,SGLESSTHAN,SGMORETHAN,
         SPACEGROUP,TRICLINIC,MONOCLINIC,ORTHORHOMBIC,TETRAGONAL,RHOMBOHEDRAL,HEXAGONAL,CUBIC,
         UNIQUE,N_ARY,NOBROKENBASIS,BASISLT,BASISGT,CHEM,PROTO,NOPARTIALOCC,ICSD_REMOVE_AND,ICSD_REMOVE_OR,REMOVEMETALS;
    bool MAKELABEL;
    bool TRI,MCL,MCLC,ORC,ORCC,ORCF,ORCI,  TET,BCT,RHL,ICSDHEX,CUB,FCC,BCC;
    uint iargv,Zi,zmin,zmax,SGmin,SGmax,SGref,SGnumber;
    int i,j,iline,ilessthan,imorethan,iscan,Nspec,Nspecref,Zformula; //,ICSDid;
    float ftmp,Vol,Volred;
    float Nbasis,Nbasismin,Nbasismax;
    float Dmin,Dmax,Dcalc;
    string sline,sword,ssub,sdata,Pearson,ChemFormula,title,sICSDid,SG,stmp;
    string::iterator siter;
    vector<bool> foundbasis,iwritelist;
    vector<string> sbuf,Zname,Znameref,ChemUnique,ChemNameRef,WyckLabel,IDref;
    vector<int> Zcount;
    vector<float> Zconc,ChemConcRef,xpos,ypos,zpos,vsof;
    vector<uint> Zref,Ztotest,Znumber,SGunique,ChemZRef;

    //  ofstream dbg;
    //  dbg.open("debug.out");

    const float TINY6=1e-6;

    Dcalc=0.0;
    iargv=0;
    Dmin=100.0;Dmax=0.00;
    ilessthan=imorethan=0;
    Nbasis=Nbasismin=Nbasismax=0.0;
    Nspec=Nspecref=0;
    //  ICSDid=0;
    iscan=0;

    // ifound=  ifound1=  ifound2=  ifound3=ifound4=
    iwrite=iread=ifoundstart=ifoundend=AtomLineFound=false;
    ALLLESSTHAN=ALLMORETHAN=BASISLT=BASISGT=CHEM=DENSLESSTHAN=DENSMORETHAN=
      ICSD=ICSD_ID=ICSD_LESSTHAN=ICSD_MORETHAN=ICSD_REMOVE_AND=ICSD_REMOVE_OR=N_ARY=NOBROKENBASIS=NOPARTIALOCC=
      PROTO=REMOVEMETALS=SGLESSTHAN=SGMORETHAN=UNIQUE=false;
    SPACEGROUP=TRICLINIC=MONOCLINIC=ORTHORHOMBIC=TETRAGONAL=RHOMBOHEDRAL=HEXAGONAL=CUBIC=false;
    TRI=MCL=MCLC=ORC=ORCC=ORCF=ORCI=TET=BCT=RHL=ICSDHEX=CUB=FCC=BCC=false;
    MAKELABEL=false;

    SGmin=231; SGmax=SGref=SGnumber=0;
    SGunique.clear(); ChemUnique.clear(); Znameref.clear();
    xpos.clear();ypos.clear();vsof.clear();
    sbuf.clear();
    for(i=1;i<(int) argv.size();i++) {
      sdata=argv.at(i);
      if(sdata=="--icsd") {ICSD=true;iargv=i+1;}
      if(sdata=="--icsd_alllessthan") {ALLLESSTHAN=true;ilessthan=i+1;}
      if(sdata=="--icsd_allmorethan") {ALLMORETHAN=true;imorethan=i+1;}
      if(sdata=="--icsd_id" || sdata=="--icsd_ID") {ICSD_ID=true;iargv=i+1;}
      if(sdata=="--icsd_lessthan") {ICSD_LESSTHAN=true;ilessthan=i+1;}
      if(sdata=="--icsd_morethan") {ICSD_MORETHAN=true;imorethan=i+1;}
      if(sdata=="--icsd_denslessthan") {DENSLESSTHAN=true;Dmax=aurostd::string2utype<float>(argv.at(i+1));}
      if(sdata=="--icsd_densmorethan") {DENSMORETHAN=true;Dmin=aurostd::string2utype<float>(argv.at(i+1));}
      if(sdata=="--icsd_sg") {SPACEGROUP=true;SGref=aurostd::string2utype<uint>(argv.at(i+1));}
      if(sdata=="--icsd_sglessthan") {SPACEGROUP=true;SGLESSTHAN=true;SGmax=aurostd::string2utype<uint>(argv.at(i+1));}
      if(sdata=="--icsd_sgmorethan") {SPACEGROUP=true;SGMORETHAN=true;SGmin=aurostd::string2utype<uint>(argv.at(i+1));}
      if(sdata=="--icsd_triclinic") {SPACEGROUP=true;TRICLINIC=true;}
      if(sdata=="--icsd_monoclinic") {SPACEGROUP=true;MONOCLINIC=true;}
      if(sdata=="--icsd_orthorhombic") {SPACEGROUP=true;ORTHORHOMBIC=true;}
      if(sdata=="--icsd_tetragonal") {SPACEGROUP=true;TETRAGONAL=true;}
      if(sdata=="--icsd_rhombohedral" || sdata=="--icsd_trigonal") {SPACEGROUP=true;RHOMBOHEDRAL=true;}
      if(sdata=="--icsd_hexagonal") {SPACEGROUP=true;HEXAGONAL=true;}
      if(sdata=="--icsd_cubic") {SPACEGROUP=true;CUBIC=true;}
      if(sdata=="--icsd_tri" || sdata=="--icsd_TRI") {SPACEGROUP=true;TRI=true;}
      if(sdata=="--icsd_mcl" || sdata=="--icsd_MCL") {SPACEGROUP=true;MCL=true;}
      if(sdata=="--icsd_mclc" || sdata=="--icsd_MCLC") {SPACEGROUP=true;MCLC=true;}
      if(sdata=="--icsd_orc" || sdata=="--icsd_ORC") {SPACEGROUP=true;ORC=true;}
      if(sdata=="--icsd_orcc" || sdata=="--icsd_ORCC") {SPACEGROUP=true;ORCC=true;}
      if(sdata=="--icsd_orcf" || sdata=="--icsd_ORCF") {SPACEGROUP=true;ORCF=true;}
      if(sdata=="--icsd_orci" || sdata=="--icsd_ORCI") {SPACEGROUP=true;ORCI=true;}
      if(sdata=="--icsd_tet" || sdata=="--icsd_TET") {SPACEGROUP=true;TET=true;}
      if(sdata=="--icsd_bct" || sdata=="--icsd_BCT") {SPACEGROUP=true;BCT=true;}
      if(sdata=="--icsd_rhl" || sdata=="--icsd_RHL") {SPACEGROUP=true;RHL=true;}
      if(sdata=="--icsd_hex" || sdata=="--icsd_HEX") {SPACEGROUP=true;ICSDHEX=true;}
      if(sdata=="--icsd_cub" || sdata=="--icsd_CUB") {SPACEGROUP=true;CUB=true;}
      if(sdata=="--icsd_fcc" || sdata=="--icsd_FCC") {SPACEGROUP=true;FCC=true;}
      if(sdata=="--icsd_bcc" || sdata=="--icsd_BCC") {SPACEGROUP=true;BCC=true;}
      if(sdata=="--icsd_unique") {UNIQUE=true;}
      if(sdata=="--icsd_n_ary") {N_ARY=true;Nspecref=aurostd::string2utype<int>(argv.at(i+1));}
      if(sdata=="--icsd_basislessthan") {BASISLT=true;Nbasismax=aurostd::string2utype<float>(argv.at(i+1));}
      if(sdata=="--icsd_basismorethan") {BASISGT=true;Nbasismin=aurostd::string2utype<float>(argv.at(i+1));}
      if(sdata=="--icsd_chem") {CHEM=true;iargv=i+1;}
      if(sdata=="--icsd_makelabel") {MAKELABEL=true;}
      if(sdata=="--icsd_proto") {PROTO=true;iargv=i+1;}
      if(sdata=="--icsd_nobrokenbasis" || sdata=="--icsd_completebasis") {NOBROKENBASIS=true;}
      if(sdata=="--icsd_nopartialocc") {NOPARTIALOCC=true;}
      if(sdata=="--icsd_remove_and") {ICSD_REMOVE_AND=true;iargv=i+1;}
      if(sdata=="--icsd_remove_or") {ICSD_REMOVE_OR=true;iargv=i+1;}
      if(sdata=="--icsd_removemetals") {REMOVEMETALS=true;}
    }

    //----------Preparation---------------
    if(ICSD || ICSD_REMOVE_AND || ICSD_REMOVE_OR) {
      Zref.clear();iwritelist.clear();
      for(i=iargv;i<(int) argv.size();i++) {
        sdata=argv.at(i);
        if(sdata[0]<'A') Zref.push_back(aurostd::string2utype<uint>(sdata));
        else Zref.push_back(GetAtomNumber(sdata));
      }    
      //sorting
      for(i=0;i<((int) Zref.size())-1;i++) {
        for(j=i+1;j<(int) Zref.size();j++) {
          if(Zref[j]<Zref[i]) {Zi=Zref[i]; Zref[i]=Zref[j]; Zref[j]=Zi;}
        }
      }
      iwritelist.resize((int)Zref.size());
    }
    if(ICSD_ID) {
      IDref.clear();
      for(i=iargv;i<(int) argv.size();i++) {
        IDref.push_back(argv.at(i));
      }
    }
    if(CHEM) {
      ChemNameRef.clear(); ChemZRef.clear(); ChemConcRef.clear();
      //copying the argument, note that white characters are removed.
      sline="";
      for(i=iargv;i<(int) argv.size();i++) sline=sline+argv.at(i);
      ParseChemicalFormula(sline,ChemNameRef,ChemConcRef);//e.g. sline=MgB2.3 -> ChemNameRef=["Mg","B"] and ChemConcRef=[1,2.3]
      for(i=0;i<(int)ChemNameRef.size();i++) ChemZRef.push_back(GetAtomNumber(ChemNameRef[i]));
      //sorting based on ascending Z
      Zi=0; ftmp=0.0; sword="";
      for(i=0;i<(int)ChemZRef.size()-1;i++) {
        for(j=i+1;j<(int)ChemZRef.size();j++) {
          if(ChemZRef[j]<ChemZRef[i]) {
            Zi=ChemZRef[i];ChemZRef[i]=ChemZRef[j];ChemZRef[j]=Zi;
            ftmp=ChemConcRef[i];ChemConcRef[i]=ChemConcRef[j];ChemConcRef[j]=ftmp;
            stmp=ChemNameRef[i];ChemNameRef[i]=ChemNameRef[j];ChemNameRef[j]=stmp;
          }
        }
      }
    }
    if(PROTO) {
      ChemConcRef.clear();
      for(i=iargv;i<(int) argv.size();i++) ChemConcRef.push_back(aurostd::string2utype<float>(argv.at(i)));
      //sorting ascending
      ChemConcRef=SortFloat(ChemConcRef,1);
    }
    zmin=0;
    if(ICSD_MORETHAN || ALLMORETHAN) {
      sdata=argv[imorethan];
      if(sdata[0]<'A') zmin=aurostd::string2utype<uint>(sdata);
      else zmin=GetAtomNumber(sdata);
    }
    zmax=0;
    if(ICSD_LESSTHAN || ALLLESSTHAN) {
      sdata=argv[ilessthan];
      if(sdata[0]<'A') zmax=aurostd::string2utype<uint>(sdata);
      else zmax=GetAtomNumber(sdata);
    }
    //---------------checking consistency among options-------------
    cerr << "Screening based on: ";
    if(ICSD) cerr << "ICSD";
    if(ALLLESSTHAN) cerr << "ALLLESSTHAN";
    if(ALLMORETHAN) cerr << "ALLMORETHAN";
    if(ICSD_ID) cerr << "ICSD_ID";
    if(ICSD_LESSTHAN) cerr << "ICSD_LESSTHAN";
    if(ICSD_MORETHAN) cerr << "ICSD_MORETHAN";
    if(MAKELABEL) cerr << "ICSD MAKE LABEL";
    if(DENSLESSTHAN) cerr << "DENSLESSTHAN";
    if(DENSMORETHAN) cerr << "DENSMORETHAN";
    if(SPACEGROUP) cerr << "SPACEGROUP";
    if(SGLESSTHAN) cerr << "SGLESSTHAN";
    if(SGMORETHAN) cerr << "SGMORETHAN";
    if(TRICLINIC) cerr << "TRICLINIC";
    if(MONOCLINIC) cerr << "MONOCLINIC";
    if(ORTHORHOMBIC) cerr << "ORTHORHOMBIC";
    if(TETRAGONAL) cerr << "TETRAGONAL";
    if(RHOMBOHEDRAL) cerr << "RHOMBOHEDRAL";
    if(HEXAGONAL) cerr << "HEXAGONAL";
    if(CUBIC) cerr << "CUBIC";
    if(TRI) cerr << "TRI";
    if(MCL) cerr << "MCL";
    if(MCLC) cerr << "MCLC";
    if(ORC) cerr << "ORC";
    if(ORCC) cerr << "ORCC";
    if(ORCF) cerr << "ORCF";
    if(ORCI) cerr << "ORCI";
    if(TET) cerr << "TET";
    if(BCT) cerr << "BCT";
    if(RHL) cerr << "RHL";
    if(ICSDHEX) cerr << "ICSDHEX";
    if(CUB) cerr << "CUB";
    if(FCC) cerr << "FCC";
    if(BCC) cerr << "BCC";
    if(UNIQUE) cerr << "UNIQUE";
    if(N_ARY) cerr << "N_ARY";
    if(BASISLT) cerr << "BASISLT";
    if(BASISGT) cerr << "BASISGT";
    if(CHEM) cerr << "CHEM";
    if(PROTO) cerr << "PROTO";
    if(NOBROKENBASIS) cerr << "NOBROKENBASIS";
    if(NOPARTIALOCC) cerr << "NOPARTIALOCC";
    if(ICSD_REMOVE_AND) cerr << "ICSD_REMOVE_AND";
    if(ICSD_REMOVE_OR) cerr << "ICSD_REMOVE_OR";
    if(REMOVEMETALS) cerr << "REMOVEMETALS";
    cerr << endl;
    //------------------------------------------------------
    iwrite=false;
    //  ifound=ifound1= ifound2=ifound3=ifound4=false;
    Zformula=0;Vol=Volred=0.0;
    ifoundstart=ifoundend=false;
    cerr << "scanning" << endl;
    while(!input.eof()) {
      getline(input,sline);

      if(sline.length()>4) sword=sline.substr(0,5);
      if(sword=="*data") {
        ifoundstart=true; ifoundend=false; sbuf.clear();
      }
      sbuf.push_back(sline);
      if(sline.length()>3) sword=sline.substr(0,4);
      if(sword=="*end") {
        ifoundend=true;
      }
      if(ifoundstart && ifoundend) {//valid block of data found
        AtomLineFound=FALSE;
        iscan=iscan+1;
        xpos.clear(); ypos.clear(); zpos.clear();
        vsof.clear(); foundbasis.clear();
        Zcount.clear(); Znumber.clear(); Zname.clear(); Zconc.clear(); WyckLabel.clear();
        //extract the key words data
        for(iline=0;iline<(int)sbuf.size();iline++) {
          //----------Coll Code---------
          sline=sbuf[iline];
          if(sline.length()>8) sword=sline.substr(0,9);
          if(sword=="Coll Code") {
            sICSDid=sline.substr(9,sline.length());
            sICSDid=RemoveCharFromString(sICSDid,' ');
            sICSDid=RemoveCharFromString(sICSDid,'\t');
            sICSDid=RemoveCharFromString(sICSDid,'\n');
            //	  ICSDid=aurostd::string2utype<int>(sICSDid);
          }
          //----------Rec Date----------
          //----------Mod Date----------
          //----------Chem Name---------
          //----------Structured--------
          //----------Sum---------------
          if(sline.length()>3) sword=sline.substr(0,4);
          if(sword=="Sum ") {
            ChemFormula=sline.substr(5,sline.length());
            ChemFormula=RemoveCharFromString(ChemFormula,' ');
            ChemFormula=RemoveCharFromString(ChemFormula,'\t');
            ChemFormula=RemoveCharFromString(ChemFormula,'\n');
            ParseChemicalFormula(ChemFormula,Zname,Zconc);//e.g. sline=MgB2.3 -> ChemNameRef=["Mg","B"] and ChemConcRef=[1,2.3]
            Nspec=Zname.size();Znumber.resize(Nspec);
            for(i=0;i<Nspec;i++) Znumber[i]=GetAtomNumber(Zname[i]);
            //sorting ascending
            for(i=0;i<Nspec-1;i++) {
              for(j=i+1;j<Nspec;j++) {
                if(Znumber[j]<Znumber[i]) {
                  Zi=Znumber[i];Znumber[i]=Znumber[j];Znumber[j]=Zi;
                  stmp=Zname[i];Zname[i]=Zname[j];Zname[j]=stmp;
                  ftmp=Zconc[i];Zconc[i]=Zconc[j];Zconc[j]=ftmp;
                }
              }
            }
          }
          //----------ANX---------------
          //----------D(calc)-----------
          if(sline.length()>6) sword=sline.substr(0,7);
          if(sword=="D(calc)") {
            ssub=sline.substr(7,sline.length());
            Dcalc=aurostd::string2utype<float>(ssub);
          }
          //----------Title-------------
          //----------Author(s)---------
          //----------Reference---------
          //----------Unit Cell---------
          //----------Vol---------------
          if(sline.length()>2) sword=sline.substr(0,3);
          if(sword=="Vol") {
            Vol=aurostd::string2utype<float>(sline.substr(3,sline.length()));
          }
          //----------Z-----------------
          if(sline.length()>1) sword=sline.substr(0,2);
          if(sword=="Z ") {
            Zformula=aurostd::string2utype<int>(sline.substr(2,sline.length()));
          }
          //----------Space Group-------
          if(sline.length()>10) sword=sline.substr(0,11);
          if(sword=="Space Group") {
            SG=sline.substr(12,sline.length());
            SG=RemoveCharFromString(SG,' ');
            SG=RemoveCharFromString(SG,'\t');
            SG=RemoveCharFromString(SG,'\n');
          }
          //----------SG Number---------
          if(sline.length()>8) sword=sline.substr(0,9);
          if(sword=="SG Number") {
            ssub=sline.substr(9,sline.length());
            SGnumber=aurostd::string2utype<uint>(ssub);
          }
          //----------Cryst Sys---------
          //----------Pearson-----------
          if(sline.length()>6) sword=sline.substr(0,7);
          if(sword=="Pearson") {
            Pearson=sline.substr(8,sline.length());
            Pearson=RemoveCharFromString(Pearson,' ');
            Pearson=RemoveCharFromString(Pearson,'\t');
            Pearson=RemoveCharFromString(Pearson,'\n');
          }
          //----------Wyckoff-----------
          //----------R Value-----------
          //----------Red Cell----------
          if(sline.length()>7) sword=sline.substr(0,8);
          if(sword=="Red Cell") {
            for(i=(int)sline.length();i>0;i--) {
              if(sline[i]==' ' || sline[i]=='\t') break;
            }
            Volred=aurostd::string2utype<float>(sline.substr(i,sline.length()));
          }
          //----------Trans Red---------
          //----------Comments----------
          //----------Atom--------------
          if(sline.length()>4) sword=sline.substr(0,5);
          //icsd 2010 version
          //if(sword=="Atom")
          //icsd 2011 version
          if(sword=="Atom " && !AtomLineFound) {
            AtomLineFound=TRUE;
            iread=TRUE;
            Zcount.resize(Nspec);foundbasis.resize(Nspec);
            for(j=0;j<Nspec;j++) {Zcount[j]=0;foundbasis[j]=false;}//clearing Zcount and foundbasis
            while(iread) {//processing one wyckoff entry
              iread=false;
              iline=iline+1;
              sline=sbuf[iline];
              ssub=StringCropAt(sline,1);
              for(j=0;j<Nspec;j++) {
                if(ssub==Zname[j]) {//if it is atom name
                  foundbasis[j]=true;//basis is found
                  iread=TRUE;
                  WyckLabel.push_back(ssub);
                  Zcount[j]=Zcount[j]+1;
                  ssub=StringCropAt(sline,9);
                  ssub=RemoveStringFromTo(ssub,'(',')');
                  //ssub=RemoveCharFromString(ssub,'(');
                  //ssub=RemoveCharFromString(ssub,')');
                  vsof.push_back(atof(ssub.data()));//save sof
                }//if Zname
              }//for i<Nspec
            }//iread
            iline=iline-1;
          }//Atom
        }//iline
        title=ChemFormula+"_ICSD_"+sICSDid;
        if((iscan%50)==49) cerr << endl;
        cerr << ".";// << title << endl;
        //-------START SCREENING--------------------------------------
        iwrite=false;
        if(ALLLESSTHAN || ALLMORETHAN) {
          iwrite=true;
          for(i=0;i<Nspec;i++) {
            Zi=GetAtomNumber(Zname[i]);
            if(ALLLESSTHAN && !ALLMORETHAN)
              if(Zi>=zmax || Zi==0) iwrite=false;
            if(!ALLLESSTHAN && ALLMORETHAN)
              if(Zi<=zmin) iwrite=false;
            if(ALLLESSTHAN && ALLMORETHAN)
              if(Zi<=zmin && Zi>=zmax) iwrite=false;
          }
        } //if alllessthan allmorethan
        if(ICSD_LESSTHAN || ICSD_MORETHAN) {
          iwrite=false;
          for(i=0;i<Nspec;i++) {
            Zi=GetAtomNumber(Zname[i]);
            if(ICSD_LESSTHAN && !ICSD_MORETHAN)
              if(Zi<zmax) iwrite=TRUE;
            if(!ICSD_LESSTHAN && ICSD_MORETHAN)
              if(Zi>zmin) iwrite=TRUE;
            if(ICSD_LESSTHAN && ICSD_MORETHAN)
              if(Zi>zmin && Zi<zmax) iwrite=TRUE;
          }//i<Nspec
        }//LESSTHAN or MORETHAN
        if(BASISLT || BASISGT) {
          Nbasis=0.0;
          for(j=0;j<(int) Zconc.size();j++) Nbasis=Nbasis+Zconc[j];
          Nbasis=Nbasis*Zformula/(floor(Vol/Volred+0.5));//number of atoms in the primitive cell
          if(BASISLT && !BASISGT) {if(Nbasis<Nbasismax) iwrite=true;else iwrite=false;}
          if(BASISGT && !BASISLT) {if(Nbasis>Nbasismin) iwrite=true;else iwrite=false;}
          if(BASISLT && BASISGT) {if(Nbasis>Nbasismin && Nbasis<Nbasismax) iwrite=true;else iwrite=false;}
        }//BASISLT or BASISGT
        if(CHEM) {
          iwrite=false;
          if((int)ChemConcRef.size()==Nspec) {
            iwrite=true;
            for(i=0;i<Nspec;i++) {
              if( (aurostd::abs(Zconc[i]-ChemConcRef[i])>TINY6) || (Znumber[i]!=ChemZRef[i])) iwrite=false;
            }
          }//if same number of species
        }//CHEM	
        if(DENSLESSTHAN || DENSMORETHAN) {
          iwrite=false;
          if(DENSMORETHAN && !DENSLESSTHAN) {if(Dcalc>Dmin) iwrite=true;}
          if(!DENSMORETHAN && DENSLESSTHAN) {if(Dcalc<Dmax) iwrite=true;}
          if(DENSMORETHAN && DENSLESSTHAN) {if(Dcalc>Dmin && Dcalc<Dmax) iwrite=true;}
        }
        if(ICSD) {
          Ztotest.clear();
          for(i=0;i<Nspec;i++) Ztotest.push_back(GetAtomNumber(Zname[i]));
          iwrite=true;
          for(i=0;i<(int)Ztotest.size();i++) {
            iwritelist[i]=false;
            for(j=0;j<Nspec;j++) {
              if(Ztotest[i]==Zref[j]) iwritelist[i]=true;
            }
            if(iwritelist[i]==false) iwrite=false;
          }

        }//if ICSD
        if(ICSD_ID) {
          iwrite=false;
          for(i=0;i<(int)IDref.size();i++) {
            if(sICSDid==IDref[i]) iwrite=true;
          }
        }
        if(N_ARY) {
          if(Nspec==Nspecref) iwrite=true;
          else iwrite=false;
        }
        if(NOBROKENBASIS) {
          iwrite=true;
          for(j=0;j<(int) foundbasis.size();j++) {
            if(foundbasis[j]==false) iwrite=false; //basis not found
          }
        }
        if(NOPARTIALOCC) {
          iwrite=true;
          for(j=0;j<(int) vsof.size();j++) {
            if(!(abs(vsof[j]-1.00)<1e-3)) {
              iwrite=false;//if partial occupancy is detected
            }
          }
        }
        if(ICSD_REMOVE_AND) {
          for(i=0;i<(int)Zref.size();i++) {
            iwritelist[i]=true;
            for(j=0;j<Nspec;j++)
              if(GetAtomNumber(Zname[j])==Zref[i]) iwritelist[i]=false;
          }
          iwrite=false;
          for(i=0;i<(int)Zref.size();i++) {
            if(iwritelist[i]) iwrite=true;
          }
        }	
        if(ICSD_REMOVE_OR) {
          iwrite=true;
          for(i=0;i<(int)Zref.size();i++) {
            for(j=0;j<Nspec;j++)
              if(GetAtomNumber(Zname[j])==Zref[i]) iwrite=false;
          }
        }
        if(REMOVEMETALS) {
          iwrite=false;
          for(i=0;i<(int)Zname.size();i++)
            if(!IsMetal(Zname[i])) iwrite=true;
        }
        if(PROTO) {
          iwrite=false;
          Zconc=SortFloat(Zconc,1);
          if(Zconc.size()==ChemConcRef.size()) {
            iwrite=true;
            for(i=0;i<(int)ChemConcRef.size();i++) {
              if( aurostd::abs(Zconc[i]-ChemConcRef[i]) > TINY6 ) iwrite=false;
            }
          }
        }//PROTO
        if(SPACEGROUP) {
          iwrite=false;
          if(TRICLINIC) {if(SGnumber>0 && SGnumber<3) iwrite=true;}
          if(MONOCLINIC) {if(SGnumber>2 && SGnumber<16) iwrite=true;}
          if(ORTHORHOMBIC) {if(SGnumber>15 && SGnumber<75) iwrite=true;}
          if(TETRAGONAL) {if(SGnumber>74 && SGnumber<143) iwrite=true;}
          if(RHOMBOHEDRAL) {if(SGnumber>142 && SGnumber<168) iwrite=true;}
          if(HEXAGONAL) {if(SGnumber>167 && SGnumber<195) iwrite=true;}
          if(CUBIC) {if(SGnumber>194 && SGnumber<231) iwrite=true;}
          if(SGnumber==SGref) iwrite=true;
          if(SGLESSTHAN && !SGMORETHAN) {if(SGnumber<SGmax) iwrite=true;}
          if(!SGLESSTHAN && SGMORETHAN) {if(SGnumber>SGmin) iwrite=true;}
          if(SGLESSTHAN && SGMORETHAN) {if(SGnumber>SGmin && SGnumber<SGmax) iwrite=true;}

          cerr << "DEBUG: calling LATTICE::SpaceGroup2Lattice from pflow::ICSD_" << endl;
          if(TRI) {if("TRI"==LATTICE::SpaceGroup2Lattice((uint)SGnumber)) iwrite=true;}
          if(MCL) {if("MCL"==LATTICE::SpaceGroup2Lattice((uint)SGnumber)) iwrite=true;}
          if(MCLC) {if("MCLC"==LATTICE::SpaceGroup2Lattice((uint)SGnumber)) iwrite=true;}
          if(ORC) {if("ORC"==LATTICE::SpaceGroup2Lattice((uint)SGnumber)) iwrite=true;}
          if(ORCC) {if("ORCC"==LATTICE::SpaceGroup2Lattice((uint)SGnumber)) iwrite=true;}
          if(ORCF) {if("ORCF"==LATTICE::SpaceGroup2Lattice((uint)SGnumber)) iwrite=true;}
          if(ORCI) {if("ORCI"==LATTICE::SpaceGroup2Lattice((uint)SGnumber)) iwrite=true;}
          if(TET) {if("TET"==LATTICE::SpaceGroup2Lattice((uint)SGnumber)) iwrite=true;}
          if(BCT) {if("BCT"==LATTICE::SpaceGroup2Lattice((uint)SGnumber)) iwrite=true;}
          if(RHL) {if("RHL"==LATTICE::SpaceGroup2Lattice((uint)SGnumber)) iwrite=true;}
          if(ICSDHEX) {if("HEX"==LATTICE::SpaceGroup2Lattice((uint)SGnumber)) iwrite=true;}
          if(CUB) {if("CUB"==LATTICE::SpaceGroup2Lattice((uint)SGnumber)) iwrite=true;}
          if(FCC) {if("FCC"==LATTICE::SpaceGroup2Lattice((uint)SGnumber)) iwrite=true;}
          if(BCC) {if("BCC"==LATTICE::SpaceGroup2Lattice((uint)SGnumber)) iwrite=true;}
        }//if SPACEGROUP
        if(UNIQUE) {
          iwrite=true;
          for(j=0;j<(int) foundbasis.size();j++) {
            if(foundbasis[j]==false) iwrite=false; //basis not found
          }
          if(SGunique.size()>0) {
            for(i=0;i<(int)SGunique.size();i++) {
              if(SGnumber==SGunique[i] && ChemFormula==ChemUnique[i]) iwrite=false;
            }
          }
          if(iwrite) {
            SGunique.push_back(SGnumber);ChemUnique.push_back(ChemFormula);
          }
        }//if UNIQUE

        //---------------------OUTPUT---------------------------------------------
        if(MAKELABEL) {	  
          iwrite=false;
          //ascending sorting of elements in the compound
          for(j=0;j<(int)Zname.size()-1;j++) {
            for(i=j+1;i<(int)Zname.size();i++) {
              if(Zname[i]<Zname[j]) {//swap
                stmp=Zname[j];Zname[j]=Zname[i];Zname[i]=stmp;
                ftmp=Zconc[j];Zconc[j]=Zconc[i];Zconc[i]=ftmp;
              }
            }
          }
          for(j=0;j<(int)Zname.size();j++) cout << Zname[j] << Zconc[j];
          cout << "_ICSD_" << sICSDid << endl;
        }//if makelabel

        if(iwrite) {
          for(i=0;i<(int)sbuf.size();i++) cout << sbuf[i] << endl;
        }
        //clearing for the next block of data
        ifoundstart=ifoundend=false;
        sword="";sline="";
      }//one valid block of data
    }//while reading input
    cerr << endl;
    //dbg.close();
  }
}

// ***************************************************************************
// pflow::ICSD_CheckRaw
// ***************************************************************************
namespace pflow {
  void ICSD_CheckRaw(vector<string> argv) {
    //Synopsis:
    //aflow --icsd_check_raw > output.dat
    //same as --icsd_check_raw 0
    //aflow --icsd_check_raw 0 > output.dat
    //generate a new list of RAW dirs and LIB dirs
    //aflow --icsd_check_raw 1 > output.dat
    //use the list of RAW and LIB dirs that have been compiled in aflow
    //i.e. XHOST.vGlobal.at(X) (ICSD_LIB);  // aflow_data_calculated.cpp
    //XHOST.vGlobal.at(X) (ICSD_RAW);  // aflow_data_calculated.cpp
    //This routine checks the validity of the init::AFLOW_Projects_Directories("ICSD")/RAW/ dirs and files
    //with respect to init::AFLOW_Projects_Directories("ICSD")/LIB/

    //output.dat will contain:
    //init::AFLOW_Projects_Directories("ICSD")/RAW/../../  OK
    //init::AFLOW_Projects_Directories("ICSD")/RAW/../../  OK
    //init::AFLOW_Projects_Directories("ICSD")/RAW/../../  ERROR file1 file2 ...
    //init::AFLOW_Projects_Directories("ICSD")/RAW/../../  NotInLIB
    //init::AFLOW_Projects_Directories("ICSD")/RAW/../../  NotInLIB
    //init::AFLOW_Projects_Directories("ICSD")/LIB/../../  NotInRAW
    //init::AFLOW_Projects_Directories("ICSD")/LIB/../../  NotInRAW
    //and so on
    //the MISSING lines mean that the structures are in LIB but not in RAW	

    if(XHOST.vext.size()!=XHOST.vcat.size()) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"XHOST.vext.size()!=XHOST.vcat.size(), aborting",_RUNTIME_ERROR_); //CO20200624
    }

    vector<string> LIBlist,RAWlist,LIBnotRAW,RAWnotLIB,LIBRAW;
    string stmp;
    bool inew=true; //whether I generate an update list of LIB and RAW
    ofstream oftmp;
    ifstream iftmp;
    if(argv.size()>2) if(argv.at(2)=="1") inew=false;
    if(inew) {
      //find LIB
      aurostd::execute("find "+init::AFLOW_Projects_Directories("ICSD")+"/LIB/ -type d | grep ICSD | sort > wLIBlist.tmp");
      //find RAW
      aurostd::execute("find "+init::AFLOW_Projects_Directories("ICSD")+"/RAW/ -type d | grep ICSD | sort > wRAWlist.tmp");
    } else {
      init::InitGlobalObject("vLIBS");
      //write LIB
      oftmp.open("wLIBlist.tmp");
      for(uint i=0;i<XHOST_Library_CALCULATED_ICSD_LIB.size();i++)
        oftmp << XHOST_Library_CALCULATED_ICSD_LIB.at(i) << endl;
      oftmp.close();
      //write RAW
      oftmp.open("wRAWlist.tmp");
      for(uint i=0;i<XHOST_Library_CALCULATED_ICSD_RAW.size();i++)
        oftmp << XHOST_Library_CALCULATED_ICSD_RAW.at(i) << endl;
      oftmp.close();		
    }
    aurostd::execute("subst \""+init::AFLOW_Projects_Directories("ICSD")+"/LIB/\" \"\" wLIBlist.tmp");
    aurostd::execute("subst \""+init::AFLOW_Projects_Directories("ICSD")+"/RAW/\" \"\" wRAWlist.tmp");

    //loading lib to a vector string
    LIBlist.clear();
    iftmp.open("wLIBlist.tmp");
    if(!iftmp) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"can not open wLIBlist.tmp, aborted",_FILE_CORRUPT_); //CO20200624
    }
    while(!iftmp.eof()) {
      iftmp >> stmp;
      LIBlist.push_back(stmp);
    }
    iftmp.close();
    //loading raw to a vector string
    RAWlist.clear();
    iftmp.open("wRAWlist.tmp");
    if(!iftmp) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"can not open wRAWlist.tmp, aborted",_FILE_CORRUPT_); //CO20200624
    }
    while(!iftmp.eof()) {
      iftmp >> stmp;
      RAWlist.push_back(stmp);
    }
    iftmp.close();
    bool ifound=false;
    int i=0,j=0;
    int Nraw=0,Nlib=0,Nrawnotlib=0,Nlibnotraw=0,Nlibraw=0;
    Nlib = LIBlist.size();
    Nraw = RAWlist.size();
    //finding LIBnotRAW (dirs in LIB but not in RAW)
    // and  LIBRAW (dirs in LIB && RAW)
    LIBnotRAW.clear();
    LIBRAW.clear();
    for(i=0;i<Nlib;i++) {
      stmp = LIBlist.at(i);
      ifound = false;
      for(j=0;j<Nraw;j++) {
        if(RAWlist[j]==stmp) {ifound=true; break;}
      }
      if(ifound) LIBRAW.push_back(stmp);
      else LIBnotRAW.push_back(stmp);
    }
    Nlibraw = LIBRAW.size();
    Nlibnotraw = LIBnotRAW.size();
    //finding RAWnotLIB
    RAWnotLIB.clear();
    for(i=0;i<Nraw;i++) {
      stmp = RAWlist.at(i);
      ifound = false;
      for(j=0;j<Nlib;j++) {
        if(LIBlist[j]==stmp) {ifound=true; break;}
      }
      if(!ifound) RAWnotLIB.push_back(stmp);
    }
    Nrawnotlib = RAWnotLIB.size();
    //output LIBnotRAW
    for(i=0;i<Nlibnotraw;i++) cout << init::AFLOW_Projects_Directories("ICSD") << "LIB/" << LIBnotRAW.at(i) << " NotInRaw" << endl;
    //output RAWnotLIB
    for(i=0;i<Nrawnotlib;i++) cout << init::AFLOW_Projects_Directories("ICSD") << "RAW/" << RAWnotLIB.at(i) << " NotInLib" << endl;	

    // processing LIBRAW
    int itmp=0;
    string directory_RAW,directory_LIB;
    for(i=0;i<Nlibraw;i++) {
      directory_LIB = init::AFLOW_Projects_Directories("ICSD")+"LIB/"+LIBRAW.at(i)+"/";
      directory_RAW = init::AFLOW_Projects_Directories("ICSD")+"RAW/"+LIBRAW.at(i)+"/";
      cerr << "Processing " << directory_RAW << endl;
      // 1. check if exist
      ifound=TRUE;
      if(ifound) ifound = (aurostd::FileExist(directory_RAW+"/POSCAR.bands") || aurostd::EFileExist(directory_RAW+"/POSCAR.bands"));
      if(ifound) ifound = (aurostd::FileExist(directory_RAW+"/KPOINTS.bands") || aurostd::EFileExist(directory_RAW+"/KPOINTS.bands"));
      if(ifound) ifound = (aurostd::FileExist(directory_RAW+"/EIGENVAL.bands") || aurostd::EFileExist(directory_RAW+"/EIGENVAL.bands"));
      if(ifound) ifound = (aurostd::FileExist(directory_RAW+"/DOSCAR.bands") || aurostd::EFileExist(directory_RAW+"/DOSCAR.bands"));
      if(!ifound) { cout << directory_RAW << "ERROR - pflow::ICSD_CheckRaw: " << "POSCAR.bands[.EXT], KPOINTS.bands[.EXT], EIGENVAL.bands[.EXT], DOSCAR.static[.EXT] in directory_RAW=" << directory_RAW << " do not exist " << endl; continue; }
      if(ifound) ifound = (aurostd::FileExist(directory_LIB+"/POSCAR.bands") || aurostd::EFileExist(directory_LIB+"/POSCAR.bands"));
      if(ifound) ifound = (aurostd::FileExist(directory_LIB+"/KPOINTS.bands") || aurostd::EFileExist(directory_LIB+"/KPOINTS.bands"));
      if(ifound) ifound = (aurostd::FileExist(directory_LIB+"/DOSCAR.bands") || aurostd::EFileExist(directory_LIB+"/DOSCAR.bands"));
      if(!ifound) { cout << directory_RAW << "ERROR - pflow::ICSD_CheckRaw: " << "POSCAR.bands[.EXT], KPOINTS.bands[.EXT], DOSCAR.static[.EXT] in directory_LIB=" << directory_LIB << " do not exist " << endl; continue; }

      // 2. diff POSCAR.bands
      for(uint iext=0;iext<XHOST.vext.size();iext++) {
        if(aurostd::FileExist(directory_LIB+"/POSCAR.bands"+XHOST.vext.at(iext))) aurostd::execute(XHOST.vcat.at(iext)+" "+directory_LIB+"/POSCAR.bands"+XHOST.vext.at(iext)+" > wLIB.tmp");
        if(aurostd::FileExist(directory_RAW+"/POSCAR.bands"+XHOST.vext.at(iext))) aurostd::execute(XHOST.vcat.at(iext)+" "+directory_RAW+"/POSCAR.bands"+XHOST.vext.at(iext)+" > wRAW.tmp");
      }
      aurostd::execute("diff -w -B wLIB.tmp wRAW.tmp | wc -l > wDIFF.tmp");      
      iftmp.open("wDIFF.tmp");
      iftmp >> itmp;
      if(itmp>0) ifound=false;
      iftmp.close();
      if(!ifound) { cout << directory_RAW << "ERROR - pflow::ICSD_CheckRaw: " << "POSCAR.bands" << endl;continue; }

      // 3. diff KPOINTS.bands
      for(uint iext=0;iext<XHOST.vext.size();iext++) {
        if(aurostd::FileExist(directory_LIB+"/KPOINTS.bands"+XHOST.vext.at(iext))) aurostd::execute(XHOST.vcat.at(iext)+" "+directory_LIB+"/KPOINTS.bands"+XHOST.vext.at(iext)+" > wLIB.tmp");
        if(aurostd::FileExist(directory_RAW+"/KPOINTS.bands"+XHOST.vext.at(iext))) aurostd::execute(XHOST.vcat.at(iext)+" "+directory_RAW+"/KPOINTS.bands"+XHOST.vext.at(iext)+" > wRAW.tmp");
      }
      aurostd::execute("diff -w -B wLIB.tmp wRAW.tmp | wc -l > wDIFF.tmp");      
      iftmp.open("wDIFF.tmp");
      iftmp >> itmp;
      if(itmp>0) ifound=false;
      iftmp.close();
      if(!ifound) {
        cout << directory_RAW << "ERROR - pflow::ICSD_CheckRaw: " << "KPOINTS.bands" << endl;
        continue;
      }
      //// 4. diff EIGENVAL.bands
      //for(uint iext=0;iext<XHOST.vext.size();iext++) {
      //if(aurostd::FileExist(directory_LIB+"/EIGENVAL.bands"+XHOST.vext.at(iext))) aurostd::execute(XHOST.vcat.at(iext)+" "+directory_LIB+"/EIGENVAL.bands"+XHOST.vext.at(iext)+" > wLIB.tmp");
      //if(aurostd::FileExist(directory_RAW+"/EIGENVAL.bands"+XHOST.vext.at(iext))) aurostd::execute(XHOST.vcat.at(iext)+" "+directory_RAW+"/EIGENVAL.bands"+XHOST.vext.at(iext)+" > wRAW.tmp");
      //}
      //aurostd::execute("diff -w -B wLIB.tmp wRAW.tmp | wc -l > wDIFF.tmp");      
      //iftmp.open("wDIFF.tmp");
      //iftmp >> itmp;
      //if(itmp>0) ifound=false;
      //iftmp.close();
      //if(!ifound) {
      //cout << directory_RAW << " ERROR " << "EIGENVAL.bands" << endl;
      //continue;
      //}

      // 5. diff DOSCAR.static
      for(uint iext=0;iext<XHOST.vext.size();iext++) {
        if(aurostd::FileExist(directory_LIB+"/DOSCAR.static"+XHOST.vext.at(iext))) aurostd::execute(XHOST.vcat.at(iext)+" "+directory_LIB+"/DOSCAR.static"+XHOST.vext.at(iext)+" > wLIB.tmp");
        if(aurostd::FileExist(directory_RAW+"/DOSCAR.static"+XHOST.vext.at(iext))) aurostd::execute(XHOST.vcat.at(iext)+" "+directory_RAW+"/DOSCAR.static"+XHOST.vext.at(iext)+" > wRAW.tmp");
      }
      aurostd::execute("diff -w -B wLIB.tmp wRAW.tmp | wc -l > wDIFF.tmp");      

      iftmp.open("wDIFF.tmp");
      iftmp >> itmp;
      if(itmp>0) ifound=false;
      iftmp.close();
      if(!ifound) {
        cout << directory_RAW << "ERROR - pflow::ICSD_CheckRaw: " << "DOSCAR.static.EXT" << endl;
        continue;
      }
      cout << directory_RAW << " OK" << endl;
    }
    //delete tmp files
    aurostd::execute("rm -f wRAWlist.tmp* wLIBlist.tmp* wLIB.tmp* wRAW.tmp* wDIFF.tmp*");
  }
}


// ***************************************************************************
// pflow::ICSD_2POSCAR
// ***************************************************************************
namespace pflow {
  void ICSD_2POSCAR(istream& input) {
    //Synopsis:	
    //aflow --icsd2poscar < input.icsd
    //Generate POSCAR from icsd library stored in input.icsd
    string sword="";
    input >> sword;
    cerr << "**** Not Implemented Yet ****" << endl;
  }
}

// ***************************************************************************
// pflow::ICSD_2PROTO
// ***************************************************************************
namespace pflow {
  void ICSD_2PROTO(istream& input) {
    //Synopsis:
    //--icsd2proto < input.icsd
    //Write output to stdout  (prototype-format) from stdin ICSD-format
    //In addition, it writes some information for each structure to the file "proto.label"

    bool ifound,iread,iwrite,AtomLineFound=false,icontinue=true;
    int j,k,
        ICSDid=0,SGnumber=0,Ncompounds=0,Nspec=0,CellChoice=1;
    float ftmp,Dcalc=0.0;
    string cstmp,stmp;
    string sline,sword,ssub,
           title,sdir,sICSDid,CompoundName,StructuredName,scommand,Pearson,SG;
    string::iterator siter;
    size_t ssize_t;
    vector<int> Zcount;
    vector<float> Zconc,unitcell(6),xpos,ypos,zpos,pocc;
    vector<string> Zname,WyckLabel;
    ofstream fout;
    //struct stat st;

    fout.open("proto.label");

    Ncompounds=0;iwrite=false;
    while(!input.eof()) {
      getline(input,sline);sline=sline+'\n';

      if(sline.length()>4) sword=sline.substr(0,5);
      if(sword=="*data") {
        Ncompounds++;title="";
      }

      if(sline.length()>8) sword=sline.substr(0,9);
      if(sword=="Coll Code") {
        ICSDid=atoi((sline.substr(9,sline.length())).data());
        cstmp=aurostd::utype2string(ICSDid);
        sICSDid=cstmp;
        title="ICSD_"+sICSDid;
        cerr << "processing " << title << endl;
      }

      if(sline.length()>9) sword=sline.substr(0,10);
      if(sword=="Structured") {
        StructuredName=sline.substr(11,sline.length());
        StructuredName=RemoveCharFromString(StructuredName,' ');
        StructuredName=RemoveCharFromString(StructuredName,'\t');
        StructuredName=RemoveCharFromString(StructuredName,'\n');
      }

      if(sline.length()>3) sword=sline.substr(0,4);
      if(sword=="Sum ") {
        ssub=sline.substr(5,sline.length());
        for(j=0;j<(int) ssub.length();j++) {
          if((ssub[j]==' ' || ssub[j]=='\t') && ssub[j+1]!=' ' && ssub[j+1]!='\t') {
            ssub=ssub.substr(j+1,ssub.length());
            break;
          }
        }
        CompoundName=RemoveCharFromString(ssub,' ');
        CompoundName=RemoveCharFromString(CompoundName,'\t');
        CompoundName=RemoveCharFromString(CompoundName,'\n');
        title=title+", "+CompoundName;

        //getting Nspec && Zname[]
        Nspec=0; ifound=false;
        for(j=0;j<(int) ssub.length();j++) {
          if(ssub[j]==' ' && j!=0) {
            ifound=TRUE;
            sword=ssub.substr(0,j);
            ssub=ssub.substr(j+1,ssub.length());
            j=0;
          }
          else if(j==(int) ssub.length()-1 && ssub!="") {
            ifound=TRUE;
            sword=ssub;
          }
          if(ifound) {
            Nspec++;
            for(k=0;k<(int) sword.length();k++) {
              if(sword[k]<'A' || sword[k]>'z') {
                Zname.push_back(sword.substr(0,k));
                Zconc.push_back(atof((sword.substr(k,sword.length())).data()));
                break;
              }
            }
            ifound=false;
          }//if ifound
        }//for j=0:ssub.length()

        //ascending sorting of elements in the compound
        for(j=0;j<(int)Zname.size()-1;j++) {
          for(k=j+1;k<(int)Zname.size();k++) {
            if(Zname[k]<Zname[j]) {//swap
              stmp=Zname[j];
              Zname[j]=Zname[k];
              Zname[k]=stmp;
              ftmp=Zconc[j];
              Zconc[j]=Zconc[k];
              Zconc[k]=ftmp;
            }
          }
        }
      }//"Sum "
      Zcount.resize(Nspec);

      //----------D(calc)-----------                                                                                                        
      if(sline.length()>6) sword=sline.substr(0,7);
      if(sword=="D(calc)") {
        ssub=sline.substr(7,sline.length());
        Dcalc=aurostd::string2utype<float>(ssub);
      }

      //Unit Cell
      if(sline.length()>8) sword=sline.substr(0,9);
      if(sword=="Unit Cell") {
        ssub=sline.substr(10,sline.length());
        for(j=0;j<(int) ssub.length();j++) {
          if((ssub[j]==' ' || ssub[j]=='\t') && ssub[j+1]!=' ' && ssub[j+1]!='\t') {
            ssub=ssub.substr(j+1,ssub.length());
            break;
          }
        }
        ssub=RemoveStringFromTo(ssub,'(',')');
        //ssub=RemoveCharFromString(ssub,'(');
        //ssub=RemoveCharFromString(ssub,')');

        k=0;
        for(j=0;j<(int) ssub.length();j++) {
          if(ssub[j]==' ') {
            unitcell[k++]=atof((ssub.substr(0,j)).data());
            ssub=ssub.substr(j+1,ssub.length());
            j=0;
          }
          else if(j==(int) ssub.length()-1 && ssub!="")
            unitcell[k++]=atof(ssub.data());
        }
      }//Unit Cell

      if(sline.length()>8) sword=sline.substr(0,9);
      if(sword=="SG Number") {
        SGnumber=atoi((sline.substr(9,sline.length())).data());      
      }

      if(sline.length()>6) sword=sline.substr(0,7);
      if(sword=="Pearson") {
        Pearson=sline.substr(8,sline.length());
        Pearson=RemoveCharFromString(Pearson,' ');
        Pearson=RemoveCharFromString(Pearson,'\t');
        Pearson=RemoveCharFromString(Pearson,'\n');
      }

      if(sline.length()>10) sword=sline.substr(0,11);
      if(sword=="Space Group") {
        SG=sline.substr(12,sline.length());
        SG=RemoveCharFromString(SG,' ');
        SG=RemoveCharFromString(SG,'\t');
        SG=RemoveCharFromString(SG,'\n');
      }

      //cell choice
      if(sline.length()>7) sword=sline.substr(0,8);
      if(sword=="Comments") {
        ssub=sline;j=0;icontinue=true;
        while(icontinue && j<1000) {
          j++;
          getline(input,sline);sline=sline+'\n';
          if(sline.length()>4) sword=sline.substr(0,5);
          if(sword=="Atom ") icontinue=false;
          else ssub=ssub+sline;
        }
        ssub=RemoveCharFromString(ssub,' ');
        ssub=RemoveCharFromString(ssub,'\t');
        ssub=RemoveCharFromString(ssub,'\n');
        ssize_t=ssub.find("cellchoice");
        if(ssize_t!=ssub.npos) {sword=ssub.at(ssize_t+10); CellChoice=atoi(sword.data());}
        else CellChoice=100;
      }

      //position of atoms
      if(sline.length()>4) sword=sline.substr(0,5);
      if(sword=="Atom " && !AtomLineFound) {
        AtomLineFound=TRUE;
        iread=TRUE;
        for(j=0;j<Nspec;j++) Zcount[j]=0;//clearing Zcount
        while(iread) {//processing one wyckoff entry
          iread=false;
          input >> sline;//atom name
          for(j=0;j<Nspec;j++) {
            if(sline==Zname[j]) {//if it is atom name
              iread=TRUE;
              WyckLabel.push_back(sline);
              Zcount[j]=Zcount[j]+1;
              input >> sline >> sline >> sline >> sline;//discarded
              input >> sline;//x
              sline=RemoveStringFromTo(sline,'(',')');
              //sline=RemoveCharFromString(sline,'(');
              //sline=RemoveCharFromString(sline,')');
              xpos.push_back(atof(sline.data()));//save x
              input >> sline;//y
              sline=RemoveStringFromTo(sline,'(',')');
              //sline=RemoveCharFromString(sline,'(');
              //sline=RemoveCharFromString(sline,')');
              ypos.push_back(atof(sline.data()));//save y
              input >> sline;//z
              sline=RemoveStringFromTo(sline,'(',')');
              //sline=RemoveCharFromString(sline,'(');
              //sline=RemoveCharFromString(sline,')');
              zpos.push_back(atof(sline.data()));//save z
              input >> sline;//sof
              sline=RemoveStringFromTo(sline,'(',')');
              //sline=RemoveCharFromString(sline,'(');
              //sline=RemoveCharFromString(sline,')');
              pocc.push_back(atof(sline.data()));//save occupancy number
              getline(input,sline);//discard the rest of the line
              sline="";
            }//if Zname
          }//Nspec
        }//iread
        //ascending sorting according to elements name
        //note that this is nore needed because it is taken care of in the outputting block
      }//atom

      if(sline.length()>3) sword=sline.substr(0,4);
      if(sword=="*end") {
        AtomLineFound=FALSE;
        iwrite=true;
        /*
        //any partial occupation?
        POF=false;
        for(j=0;j<(int) pocc.size();j++) {
        if(pocc[j]<1.00) POF=true;
        }
        //outputting prototype library
        iwrite=false; iwarning=false;
        if(POF) {//partial occupation detected
        iwrite=false;
        cerr << "  Partial occupation deteced, skipping (not written)." << endl;
        }
        else {//no partial occupation
        iwrite=true;
        }
        */
        if(iwrite) {
          for(j=0;j<Nspec;j++) fout << Zname[j] << Zconc[j];
          fout << "_ICSD_" << ICSDid << " "
            << StructuredName
            << " Pearson " << Pearson
            << " SG " << SGnumber << " " << SG
            << " rhomass " << Dcalc << endl;

          cout << "// ***************************************************************************" << endl
            << "// ";
          for(j=0;j<Nspec;j++) cout << Zname[j] << Zconc[j];
          cout << " ICSD_" << ICSDid << " ";
          for(j=0;j<Nspec;j++) cout << Zname[j] << Zconc[j];
          cout << "_ICSD_" << ICSDid << endl
            << "// ***************************************************************************" << endl;

          cout << "// "; for(j=0;j<Nspec;j++) cout << Zname[j] << Zconc[j];
          cout << " " << Pearson << " " << SG << " " << StructuredName << " " << SGnumber << " ICSD_" << ICSDid << " ";
          for(j=0;j<Nspec;j++) cout << Zname[j] << Zconc[j];
          cout << "_ICSD_" << ICSDid << endl;

          cout << "STRUCTURE ";
          for(j=0;j<Nspec;j++) cout << Zname[j] << Zconc[j];
          cout << " ICSD_" << ICSDid << " ";
          for(j=0;j<Nspec;j++) cout << Zname[j] << Zconc[j];
          cout << "_ICSD_" << ICSDid << endl;

          cout << "PROTOTYPE "; for(j=0;j<Nspec;j++) cout << Zname[j] << Zconc[j];
          cout << " [";for(j=0;j<Nspec;j++) cout << Zname[j] << Zconc[j];cout << "]"
            << " " << Pearson << " " << SG << " " << StructuredName << " " << SGnumber << " ";
          for(j=0;j<Nspec;j++) cout << Zname[j] << Zconc[j];
          cout << "_ICSD_" << ICSDid << " ICSD_" << ICSDid << " " << endl;

          cout << "INFO ICSD_" << ICSDid << " ";
          for(j=0;j<Nspec;j++) cout << Zname[j] << Zconc[j];
          cout << "_ICSD_" << ICSDid << endl;

          cout << "MODE WYC_ICSD PRIM CONVENTIONAL POCC" << endl;
          for(j=0;j<6;j++) cout << " " << unitcell[j];
          cout << " " << SGnumber << " ";
          if(CellChoice!=100) cout << CellChoice;
          cout << endl;

          for(j=0;j<Nspec;j++) {
            for(k=0;k<(int) xpos.size();k++) {
              if(WyckLabel[k]==Zname[j]) {
                cout << " " << setw(10) << xpos[k]
                  << " " << setw(10) << ypos[k]
                  << " " << setw(10) << zpos[k]
                  << " " << setw(4) << WyckLabel[k];
                if(pocc[k]<1.0)
                  cout << "  " << pocc[k];
                else
                  cout << "  -";
                cout << endl;
              }
            }
          }
        }//iwrite

        //clearing vectors for the next compounds
        xpos.clear(); ypos.clear(); zpos.clear(); pocc.clear();
        Zcount.clear(); Zname.clear(); Zconc.clear(); WyckLabel.clear();
        sword="";
      }//"*end"
    }//input.eof
    fout.close();
  }
}

// ***************************************************************************
// pflow::ICSD_2WYCK
// ***************************************************************************
namespace pflow {
  void ICSD_2WYCK(istream& input,bool SOF) {
    //   Write output to file WYCKCAR  (WYCKCAR-format) from stdin ICSD-format
    //   aflow --icsd2wyck < input.icsd
    //   If the input contains multiple compounds, a subfolder will
    //   be created for each compound using the following folder name template:
    //   CompoundFormula_ICSD_ICSDid/
    //   The ICSDid is the compound id in the original ICSD database. Only compounds
    //   with elements having full occupation (sof=1) will be processed.

    //   aflow --icsd2wyck --sof < input.icsd
    //   All compounds will be processed. If partial occupation is detected, the sof
    //   will be written as part of the label of the element and WARNING will be
    //   written in the title and cerr.

    //   2009, wahyu@alumni.duke.edu

    bool ifound,iread,iwrite,iwarning,POF,AtomLineFound=false;
    int j,k,ICSDid=0,SGnumber=0,Ncompounds=0,Nspec=0;
    int itmp=0; if(itmp) {;} // dummy load
    string cstmp;
    string sline,sword,ssub,
           title,sdir,sICSDid,CompoundName,scommand;
    string::iterator siter;
    vector<int> Zcount;
    vector<float> Zconc,unitcell(6),xpos,ypos,zpos,pocc;
    vector<string> Zname,WyckLabel;
    ofstream fout;
    struct stat st;

    Ncompounds=0;iwrite=false;
    while(!input.eof()) {
      getline(input,sline);sline=sline+'\n';

      if(sline.length()>4) sword=sline.substr(0,5);
      if(sword=="*data") {
        Ncompounds++;title="";
      }

      if(sline.length()>8) sword=sline.substr(0,9);
      if(sword=="Coll Code") {
        ICSDid=atoi((sline.substr(9,sline.length())).data());
        cstmp=aurostd::utype2string(ICSDid);
        sICSDid=cstmp;
        title="ICSD_"+sICSDid;
        cerr << "processing " << title << endl;
      }

      if(sline.length()>3) sword=sline.substr(0,4);
      if(sword=="Sum ") {
        ssub=sline.substr(5,sline.length());
        for(j=0;j<(int) ssub.length();j++) {
          if((ssub[j]==' ' || ssub[j]=='\t') && ssub[j+1]!=' ' && ssub[j+1]!='\t') {
            ssub=ssub.substr(j+1,ssub.length());
            break;
          }
        }
        CompoundName=ssub;
        CompoundName=RemoveCharFromString(CompoundName,' ');
        CompoundName=RemoveCharFromString(CompoundName,'\n');
        title=title+", "+CompoundName;

        //getting Nspec && Zname[]
        Nspec=0; ifound=false;
        for(j=0;j<(int) ssub.length();j++) {
          if(ssub[j]==' ' && j!=0) {
            ifound=TRUE;
            sword=ssub.substr(0,j);
            ssub=ssub.substr(j+1,ssub.length());
            j=0;
          }
          else if(j==(int) ssub.length()-1 && ssub!="") {
            ifound=TRUE;
            sword=ssub;
          }
          if(ifound) {
            Nspec++;
            for(k=0;k<(int) sword.length();k++) {
              if(sword[k]<'A' || sword[k]>'z') {
                Zname.push_back(sword.substr(0,k));
                Zconc.push_back(atof((sword.substr(k,sword.length())).data()));
                break;
              }
            }
            ifound=false;
          }//if ifound
        }//for j=0:ssub.length()
      }//"Sum "
      Zcount.resize(Nspec);

      //Unit Cell
      if(sline.length()>8) sword=sline.substr(0,9);
      if(sword=="Unit Cell") {
        ssub=sline.substr(10,sline.length());
        for(j=0;j<(int) ssub.length();j++) {
          if((ssub[j]==' ' || ssub[j]=='\t') && ssub[j+1]!=' ' && ssub[j+1]!='\t') {
            ssub=ssub.substr(j+1,ssub.length());
            break;
          }
        }
        ssub=RemoveStringFromTo(ssub,'(',')');
        //ssub=RemoveCharFromString(ssub,'(');
        //ssub=RemoveCharFromString(ssub,')');

        k=0;
        for(j=0;j<(int) ssub.length();j++) {
          if(ssub[j]==' ') {
            unitcell[k++]=atof((ssub.substr(0,j)).data());
            ssub=ssub.substr(j+1,ssub.length());
            j=0;
          }
          else if(j==(int) ssub.length()-1 && ssub!="")
            unitcell[k++]=atof(ssub.data());
        }
      }//Unit Cell

      if(sline.length()>8) sword=sline.substr(0,9);
      if(sword=="SG Number") {
        SGnumber=atoi((sline.substr(9,sline.length())).data());      
      }

      //position of atoms
      if(sline.length()>4) sword=sline.substr(0,5);
      if(sword=="Atom " && !AtomLineFound) {
        AtomLineFound=TRUE;
        iread=TRUE;
        for(j=0;j<Nspec;j++) Zcount[j]=0;//clearing Zcount
        while(iread) {//processing one wyckoff entry
          iread=false;
          input >> sline;//atom name
          for(j=0;j<Nspec;j++) {
            if(sline==Zname[j]) {//if it is atom name
              iread=TRUE;
              WyckLabel.push_back(sline);
              Zcount[j]=Zcount[j]+1;
              input >> sline >> sline >> sline >> sline;//discarded
              input >> sline;//x
              sline=RemoveStringFromTo(sline,'(',')');
              //sline=RemoveCharFromString(sline,'(');
              //sline=RemoveCharFromString(sline,')');
              xpos.push_back(atof(sline.data()));//save x
              input >> sline;//y
              sline=RemoveStringFromTo(sline,'(',')');
              //sline=RemoveCharFromString(sline,'(');
              //sline=RemoveCharFromString(sline,')');
              ypos.push_back(atof(sline.data()));//save y
              input >> sline;//z
              sline=RemoveStringFromTo(sline,'(',')');
              //sline=RemoveCharFromString(sline,'(');
              //sline=RemoveCharFromString(sline,')');
              zpos.push_back(atof(sline.data()));//save z
              input >> sline;//sof
              sline=RemoveStringFromTo(sline,'(',')');
              //sline=RemoveCharFromString(sline,'(');
              //sline=RemoveCharFromString(sline,')');
              pocc.push_back(atof(sline.data()));//save sof
              getline(input,sline);//discard the rest of the line
              sline="";
            }//if Zname
          }//Nspec
        }//iread
      }//atom

      if(sline.length()>3) sword=sline.substr(0,4);
      if(sword=="*end") {
        AtomLineFound=false;
        //any partial occupation?
        POF=false;
        for(j=0;j<(int) pocc.size();j++) {
          if(pocc[j]<1.00) POF=true;
        }
        //writing WYCKCAR
        iwrite=false; iwarning=false;
        if(POF) {//partial occupation detected
          iwarning=true;
          if(SOF) iwrite=true;
          else {iwrite=false; cerr << "  Partial occupation deteced, skipping (not written)." << endl;}
        }
        else {//no partial occupation
          iwarning=false;iwrite=true;
        }
        if(iwrite) {
          fout.open("WYCKCAR");
          if(iwarning) {
            cerr << "  WARNING: partial occupation detected but written." << endl;
            fout << "WARNING(partial occup.) ";
          }
          fout << title << endl;
          fout << " 1.0000" << endl;
          for(j=0;j<6;j++) fout << " " << unitcell[j];
          fout << " " << SGnumber << endl;
          for(j=0;j<Nspec;j++) fout << " " << Zcount[j];
          fout << "\nDirect" << endl;
          itmp=0;
          for(j=0;j<Nspec;j++) {
            for(k=0;k<(int) xpos.size();k++) {
              if(WyckLabel[k]==Zname[j]) {
                fout << " " << setw(10) << xpos[k]
                  << " " << setw(10) << ypos[k]
                  << " " << setw(10) << zpos[k]
                  << " " << WyckLabel[k];
                if(iwarning) fout << "_sof_" << pocc[k];
                fout << endl;
              }
            }
          }
          fout.close();

          //saving in subdirectory
          sdir=CompoundName+"_ICSD_"+sICSDid+"/";
          scommand="mkdir "+sdir;
          if(stat(sdir.data(),&st)!=0) system(scommand.data()); //create subdirectory if it does not exist
          scommand="cp WYCKCAR "+sdir;
          system(scommand.data());
        }//iwrite

        //clearing vectors for the next compounds
        xpos.clear(); ypos.clear(); zpos.clear(); pocc.clear();
        Zcount.clear(); Zname.clear(); Zconc.clear(); WyckLabel.clear();
        sword="";
      }//"*end"
    }//input.eof
    if(iwrite) {
      if(Ncompounds==1) {
        scommand="cp "+sdir+"WYCKCAR . && rm -rf "+sdir;
        system(scommand.data());
      }
      if(Ncompounds>1) system("rm -f WYCKCAR");//clean up
    }
  }
}

// ***************************************************************************
// pflow::ICSD_ListMetals
// ***************************************************************************
namespace pflow {
  void ICSD_ListMetals() {
    //   Print out to stdout a list of metalic elements.
    vector<string> M;
    M=GetMetalsList(true);
  }
}

// ***************************************************************************
// GetMetalsList
// ***************************************************************************
vector<string> GetMetalsList(bool v) {
  //   Return a vector<string> of a list of metalic elements. If v=true, it will
  //   also print it out to stdout
  int i,j;
  vector<string> vmetal;
  vmetal.clear();
  //Alkali metals
  vmetal.push_back("Li");vmetal.push_back("Na");vmetal.push_back("K");vmetal.push_back("Rb");
  vmetal.push_back("Cs");vmetal.push_back("Fr");
  j=0;
  if(v) {
    cout << "Alkali Metals:";
    for(i=j;i<(int)vmetal.size();i++) cout << " " << vmetal.at(i);
    cout << '\n';
  }
  j=(int)vmetal.size();
  //Alkaline earth metals
  vmetal.push_back("Be");vmetal.push_back("Mg");vmetal.push_back("Ca");vmetal.push_back("Sr");
  vmetal.push_back("Ba");vmetal.push_back("Ra");
  if(v) {
    cout << "Alkaline Earth Metals:";
    for(i=j;i<(int)vmetal.size();i++) cout << " " << vmetal.at(i);
    cout << '\n';
  }
  j=(int)vmetal.size();
  //Transition metals
  vmetal.push_back("Sc");vmetal.push_back("Ti");vmetal.push_back("V");vmetal.push_back("Cr");
  vmetal.push_back("Mn");vmetal.push_back("Fe");vmetal.push_back("Co");vmetal.push_back("Ni");
  vmetal.push_back("Cu");vmetal.push_back("Zn");vmetal.push_back("Y");vmetal.push_back("Zr");
  vmetal.push_back("Nb");vmetal.push_back("Mo");vmetal.push_back("Tc");vmetal.push_back("Ru");
  vmetal.push_back("Rh");vmetal.push_back("Pd");vmetal.push_back("Ag");vmetal.push_back("Cd");
  vmetal.push_back("Hf");vmetal.push_back("Ta");vmetal.push_back("W");vmetal.push_back("Re");
  vmetal.push_back("Os");vmetal.push_back("Ir");vmetal.push_back("Pt");vmetal.push_back("Au");
  vmetal.push_back("Hg");
  vmetal.push_back("Rf");vmetal.push_back("Db");vmetal.push_back("Sg");vmetal.push_back("Bh");
  vmetal.push_back("Hs");vmetal.push_back("Mt");vmetal.push_back("Ds");vmetal.push_back("Rg");vmetal.push_back("Uub");
  if(v) {
    cout << "Transition Metals:";
    for(i=j;i<(int)vmetal.size();i++) cout << " " << vmetal.at(i);
    cout << '\n';
  }
  j=(int)vmetal.size();
  //Lanthanoids
  vmetal.push_back("La");vmetal.push_back("Ce");vmetal.push_back("Pr");vmetal.push_back("Nd");
  vmetal.push_back("Pm");vmetal.push_back("Sm");vmetal.push_back("Eu");vmetal.push_back("Gd");
  vmetal.push_back("Tb");vmetal.push_back("Dy");vmetal.push_back("Ho");vmetal.push_back("Er");
  vmetal.push_back("Tm");vmetal.push_back("Yb");vmetal.push_back("Lu");
  if(v) {
    cout << "Lanthanoids:";
    for(i=j;i<(int)vmetal.size();i++) cout << " " << vmetal.at(i);
    cout << '\n';
  }
  j=(int)vmetal.size();
  //Actinoids
  vmetal.push_back("Ac");vmetal.push_back("Th");vmetal.push_back("Pa");vmetal.push_back("U");
  vmetal.push_back("Np");vmetal.push_back("Pu");vmetal.push_back("Am");vmetal.push_back("Cm");
  vmetal.push_back("Bk");vmetal.push_back("Cf");vmetal.push_back("Es");vmetal.push_back("Fm");
  vmetal.push_back("Md");vmetal.push_back("No");vmetal.push_back("Lr");
  if(v) {
    cout << "Actinoids:";
    for(i=j;i<(int)vmetal.size();i++) cout << " " << vmetal.at(i);
    cout << '\n';
  }
  j=(int)vmetal.size();
  //Other metals
  vmetal.push_back("Al");vmetal.push_back("Ga");vmetal.push_back("In");vmetal.push_back("Sn");
  vmetal.push_back("Tl");vmetal.push_back("Pb");vmetal.push_back("Bi");
  if(v) {
    cout << "Other Metals:";
    for(i=j;i<(int)vmetal.size();i++) cout << " " << vmetal.at(i);
    cout << '\n';
  }
  return vmetal;
}

// ***************************************************************************
// GetBandgap
// ***************************************************************************
float GetBandGap_WAHYU(stringstream& ein,float Efermi,char& gaptype) {
  //Calculate bandgap Egap from stringstream EIGENVAL.bands (ein)
  //using Efermi to determine the occupied/unocc bands.
  //Egap = CBM-VBM
  //Set Egap = -1.0 if not found.
  //gaptype is either 'D' or 'I'

  float metal_gap_tol = DEFAULT_METAL_GAP_TOLERANCE;
  int Nk,Nbands,i,ik,ib,ispin=1,count,itmp;
  float ftmp;
  string stmp;

  //cerr << "Efermi: " << Efermi << endl;

  for(i=0;i<5;i++) getline(ein,stmp);
  ein >> itmp >> Nk >> Nbands;
  getline(ein,stmp);

  //first kpoint
  getline(ein,stmp);//empty line
  getline(ein,stmp);//kpoint coordinate
  getline(ein,stmp);//first entry line
  vector<string> sword;
  count=StringCrop(stmp,sword);
  if(count==2) ispin=1;
  if(count==3) ispin=2;
  if(ispin==0) {
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ispin = 0, aborted",_INPUT_ILLEGAL_); //CO20200624
  }
  vector<vector<float> > data(Nk);
  for(i=0;i<Nk;i++) {
    data[i].resize(Nbands*ispin);
  }
  ik=0;
  if(ispin==1) {
    ib=0;
    data[ik][ib]=aurostd::string2utype<float>(sword[1]);
    for(i=1;i<Nbands;i++) {
      ein >> ftmp >> ftmp;
      ib=ib+1;
      data[ik][ib]=ftmp;
    }
    getline(ein,stmp);
  }
  ik=0;
  if(ispin==2) {
    ib=0;
    data[ik][ib]=aurostd::string2utype<float>(sword[1]);
    data[ik][ib+Nbands]=aurostd::string2utype<float>(sword[2]);
    for(i=1;i<Nbands;i++) {
      ein >> ftmp;
      ein >> ftmp;//up
      ib=ib+1;
      data[ik][ib]=ftmp;
      ein >> ftmp;//down
      data[ik][ib+Nbands]=ftmp;
    }
    getline(ein,stmp);
  }
  //--now read all other kpoints
  if(ispin==1) {
    for(ik=1;ik<Nk;ik++) {
      getline(ein,stmp);//empty line
      getline(ein,stmp);//kpoint coordinate
      ib=-1;
      for(i=0;i<Nbands;i++) {
        ein >> ftmp >> ftmp;
        ib=ib+1;
        data[ik][ib]=ftmp;
      }
      getline(ein,stmp);
    }
  }
  if(ispin==2) {
    for(ik=1;ik<Nk;ik++) {
      getline(ein,stmp);//empty line
      getline(ein,stmp);//kpoint coordinate
      ib=-1;
      for(i=0;i<Nbands;i++) {
        ein >> ftmp;
        ein >> ftmp;//up
        ib=ib+1;
        data[ik][ib]=ftmp;
        ein >> ftmp;//down
        data[ik][ib+Nbands]=ftmp;
      }
      getline(ein,stmp);
    }
  }
  //cerr << "reading done" << endl;

  //finding ivb and icb for each up and down spin separately
  //correspoding to index of vb and cb, we need to do separately
  //because up and down vb and cb don't usually occur at the same band index.
  int ivbup,icbup,ivbdw,icbdw;
  float Egap,vbm,vbmup,vbmdw=-1e6,cbm,cbmup,cbmdw=-1e6;
  float maxtest,mintest,vbmtestup,vbmtestdw,cbmtestup,cbmtestdw;

  //KY FIXES THIS BUG
  bool FLAG_HALF_BANDS_UP = false;
  bool FLAG_HALF_BANDS_DN = false;

  ivbup=0;icbup=0;
  ivbdw=0;icbdw=0;
  vbmup=Efermi; cbmup=Efermi; //initialization
  vector<float> vband,cband;
  vband.clear(); vband.resize(Nk);
  cband.clear(); cband.resize(Nk);
  for(ib=0;ib<Nbands-1;ib++) {
    for(i=0;i<Nk;i++) vband[i]=data[i][ib];
    vbmtestup=max(vband);
    for(i=0;i<Nk;i++) vband[i]=data[i][ib+1];
    cbmtestup=min(vband);
    if(vbmtestup<Efermi && cbmtestup>Efermi) {
      vbmup=vbmtestup;
      cbmup=cbmtestup;
      ivbup=ib;icbup=ib+1;
      break;
    }
  }
  //if vbmup and cmbup are both above or below Efermi
  //KY FIXES THIS BUG
  //This includes the case:  vbmdw == cbmdwa && vbmdw == Efermi
  if(((vbmup-Efermi)<1E-8 && (cbmup-Efermi)<1E-8)||
      ((vbmup-Efermi)>1E-8 && (cbmup-Efermi)>1E-8)) {
    FLAG_HALF_BANDS_UP = true;
  }

  if(ispin==2) {
    vbmdw=Efermi; cbmdw=Efermi; //initialization
    ivbdw=0;icbdw=0;
    for(ib=0;ib<Nbands-1;ib++) {
      for(i=0;i<Nk;i++) vband[i]=data[i][ib+Nbands];
      vbmtestdw=max(vband);
      for(i=0;i<Nk;i++) vband[i]=data[i][ib+1+Nbands];
      cbmtestdw=min(vband);
      if(vbmtestdw<Efermi && cbmtestdw>Efermi) {
        vbmdw=vbmtestdw;
        cbmdw=cbmtestdw;
        ivbdw=ib;icbdw=ib+1;
        break;
      }
    }
    //if vbmdw and cmbdw are both above or below Efermi
    //KY FIXES THIS BUG
    //This includes the case:  vbmdw == cbmdwa && vbmdw == Efermi
    if(((vbmdw-Efermi)<1E-8 && (cbmdw-Efermi)<1E-8)||
        ((vbmdw-Efermi)>1E-8 && (cbmdw-Efermi)>1E-8)) {
      FLAG_HALF_BANDS_DN = true;
    }
  }


  //cerr << "ivb icb: " << ivb << " " << icb << endl;
  //now determining the overal vbm and cbm
  bool vbmfromdw,cbmfromdw;
  vbm=vbmup;  cbm=cbmup;//for spin off
  vbmfromdw=false; cbmfromdw=false;
  if(ispin==2) {//for spin on, we need to check further
    //KY FIXES THIS BUG
    if(!FLAG_HALF_BANDS_UP && FLAG_HALF_BANDS_DN) {
      vbm=vbmup;
      cbm=cbmup;
    }
    else if(FLAG_HALF_BANDS_UP && !FLAG_HALF_BANDS_DN) {
      vbm=vbmdw;
      cbm=cbmdw;
    }
    else {
      if(vbmdw>vbm) {vbmfromdw=true; vbm=vbmdw;}
      if(cbmdw<cbm) {cbmfromdw=true; cbm=cbmdw;}
    }
  }
  int kvbm=0,kcbm=0;
  if(aurostd::abs(cbm-vbm) < metal_gap_tol) {
    //   cout << "ERROR [GetBandGap] VBM and CBM can not be found. Setting VBM and CBM to Efermi." << endl;
    cout << "WARNING [GetBandGap] VBM and CBM can not be found. Setting VBM and CBM to Efermi." << endl;
    gaptype='I';
    vbm=Efermi; cbm=Efermi; Egap=0;
    return(Egap);
  }
  else {
    Egap=cbm-vbm;
    //checking I/D gaptype
    kvbm=0;kcbm=0;
    maxtest=vbm-10;
    mintest=cbm+10;
    if(vbmfromdw==false) {for(i=0;i<Nk;i++) vband[i]=data[i][ivbup];}
    else {for(i=0;i<Nk;i++) vband[i]=data[i][ivbdw+Nbands];}
    if(cbmfromdw==false) {for(i=0;i<Nk;i++) cband[i]=data[i][icbup];}
    else {for(i=0;i<Nk;i++) cband[i]=data[i][icbdw+Nbands];}
    for(ik=0;ik<Nk;ik++) {
      if(vband[ik]>maxtest) {maxtest=vband[ik];kvbm=ik;}
      if(cband[ik]<mintest) {mintest=cband[ik];kcbm=ik;}
    }
    //cerr << "kvbm kcbm: " << kvbm << " " << kcbm << endl;
    gaptype='I';
    if(kvbm==kcbm) gaptype='D';
    return(Egap);
  }
}
// ***************************************************************************
// GetBandgapFromDOS
// ***************************************************************************
float GetBandgapFromDOS(ifstream& dos) {
  //   Return band gap calculated from DOSCAR file.
  bool isMetal=false;
  int i,iclosest,Ngrids;
  float Emin,Emax,EF,bandgap,Ec,Ev,tol,delE,Etmp,Dtmp,delEE,Egrid;
  string sbuf;

  // aurostd::execute(command);
  for(i=0;i<5;i++)  getline(dos,sbuf);
  dos >> Emax >> Emin >> Ngrids >> EF; getline(dos,sbuf);
  Egrid=(Emax-Emin)/Ngrids;
  //note that the EF here is not the actual fermi energy
  vector<float> E(Ngrids),D(Ngrids);
  tol=5e-3;
  delE=1e6;delEE=delE;
  iclosest=0;
  for(i=0;i<Ngrids;i++) {
    dos >> E[i] >> D[i]; getline(dos,sbuf);
  }
  for(i=0;i<Ngrids;i++) {
    if(aurostd::abs(EF-E[i])<Egrid) {
      if(D[i]>tol && D[i+1]>tol) {isMetal=true;break;}//metal
    }
  }
  if(isMetal)  bandgap=-1;
  else {
    for(i=0;i<Ngrids;i++) {
      Dtmp=D[i];Etmp=E[i];
      if(Dtmp<tol) delEE=aurostd::abs(EF-Etmp);
      if(delEE<delE) {
        iclosest=i; delE=delEE;
      }
    }
    //    cerr << "iclosest E: " << iclosest << " " << E[iclosest] << endl;
    for(i=iclosest-1;iclosest>-1;i--) {
      if(D[i]>tol) break;
    }
    Ev=E[i];
    for(i=iclosest+1;i<Ngrids;i++) {
      if(D[i]>tol) break;
    }
    Ec=E[i];
    EF=(Ev+Ec)/2;//fix the EF
    bandgap=Ec-Ev;
    if(aurostd::abs(Ev-Emin)<tol || aurostd::abs(Emax-Ec)<tol) bandgap=-1.0;
  }
  E.clear(); D.clear();
  return bandgap;
}
// ***************************************************************************
// GetBandgapFromDOS
// ***************************************************************************
float GetBandgapFromDOS(istream& dos) {
  //   Return band gap calculated from DOSCAR as input stream.
  bool isMetal=false;
  int i,iclosest,Ngrids;
  float Emin,Emax,EF,bandgap,Ec,Ev,tol,delE,Etmp,Dtmp,delEE,Egrid;
  string sbuf;

  for(i=0;i<5;i++)  getline(dos,sbuf);
  dos >> Emax >> Emin >> Ngrids >> EF; getline(dos,sbuf);
  Egrid=(Emax-Emin)/Ngrids;
  //note that the EF here is not the actual fermi energy
  vector<float> E(Ngrids),D(Ngrids);
  tol=5e-3;
  delE=1e6;delEE=delE;
  iclosest=0;
  for(i=0;i<Ngrids;i++) {
    dos >> E[i] >> D[i]; getline(dos,sbuf);
  }
  for(i=0;i<Ngrids;i++) {
    if(aurostd::abs(EF-E[i])<Egrid) {
      if(D[i]>tol && D[i+1]>tol) {isMetal=true;break;}//metal
    }
  }
  if(isMetal)  bandgap=-1;
  else {
    for(i=0;i<Ngrids;i++) {
      Dtmp=D[i];Etmp=E[i];
      if(Dtmp<tol) delEE=aurostd::abs(EF-Etmp);
      if(delEE<delE) {
        iclosest=i; delE=delEE;
      }
    }
    //    cerr << "iclosest E: " << iclosest << " " << E[iclosest] << endl;
    for(i=iclosest-1;iclosest>-1;i--) {
      if(D[i]>tol) break;
    }
    Ev=E[i];
    for(i=iclosest+1;i<Ngrids;i++) {
      if(D[i]>tol) break;
    }
    Ec=E[i];
    EF=(Ev+Ec)/2;//fix the EF
    bandgap=Ec-Ev;
    if(aurostd::abs(Ev-Emin)<tol || aurostd::abs(Emax-Ec)<tol) bandgap=-1.0;
  }
  E.clear(); D.clear();
  return bandgap;
}
// ***************************************************************************
// IsMetal
// ***************************************************************************
bool IsMetal(const string s) {
  //   Check if an element name in string s is a metal.
  bool ibool=false;
  vector<string> vmetal;
  vmetal=GetMetalsList(false);
  for(int i=0;i<(int)vmetal.size();i++)
    if(s==vmetal[i]) ibool=true;
  return(ibool);
}
// ***************************************************************************
// ParseChemicalFormula
// ***************************************************************************
void ParseChemicalFormula(string Formula,vector<string>& Zname,vector<float>& Zconc) {
  //   Parsing elements' name and concentration from a chemical formula
  //   example:
  //    Formula is MgB2
  //    output: Zname=[Mg,B], Zconc=[1,2].
  int i,j,k,n,p;
  vector<string> sbuf;
  vector<int> CapPos;

  //removing white characters from s
  Formula=RemoveCharFromString(Formula,' ');
  Formula=RemoveCharFromString(Formula,'\t');
  Formula=RemoveCharFromString(Formula,'\n');
  //splitting each atom name and its concentration pair
  CapPos.clear();
  for(i=0;i<(int)Formula.length();i++) {
    if(Formula[i]>='A' && Formula[i]<='Z') CapPos.push_back(i);
  }
  n=(int)CapPos.size();
  sbuf.clear();
  for(i=0;i<n-1;i++) {
    j=CapPos[i+1]-CapPos[i];
    sbuf.push_back(Formula.substr(CapPos[i],j));
  }
  sbuf.push_back(Formula.substr(CapPos[n-1],Formula.length()-CapPos[n-1]));
  //processing each pair
  Zname.clear();
  Zconc.clear();
  k=0;
  for(i=0;i<n;i++) {
    p=(int)sbuf[i].length();
    for(j=p;j>0;j--)
      if(sbuf[i][j]<'A') k=j;
    Zname.push_back(sbuf[i].substr(0,k));
    if(k==p) Zconc.push_back(1);
    else Zconc.push_back(aurostd::string2utype<float>(sbuf[i].substr(k,p)));
  }
}
// ***************************************************************************
// RemoveCharFromString
// ***************************************************************************
string RemoveCharFromString(const string s, char c) {
  // Remove a character from a string and return the resultant string
  string ss;
  string::iterator siter;
  ss=s;
  int i;
  for(i=0;i<(int) ss.length();i++) {
    if(ss[i]==c) {
      siter=ss.begin()+i; ss.erase(siter);i--;
    }
  }
  return ss;
}
// ***************************************************************************
// RemoveStringFromTo
// ***************************************************************************
string RemoveStringFromTo(const string s, char cstart, char cstop) {
  // Remove substring from character cstart to cstop
  // example: s="0.34(21)65"
  // ss=RemoveStringFromTo(s,'(',')');
  // ss will be s="0.3465"
  string ss=s;
  int i,istart=0,istop=0;
  bool fstart=false,fstop=false;
  for(i=0;i<(int) ss.length();i++) {
    if(ss[i]==cstart) {fstart=true;istart=i;}
    if(ss[i]==cstop && fstart==true) {fstop=true;istop=i;break;}
  }
  if(fstop) {ss.erase(istart,istop-istart+1);ss=RemoveStringFromTo(ss,cstart,cstop);}
  return ss;
}
// ***************************************************************************
// StringCrop
// ***************************************************************************
int StringCrop(string str,vector<string>& vstr) {
  //   given for example a string str=" this is an example 1 2 3 4"
  //   return the number of words using white characters as separators,
  //   and put the words in the vector<string> vstr
  //   ignore any leading white characters.
  string s,stmp;
  int i,Ns,Nw;
  //preadd and postadd a space and replacing any white character with space character
  s = " "+str+" ";
  Ns = strlen(s.c_str());
  for(i=0;i<Ns;i++) if(s[i]=='\t') s[i]=' ';
  //calculating the number of words
  Nw = 0;
  bool fbeg,fend;
  fbeg=fend=false;
  for(i=0;i<Ns;i++) {
    if(s[i]==' ' && s[i+1]!=' ') fbeg=true;
    if(s[i]!=' ' && s[i+1]==' ') fend=true;
    if(fbeg && fend) {
      Nw++;
      fbeg=fend=false;
    }
  }
  vstr.clear();
  if(Nw>0) {
    for(i=0;i<Nw;i++) {
      stmp=StringCropAt(str,i+1);
      vstr.push_back(stmp);
    }
  }
  return Nw;
}
// ***************************************************************************
// StringCropAt
// ***************************************************************************
string StringCropAt(const string input,int icrop) {
  //   Crop a word frong a string.
  //   given for example a string: "this is an example 1 2 3 4"
  //   icrop=1 : returns "this"
  //   icrop=2 : returns "is"
  //   and so on.
  bool istartfound=false,istopfound=false;
  int i=0,istart=0,ilong=0,iicrop=0;

  if(istartfound) {;} // dummy load

  istartfound=istopfound=false;
  i=istart=ilong=iicrop=0;
  if(input.size()==0) return("");
  else {
    //finding the istart to crop
    for(i=0;i<(int)input.size();i++) {
      if(i==0) {
        if(input[i]!=' ' && input[i]!='\t' && input[i]!='\n') {
          istart=i; iicrop++;
        }
      }
      else {
        if( (input[i]!=' ' && input[i]!='\t' && input[i]!='\n') &&
            (input[i-1]==' ' || input[i-1]=='\t' || input[i-1]=='\n') ) {
          istart=i; iicrop++;
        }
      }
      if(iicrop==icrop) {i=(int)input.size();}
    }
    if(iicrop==icrop) {
      //finding the length of the word to crop
      for(i=istart;i<(int)input.size();i++) {
        if(i==(int)input.size()-1) {
          if(input[i]!=' ' && input[i]!='\t' && input[i]!='\n') {
            ilong=i+1-istart;istopfound=true;i=(int)input.size();
          }
        }
        else {
          if( (input[i]!=' ' && input[i]!='\t' && input[i]!='\n') &&
              (input[i+1]==' ' || input[i+1]=='\t' || input[i+1]=='\n') ) {
            ilong=i+1-istart;istopfound=true;i=(int)input.size();
          }
        }
      }//for
      if(istopfound) return(input.substr(istart,ilong));
    }//iicrop==icrop
    return("");
  }//else
}
// ***************************************************************************
// SortFloat
// ***************************************************************************
vector<float> SortFloat(vector<float> va, int mode) {
  //   Sort vector<float>.
  //   mode = 1 (ascending)
  //   mode = -1 (descending)


  int i,j;
  float ftmp;
  vector<float> vb;
  vb=va;

  //sorting descending
  if(mode==-1) {
    ftmp=0.0;
    for(i=0;i<(int)vb.size()-1;i++) {
      for(j=i+1;j<(int)vb.size();j++) {
        if(vb[j]>vb[i]) {
          ftmp=vb[i];vb[i]=vb[j];vb[j]=ftmp;
        }
      }
    }
  }
  //sorting ascending
  else {
    ftmp=0.0;
    for(i=0;i<(int)vb.size()-1;i++) {
      for(j=i+1;j<(int)vb.size();j++) {
        if(vb[j]<vb[i]) {
          ftmp=vb[i];vb[i]=vb[j];vb[j]=ftmp;
        }
      }
    }
  }
  return vb;
}

#endif

// ***************************************************************************
// *                                                                         *
// *              AFlow WAHYU SETYAWAN - Duke University 2007-2011           *
// *                                                                         *
// ***************************************************************************
