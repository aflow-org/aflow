// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// aflow_ifrozsl.cpp
// This file contains the routines to prepare FROZSL input files (standard version).
// frozsl and frozsl_init must be present as binaries 
// written by Stefano Curtarolo and Kesong Yang

#ifndef _AFLOW_IFROZSL_CPP
#define _AFLOW_IFROZSL_CPP

#include "aflow.h"

using aurostd::StringSubst;

const string FROZSL_RELATIVE_ENERGY_DATA = "Frozsl_phonos_relative_energy.e";
const string FROZSL_DIELECTRIC_DATA = "Frozsl_DIELECTRIC.dat";
//const string FROZSL_RELATIVE_ENERGY_DATA = DEFAULT_AFLOW_FROZSL_MODES_OUT;

// ***************************************************************************
// [OBSOLETE] string STOKES_LIBRARY_DAT;
// [OBSOLETE] string STOKES_LIBRARY_DAT_file;
// ***************************************************************************
namespace FROZSL {
  bool Setup_frozsl_init_input(const string& AflowIn,ofstream &FileMESSAGE,stringstream &input_file,_aflags &aflags,_kflags &kflags);
  bool input_TO_FROZSL_poscar(ofstream &FileMESSAGE,stringstream &input_file,_aflags &aflags,_kflags &kflags);
  bool Already_Calculated_Input(const string& AflowIn) {
    return aurostd::substring2bool(AflowIn,"Harold");
  }
}

namespace KBIN {
  void VASP_RunPhonons_FROZSL( _xvasp&  xvasp,
      string  AflowIn,_aflags& aflags,_kflags& kflags,_vflags& vflags, 
      ofstream& messageFile) {
    // Test
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    if( !kflags.KBIN_PHONONS_CALCULATION_FROZSL ) return;

    FROZSL::WGET_INPUT(messageFile,AflowIn,aflags,kflags);

    if(LDEBUG) cerr << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size() << endl;
    string input_file=kflags.KBIN_PHONONS_CALCULATION_FROZSL_poscars;
    aurostd::substring2strings(input_file,vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING,_VASP_POSCAR_MODE_EXPLICIT_START_P_);  //CO20200624
    if(LDEBUG) cerr << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size() << endl;
    // some verbose
    for(uint i=0;i<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++)
      if(LDEBUG) cerr << "DEBUG= " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i) << endl;
    // load up the structures
    for(uint i=0;i<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++) {
      string START=_VASP_POSCAR_MODE_EXPLICIT_START_P_+vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i); //CO20200624
      string STOP=_VASP_POSCAR_MODE_EXPLICIT_STOP_P_+vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i); //CO20200624
      stringstream POSCAR;POSCAR.clear();POSCAR.str(std::string());
      POSCAR.str(aurostd::substring2string(input_file,START,STOP,-1));
      //[SD20220520 - OBSOLETE]if(aurostd::substring2bool(input_file,START) && aurostd::substring2bool(input_file,STOP))
        //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(input_file,POSCAR,START,STOP);
      if(!POSCAR.str().empty()) {vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.push_back(xstructure(POSCAR,IOVASP_AUTO));}
    }

    xvasp.str.species.clear();
    xvasp.str.species_pp.clear();
    vector<string> species,species_pp;

    KBIN::VASP_Produce_POTCAR(xvasp,AflowIn,messageFile,aflags,kflags,vflags);

    string potcar=aurostd::TmpFileCreate("potcar");
    aurostd::stringstream2file(xvasp.POTCAR,potcar);
    vector<string> tokens,vtitel;
    aurostd::string2vectorstring(aurostd::execute2string("cat "+potcar+" | grep TITEL"),vtitel);
    aurostd::RemoveFile(potcar);
    for(uint i=0;i<vtitel.size();i++) {
      aurostd::string2tokens(vtitel.at(i),tokens," ");
      if(tokens.size()!=4 && tokens.size()!=5) {
        string message = "POTCAR KBIN_VASP_RunPhonons_FROZSL " + aurostd::joinWDelimiter(vtitel, ", ");
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      species_pp.push_back(tokens.at(3));
      species.push_back(KBIN::VASP_PseudoPotential_CleanName(tokens.at(3)));
    }

    vector<_xvasp> vaspRuns;
    for(uint i=0;i<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++) {
      _xvasp xvasp_aus;
      xvasp_aus=xvasp;
      xvasp_aus.str=vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.at(i);
      xvasp_aus.str.species.clear(); for(uint j=0;j<species.size();j++) xvasp_aus.str.species.push_back(species.at(j));
      xvasp_aus.str.species_pp.clear();for(uint j=0;j<species_pp.size();j++) xvasp_aus.str.species_pp.push_back(species_pp.at(j));
      xvasp_aus.Directory=aflags.Directory+"/ARUN."+vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i);
      aurostd::StringSubst(xvasp_aus.Directory,"//","/");
      for(uint l=0,j=0;j<species.size();j++)
        for(uint k=0;k<(uint) xvasp_aus.str.num_each_type.at(j);k++) {
          xvasp_aus.str.atoms.at(l).name_is_given=TRUE;
          xvasp_aus.str.atoms.at(l).name=species.at(j);
          l++;
        }     
      // xvasp.str.species_pp.at(isp)=AVASP_Get_PseudoPotential_PAW_PBE(KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species.at(isp)));
      xvasp_aus.AVASP_flag_RUN_RELAX=FALSE;
      xvasp_aus.AVASP_flag_RUN_RELAX_STATIC=FALSE;
      xvasp_aus.AVASP_flag_RUN_STATIC=TRUE;
      xvasp_aus.AVASP_prototype_mode=LIBRARY_MODE_XSTRUCTURE;
      xvasp_aus.aopts.flag("AVASP_flag_RELAX_FORCES",TRUE);
      xvasp_aus.aopts.flag("FLAG::AVASP_SPIN",vflags.KBIN_VASP_FORCE_OPTION_SPIN.option);
      xvasp_aus.aopts.flag("FLAG::AVASP_CHGCAR",vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.option);
      xvasp_aus.aopts.flag("FLAG::AVASP_WAVECAR",vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.option);
      xvasp_aus.AVASP_flag_PRECISION_scheme="ACCURATE";
      xvasp_aus.aopts.flag("FLAG::PRECISION_SET",TRUE);
      xvasp_aus.AVASP_flag_ALGO_scheme=vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme;
      xvasp_aus.aopts.flag("FLAG::ALGO_SET",TRUE);
      xvasp_aus.AVASP_flag_TYPE=vflags.KBIN_VASP_FORCE_OPTION_TYPE;
      xvasp_aus.aopts.flag("FLAG::AVASP_SYMMETRY=OFF",TRUE);
      xvasp_aus.aopts.flag("FLAG::AVASP_NEIGHBORS=OFF",TRUE);
      xvasp_aus.aopts.flag("FLAG::AVASP_APL=OFF",TRUE);
      // xvasp_aus.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_STANDARD_PRIMITIVE",FALSE);
      // xvasp_aus.aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_MINKOWSKI",FALSE);
      xvasp_aus.AVASP_value_KPPRA=vflags.KBIN_VASP_KPOINTS_KPPRA.content_int;
      vaspRuns.push_back(xvasp_aus);
      // cerr << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i) << endl;
    }

    for(uint i=0;i<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++) 
      AVASP_MakeSingleAFLOWIN(vaspRuns[i]);

  }
}

//**************************************************New Code by K.S.Y****************************************************
// Generating inputfile of frozsl_init
namespace FROZSL {
  string Generate_Input_file(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) {
    ostringstream aus;
    aus << "00000  MESSAGE FROZSL running GENERATE INPUT files " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Writing original input data into a file
    string aflow_frozsl_ori_input = aflags.Directory+"/"+ "aflow.frozsl_ori_input."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".in";  //CO20200502 - threadID
    ofstream orifin;
    orifin.open(aflow_frozsl_ori_input.c_str());
    orifin << kflags.KBIN_FROZSL_STRUCTURE_STRING << endl;
    orifin.close();

    // Open the original file
    // Create a new input file for frozsl_init
    ifstream orifin_new;
    orifin_new.open(aflow_frozsl_ori_input.c_str());
    string aflow_frozsl_input = aflags.Directory+"/"+"aflow.frozsl_input."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".in"; //CO20200502 - threadID
    ofstream curfin;
    curfin.open(aflow_frozsl_input.c_str());

    // Formating the data, and adding some string
    string stmp;
    for (int i=0; i<7; i++) {
      getline(orifin_new, stmp);
      curfin << stmp << endl;
    }

    // string FROZSL_ENERGY_STRING = "! name of file containing energies";
    string FROZSL_ENERGY_STRING = FROZSL_RELATIVE_ENERGY_DATA;
    string FROZSL_EIGENVECTOR_STRING = " ! name of file containing eigenvectors (none)";
    string FROZSL_SUBGROUP = " ! name of file containing subgroup information (none)";

    string FROZSL_DIELECTRIC_STRING;
    if(kflags.KBIN_FROZSL_DIELECTRIC_STRING.length()==0) {
      FROZSL_DIELECTRIC_STRING = " ! name of file containing effective charge tensors";
    } else {
      FROZSL_DIELECTRIC_STRING = FROZSL_DIELECTRIC_DATA;
      aurostd::string2file(kflags.KBIN_FROZSL_DIELECTRIC_STRING, FROZSL_DIELECTRIC_STRING);
    }

    curfin << FROZSL_ENERGY_STRING << endl;
    curfin << FROZSL_EIGENVECTOR_STRING << endl;
    curfin << FROZSL_SUBGROUP << endl;
    curfin << FROZSL_DIELECTRIC_STRING << endl;

    while (getline(orifin_new, stmp)) {
      curfin << stmp << endl;
    }

    curfin.close();
    orifin_new.close();

    string command;
    command.clear();
    aurostd::RemoveFile(aflow_frozsl_ori_input);
    return aflow_frozsl_input;
  }
}

//**************************************************New Code by K.S.Y****************************************************
namespace FROZSL {
  bool WGET_INPUT(ofstream &FileMESSAGE,string AflowIn,_aflags &aflags,_kflags &kflags) {
    ostringstream aus;
    aus << "00000  MESSAGE FROZSL " << _AFLOWIN_ << " self-modification for input files " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // [OBSOLETE] string aflowin=aflags.Directory+"/"+_AFLOWIN_;
    // [OBSOLETE] if(!aurostd::FileExist(aflowin)) {FileMESSAGE << "ERROR" << ": file not found " << aflowin << endl;return FALSE;}
    // [OBSOLETE] string AflowIn;aurostd::file2string(aflowin,AflowIn);

    kflags.KBIN_PHONONS_CALCULATION_FROZSL_output="";
    kflags.KBIN_PHONONS_CALCULATION_FROZSL_poscars="";

    // CHECK FOR INSIDE STUFF
    if(!FROZSL::Already_Calculated_Input(AflowIn)) {
      stringstream input_file; input_file.clear();input_file.str(std::string());
      stringstream input_file_aus; input_file_aus.clear();input_file_aus.str(std::string());

      FROZSL::Setup_frozsl_init_input(AflowIn,FileMESSAGE,input_file,aflags,kflags); // must know about stuff before

      // MAKE FROZSL INPUT
      string FROZSL_INPUT= FROZSL::Generate_Input_file(FileMESSAGE,aflags,kflags);   
      aus << "00000  MESSAGE FROZSL loading data_XXXX.txt" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      FROZSL::Write("data_space.txt;data_wyckoff.txt;data_images.txt;data_irreps.txt;data_little.txt;data_isotropy.txt;symmetry2.dat;const.dat",aflags.Directory);

      aus << "00000  MESSAGE FROZSL running frozsl_init" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      kflags.KBIN_PHONONS_CALCULATION_FROZSL_output=aurostd::execute2string("cat "+FROZSL_INPUT+" | "+XHOST.command("frozsl_init"));
      aurostd::string2file("[AFLOW_FROZSL]CALC.START\n",aflags.Directory+"/"+_AFLOWIN_,"POST");
      aurostd::string2file(kflags.KBIN_PHONONS_CALCULATION_FROZSL_output+"\n",aflags.Directory+"/"+_AFLOWIN_,"POST");
      aurostd::string2file("[AFLOW_FROZSL]CALC.STOP\n",aflags.Directory+"/"+_AFLOWIN_,"POST");
      //    cout << kflags.KBIN_PHONONS_CALCULATION_FROZSL_output << endl;

      // make the POSCARS
      string AflowIn2=AflowIn;
      string file_frozsl_input=aflags.Directory+"/aflow.frozsl_input."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";  //CO20200502 - threadID
      string file_frozsl_perl=aflags.Directory+"/aflow.frozsl_perl."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";  //CO20200502 - threadID
      aurostd::string2file(kflags.KBIN_PHONONS_CALCULATION_FROZSL_output+"\n",file_frozsl_input);

      //postprocess    
      if(aurostd::FileExist(FROZSL_INPUT)) aurostd::RemoveFile(FROZSL_INPUT);
      if(aurostd::FileExist(FROZSL_DIELECTRIC_DATA)) aurostd::RemoveFile(FROZSL_DIELECTRIC_DATA);

      FROZSL::Delete("data_space.txt;data_wyckoff.txt;data_images.txt;data_irreps.txt;data_little.txt;data_isotropy.txt;symmetry2.dat;const.dat",aflags.Directory);
      if(XHOST_FROZSL_phvaspsetup_AFLOW.length()==0) init::InitGlobalObject("FROZSL_phvaspsetup_AFLOW");
      aurostd::string2file(XHOST_FROZSL_phvaspsetup_AFLOW,file_frozsl_perl);

      aurostd::ChmodFile("755",file_frozsl_perl);

      kflags.KBIN_PHONONS_CALCULATION_FROZSL_poscars=aurostd::execute2string("cat "+file_frozsl_input+" | "+file_frozsl_perl);
      aurostd::string2file("[AFLOW_FROZSL]STRUCTURES\n",aflags.Directory+"/"+_AFLOWIN_,"POST");
      aurostd::string2file("[AFLOW_MODE_POSTSCRIPT]\n",aflags.Directory+"/"+_AFLOWIN_,"POST");
      aurostd::string2file(kflags.KBIN_PHONONS_CALCULATION_FROZSL_poscars+"\n",aflags.Directory+"/"+_AFLOWIN_,"POST");

      aurostd::RemoveFile(file_frozsl_input);
      aurostd::RemoveFile(file_frozsl_perl);

      aus << "00000  MESSAGE FROZSL Clean up \" aflow --clean -D ./\" " << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "00000  MESSAGE FROZSL Restart. " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      return TRUE;
    } else {
      aus << "00000  MESSAGE FROZSL input file already created " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return FALSE;
    }
    return FALSE;
  }
}

namespace FROZSL {
  string Relative_Energies(string input);
}

namespace FROZSL {
  bool WGET_OUTPUT(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) {
    string message = "";
    ostringstream aus;
    aus << "00000  MESSAGE FROZSL running WGET OUTPUT files " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    string aflowin=aflags.Directory+"/"+_AFLOWIN_;

    if(!aurostd::FileExist(aflowin)) {FileMESSAGE << "ERROR" << ": file not found " << aflowin << endl;return FALSE;}
    string AflowIn;aurostd::file2string(aflowin,AflowIn);

    if(XHOST.vext.size()!=XHOST.vcat.size()) {
      message = "XHOST.vext.size()!=XHOST.vcat.size().";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }

    // CHECK FOR INSIDE STUFF
    if(FROZSL::Already_Calculated_Input(AflowIn)) {
      stringstream input_file; input_file.clear();input_file.str(std::string());
      vector<string> vdirectories,vfiles,vcommands;

      //Extract Structure and Dielectric data
      FROZSL::Setup_frozsl_init_input(AflowIn,FileMESSAGE,input_file,aflags,kflags); // must know about stuff before

      string command,aflow_frozsl_out;

      command="grep \"\\[VASP_POSCAR_MODE_EXPLICIT\\]START\\.\" "+aflags.Directory+"/"+_AFLOWIN_;
      aflow_frozsl_out=aurostd::execute2string(command);
      aurostd::StringSubst(aflow_frozsl_out,_VASP_POSCAR_MODE_EXPLICIT_START_P_,""); //CO20200624
      aurostd::string2vectorstring(aflow_frozsl_out,vdirectories);

      for(uint i=0;i<vdirectories.size();i++) {
        // cout << vdirectories.at(i) << "hihi" << endl;
        for(uint iext=0;iext<XHOST.vext.size();iext++) {
          if(aurostd::FileExist(aflags.Directory+"/ARUN."+vdirectories.at(i)+"/"+DEFAULT_AFLOW_QMVASP_OUT+XHOST.vext.at(iext)+"")) {
            vfiles.push_back("./ARUN."+vdirectories.at(i)+"/"+DEFAULT_AFLOW_QMVASP_OUT+XHOST.vext.at(iext)+"");
            vcommands.push_back(XHOST.vcat.at(iext)+" "+aflags.Directory+"/ARUN."+vdirectories.at(i)+"/"+DEFAULT_AFLOW_QMVASP_OUT+XHOST.vext.at(iext)+""+" | grep \"H=\"");// >> aflow.frozsl.out");
          }
        }
      }
      for(uint i=0;i<vfiles.size();i++) {
        if(aurostd::FileExist(vfiles.at(i))) {
          aus << "00000  MESSAGE FROZSL file OK: " << vfiles.at(i) << "" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          message = "file not found " + vfiles[i];
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
        }
      }
      // ENERGIES
      aus << "00000  MESSAGE FROZSL Frozen phonons energies: " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      stringstream aflow_frozsl_out_stringstream;
      aflow_frozsl_out_stringstream.str(std::string());
      for(uint i=0;i<vfiles.size();i++)
        aflow_frozsl_out_stringstream << vdirectories.at(i) << " : " << aurostd::execute2string(vcommands.at(i));// << endl;
      aus << aflow_frozsl_out_stringstream.str();
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // RELATIVE ENERGIES
      aus << "00000  MESSAGE FROZSL Frozen phonons relative energies: " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      stringstream aflow_frozsl_modes_stringstream;
      aflow_frozsl_modes_stringstream.str(std::string());
      aflow_frozsl_modes_stringstream << FROZSL::Relative_Energies(aflow_frozsl_out_stringstream.str());
      aus << aflow_frozsl_modes_stringstream.str();
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // EIGENVALUES
      aus << "00000  MESSAGE FROZSL Frozen phonons eigenvalues: " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      //**************************************************New Code by K.S.Y****************************************************
      string _dielectric=kflags.KBIN_FROZSL_DIELECTRIC_STRING;  //Dielectric data
      //   string _energies=aflow_frozsl_modes_stringstream.str();   //Frozen phonos relative energies
      string _energies; 
      if(aurostd::FileExist(DEFAULT_AFLOW_FROZSL_MODES_OUT)) aurostd::file2string(DEFAULT_AFLOW_FROZSL_MODES_OUT,_energies);
      if(aurostd::EFileExist(DEFAULT_AFLOW_FROZSL_MODES_OUT)) aurostd::efile2string(DEFAULT_AFLOW_FROZSL_MODES_OUT,_energies);

      string FROZSL_RELATIVE_ENERGY = aflags.Directory+"/"+FROZSL_RELATIVE_ENERGY_DATA;
      aurostd::string2file(_energies, FROZSL_RELATIVE_ENERGY);

      string FROZSL_INPUT= FROZSL::Generate_Input_file(FileMESSAGE,aflags,kflags);

      FROZSL::Write("data_space.txt;data_wyckoff.txt;data_images.txt;data_irreps.txt;data_little.txt;data_isotropy.txt;symmetry2.dat;const.dat",aflags.Directory);

      stringstream aflow_frozsl_eigen_stringstream;
      command.clear();
      command = XHOST.command("frozsl_init")+" < "+FROZSL_INPUT+" | frozsl ";
      aflow_frozsl_eigen_stringstream << aurostd::execute2string(command);

      //Output the calculated results
      aurostd::stringstream2file(aflow_frozsl_eigen_stringstream,aflags.Directory+"/"+DEFAULT_AFLOW_FROZSL_EIGEN_OUT);
      aus << aflow_frozsl_eigen_stringstream.str();
      cout << aflow_frozsl_eigen_stringstream.str() << endl;

      FROZSL::Delete("data_space.txt;data_wyckoff.txt;data_images.txt;data_irreps.txt;data_little.txt;data_isotropy.txt;symmetry2.dat;const.dat",aflags.Directory);

      if(aurostd::FileExist(FROZSL_INPUT)) aurostd::RemoveFile(FROZSL_INPUT);
      if(aurostd::FileExist(FROZSL_DIELECTRIC_DATA))  aurostd::RemoveFile(FROZSL_DIELECTRIC_DATA);
      if(aurostd::FileExist(FROZSL_RELATIVE_ENERGY_DATA))  aurostd::RemoveFile(FROZSL_RELATIVE_ENERGY_DATA);
      //**************************************************New Code by K.S.Y****************************************************

      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aus << "00000  MESSAGE FROZSL calculation finished " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aus << "00000  MESSAGE FROZSL you have to generate the input and run it " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    return FALSE;
  }
}


// ---------------------------------------------------------------------------
namespace FROZSL {
  string Relative_Energies(string input) {
    ostringstream oss;
    oss.clear();
    oss.setf(std::ios::fixed,std::ios::floatfield);
    uint _precision_=10; //was 16 SC 10 DM
    oss.precision(_precision_);

    vector<string> vinput,tokens,Mrefs;
    vector<double> Erefs;
    aurostd::string2vectorstring(input,tokens);
    // take the good ones
    for(uint i=0;i<tokens.size();i++)
      if(aurostd::substring2bool(tokens.at(i),":"))
        vinput.push_back(tokens.at(i));
    // find m0a0
    for(uint i=0;i<vinput.size();i++) {
      if(aurostd::substring2bool(vinput.at(i),"m0a0")) {
        aurostd::string2tokens(vinput.at(i),tokens,"m0a0");
        Mrefs.push_back(tokens.at(0)+"m");
        aurostd::string2tokens(vinput.at(i),tokens,"=");
        Erefs.push_back(aurostd::string2utype<double>(tokens.at(tokens.size()-1)));
      }
    }
    //  for(uint i=0;i<Mrefs.size();i++) cout << Mrefs.at(i) << " " << Erefs.at(i) << endl;
    // now print
    double frozsl_eps=1.0e-6;
    double ev2hartree=1.0/27.211383;
    for(uint i=0;i<vinput.size();i++) {
      for(uint j=0;j<Mrefs.size();j++) {
        if(aurostd::substring2bool(vinput.at(i),Mrefs.at(j))) {
          //cout << Mrefs.at(j) << " " << Erefs.at(j) << endl;
          aurostd::string2tokens(vinput.at(i),tokens,"=");
          if(abs(aurostd::string2utype<double>(tokens.at(tokens.size()-1))-Erefs.at(j))>frozsl_eps)
          {
            oss << aurostd::PaddedPOST(aurostd::utype2string(1000*ev2hartree*(aurostd::string2utype<double>(tokens.at(tokens.size()-1))-Erefs.at(j))),23," ") << " ! ";
            aurostd::string2tokens(vinput.at(i),tokens);
            oss << tokens.at(0) << " - " << Mrefs.at(j) << "0a0" << endl;
          }
        }
      }
    }
    return oss.str();
  }
}

// ---------------------------------------------------------------------------
namespace FROZSL {
  bool Extract_INPUT(const string& AflowIn,ofstream &FileMESSAGE,stringstream &input_file,_aflags &aflags,_kflags &kflags) {
    ostringstream aus;
    bool Krun=TRUE;

    input_file.clear();input_file.str(std::string());
    // aus << "00000  MESSAGE FROZSL from [AFLOW_FROZSL]CALC " << Message(__AFLOW_FILE__,aflags) << endl;
    // aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // CHECK FOR INSIDE STUFF
    if(!FROZSL::Already_Calculated_Input(AflowIn)) {
      FROZSL::WGET_INPUT(FileMESSAGE,AflowIn,aflags,kflags);
      return Krun;
    }

    aurostd::ExtractJustAfterToStringstreamEXPLICIT(AflowIn,input_file,"[AFLOW_FROZSL]STRUCTURES");
    aurostd::stringstream2file(input_file,string(aflags.Directory+"/"+DEFAULT_AFLOW_FROZSL_POSCAR_OUT));

    // CREATION OF FOUTPUT
    if(FROZSL::Already_Calculated_Input(AflowIn)) {
      stringstream foutput;
      foutput.str(std::string());
    }
    return Krun;
  }
}

namespace FROZSL {
  bool input_TO_FROZSL_poscar(ofstream &FileMESSAGE,stringstream &input_file,_aflags &aflags,_kflags &kflags) {
    string aflowinfake="";
    aflowinfake=aflowinfake+"[AFLOW_FROZSL]CALC "+"\n"+input_file.str();
    input_file.clear();input_file.str(std::string());
    return FROZSL::Extract_INPUT(aflowinfake,FileMESSAGE,input_file,aflags,kflags);
  }
}
// ---------------------------------------------------------------------------
namespace FROZSL {
  bool Setup_frozsl_init_input(const string& AflowIn,ofstream &FileMESSAGE,stringstream &input_file,_aflags &aflags,_kflags &kflags) {
    ostringstream aus;
    bool Krun=TRUE;
    stringstream ssfrozslSTRUCTURE,ssfrozslENERGY,ssfrozslDIELECTRIC;
    ssfrozslSTRUCTURE.clear();ssfrozslENERGY.clear();ssfrozslDIELECTRIC.clear();

    //bool flagSTRUCTURE=FALSE,flagENERGY=FALSE,flagDIELECTRIC=FALSE;
    input_file.clear();input_file.str(std::string());
    kflags.KBIN_FROZSL_STRUCTURE_STRING="";
    kflags.KBIN_FROZSL_DIELECTRIC_STRING="";
    kflags.KBIN_FROZSL_DIELECTRIC_ZEFF=FALSE;

    // SOME VERBOSE
    aus << "00000  MESSAGE FROZSL from [AFLOW_FROZSL]CALC " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // GET STRUCTURE
    ssfrozslSTRUCTURE.str(aurostd::substring2string(AflowIn,"[FROZSL_STRUCTURE_FILE]",0));
    kflags.KBIN_FROZSL_STRUCTURE_STRING=ssfrozslSTRUCTURE.str();
    kflags.KBIN_FROZSL_STRUCTURE_MODE_FILE=!kflags.KBIN_FROZSL_STRUCTURE_STRING.empty();
    if(kflags.KBIN_FROZSL_STRUCTURE_MODE_FILE) {
      aus << "00000  MESSAGE FROZSL found FROZSL_STRUCTURE_MODE_FILE " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aus << "00000  MESSAGE FROZSL_STRUCTURE generation EXPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn,ssfrozslSTRUCTURE,"[FROZSL_STRUCTURE_FILE]");
      //    aus << kflags.KBIN_FROZSL_STRUCTURE_STRING;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //flagSTRUCTURE=TRUE;
    }
    ssfrozslSTRUCTURE.str(aurostd::substring2string(AflowIn,"[FROZSL_MODE_EXPLICIT]START.FROZSL_STRUCTURE","[FROZSL_MODE_EXPLICIT]STOP.FROZSL_STRUCTURE",0));
    kflags.KBIN_FROZSL_STRUCTURE_STRING=ssfrozslSTRUCTURE.str();
    kflags.KBIN_FROZSL_STRUCTURE_MODE_EXPLICIT_START_STOP=!kflags.KBIN_FROZSL_STRUCTURE_STRING.empty();
    if(kflags.KBIN_FROZSL_STRUCTURE_MODE_EXPLICIT_START_STOP) {
      aus << "00000  MESSAGE FROZSL found FROZSL_STRUCTURE_MODE_EXPLICIT_START_STOP " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aus << "00000  MESSAGE FROZSL_STRUCTURE generation EXPLICIT file from " << _AFLOWIN_ << " with START/STOP  " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn,ssfrozslSTRUCTURE,"[FROZSL_MODE_EXPLICIT]START.FROZSL_STRUCTURE","[FROZSL_MODE_EXPLICIT]STOP.FROZSL_STRUCTURE");
      //   aus << kflags.KBIN_FROZSL_STRUCTURE_STRING;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //flagSTRUCTURE=TRUE;
    }
    // NO STRUCTURE
    if(kflags.KBIN_FROZSL_STRUCTURE_MODE_FILE==FALSE && kflags.KBIN_FROZSL_STRUCTURE_MODE_EXPLICIT_START_STOP==FALSE) {
      aus << "EEEEE  [FROZSL_MODE_EXPLICIT] do not confuse aflow !! " << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  [FROZSL_MODE_EXPLICIT] Possible modes " << Message(__AFLOW_FILE__,aflags) << endl;
      //     aus << "----------------------------------------------------------------------------------------------------" << endl;
      //     aus << "[AFLOW] FROSZL EXPLICIT MODE without START/STOP" << endl;
      //     aus << "[FROZSL_STRUCTURE_FILE]C2Ca2O3_calcite_real  ! title line" << endl;
      //     aus << "[FROZSL_STRUCTURE_FILE]167 ! space group number (1-230)" << endl;
      //     aus << "[FROZSL_STRUCTURE_FILE]9.3910974581 9.3910974581 30.9264884447 90.000 90.000 120.000   ! lattice " << endl;
      //     aus << "[FROZSL_STRUCTURE_FILE]1 ! number of displacement amplitudes for frozen phonons" << endl;
      //     aus << "[FROZSL_STRUCTURE_FILE]0.01 ! displacement amplitudes" << endl;
      //     aus << "[FROZSL_STRUCTURE_FILE]1 ! number of terms in fitting polynomial" << endl;
      //     aus << "[FROZSL_STRUCTURE_FILE]2 ! powers of polynomial terms" << endl;
      //     aus << "[FROZSL_STRUCTURE_FILE]3  ! number of Wyckoff positions." << endl;
      //     aus << "[FROZSL_STRUCTURE_FILE]C 12.011  a  0.00000000000000  0.00000000000000  0.25000000000000" << endl;
      //     aus << "[FROZSL_STRUCTURE_FILE]Ca 40.08  b  0.00000000000000  0.00000000000000  0.00000000000000" << endl;
      //     aus << "[FROZSL_STRUCTURE_FILE]O 15.9994 e  0.25879829839393  0.00000000000000  0.25000000000000" << endl;
      //     aus << "[FROZSL_STRUCTURE_FILE]1  ! number of k vectors: symbol and parameters a,b,c (if needed)" << endl;
      //     aus << "[FROZSL_STRUCTURE_FILE]GM" << endl;
      //     aus << "[FROZSL_STRUCTURE_FILE]0 " << endl;
      //     aus << "[AFLOW]" << endl;
      aus << "----------------------------------------------------------------------------------------------------" << endl;
      aus << "[AFLOW] POSCAR EXPLICIT MODE with START/STOP defaulkt" << endl;
      aus << "[FROZSL_MODE_EXPLICIT]START.FROZSL_STRUCTURE" << endl;
      aus << "C2Ca2O3_calcite_real  ! title line" << endl;
      aus << "167 ! space group number (1-230)" << endl;
      aus << "9.3910974581 9.3910974581 30.9264884447 90.000 90.000 120.000   ! lattice " << endl;
      aus << "1 ! number of displacement amplitudes for frozen phonons" << endl;
      aus << "0.01 ! displacement amplitudes" << endl;
      aus << "1 ! number of terms in fitting polynomial" << endl;
      aus << "2 ! powers of polynomial terms" << endl;
      aus << "3  ! number of Wyckoff positions." << endl;
      aus << "C 12.011  a  0.00000000000000  0.00000000000000  0.25000000000000" << endl;
      aus << "Ca 40.08  b  0.00000000000000  0.00000000000000  0.00000000000000" << endl;
      aus << "O 15.9994 e  0.25879829839393  0.00000000000000  0.25000000000000" << endl;
      aus << "1  ! number of k vectors: symbol and parameters a,b,c (if needed)" << endl;
      aus << "GM" << endl;
      aus << "0 " << endl;
      aus << "[FROZSL_MODE_EXPLICIT]STOP.FROZSL_STRUCTURE" << endl;
      aus << "----------------------------------------------------------------------------------------------------" << endl;
      aus << "EEEEE  [FROZSL_MODE_EXPLICIT] Note " << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  [FROZSL_MODE_EXPLICIT]START.FROZSL_STRUCTURE must be present and no [FROZSL_FILE] " << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  [FROZSL_MODE_EXPLICIT]STOP.FROZSL_STRUCTURE  must be present and no [FROZSL_FILE] " << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  or [FROZSL_FILE] present and NO START/STOP " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }

    // GET DIELECTRIC
    ssfrozslDIELECTRIC.str(aurostd::substring2string(AflowIn,"[FROZSL_DIELECTRIC_FILE]",0));
    kflags.KBIN_FROZSL_DIELECTRIC_STRING=ssfrozslDIELECTRIC.str();
    kflags.KBIN_FROZSL_DIELECTRIC_MODE_FILE=!kflags.KBIN_FROZSL_DIELECTRIC_STRING.empty();
    if(kflags.KBIN_FROZSL_DIELECTRIC_MODE_FILE) {
      aus << "00000  MESSAGE FROZSL found FROZSL_DIELECTRIC_MODE_FILE " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aus << "00000  MESSAGE FROZSL_DIELECTRIC generation EXPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      kflags.KBIN_FROZSL_DIELECTRIC_ZEFF=TRUE;
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn,ssfrozslDIELECTRIC,"[FROZSL_DIELECTRIC_FILE]");
      //    aus << kflags.KBIN_FROZSL_DIELECTRIC_STRING;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //flagDIELECTRIC=TRUE;
    }
    ssfrozslDIELECTRIC.str(aurostd::substring2string(AflowIn,"[FROZSL_MODE_EXPLICIT]START.FROZSL_DIELECTRIC","[FROZSL_MODE_EXPLICIT]STOP.FROZSL_DIELECTRIC",0));
    kflags.KBIN_FROZSL_DIELECTRIC_STRING=ssfrozslDIELECTRIC.str();
    kflags.KBIN_FROZSL_DIELECTRIC_MODE_EXPLICIT_START_STOP=!kflags.KBIN_FROZSL_DIELECTRIC_STRING.empty();
    if(kflags.KBIN_FROZSL_DIELECTRIC_MODE_EXPLICIT_START_STOP) {
      aus << "00000  MESSAGE FROZSL found FROZSL_DIELECTRIC_MODE_EXPLICIT_START_STOP " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aus << "00000  MESSAGE FROZSL_DIELECTRIC generation EXPLICIT file from " << _AFLOWIN_ << " with START/STOP " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      kflags.KBIN_FROZSL_DIELECTRIC_ZEFF=TRUE;
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn,ssfrozslDIELECTRIC,"[FROZSL_MODE_EXPLICIT]START.FROZSL_DIELECTRIC","[FROZSL_MODE_EXPLICIT]STOP.FROZSL_DIELECTRIC");
      //   aus << kflags.KBIN_FROZSL_DIELECTRIC_STRING;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //flagDIELECTRIC=TRUE;
    }
    // NO DIELECTRIC
    if(kflags.KBIN_FROZSL_DIELECTRIC_MODE_FILE==FALSE && kflags.KBIN_FROZSL_DIELECTRIC_MODE_EXPLICIT_START_STOP==FALSE) {
      aus << "00000  MESSAGE FROZSL No DIELECTRIC found " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      kflags.KBIN_FROZSL_DIELECTRIC_ZEFF=FALSE;
    }

    //   Krun=Krun && FROZSL_Download_FROZSL(input_file,ssfrozslSTRUCTURE,flagENERGY,ssfrozslENERGY,flagDIELECTRIC,ssfrozslDIELECTRIC);
    //   aurostd::stringstream2file(input_file,string(aflags.Directory+"/"+DEFAULT_AFLOW_FROZSL_INPUT_OUT));
    //   Krun=Krun && FROZSL::input_TO_FROZSL_poscar(FileMESSAGE,input_file,aflags,kflags);
    //   //  aurostd::stringstream2file(input_file,string(aflags.Directory+"/"+DEFAULT_AFLOW_FROZSL_POSCAR_OUT));
    //   cerr << "Exit to make some scripts" << endl;
    return Krun;
  }
}
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
namespace FROZSL {
  bool File_INPUT(const string& AflowIn,ofstream &FileMESSAGE,stringstream &input_file,_aflags &aflags,_kflags &kflags) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    bool Krun=TRUE;
    input_file.clear();input_file.str(std::string());

    if(aurostd::substring2bool(AflowIn,"[AFLOW_FROZSL]FILE=",TRUE)) {
      kflags.KBIN_FROZSL_FILE_NAME=aurostd::substring2string(AflowIn,"[AFLOW_FROZSL]FILE=",1,FALSE);
    } else {
      kflags.KBIN_FROZSL_FILE_NAME=DEFAULT_AFLOW_FROZSL_INPUT_OUT;
    }
    aus << "00000  MESSAGE FROZSL_FILE_NAME= " << kflags.KBIN_FROZSL_FILE_NAME << " " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    Krun=Krun && aurostd::FileExist(string(aflags.Directory+"/"+kflags.KBIN_FROZSL_FILE_NAME));
    if(Krun) {
      aurostd::file2stringstream(aflags.Directory+"/"+kflags.KBIN_FROZSL_FILE_NAME,input_file);
      if(LDEBUG) aurostd::stringstream2file(input_file,string(aflags.Directory+"/test_input"));
      Krun=Krun && FROZSL::input_TO_FROZSL_poscar(FileMESSAGE,input_file,aflags,kflags);
      if(LDEBUG) aurostd::stringstream2file(input_file,string(aflags.Directory+"/test_poscar"));
      //  aurostd::stringstream2file(input_file,string(aflags.Directory+"/"+DEFAULT_AFLOW_FROZSL_POSCAR_OUT));
    }
    return Krun;
  }
}

// ***************************************************************************
namespace FINDSYM {
  bool Write(string data,string directory) {
    if(data=="data_space.txt") {
      init::InitGlobalObject("FINDSYM_data_space_txt");
      if(!aurostd::FileExist(directory+"./data_space.txt")) aurostd::string2file(XHOST_FINDSYM_data_space_txt,directory+"./data_space.txt");
      return TRUE;}
    if(data=="data_wyckoff.txt") {
      init::InitGlobalObject("FINDSYM_data_wyckoff_txt");
      if(!aurostd::FileExist(directory+"./data_wyckoff.txt")) aurostd::string2file(XHOST_FINDSYM_data_wyckoff_txt,directory+"./data_wyckoff.txt");
      return TRUE;}
    return FALSE;
  }
}

namespace FROZSL {
  bool Write(string data,string directory) {
    if(aurostd::substring2bool(data,";")) {
      vector<string> vdata;
      aurostd::string2tokens(data,vdata,";");
      bool out=TRUE;
      for(uint i=0;i<vdata.size();i++) out=out && FROZSL::Write(vdata.at(i),directory);
      return out;
    }
    if(data=="data_space.txt") {
      init::InitGlobalObject("FROZSL_data_space_txt");
      if(!aurostd::FileExist(directory+"./data_space.txt")) aurostd::string2file(XHOST_FROZSL_data_space_txt,directory+"./data_space.txt");
      return TRUE;}
    if(data=="data_wyckoff.txt") {
      init::InitGlobalObject("FROZSL_data_wyckoff_txt");
      if(!aurostd::FileExist(directory+"./data_wyckoff.txt")) aurostd::string2file(XHOST_FROZSL_data_wyckoff_txt,directory+"./data_wyckoff.txt");
      return TRUE; }
    if(data=="data_images.txt") {
      init::InitGlobalObject("FROZSL_data_images_txt");
      if(!aurostd::FileExist(directory+"./data_images.txt")) aurostd::string2file(XHOST_FROZSL_data_images_txt,directory+"./data_images.txt");
      return TRUE;}
    if(data=="data_irreps.txt") {
      init::InitGlobalObject("FROZSL_data_irreps_txt");
      if(!aurostd::FileExist(directory+"./data_irreps.txt")) aurostd::string2file(XHOST_FROZSL_data_irreps_txt,directory+"./data_irreps.txt");
      return TRUE;}
    if(data=="data_isotropy.txt") {
      init::InitGlobalObject("FROZSL_data_isotropy_txt");
      if(!aurostd::FileExist(directory+"./data_isotropy.txt")) aurostd::string2file(XHOST_FROZSL_data_isotropy_txt,directory+"./data_isotropy.txt");
      return TRUE;}
    if(data=="data_little.txt") {
      init::InitGlobalObject("FROZSL_data_little_txt");
      if(!aurostd::FileExist(directory+"./data_little.txt")) aurostd::string2file(XHOST_FROZSL_data_little_txt,directory+"./data_little.txt");
      return TRUE;}
    if(data=="symmetry2.dat") {
      init::InitGlobalObject("FROZSL_symmetry2_dat");
      if(!aurostd::FileExist(directory+"./symmetry2.dat")) aurostd::string2file(XHOST_FROZSL_symmetry2_dat,directory+"./symmetry2.dat");
      return TRUE;}
    if(data=="const.dat") {
      init::InitGlobalObject("FROZSL_const_dat");
      if(!aurostd::FileExist(directory+"./const.dat")) aurostd::string2file(XHOST_FROZSL_const_dat,directory+"./const.dat");
      return TRUE;}
    if(data=="phvaspsetup_AFLOW") {
      init::InitGlobalObject("FROZSL_phvaspsetup_AFLOW");
      if(!aurostd::FileExist(directory+"./phvaspsetup_AFLOW")) aurostd::string2file(XHOST_FROZSL_phvaspsetup_AFLOW,directory+"./phvaspsetup_AFLOW");
      return TRUE;}
    if(data=="phvaspsetup_POSCAR") {
      init::InitGlobalObject("FROZSL_phvaspsetup_POSCAR");
      if(!aurostd::FileExist(directory+"./phvaspsetup_POSCAR")) aurostd::string2file(XHOST_FROZSL_phvaspsetup_POSCAR,directory+"./phvaspsetup_POSCAR");
      return TRUE;}
    return FALSE;
  }
}

namespace FROZSL {
  bool Delete(string data,string directory) {
    if(aurostd::substring2bool(data,";")) {
      vector<string> vdata;
      aurostd::string2tokens(data,vdata,";");
      bool out=TRUE;
      for(uint i=0;i<vdata.size();i++) out=out && FROZSL::Delete(vdata.at(i),directory);
      return out;
    }
    if(data=="data_space.txt") {
      if(aurostd::FileExist(directory+"./data_space.txt")) aurostd::RemoveFile(directory+"./data_space.txt");
      return TRUE;}
    if(data=="data_wyckoff.txt") {
      if(aurostd::FileExist(directory+"./data_wyckoff.txt")) aurostd::RemoveFile(directory+"./data_wyckoff.txt");
      return TRUE; }
    if(data=="data_images.txt") {
      if(aurostd::FileExist(directory+"./data_images.txt")) aurostd::RemoveFile(directory+"./data_images.txt");
      return TRUE;}
    if(data=="data_irreps.txt") {
      if(aurostd::FileExist(directory+"./data_irreps.txt")) aurostd::RemoveFile(directory+"./data_irreps.txt");
      return TRUE;}
    if(data=="data_isotropy.txt") {
      if(aurostd::FileExist(directory+"./data_isotropy.txt")) aurostd::RemoveFile(directory+"./data_isotropy.txt");
      return TRUE;}
    if(data=="data_little.txt") {
      if(aurostd::FileExist(directory+"./data_little.txt")) aurostd::RemoveFile(directory+"./data_little.txt");
      return TRUE;}
    if(data=="symmetry2.dat") {
      if(aurostd::FileExist(directory+"./symmetry2.dat")) aurostd::RemoveFile(directory+"./symmetry2.dat");
      return TRUE;}
    if(data=="const.dat") {
      if(aurostd::FileExist(directory+"./const.dat")) aurostd::RemoveFile(directory+"./const.dat");
      return TRUE;}
    if(data=="phvaspsetup_AFLOW") {
      if(aurostd::FileExist(directory+"./phvaspsetup_AFLOW")) aurostd::RemoveFile(directory+"./phvaspsetup_AFLOW");
      return TRUE;}
    if(data=="phvaspsetup_POSCAR") {
      if(aurostd::FileExist(directory+"./phvaspsetup_POSCAR")) aurostd::RemoveFile(directory+"./phvaspsetup_POSCAR");
      return TRUE;}
    return FALSE;
  }
}

#endif     // _AFLOW_IFROZSL_CPP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

//   nohup aflow --DIRECTORY=./ &
//   sleep 1
//   mv LOCK LOCK.0

//   nohup aflow --DIRECTORY=./ &
//   sleep 2
//   mv LOCK LOCK.1

//   nohup aflow --DIRECTORY=./ &
//   sleep 3
//   mv LOCK LOCK.2

//   nohup aflow --DIRECTORY=./ &
//   sleep 4
//   mv LOCK LOCK.3

//   nohup aflow --DIRECTORY=./ &
//   sleep 5
//   mv LOCK LOCK.4


