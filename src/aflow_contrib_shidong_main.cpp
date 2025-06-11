// ***************************************************************************
// *                                                                         *
// *               AFlow SHIDONG WANG - Duke University 2010-2011            *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_shidong_main.cpp
// functions written by
// 2010-2011: shidong.wang@duke.edu

#ifndef _AFLOW_CONTRIB_SHIDONG_MAIN_CPP_
#define _AFLOW_CONTRIB_SHIDONG_MAIN_CPP_

// ***************************************************************************
#include "aflow_contrib_shidong_main.h"
#include "aflow_pflow.h"

#include "aflow_contrib_shidong.h"
#include "aflow_contrib_shidong_auxiliary.h"

// ***************************************************************************
// Main program for cluster expansion
// pflow::AClusterExpansionMethodMain
// ***************************************************************************
namespace pflow {
  void AClusterExpansionMethodMain(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=5) {
      init::ErrorOption(options,__AFLOW_FUNC__,aurostd::liststring2string("aflow --cluster-expansion=... | --ce=structure_type,A,B,EA,EB"));
    } 

    string structure_type;
    string AlloyName; // alloy name

    cerr << "***********************************************************" << endl;
    cerr << "*              Cluster Expansion Method                   *" << endl;
    cerr << "***********************************************************" << endl;
    cerr << aflow::Banner("BANNER_TINY") << endl;    
    cerr << endl;

    vector<cestructure> str_list;
    int num_input_str = 0;

    structure_type = tokens.at(0);

    for (uint i = 0; i < structure_type.size(); i++) {
      structure_type[i] = tolower(structure_type[i]);
    }
    bool ACEFlag = ACEGetStructureType(structure_type);

    cerr << "***********************************************************" << endl;
    cerr << "*       Check the existence of all input files            *" << endl;
    cerr << "***********************************************************" << endl;

    CheckAllInputFileExistence(structure_type);


    cerr << "***********************************************************" << endl;
    cerr << "*                 Get the input structures                *" << endl;
    cerr << "***********************************************************" << endl;
    cerr << endl;


    if(ACEFlag) {
      str_list =  ReadInFitStructure(cin, structure_type);
      num_input_str = str_list.size();
    } else {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No structure",_INPUT_MISSING_); //CO20200624
    }


    if(num_input_str==0) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No input structure matches "+structure_type,_FILE_CORRUPT_); //CO20200624
    }

    vector<_ceatom> atom_species;
    for (uint i=3; i<5; i++) { // only for binary alloy
      _ceatom species_tmp;
      species_tmp.name = tokens.at(i-3+1);

      for (uint j=1; j<vatom_symbol.size(); j++) {
        if(species_tmp.name==vatom_symbol.at(j)) {
          species_tmp.volume = vatom_volume.at(j);
          break;
        }
      }

      atom_species.push_back(species_tmp);
    }

    for (uint i=0; i< atom_species.size(); i++) {
      AlloyName += atom_species.at(i).name;
    }
    cerr << "AlloyName " << AlloyName << endl;
    cerr << "Number of input structures is " << str_list.size() << endl;
    cerr << endl;

    // needed to calcualte formation energy of SL's from Ab-initio results
    double energy_A, energy_B;
    energy_A = aurostd::string2utype<double>(tokens.at(tokens.size()-2));
    energy_B = aurostd::string2utype<double>(tokens.at(tokens.size()-1));


    cerr << "***********************************************************" << endl;
    cerr << "*       Generate all clusters in a base structure         *" << endl;
    cerr << "***********************************************************" << endl;

    ceallclusters allcluster(structure_type);

    string filename;
    filename = structure_type+_FILENAME;

    ifstream myfile;
    myfile.open(filename.c_str());

    allcluster.SetCluster(filename);

    cerr << allcluster << endl;
    allcluster.PrintRepCluster();
    allcluster.PrintAllCluster();
    cerr << endl;

    cerr << "size of cecluster " << sizeof(cecluster) << endl;

    //////////////////////////////////////////////////////////////
    // Set ECI clusters
    //////////////////////////////////////////////////////////////

    ceECIcluster ECIcluster;
    ceECIcluster ECIcluster_orig;

    ECIcluster.SetUp(allcluster);
    ECIcluster_orig.SetUp(allcluster);

    bool flag_stop = false;
    int cal_count = 0;

    bool cal_eci = false;

    GenerateGNUplotScript();

    //////////////////////////////////////////////////////////////
    // iteration to find "true" ground states
    //////////////////////////////////////////////////////////////
    while(!flag_stop) {


      if(cal_count > 0) {
        str_list.clear();

        ifstream myfilein;
        myfilein.open(_FITSTRUCTUREFILE.c_str());

        str_list = ReadInFitStructure(myfilein, structure_type);

        myfilein.close();

        //RenameFiles(cal_count-1);
      } else {
        ofstream fit_file;
        fit_file.open(_FITSTRUCTUREFILE.c_str());
        WriteOutFitStructure(fit_file, AlloyName, str_list);
        fit_file.close();
      }

      cerr << "***********************************************************" << endl;
      cerr << "*              Get the correlation functions              *" << endl;
      cerr << "***********************************************************" << endl;

      cerr << "str list size " << str_list.size() << endl;

      ifstream corfilein;
      string corfilename = structure_type+_CORFILENAME;
      corfilein.open(corfilename.c_str());

      vector<cestructure>::iterator str_list_itr;
      if(!corfilein.is_open()) {
        // calcuate correlations and output them
        ofstream corfile;
        string corfilename = structure_type+_CORFILENAME;
        corfile.open(corfilename.c_str(), ios_base::out | ios_base::app);

        for (str_list_itr = str_list.begin();
            str_list_itr < str_list.end(); str_list_itr++) {
          //(*str_list_itr).SetUp(allcluster, ECIcluster_orig);
          (*str_list_itr).SetUp(allcluster, ECIcluster);
          (*str_list_itr).WriteFile(corfile);
        }

        corfile.close();
      } else {

        // read calculated correlations

        for (str_list_itr = str_list.begin();
            str_list_itr < str_list.end(); str_list_itr++) {
          (*str_list_itr).SetUp(allcluster, ECIcluster, corfilein);
        }
        corfilein.close();
      }

      /////////////////////////////////////////////////////////////////
      // get the lowest energies of input file
      vector<vector<double> > alist;
      vector<double> blist;
      srand(time(0));
      for (str_list_itr = str_list.begin();
          str_list_itr < str_list.end(); str_list_itr++) {
        blist.clear();
        blist.push_back((*str_list_itr).StoichB());
        if(cal_count==0) {
          blist.push_back((*str_list_itr).EnergyIn());
        } else {
          blist.push_back((*str_list_itr).Energy());
        }

        alist.push_back(blist);
      }

      vector<int> ground_state_candidates;
      ground_state_candidates = GroundStateCandidate(alist);
      string option = "min";
      vector<int> hull_indices;
      //hull_indices = CEConvexHull(alist, option);
      hull_indices = CEConvexHull(alist);

      // get the name of the ground state candidates

      vector<string> ground_state_names;
      vector<string> hull_state_names;

      vector<int>::iterator int_vec_itr;

      for (int_vec_itr = ground_state_candidates.begin();
          int_vec_itr < ground_state_candidates.end();
          int_vec_itr++) {

        int index = *int_vec_itr;
        double x_b = alist.at(index).at(0);
        double energy = alist.at(index).at(1);

        for (uint j=0; j<str_list.size(); j++) {
          if(str_list.at(j).StoichB()==x_b) {
            if(str_list.at(j).EnergyIn()==energy) {
              ground_state_names.push_back(str_list.at(j).Name());
              break;
            }
          } else {
            continue;
          }
        }
      }

      for (int_vec_itr=hull_indices.begin();
          int_vec_itr < hull_indices.end(); int_vec_itr++) {
        //int index = hull_indices.at(i);
        int index = *int_vec_itr;
        double x_b = alist.at(index).at(0);
        double energy = alist.at(index).at(1);

        for (uint j=0; j<str_list.size(); j++) {
          if(str_list.at(j).StoichB()==x_b) {
            if(str_list.at(j).EnergyIn()==energy) {
              hull_state_names.push_back(str_list.at(j).Name());
              break;
            }
          } else {
            continue;
          }
        }
      }

      vector<string>::iterator str_it;

      cerr << "ground state points \n";
      for (uint i=0; i<ground_state_names.size(); i++) {
        int index = ground_state_candidates.at(i);
        cerr << setw(12)
          << ground_state_names.at(i)
          << setw(12)
          << alist.at(index).at(0)
          << setw(12)
          << alist.at(index).at(1)
          << endl;
      }

      cerr << "hull points \n";
      for (uint i=0; i<hull_state_names.size(); i++) {
        int index = hull_indices.at(i);
        cerr << setw(12)
          << hull_state_names.at(i)
          << setw(12)
          << alist.at(index).at(0)
          << setw(12)
          << alist.at(index).at(1)
          << endl;
      }

      /////////////////////////////////////////////////////////////////

      cerr << "***********************************************************" << endl;
      cerr << "*                       Get ECI                           *" << endl;
      cerr << "***********************************************************" << endl;

      cerr << "size of str_list " << str_list.size() << endl;

      // calculate eci or read them from the input file

      ifstream ecifilein;
      string ecifilename = structure_type+_ECIFILENAME;
      ecifilein.open(ecifilename.c_str());

      if(cal_count==0) {
        cal_eci = (! ecifilein.is_open());
      }

      ecifilein.close();

      if(cal_eci) {
        ECIcluster = ECIcluster_orig;
        ACEOptimalECIClusters( str_list, ECIcluster);

        cerr << " Optimal ECI size " << ECIcluster.ECIValue().size() << endl;
        for (uint i=0; i<ECIcluster.ECIValue().size(); i++) {
          cerr << ECIcluster.ECIValue().at(i) << endl;
        }

        BackUpFile(ecifilename, "mv", cal_count);

        ofstream ecifile;
        ecifile.open(ecifilename.c_str());

        ECIcluster.WriteFile(ecifile);

        ecifile.close();
      } else {

        BackUpFile(ecifilename, "cp", cal_count);

        ecifilein.open(ecifilename.c_str());

        ECIcluster.ReadIn(ecifilein);

        ecifilein.close();
      }

      ofstream fout_eciresult;
      fout_eciresult.open(_ECIFILERESULT.c_str());
      ECIcluster.PrintECI(cout);
      ECIcluster.PrintECI(fout_eciresult);
      fout_eciresult.close();


      cerr << "***********************************************************" << endl;
      cerr << "*                    Get Energy by ECI                    *" << endl;
      cerr << "***********************************************************" << endl;

      // show the predict (formation) energy and entropy of each input strucutre

      for (str_list_itr = str_list.begin();
          str_list_itr < str_list.end(); str_list_itr++) {

        (*str_list_itr).SetECICluster(ECIcluster); // get the optimal ECI_cluster
        (*str_list_itr).GetECICorrelation();
        (*str_list_itr).SetECI(ECIcluster.ECIValue()); // store ECI to cestructure object
        (*str_list_itr).GetEnergy();
        cout << *str_list_itr << endl;

      }

      ofstream fout_comparison;
      fout_comparison.open(_FITCOMPARISONFILE.c_str());

      for (str_list_itr = str_list.begin();
          str_list_itr < str_list.end(); str_list_itr++) {
        (*str_list_itr).PrintOutComparison(fout_comparison);
      }

      fout_comparison.close();


      //////////////////////////////////////////////////////////////
      // random alloy
      //////////////////////////////////////////////////////////////

      // after obtain the ECI's, calculate energies of random alloy
      // and Superllatice

      ceralloy ralloy;
      double stoich_b_ra;

      int nstep = 50;
      double ds = 1.0/double(nstep);

      ofstream fout_ralloy;
      fout_ralloy.open(_RALLOYRESULT.c_str());

      for (int i=0; i<nstep+1; i++) {
        stoich_b_ra = ds * i;
        ralloy = ceralloy(structure_type, stoich_b_ra);

        cerr << endl << ralloy << endl;

        ralloy.SetUp(allcluster, ECIcluster);
        ralloy.GetECICorrelation();
        ralloy.GetEnergy();


        fout_ralloy << ralloy << endl;

      }

      fout_ralloy.close();


      //////////////////////////////////////////////////////////////
      // superlattice
      //////////////////////////////////////////////////////////////

      ceSL sl1;
      sl1 = ceSL(structure_type);

      cerr << "superlattice structre type " << structure_type << endl;

      ifstream mySLfilein;
      mySLfilein.open(_SLFILENAME.c_str());
      ACESLProperties_Readin_Corfile(mySLfilein, structure_type,  allcluster,
          ECIcluster);
      mySLfilein.close();

      string SL_name;

      // read in the SL results and output the convexhull and ground states

      mySLfilein.open(_SLFILERESULT.c_str());

      ACEGroundStateConvexHull(mySLfilein, str_list);

      mySLfilein.close();

      string gnuplot_command;
      gnuplot_command = XHOST.command("gnuplot")+" " + _GNUPLOTPTFILE;
      aurostd::execute(gnuplot_command);

      vector<string> new_hull_state, new_gnd_state;
      new_hull_state = GetNewHullState(str_list);
      new_gnd_state = GetNewGroundState(str_list);
      flag_stop =(new_hull_state.size()==0 );


      if(new_hull_state.size()==0) {
        // no new ground state structure found
        flag_stop = true;

        aurostd::execute("cp " + _FITSTRUCTUREFILE + " " + _FITSTRUCTUREFILE + "tmp");

        cerr << "no new hull state is found \n";
      } else {

        cerr << "new hull state\n";
        for (str_it = new_hull_state.begin();
            str_it < new_hull_state.end(); str_it++) {
          cerr << *str_it << endl;
        }

        for (uint i=0; i<1; i++) {
          // calculate the correlation functions
          // and insert into correlation input file

          // calculate only one new ground state

          SL_name = new_hull_state.at(i);
          sl1.SetUp(SL_name);

          sl1.SetStrCluster(allcluster);
          sl1.SetAllECICluster(ECIcluster);

          sl1.GetAllECICorrelation();
          sl1.PrintOutCorrelation();

          ofstream corfile;
          corfile.open(corfilename.c_str(), ios::app);
          sl1.WriteFile(corfile);
          corfile.close();

          // calculate the fit_quality by other method
          // here to calculate the formation energy per atom
          // by VASP
          bool mpi_flag = false;
          GenerateAflowInputFile(structure_type, AlloyName,
              SL_name, atom_species, mpi_flag);

          // backup result files
          RenameFiles(cal_count);
          MoveFiles(cal_count);

          // create a working directory
          // and insert the calculated resulst to fit structure input file
          GenerateNewTrainingItem(SL_name);

          // output the new fit structure;
          double energy, fenergy;

          energy = GetResultFromAFLOW();

          fenergy = energy - (1.0-sl1.StoichB()) * energy_A
            - sl1.StoichB() * energy_B;

          ofstream fout;
          fout.open(_FITSTRUCTUREFILE.c_str(), ios::app);

          fout.precision(8);
          fout.setf(ios_base::fixed, ios_base::floatfield);

          fout
            << "(" << AlloyName << ") "
            << setw(6)
            << SL_name
            << setw(12)
            << sl1.StoichB()
            << setw(12)
            << fenergy
            << endl;
          fout.close();

          cal_eci = true;
        }

      }

      ++cal_count;

    }

    CleanUp();
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " END" << endl;
  }
} // namespace pflow

namespace pflow {
  void SQS(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=7 && tokens.size()!=3) {
      init::ErrorOption(options,__AFLOW_FUNC__,aurostd::liststring2string("aflow --sqs=structure_type,atom_num,neighbor_num,sl_num_min,sl_num_max,A,B | --special-quasirandom-structure=...","aflow --sqs=structure_type n1 n2 < POSCAR | --special-quasirandom-structure=..."));
    } 

    string structure_type;
    structure_type = tokens.at(0);
    for (uint i = 0; i < structure_type.size(); i++) {
      structure_type[i] = tolower(structure_type[i]);
    }
    vector<int> atom_config;

    bool ACEFlag = ACEGetStructureType(structure_type);

    if(!ACEFlag) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"structure types are fcc, bcc, hcp",_INPUT_ILLEGAL_); //CO20200624
    }

    if((tokens.size() != 3) && (tokens.size() != 7)) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"incorrect argument numbers",_INPUT_NUMBER_); //CO20200624
    }

    //////////////////////////////////////////////////////////////
    // check the existences of all input files
    //////////////////////////////////////////////////////////////
    CheckAllInputFileExistence(structure_type);


    // read the clusters
    ceallclusters allcluster(structure_type);

    string filename;
    filename = structure_type+_FILENAME;

    allcluster.SetCluster(filename);

    //////////////////////////////////////////////////////////////
    // Set ECI clusters
    //////////////////////////////////////////////////////////////

    ceECIcluster ECIcluster;

    ECIcluster.SetUp(allcluster);

    int site_num = aurostd::string2utype<int>(tokens.at(1));
    int NNNum =  aurostd::string2utype<int>(tokens.at(2));

    if(tokens.size()==3) {
      // read from POSCAR
      xstructure a(cin,IOAFLOW_AUTO);

      vector<_ceatom> atom_species;
      for (uint i=0; i<a.species.size(); i++) { // only for binary alloy
        _ceatom species_tmp;
        species_tmp.name = a.species.at(i);

        for (uint j=1; j<vatom_symbol.size(); j++) {
          if(species_tmp.name==vatom_symbol.at(j)) {
            species_tmp.volume = vatom_volume.at(j);
            break;
          }
        }

        atom_species.push_back(species_tmp);
      }

      ceSL SLtmp;
      string str_type = structure_type;
      SLtmp = ceSL(str_type, a);
      //SLtmp = ceSL(str_type);

      if(a.species.size() > 2) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"only binary alloy is implemented",_INPUT_ILLEGAL_); //CO20200624
      }

      SLtmp.SetStrCluster(allcluster);
      SLtmp.SetAllECICluster(ECIcluster);
      SLtmp.GetAllECICorrelation();

      cout << setw(20)
        << "weight "
        << setw(20)
        << SLtmp.IsSQS(site_num, NNNum)
        << endl;
      SLtmp.OutputSQS(cout);
      SLtmp.PrintStructure(cout, atom_species);

    } else {
      // search from all SL's

      int min_SLcell_nr = aurostd::string2utype<int>(tokens.at(3));
      int max_SLcell_nr = aurostd::string2utype<int>(tokens.at(4));

      vector<_ceatom> atom_species;
      for (uint i=7; i<9; i++) { // only for binary alloy
        _ceatom species_tmp;
        species_tmp.name = tokens.at(i-2);

        for (uint j=1; j<vatom_symbol.size(); j++) {
          if(species_tmp.name==vatom_symbol.at(j)) {
            species_tmp.volume = vatom_volume.at(j);
            break;
          }
        }

        atom_species.push_back(species_tmp);
      }

      ifstream mySLfilein;

      cerr << "superlattice structre type " << structure_type << endl;

      mySLfilein.open(_SLFILENAME.c_str());
      ACESLPropertiesSQS_Readin_Corfile(mySLfilein, structure_type,  allcluster,
          ECIcluster, site_num, NNNum,
          min_SLcell_nr, max_SLcell_nr, atom_species);
      mySLfilein.close();

    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " END" << endl;
  }
}

namespace pflow {
  void Superlattice(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=3 && tokens.size()!=4) {
      init::ErrorOption(options,__AFLOW_FUNC__,aurostd::liststring2string("aflow --superlattice=structure_type,n_min,n_max < POSCAR","--superlattice=VASP,structure_type,A,B < superlattice_name"));
    }      

    // calculate all superlattices with atom number in a unit cell
    // from a given range
    // SL are stored in files
    //_SLFILE, _SLFILENAME, and _SLFILECOR (see aflow_contrib_shidong_auxiliary.h)

    string structure_type;
    string SL_name;

    if(tokens.at(0)=="VASP") {
      // Output VASP POSCAR from the given superlattice name

      cin >> SL_name;

      string affix = "_";
      size_t pos;
      pos = SL_name.find(affix);

      structure_type = SL_name.substr(0, pos);

    } else {

      structure_type = tokens.at(0);
      cout << structure_type << endl;
      for (uint i = 0; i < structure_type.size(); i++) {
        structure_type[i] = tolower(structure_type[i]);
      }

    }

    bool ACEFlag = ACEGetStructureType(structure_type);
    if(!ACEFlag) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No structure",_INPUT_MISSING_); //CO20200624
    }

    if(tokens.at(0)=="VASP") {
      vector<_ceatom> atom_species;
      for (uint i=3; i<5; i++) { // only for binary alloy
        _ceatom species_tmp;
        species_tmp.name = tokens.at(i-2);

        for (uint j=1; j<vatom_symbol.size(); j++) {
          if(species_tmp.name==vatom_symbol.at(j)) {
            species_tmp.volume = vatom_volume.at(j);
            break;
          }
        }

        atom_species.push_back(species_tmp);
      }

      ceSL sl1;
      sl1 = ceSL(structure_type);
      sl1.SetUp(SL_name);
      sl1.PrintStructure(cout, atom_species);

    } else {

      int min_SLcell_nr = aurostd::string2utype<int>(tokens.at(1));
      int max_SLcell_nr = aurostd::string2utype<int>(tokens.at(2));

      string filename;
      filename = structure_type+_FILENAME;

      ifstream myfile;
      myfile.open(filename.c_str());

      if(!myfile.is_open()) {

        message << "The cluster data file is missing" << endl;
        message << "Please use the command" << endl;
        message << "    aflow --cluster=structure_type,minimun_site_num,maximum_site_num,minimum_nearest_neighbor,maximum_nearest_neighbor" << endl;
        message << "to generate it";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_FILE_CORRUPT_); //CO20200624

      } else {
        myfile.close();
      }

      // read the clusters
      ceallclusters allcluster(structure_type);

      allcluster.SetCluster(filename);

      //////////////////////////////////////////////////////////////
      // Set ECI clusters
      //////////////////////////////////////////////////////////////

      ceECIcluster ECIcluster;

      ECIcluster.SetUp(allcluster);

      // get superlattice names
      ACESLProperties(cerr, structure_type, min_SLcell_nr, max_SLcell_nr, allcluster, ECIcluster);
    }

    if(LDEBUG) cerr << __AFLOW_FUNC__ << " END" << endl;
  }
}


// ***************************************************************************
// pflow::Cluster
// ***************************************************************************
namespace pflow {
  void Cluster(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=5) {
      init::ErrorOption(options,__AFLOW_FUNC__,"aflow --cluster=structure_type,n_min,n_max,m_min,m_max");
    } 
    // generate the clusters

    string structure_type;

    structure_type = tokens.at(0);

    for (uint i = 0; i < structure_type.size(); i++) {
      structure_type[i] = tolower(structure_type[i]);
    }
    bool ACEFlag = ACEGetStructureType(structure_type);
    if(!ACEFlag) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No structure",_INPUT_MISSING_); //CO20200624
    }

    ceallclusters allcluster(structure_type);

    string open_opt;
    open_opt = "new";

    string filename;
    filename = structure_type+_FILENAME;

    int site_num_low=0, site_num_up=0, NNNum_low=0, NNNum_up=0;

    site_num_low = aurostd::string2utype<int>(tokens.at(1));
    site_num_up = aurostd::string2utype<int>(tokens.at(2));
    NNNum_low = aurostd::string2utype<int>(tokens.at(3));
    NNNum_up = aurostd::string2utype<int>(tokens.at(4));
    allcluster.SetCluster(site_num_low, site_num_up, NNNum_low, NNNum_up);
    allcluster.WriteFile(filename, open_opt);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " END" << endl;
  }
}

// ***************************************************************************
// Create animation gif by using Jmol and ImageMagick
// ***************************************************************************
namespace pflow {
  void JMOLAnimation(istream& input, vector<string> argv) {
    string cif_tmp=aurostd::TmpFileCreate("pflow::JmolAnimation");
    aurostd::RemoveFile(cif_tmp);
    string jmol_script_tmp=aurostd::TmpFileCreate("pflow::JmolAnimation");
    aurostd::RemoveFile(jmol_script_tmp);
    // generate cif file
    xstructure str_in(input,IOAFLOW_AUTO);
    ofstream FoutCIF;
    FoutCIF.open(cif_tmp.c_str());
    pflow::PrintCIF(FoutCIF,str_in);
    FoutCIF.close();

    // generate script file to get gifs of each steps
    string movieScript = "\
                          save orientation; load \"\" {444 666 1} ; restore orientation;\n \
                          unitcell on; display cell=555; center visible; zoom 280;\n \
                          color background [xFFFFFF]\n \
                          /* \n \
                           *   add any other rendering commands \n \
                           *     */ \n \
                          for (var i=0; i<36; i=i+1) \n \
                            write image 1000 1000 @{\"" + cif_tmp + "_movie\" + (\"0000\" + i)[-3][0] + \".gif\"} \n \
                              /* 200 and 200 are width and height */ \n \
                              rotate axisangle {1 1 0} 10 \n \
                              /* axis is defined by X Y Z lengths between braces; this one is at \n \
                               *  * 45 degrees \n \
                               *   and 10 (degrees) is angle of rotation, so the 36-loop gives \n \
                               *    a full turn \n \
                               *     */ \n \
                              end for \n";

    ofstream FoutScript;
    FoutScript.open(jmol_script_tmp.c_str());
    FoutScript << movieScript;
    FoutScript.close();

    string jmol_command = XHOST.command("jmol")+" " + cif_tmp + " -s " + jmol_script_tmp + " -x -i -L ";

    string gif_command = XHOST.command("convert")+" -delay 20 -loop 0 " + cif_tmp + "_movie00* ";
    string gif_file = "animation_" + XHOST.ostrPID.str() + "_" + XHOST.ostrTID.str() + ".gif";  //CO20200502 - threadID

    if(argv.size()==3) { // give output file name
      gif_file = argv.at(2);
    }

    // check if this file is writable
    ofstream foutAnime;
    foutAnime.open(gif_file.c_str());
    if(!foutAnime.is_open()) { // given file is not writable
      gif_file = "/tmp/animation_" + XHOST.ostrPID.str() + "_" + XHOST.ostrTID.str() + ".gif";  //CO20200502 - threadID
      cerr << "Cannot write to " << argv.at(2) << ". ";
    }
    foutAnime.close();
    cerr << "Gif file is written to " << gif_file << endl;

    gif_command += gif_file;

    ostringstream oss;
    oss << jmol_command << " && " << gif_command;
    aurostd::execute(oss);

    aurostd::RemoveFile(cif_tmp);
    aurostd::RemoveFile(jmol_script_tmp);
    aurostd::execute("rm -f " + cif_tmp + "_movie00*");
  }
}


#endif
