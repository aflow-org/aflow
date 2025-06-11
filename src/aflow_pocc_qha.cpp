//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow Andriy Smolyanyuk - Duke University 2021-2021           *
// *                                                                         *
//****************************************************************************
// Written by Andriy Smolyanyuk, 2021.
//
//    This file is part of Aflow.
//
//    Aflow is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Aflow is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Aflow.  If not, see <https://www.gnu.org/licenses/>.


// This file provides a framework to calculate thermal properties for
// disordered materials modeled using the POCC + QHA methodology.

#ifndef _AFLOW_POCC_QHA_CPP_
#define _AFLOW_POCC_QHA_CPP_

#include "aflow.h"
#include "aflow_pocc.h"

#define _DEBUG_POCC_QHA_ false

#define SW 5  // width of columns with blank space separator
#define TW 15 // width of columns containing label/number
#define PRECISION 10

namespace pocc {
  void writeQHAdatablock(stringstream &output, vector<vector<double> > data,
      string block_name)
  {
    string blockname = "[" + block_name + "]";

    output << AFLOWIN_SEPARATION_LINE << std::endl;
    output << blockname + "START" << std::endl;

    // print header
    output << setw(5)  << "#T[K]"  << setw(SW) << ' ' <<
      setw(TW) << "V[A^3/atom]"          << setw(SW) << ' ' <<
      setw(TW) << "F(V)[eV/atom]"        << setw(SW) << ' ' <<
      setw(TW) << "B[GPa]"               << setw(SW) << ' ' <<
      setw(TW) << "beta[10^-5/K]"        << setw(SW) << ' ' <<
      setw(TW) << "Cv(V)[kB/atom]"       << setw(SW) << ' ' <<
      setw(TW) << "Cp(V)[kB/atom]"       << setw(SW) << ' ' <<
      setw(TW) << "gamma"                << setw(SW) << ' ' <<
      setw(TW) << "Bprime";
    output << std::endl;

    // print data
    for (uint row=0; row<data.size(); row++){
      // special formatting for temperature
      output.unsetf(ios_base::floatfield);
      output << setw(5) << data[row][0] << setw(SW) << ' '
             << std::fixed << std::setprecision(PRECISION);
      // output other properties
      for (uint col=1; col<data[row].size(); col++){
        output << setw(TW) << data[row][col] << setw(SW) << ' ';
      }
      output << std::endl;
    }

    // print footer
    output << blockname + "STOP" << std::endl;
    output << AFLOWIN_SEPARATION_LINE << std::endl;
  }

  /// Calcualates POCC-average of QHA-related properties.
  void POccCalculator::calculateQHAPropertiesAVG(const vector<double>& v_temperatures) {
    string msg = "", filename = "";

    bool LDEBUG = false || _DEBUG_POCC_QHA_ || XHOST.DEBUG;

    msg = "Performing POCC+QHA post-processing step.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, msg, m_aflags.Directory,
      *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    apl::QHAmethod qha_method = apl::QHA_CALC;
    vector<string> dirs2ignore;

    // check if there are unstable structures we want to ignore
    for (uint i=0; i<m_ARUN_directories.size(); i++){
      filename = m_aflags.Directory + "/" + m_ARUN_directories[i] + "/";
      filename += DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_IMAG_FILE;

      if (apl::hasImaginary(filename, apl::QHAmethod2label(qha_method))){
        dirs2ignore.push_back(m_ARUN_directories[i]);
      }
    }

    // if there are: ignore and reload POCC data
    if (dirs2ignore.size()){
      XHOST.vflag_control.flag("ARUNS2SKIP", true);
      XHOST.vflag_control.addattachedscheme("ARUNS2SKIP",
           aurostd::joinWDelimiter(dirs2ignore, ","), true);

      loadDataIntoCalculator();
      setDFTEnergies();
    }

    // read properties
    vector<vector<vector<double> > > pocc_qha_thermo_properties;
    vector<double> line;
    uint ncols = 0;
    string data;
    for (uint i=0; i<m_ARUN_directories.size(); i++){
      data = "";

      filename = m_aflags.Directory + "/" + m_ARUN_directories[i] + "/";
      filename += DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_THERMO_FILE;
      if (LDEBUG) cerr << __AFLOW_FUNC__ << " file: " << filename << endl;

      if (aurostd::EFileExist(filename)){
         data = aurostd::substring2string(aurostd::efile2string(filename),
          "[QHA_SJ_THERMO]START", "[QHA_SJ_THERMO]STOP", 0);

        vector<string> lines = aurostd::string2vectorstring(data);
        if (LDEBUG) cerr << __AFLOW_FUNC__ << " datasize: " <<  lines.size() << endl;
        if (lines.size() > 1){
          lines.erase(lines.begin()); // remove the first line containing properties labels

          vector<vector<double> > thermo_properties;
          if (i==0){
            aurostd::string2tokens(lines[0], line);
            ncols = line.size();
          }

          for (uint j=0; j<lines.size(); j++){
            aurostd::string2tokens(lines[j], line);

            if (ncols==line.size()){
              thermo_properties.push_back(line);
              ncols = line.size();
            }
            else{
              msg="Inconsistent amount of properties (columns) among different";
              msg+="POCC-QHA calculations.";
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, msg,
                  _INDEX_MISMATCH_);
            }
          }
          pocc_qha_thermo_properties.push_back(thermo_properties);
        }
        else{
          msg = "The " + filename + " file is empty.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, msg, _FILE_ERROR_);
        }
      }
      else{
        msg = "The " + filename + " file is missing.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, msg, _FILE_NOT_FOUND_);
      }
    }

    if (!ncols) {
      msg = "POCC+QHA was not able to extract any data.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_CORRUPT_);
    }

    if (LDEBUG){
      cerr << __AFLOW_FUNC__ << " ncols: " << ncols << endl;
      cerr << __AFLOW_FUNC__ << " pocc_qha size: " << pocc_qha_thermo_properties.size();
      cerr << endl;
    }

    uint n = pocc_qha_thermo_properties.size();
    uint minsize = pocc_qha_thermo_properties[0].size();
    for (uint i=1; i<n; i++){
      if (pocc_qha_thermo_properties[i].size() < minsize)
        minsize = pocc_qha_thermo_properties[i].size();
    }

    for (uint i=0; i<n; i++){
      if (pocc_qha_thermo_properties[i].size()){
        pocc_qha_thermo_properties[i].resize(minsize);
      }
      else{
        pocc_qha_thermo_properties[i] = vector<vector<double> >(minsize,
          vector<double> (ncols));
      }
    }

    // check that temperatures are the same for all data blocks
    for (uint i=0; i<n-1; i++){
      for (uint row=0; row<minsize; row++){
        if (!aurostd::isequal(pocc_qha_thermo_properties[i][row][0],
                              pocc_qha_thermo_properties[i+1][row][0])){
          msg="Inconsistent list of temperatures among different";
          msg+="POCC-QHA calculations.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, msg,
              _VALUE_ILLEGAL_);
        }
      }
    }

    // collect degeneracies
    vector<int> degeneracies;
    int n_total = 0;
    unsigned long long int isupercell=0;
    for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin();
         it != l_supercell_sets.end(); ++it)
    {
      degeneracies.push_back((*it).getDegeneracy());
      n_total += degeneracies.back();
    }

    // average for T_pocc=inf
    vector<vector<double> > averaged_data(minsize, vector<double> (ncols));
    for (uint i=0; i<n; i++){
      for (uint row=0; row<minsize; row++){
        for (uint col=0; col<ncols; col++){
          averaged_data[row][col] += degeneracies[i]*pocc_qha_thermo_properties[i][row][col]/n_total;
        }
      }
    }

    // output averaged data
    stringstream file;
    writeQHAdatablock(file, averaged_data, "POCC_QHA_SJ_THERMO_T=INFTY");

    // average for const Tpocc, where Tpocc is a given (sintering) temperature
    double T = 0;
    for (uint i=0; i<v_temperatures.size(); i++){
      T = v_temperatures[i];
      setPOccStructureProbabilities(T);

      vector<vector<double> > averaged_qha_data(minsize, vector<double> (ncols));
      for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin();
           it != l_supercell_sets.end(); ++it)
      {
        isupercell=std::distance(l_supercell_sets.begin(),it);
        for (uint row=0; row<minsize; row++){
          for (uint col=0; col<ncols; col++){
            averaged_qha_data[row][col] += (*it).m_probability *
              pocc_qha_thermo_properties[isupercell][row][col];
          }
        }
      }
      writeQHAdatablock(file, averaged_qha_data, "POCC_QHA_SJ_THERMO_T=" +
          aurostd::utype2string(T));
    }

    // average for Tpocc = T using the DFT energies to calculate the probabilities
    vector<vector<double> > averaged_qha_data_T(minsize, vector<double> (ncols));
    for (uint row=0; row<minsize; row++){
      T = pocc_qha_thermo_properties[0][row][0];
      setPOccStructureProbabilities(T);

      for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin();
           it != l_supercell_sets.end(); ++it)
      {
        isupercell=std::distance(l_supercell_sets.begin(),it);
        for (uint col=0; col<ncols; col++){
          averaged_qha_data_T[row][col] += (*it).m_probability *
                  pocc_qha_thermo_properties[isupercell][row][col];
        }
      }
    }
    writeQHAdatablock(file, averaged_qha_data_T, "POCC_QHA_SJ_THERMO_E_T=T");

    // average for Tpocc = T using the free energies to calculate the probabilities
    vector<vector<double> > averaged_qha_data_F_T(minsize, vector<double> (ncols));
    for (uint row=0; row<minsize; row++){
      m_energy_dft_ground = AUROSTD_MAX_DOUBLE;
      for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin();
           it != l_supercell_sets.end(); ++it)
      {
        isupercell=std::distance(l_supercell_sets.begin(),it);
        (*it).m_energy_dft = pocc_qha_thermo_properties[isupercell][row][2];
        m_energy_dft_ground = std::min(m_energy_dft_ground, (*it).m_energy_dft);
      }

      T = pocc_qha_thermo_properties[0][row][0];
      setPOccStructureProbabilities(T);

      for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin();
           it != l_supercell_sets.end(); ++it)
      {
        isupercell=std::distance(l_supercell_sets.begin(),it);
        for (uint col=0; col<ncols; col++){
          averaged_qha_data_F_T[row][col] += (*it).m_probability *
                  pocc_qha_thermo_properties[isupercell][row][col];
        }
      }
    }
    writeQHAdatablock(file, averaged_qha_data_F_T, "POCC_QHA_SJ_THERMO_F_T=T");

    // calculate the thermal expansion as a logarithmic derivative of volume
    double dT = (pocc_qha_thermo_properties[0][minsize-1][0] -
                 pocc_qha_thermo_properties[0][0][0])/(minsize - 1);

    xvector<double> logV(minsize);
    for (uint row=0; row<minsize; row++){
      logV[row+1] = log(averaged_qha_data_T[row][1]);
    }
    xvector<double> beta = aurostd::diffSG(logV, dT);
    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << "[POCC_QHA_SJ_BETA_E]START" << std::endl;
    file << setw(5)  << "#T[K]"         << setw(SW) << ' ' <<
            setw(TW) << "beta[10^-5/K]" << std::endl;
    for (uint row=0; row<minsize; row++){
      file.unsetf(ios_base::floatfield);
      file << setw(5)    << averaged_qha_data_T[row][0] << setw(SW) << ' '
           << std::fixed << std::setprecision(PRECISION)
           << setw(TW)   << beta[row+1]*1e5
           << std::endl;
    }
    file << "[POCC_QHA_SJ_BETA_E]STOP" << std::endl;
    file << AFLOWIN_SEPARATION_LINE << std::endl;

    for (uint row=0; row<minsize; row++){
      logV[row+1] = log(averaged_qha_data_F_T[row][1]);
    }
    beta = aurostd::diffSG(logV, dT);
    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << "[POCC_QHA_SJ_BETA_F]START" << std::endl;
    file << setw(5)  << "#T[K]"         << setw(SW) << ' ' <<
            setw(TW) << "beta[10^-5/K]" << std::endl;
    for (uint row=0; row<minsize; row++){
      file.unsetf(ios_base::floatfield);
      file << setw(5)    << averaged_qha_data_F_T[row][0] << setw(SW) << ' '
           << std::fixed << std::setprecision(PRECISION)
           << setw(TW)   << beta[row+1]*1e5
           << std::endl;
    }
    file << "[POCC_QHA_SJ_BETA_F]STOP" << std::endl;
    file << AFLOWIN_SEPARATION_LINE << std::endl;

    string output_file = POCC_FILE_PREFIX + DEFAULT_POCC_QHA_AVGTHERMO_FILE;
    msg = "Writing the averaged POCC+QHA data to " + output_file+" file.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, msg, m_aflags.Directory,
        *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    output_file = m_aflags.Directory + "/" + output_file;

    if (!aurostd::stringstream2file(file, output_file)){
      msg = "Error writing to " + output_file + " file.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
    }
  }
}

#endif
