// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu
#ifndef _AFLOW_AGL_ELECTRONIC_CPP
#define _AFLOW_AGL_ELECTRONIC_CPP
#include "aflow.h"
#include "aflow_agl_debye.h"


// ###############################################################################
//                  AFLOW Automatic GIBBS Library (AGL) (2013-2021)
// ###############################################################################
//
// Uses quasi-harmonic Debye model to obtain thermodynamic properties of materials
// Based on original Fortran program written by M. A. Blanco et al.
// See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details of original GIBBS program
// See C. Toher et al., Phys. Rev. B 90, 174107 (2014), Phys. Rev. 1, 015401 (2017) and references therein for description of this AGL implementation
// Please cite these works in addition to the general AFLOW papers if you use results generated using AGL
//

// **********************************************************************************************
// The following functions are for extracting electronic properties as a function of pressure
// **********************************************************************************************

// ***************************************************************************
// AGL_functions::extractedos
// ***************************************************************************
namespace AGL_functions {
  //
  // extractedos: Extract electronic density of states from the completed VASP calculations
  // Adapted from section of AFLOW APL function DirectMethodPC::runVASPCalculations()
  //
  uint extractedos(vector<_xvasp>& vaspRuns, _AGL_data& AGL_data, vector<string>& dirrunname, ofstream& FileMESSAGE) {
    bool LVERBOSE=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    vector<string> vfile, dfile, vofile, dofile;
    aurostd::string2tokens(string("DOSCAR.static.bz2,DOSCAR.static.gz,DOSCAR.static.xz,DOSCAR.static"),vfile,",");
    aurostd::string2tokens(string("OUTCAR.static.bz2,OUTCAR.static.gz,OUTCAR.static.xz,OUTCAR.static"),vofile,",");    
    xDOSCAR doscar;
    xOUTCAR outcar;
    string vfilename, dfilename, ffilename;
    bool skipdir = false;
    // [OBSOLETE] aurostd::xmatrix<double> stress_tensor(3, 3);
    AGL_dos_pressures AGL_dos_pressure;
    // [OBSOLETE] vector<AGL_dos_pressures> AGL_dos_pressure_list;
    double sumdos, sumidos;
    for(uint idVaspRun = 0; idVaspRun < vaspRuns.size(); idVaspRun++) {
      skipdir = false;
      aurostd::StringstreamClean(aus);
      // Print out total energy
      aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aurostd::string2tokens(vaspRuns.at(idVaspRun).Directory, dfile, "/");
      dfilename = dfile.at(dfile.size()-1);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Directory name = " << dfilename << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // Check if structure is on list of failed runs to be skipped
      // If so, then skip reading and continue to next structure
      for (uint ij = 0; ij < AGL_data.failed_arun_list.size(); ij++) {
        ffilename = AGL_data.failed_arun_list.at(ij);
        if(LVERBOSE) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "dfilename = " << dfilename << endl;
          aus << _AGLSTR_MESSAGE_ + "ffilename = " << ffilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(aurostd::substring2bool(dfilename,ffilename,TRUE)) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Found directory in to-skip list: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          skipdir = true;
        }
      }
      if(skipdir) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        skipdir = false;
        continue;
      }

      // If tarred and compressed directory exists...
      // [OBSOLETE]      string tarfilename = vaspRuns.at(idVaspRun).Directory + ".tar.bz2";
      // [OBSOLETE]    if( aurostd::FileExist(tarfilename) ) {
      // [OBSOLETE]   	// Extract all...
      // [OBSOLETE]   	aurostd::execute( string("tar -xf ") + tarfilename );
      // [OBSOLETE]     }
      if( aurostd::FileExist(vaspRuns.at(idVaspRun).Directory + ".tar.bz2") ) {
        aurostd::execute( string("tar -xf ") + vaspRuns.at(idVaspRun).Directory + ".tar.bz2" );
      } else if( aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.bz2") ) {
        aurostd::execute( string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.bz2" );
      } // Extract all...
      if( aurostd::FileExist(vaspRuns.at(idVaspRun).Directory + ".tar.gz") ) {
        aurostd::execute( string("tar -xf ") + vaspRuns.at(idVaspRun).Directory + ".tar.gz" );
      } else if( aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.gz") ) {
        aurostd::execute( string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.gz" );
      } // Extract all...
      if( aurostd::FileExist(vaspRuns.at(idVaspRun).Directory + ".tar.xz") ) {
        aurostd::execute( string("tar -xf ") + vaspRuns.at(idVaspRun).Directory + ".tar.xz" );
      } else if( aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.xz") ) {
        aurostd::execute( string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.xz" );
      } // Extract all...

      // If the LOCK file is missing, then it is probably a corrupted run
      // Do not accept it and wait for the new run
      if( !aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/"+_AFLOWLOCK_ )  && !aurostd::FileExist( dirrunname.at(idVaspRun) + "/"+_AFLOWLOCK_ ) &&
          !(((XHOST.ARUN_POSTPROCESS || AGL_data.postprocess) &&
              ((aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/agl.LOCK")) || (aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/LOCK")) ||
               (aurostd::FileExist( dirrunname.at(idVaspRun) + "/agl.LOCK")) || (aurostd::FileExist( dirrunname.at(idVaspRun) + "/LOCK")))))) {
        aurostd::StringstreamClean(aus);
        aus <<  _AGLSTR_WARNING_ + "The " << _AFLOWLOCK_ << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        if(AGL_data.autoskipfailedaruns) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          continue;
        } else {
          throw AGLStageBreak();
        }
      }

      //for(uint i=0;i<vfile.size()&&(outcar.outcar=="");i++)
      for(uint i=0;i<vfile.size()&&(doscar.content=="");i++) {
        if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(i))) {
          doscar.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(i));
        } else if(aurostd::FileExist(dirrunname.at(idVaspRun)+"/"+vfile.at(i))) {
          doscar.GetPropertiesFile(dirrunname.at(idVaspRun)+"/"+vfile.at(i));
        }
      }
      //if(outcar.outcar=="")
      if(doscar.content=="") {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "The DOSCAR.static file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        if(AGL_data.autoskipfailedaruns) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          continue;
        } else {
          throw AGLStageBreak();
        }
      }

      //for(uint i=0;i<vfile.size()&&(outcar.outcar=="");i++)
      for(uint i=0;i<vofile.size()&&(outcar.content=="");i++) {
        if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vofile.at(i))) {
          outcar.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vofile.at(i));
        } else if(aurostd::FileExist(dirrunname.at(idVaspRun)+"/"+vofile.at(i))) {
          outcar.GetPropertiesFile(dirrunname.at(idVaspRun)+"/"+vofile.at(i));
        }
      }
      //if(outcar.outcar=="")
      if(outcar.content=="") {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "The OUTCAR.static file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        if(AGL_data.autoskipfailedaruns) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          continue;
        } else {
          throw AGLStageBreak();
        }
      }      

      // Store pressure in units of GPa
      AGL_dos_pressure.pressure_external = outcar.pressure_residual / 10.0;
      // Store Fermi energy in units of eV
      AGL_dos_pressure.Fermi_energy = doscar.Efermi;
      // Store DOS and IDOS
      AGL_dos_pressure.energy.clear();
      AGL_dos_pressure.dosval.clear();
      AGL_dos_pressure.idosval.clear();
      for (uint i = 0; i < doscar.vDOS.at(0).at(0).size(); i++) { //CT20190626 - check vector size is correct for all spin states
        if (doscar.vDOS.at(0).at(0).at(i).size() < doscar.venergy.size()) {  //ME20190614 - new vDOS format
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "Mismatch between number of DOS values and number of energy values" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          continue;
        }
      }
      for (uint i = 0; i < doscar.viDOS.size(); i++) { //CT20190626 - check vector size is correct for all spin states
        if (doscar.viDOS.at(i).size() < doscar.venergy.size()) {  //ME20190614 - new vDOS format
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "Mismatch between number of DOS values and number of energy values" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          continue;
        }
      }
      for (uint i = 0; i < doscar.venergy.size(); i++) {
        AGL_dos_pressure.energy.push_back(doscar.venergy.at(i));
        sumdos = 0.0;
        sumidos = 0.0;
        for (uint j = 0; j < doscar.spin + 1; j++) {  //ME20190614 - new vDOS format
          sumdos = sumdos + doscar.vDOS.at(0).at(0).at(j).at(i);  //ME20190614 - new vDOS format
          sumidos = sumidos + doscar.viDOS.at(j).at(i);  //ME20190614 - new viDOS format
        }
        AGL_dos_pressure.dosval.push_back(sumdos);
        AGL_dos_pressure.idosval.push_back(sumidos);
      }
      AGL_data.AGL_edos_properties.push_back(AGL_dos_pressure);
      // [OBSOLETE] stress_tensor = -outcar.stress;
      // [OBSOLETE] AGL_data.stresscalculated.push_back(stress_tensor);
      // [OBSOLETE] AGL_data.structurecalculated.push_back(idVaspRun);

      // [OBSOLETE] aurostd::StringstreamClean(aus);
      // Print out total energy, volume and calculated residual pressure
      // [OBSOLETE] aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << ", Energy (eV) = " << AGL_data.energyinput.at(idVaspRun) << endl;
      // [OBSOLETE] aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << ", Energy (eV) = " << AGL_data.energyinput.at(AGL_data.energyinput.size()-1) << endl;
      // [OBSOLETE] aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << ", Volume (Ang^3) = " << AGL_data.volumeinput.at(AGL_data.volumeinput.size()-1) << endl;
      // [OBSOLETE] aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << ", Pressure (kB) = " << AGL_data.pressurecalculated.at(AGL_data.pressurecalculated.size()-1) << endl;      
      // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      outcar.clear();
      doscar.clear();
    }
    return 0;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::edosbandgap
// ***************************************************************************
namespace AGL_functions {
  //
  // edosbandgap: Calculate electronic band gap from density of states data
  //
  uint edosbandgap(_AGL_data& AGL_data, ofstream& FileMESSAGE) {
    bool LVERBOSE=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    double egap_tol = 1.0e-6;
    // [OBSOLETE] double energ_prev;
    // [OBSOLETE] double egapval;
    uint below_EF = 0, above_EF = 0, bottom_EF = 0, top_EF = 0, min_EF = 0;
    uint valbandedge = 0, condbandedge = 0, valregionedge = 0, condregionedge = 0;
    int j_min_EF = 0, j_below_EF = 0, j_above_EF = 0;
    // [OBSOLETE] bool energy_ordered = true;
    // [OBSOLETE] bool metallic_dos = false;
    bool ingap = false;
    bool polyfitfail = false;
    uint vallower = 0;
    uint valupper = 0;
    uint condlower = 0;
    uint condupper = 0;
    uint nminlimit, nmaxlimit;
    uint aglerror = 0;
    uint valinit = 0;
    uint condinit = 0;
    uint ienergdiffmin = 0, iter = 0, itermax = 0;
    double dvallower = 0.0, dvalupper = 0.0, dcondlower = 0.0, dcondupper = 0.0;
    double valbandmax = 0.0, condbandmin = 0.0;
    double dosminEF = 0.0, energdiffmin = 0.0;
    double valdosmin = 0.0, condosmin = 0.0, valpolymin = 0.0, condpolymin = 0.0;
    vector<double> energtofit, dostofit, energtoeval, dospolyeval, pressure_list, edosgap_list, dosvalEF_list, edosgapfit_list, dosvalEFfit_list;
    pressure_list.clear();
    edosgap_list.clear();
    for (uint i = 0; i < AGL_data.AGL_edos_properties.size(); i++) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: pressure = " << AGL_data.AGL_edos_properties.at(i).pressure_external << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // Initialize band gap and value of DOS at Fermi energy to zero
      AGL_data.AGL_edos_properties.at(i).edos_band_gap = 0.0;
      AGL_data.AGL_edos_properties.at(i).dosval_EF = 0.0;
      // Initialize record of whether band gap and DOS(EF) values have been calculated to false
      AGL_data.AGL_edos_properties.at(i).edos_gap_set = false;
      AGL_data.AGL_edos_properties.at(i).dosval_EF_set = false;
      AGL_data.AGL_edos_properties.at(i).edos_gap_poly_set = false;
      AGL_data.AGL_edos_properties.at(i).dosval_EF_poly_set = false;
      polyfitfail = false;
      // Check that energy value is monotonically increasing
      // [OBSOLETE] energ_prev = AGL_data.AGL_edos_properties.at(i).energy.at(0);
      // [OBSOLETE] energy_ordered = true;
      // [OBSOLETE] metallic_dos = false;
      energtofit.clear();
      dostofit.clear();
      // Calls qcksortev to sort in order of increasing energy
      aglerror = AGL_functions::qcksortevt(AGL_data.AGL_edos_properties.at(i).energy, AGL_data.AGL_edos_properties.at(i).dosval, AGL_data.AGL_edos_properties.at(i).idosval, FileMESSAGE);
      if (aglerror != 0) {
        return aglerror;
      }
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Sort successful " << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 	
      // [OBSOLETE] for (uint j = 0; j < AGL_data.AGL_edos_properties.at(i).energy.size()) {
      // [OBSOLETE]   if (AGL_data.AGL_edos_properties.at(i).energy.at(j) < energ_prev) {
      // [OBSOLETE]     aurostd::StringstreamClean(aus);
      // [OBSOLETE]     aus << _AGLSTR_WARNING_ + "The DOS for j =  " << j << " is not in monotonically ascending values of the energy." << endl;
      // [OBSOLETE]     aus << _AGLSTR_WARNING_ + "Skipping band gap calculation for this system." << endl;	  
      // [OBSOLETE]     aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // [OBSOLETE]     energy_ordered = false;
      // [OBSOLETE]   } else {
      // [OBSOLETE]     energ_prev = AGL_data.AGL_edos_properties.at(i).energy.at(j);
      // [OBSOLETE]   }
      // [OBSOLETE] }
      // [OBSOLETE] if (energy_ordered)
      // Find energy points immediately above and below Fermi energy
      below_EF = 0;
      above_EF = 1;
      for (uint j = 0; j < AGL_data.AGL_edos_properties.at(i).energy.size(); j++) {	
        if (AGL_data.AGL_edos_properties.at(i).energy.at(j) < AGL_data.AGL_edos_properties.at(i).Fermi_energy) {
          below_EF = j;
          above_EF = j + 1;
        } else if (AGL_data.AGL_edos_properties.at(i).energy.at(j) > AGL_data.AGL_edos_properties.at(i).Fermi_energy) {
          break;
        }
      }
      bottom_EF = 0;
      // Determine the energy region around the Fermi level where the capacity of the DOS changes by less than 1 electron
      // If the DOS is less than the tolerance at any point within this region, then there is a band gap
      // Otherwise, the material is a metal, since there are sufficient states to fit an electron or hole at the Fermi level
      // This procedure is necessary since VASP doesn't always put the Fermi level inside the gap (usually it is just before the start of the gap)
      for (uint j = 0; j < above_EF; j++) {
        if ((AGL_data.AGL_edos_properties.at(i).idosval.at(below_EF) - AGL_data.AGL_edos_properties.at(i).idosval.at(j)) < 1.0) {
          bottom_EF = j;
          break;
        } else if (j == below_EF) {
          bottom_EF = below_EF;
        }
      }
      top_EF = above_EF;
      for (uint j = above_EF; j < AGL_data.AGL_edos_properties.at(i).energy.size(); j++) {
        if ((AGL_data.AGL_edos_properties.at(i).idosval.at(j) - AGL_data.AGL_edos_properties.at(i).idosval.at(above_EF)) > 1.0) {
          top_EF = j;
          break;
        } else if (j == (AGL_data.AGL_edos_properties.at(i).energy.size() - 1)) {
          top_EF = j;
        }
      }
      dosminEF = AGL_data.AGL_edos_properties.at(i).dosval.at(bottom_EF);
      min_EF = bottom_EF;
      // Find point minimum point of DOS in region around Fermi level
      for (uint j = bottom_EF; j <= top_EF; j++) {
        if (AGL_data.AGL_edos_properties.at(i).dosval.at(j) < dosminEF) {
          dosminEF = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
          min_EF = j;
        }
      }
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Minimum in DOS is at energy = " << AGL_data.AGL_edos_properties.at(i).energy.at(min_EF) << endl;
      aus << _AGLSTR_MESSAGE_ + "Lower edge band gap region is at energy = " << AGL_data.AGL_edos_properties.at(i).energy.at(bottom_EF) << endl;
      aus << _AGLSTR_MESSAGE_ + "Upper edge band gap region is at energy = " << AGL_data.AGL_edos_properties.at(i).energy.at(top_EF) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 	
      // Check if the minimum point is greater than or equal to the DOS tolerance for a band gap
      if (AGL_data.AGL_edos_properties.at(i).dosval.at(min_EF) >= egap_tol) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Minimum in DOS is = " << AGL_data.AGL_edos_properties.at(i).dosval.at(min_EF) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 	
        // If it is one of the endpoints of the region around the Fermi level, then the structure is metallic
        if ((min_EF == bottom_EF) || (min_EF == top_EF)) {
          AGL_data.AGL_edos_properties.at(i).edos_band_gap = 0.0;
          AGL_data.AGL_edos_properties.at(i).edos_gap_set = true;
          // Interpolate DOS to get value at Fermi energy
          // Find region around Fermi energy where DOS is monotonic
          if (AGL_data.AGL_edos_properties.at(i).dosval.at(below_EF) < AGL_data.AGL_edos_properties.at(i).dosval.at(above_EF)) {
            dvallower = AGL_data.AGL_edos_properties.at(i).dosval.at(below_EF);
            vallower = below_EF;
            uint j = below_EF - 1;
            while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) < dvallower) && (j > 0)) {
              vallower = j;
              dvallower = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
              j--;
            }
            dvalupper = AGL_data.AGL_edos_properties.at(i).dosval.at(above_EF);
            valupper = above_EF;
            j = above_EF - 1;
            while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) > dvalupper) && (j > 0)) {
              valupper = j;
              dvalupper = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
              j--;
            }
          } else {
            dvallower = AGL_data.AGL_edos_properties.at(i).dosval.at(below_EF);
            vallower = below_EF;
            uint j = below_EF - 1;
            while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) > dvallower) && (j > 0)) {
              vallower = j;
              dvallower = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
              j--;
            }
            dvalupper = AGL_data.AGL_edos_properties.at(i).dosval.at(above_EF);
            valupper = above_EF;
            j = above_EF - 1;
            while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) < dvalupper) && (j > 0)) {
              valupper = j;
              dvalupper = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
              j--;
            }
          }
          energtofit.clear();
          dostofit.clear();
          // Range of values to fit
          for (uint j = vallower; j <= valupper; j++) {
            // Only fit values within two DOS points or within factor 10 of points on each side of Fermi level
            // This prioritizes fitting the bottom of the DOS
            j_below_EF = below_EF - j;
            j_above_EF = above_EF - j;
            if ((aurostd::abs(AGL_data.AGL_edos_properties.at(i).dosval.at(j) / AGL_data.AGL_edos_properties.at(i).dosval.at(below_EF)) < 10.0) || (aurostd::abs(AGL_data.AGL_edos_properties.at(i).dosval.at(j) / AGL_data.AGL_edos_properties.at(i).dosval.at(above_EF)) < 10.0) || (std::abs(j_below_EF) <= 2) || (std::abs(j_above_EF) <= 2)) {
              energtofit.push_back(AGL_data.AGL_edos_properties.at(i).energy.at(j));
              dostofit.push_back(AGL_data.AGL_edos_properties.at(i).dosval.at(j));
            }
          }
          // Fit by polynomial
          nminlimit = 0;
          nmaxlimit = 0;
          if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "Calling edosfit" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          aglerror = AGL_functions::edosfit(energtofit, dostofit, energtoeval, dospolyeval, nminlimit, nmaxlimit, AGL_data.gaussxm_debug, FileMESSAGE);
          if (aglerror != 0) {
            if (aglerror == 2) {
              polyfitfail = true;
            } else {
              return aglerror;
            }
          }
          if (polyfitfail) {
            if (LVERBOSE) {	  
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "edosfit unsuccessful, using DOS values from file" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            // Evaluate DOS at Fermi energy
            energdiffmin = aurostd::abs(energtofit.at(0) - AGL_data.AGL_edos_properties.at(i).Fermi_energy);
            ienergdiffmin = 0;
            for (uint j = 0; j < energtofit.size(); j++) {
              if (aurostd::abs(energtofit.at(j) - AGL_data.AGL_edos_properties.at(i).Fermi_energy) < energdiffmin) {
                energdiffmin = aurostd::abs(energtofit.at(j) - AGL_data.AGL_edos_properties.at(i).Fermi_energy);
                ienergdiffmin = j;
              }
            }
            AGL_data.AGL_edos_properties.at(i).dosval_EF = dostofit.at(ienergdiffmin);
            if (AGL_data.AGL_edos_properties.at(i).dosval_EF < 0.0) {
              AGL_data.AGL_edos_properties.at(i).dosval_EF = 0.0;
            }
            AGL_data.AGL_edos_properties.at(i).dosval_EF_set = true;
            // Record pressure, band gap, and DOS value at Fermi energy
            pressure_list.push_back(AGL_data.AGL_edos_properties.at(i).pressure_external);
            edosgap_list.push_back(AGL_data.AGL_edos_properties.at(i).edos_band_gap);
            dosvalEF_list.push_back(AGL_data.AGL_edos_properties.at(i).dosval_EF);
            continue;	    
          } else {
            if (LVERBOSE) {	  
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "edosfit successful" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            // Evaluate DOS at Fermi energy
            energdiffmin = aurostd::abs(energtoeval.at(0) - AGL_data.AGL_edos_properties.at(i).Fermi_energy);
            ienergdiffmin = 0;
            for (uint j = 0; j < energtoeval.size(); j++) {
              if (aurostd::abs(energtoeval.at(j) - AGL_data.AGL_edos_properties.at(i).Fermi_energy) < energdiffmin) {
                energdiffmin = aurostd::abs(energtoeval.at(j) - AGL_data.AGL_edos_properties.at(i).Fermi_energy);
                ienergdiffmin = j;
              }
            }
            AGL_data.AGL_edos_properties.at(i).dosval_EF = dospolyeval.at(ienergdiffmin);
            if (AGL_data.AGL_edos_properties.at(i).dosval_EF < 0.0) {
              AGL_data.AGL_edos_properties.at(i).dosval_EF = 0.0;
            }
            AGL_data.AGL_edos_properties.at(i).dosval_EF_set = true;
            // Record pressure, band gap, and DOS value at Fermi energy
            pressure_list.push_back(AGL_data.AGL_edos_properties.at(i).pressure_external);
            edosgap_list.push_back(AGL_data.AGL_edos_properties.at(i).edos_band_gap);
            dosvalEF_list.push_back(AGL_data.AGL_edos_properties.at(i).dosval_EF);
            continue;
          }
          // Otherwise, the lowest DOS value near the Fermi energy is a local minimum
          // Take monotonically increasing and decreasing points immediately above and below to define fitting region
        } else {
          dvallower = AGL_data.AGL_edos_properties.at(i).dosval.at(min_EF);
          vallower = min_EF;
          uint j = min_EF - 1;
          while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) > dvallower) && (j > 0)) {
            vallower = j;
            dvallower = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
            j--;
          }
          dvalupper = AGL_data.AGL_edos_properties.at(i).dosval.at(min_EF);
          valupper = min_EF;
          j = min_EF + 1;
          while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) < dvalupper) && (j < AGL_data.AGL_edos_properties.at(i).dosval.size())) {
            valupper = j;
            dvalupper = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
            j++;
          }
          dcondlower = AGL_data.AGL_edos_properties.at(i).dosval.at(min_EF);
          condlower = min_EF;
          j = min_EF - 1;
          while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) < dcondlower) && (j > 0)) {
            condlower = j;
            dcondlower = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
            j--;
          }	  
          dcondupper = AGL_data.AGL_edos_properties.at(i).dosval.at(min_EF);
          condupper = min_EF;
          j = min_EF + 1;
          while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) > dcondupper) && (j < AGL_data.AGL_edos_properties.at(i).dosval.size())) {
            condupper = j;
            dcondupper = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
            j++;
          }
          // [OBSOLETE] valbandedge = vallower;
          // [OBSOLETE] condbandedge = condupper;
          valbandedge = min_EF;
          condbandedge = min_EF;
          valregionedge = vallower;
          condregionedge = condupper;
        }
      } else {
        // Lowest point is less than the DOS tolerance for a band gap
        // Therefore, perform a search starting around the Fermi level for the closest point that goes to zero
        // Check if points around Fermi level are less than DOS band gap tolerance
        if (AGL_data.AGL_edos_properties.at(i).dosval.at(below_EF) < egap_tol) {
          // Search below Fermi level for valence band edge
          valbandedge = below_EF;
          uint j = below_EF - 1;
          while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) < egap_tol) && (j > 0)) {
            valbandedge = j;
            j--;
          }
          condbandedge = below_EF;
          j = above_EF;
          while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) < egap_tol) && (j < AGL_data.AGL_edos_properties.at(i).dosval.size())) {
            condbandedge = j;
            j++;
          }
        } else if ((AGL_data.AGL_edos_properties.at(i).dosval.at(above_EF) < egap_tol)) {
          // Search above Fermi level for valence band edge
          valbandedge = above_EF;
          uint j = below_EF;
          while (AGL_data.AGL_edos_properties.at(i).dosval.at(j) < egap_tol && (j > 0)) {
            valbandedge = j;
            j--;
          }
          condbandedge = above_EF;
          j = above_EF + 1;
          while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) < egap_tol) && (j < AGL_data.AGL_edos_properties.at(i).dosval.size())) {
            condbandedge = j;
            j++;
          }
        } else {
          // Points on either side of Fermi level are non-zero
          // Therefore, perform a search starting around the Fermi level for the closest point that goes to zero
          // First search below Fermi level
          condinit = 0;
          for (uint j = below_EF; j > 0; j--) {
            if (AGL_data.AGL_edos_properties.at(i).dosval.at(j) < egap_tol) {
              condinit = j;
              break;
            } 
          }
          // Search above Fermi level
          valinit = AGL_data.AGL_edos_properties.at(i).dosval.size() - 1;
          for (uint j = above_EF; j < AGL_data.AGL_edos_properties.at(i).dosval.size(); j++) {
            if (AGL_data.AGL_edos_properties.at(i).dosval.at(j) < egap_tol) {
              valinit = j;
              break;
            } 
          }
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valinit = " << valinit << endl;
          aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: condinit = " << condinit << endl;
          aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valinit energy = " << AGL_data.AGL_edos_properties.at(i).energy.at(valinit) << endl;
          aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: condinit energy = " << AGL_data.AGL_edos_properties.at(i).energy.at(condinit) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 	
          // Check which is closer to Fermi level
          if ((aurostd::abs(AGL_data.AGL_edos_properties.at(i).energy.at(below_EF) - AGL_data.AGL_edos_properties.at(i).energy.at(condinit))) < (aurostd::abs(AGL_data.AGL_edos_properties.at(i).energy.at(valinit) - AGL_data.AGL_edos_properties.at(i).energy.at(above_EF)))) {
            // Band gap is below Fermi level
            // condinit is just before the conduction band edge
            condbandedge = condinit;
            valbandedge = 0;
            // Find valence band edge
            for (uint j = condinit; j > 0; j--) {
              if (AGL_data.AGL_edos_properties.at(i).dosval.at(j) > egap_tol) {
                valbandedge = j + 1;
                break;
              }
            }
          } else if ((aurostd::abs(AGL_data.AGL_edos_properties.at(i).energy.at(below_EF) - AGL_data.AGL_edos_properties.at(i).energy.at(condinit))) > (aurostd::abs(AGL_data.AGL_edos_properties.at(i).energy.at(valinit) - AGL_data.AGL_edos_properties.at(i).energy.at(above_EF)))) {
            // Band gap is above Fermi level
            // valinit is just after valence band edge
            valbandedge = valinit;
            condbandedge = AGL_data.AGL_edos_properties.at(i).dosval.size() - 1;
            for (uint j = valinit; j < AGL_data.AGL_edos_properties.at(i).dosval.size(); j++) {
              if (AGL_data.AGL_edos_properties.at(i).dosval.at(j) > egap_tol) {
                condbandedge = j - 1;
                break;
              }
            }
          } else {
            // Points are same distance from Fermi level
            // Check whether DOS minimum is above or below Fermi level
            if (AGL_data.AGL_edos_properties.at(i).energy.at(min_EF) < AGL_data.AGL_edos_properties.at(i).Fermi_energy) {
              // Band gap is below Fermi level
              // condinit is just before the conduction band edge
              condbandedge = condinit;
              valbandedge = 0;
              // Find valence band edge
              for (uint j = condinit; j > 0; j--) {
                if (AGL_data.AGL_edos_properties.at(i).dosval.at(j) > egap_tol) {
                  valbandedge = j + 1;
                  break;
                }
              }
            } else {
              // Band gap is above Fermi level
              // valinit is just after valence band edge
              valbandedge = valinit;
              condbandedge = AGL_data.AGL_edos_properties.at(i).dosval.size() - 1;
              for (uint j = valinit; j < AGL_data.AGL_edos_properties.at(i).dosval.size(); j++) {
                if (AGL_data.AGL_edos_properties.at(i).dosval.at(j) > egap_tol) {
                  condbandedge = j - 1;
                  break;
                }
              }
            }
          }
        }
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence edge = " << AGL_data.AGL_edos_properties.at(i).dosval.at(valbandedge) << endl;
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence edge energy = " << AGL_data.AGL_edos_properties.at(i).energy.at(valbandedge) << endl;	
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction edge = " << AGL_data.AGL_edos_properties.at(i).dosval.at(condbandedge) << endl;
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction edge energy = " << AGL_data.AGL_edos_properties.at(i).energy.at(condbandedge) << endl;	
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 	
        // Find monotonic region around valence band edge
        dvallower = AGL_data.AGL_edos_properties.at(i).dosval.at(valbandedge);
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence band region edge = " << dvallower << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 
        vallower = valbandedge;
        uint j = valbandedge - 1;
        while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) > dvallower) && (j > 0)) {
          vallower = j;
          dvallower = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
          j--;
        }
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence band region edge = " << dvallower << endl;
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence band region edge energy = " << AGL_data.AGL_edos_properties.at(i).energy.at(vallower) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 
        dvalupper = AGL_data.AGL_edos_properties.at(i).dosval.at(valbandedge);
        valupper = valbandedge;
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence band region edge = " << dvalupper << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        j = valbandedge + 1;
        while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) < dvalupper) && (j < AGL_data.AGL_edos_properties.at(i).dosval.size())) {
          valupper = j;
          dvalupper = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
          j++;
        }
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence band region edge = " << dvalupper << endl;
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence band region edge energy = " << AGL_data.AGL_edos_properties.at(i).energy.at(valupper) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 	
        // Find monotonic region around conduction band edge
        dcondlower = AGL_data.AGL_edos_properties.at(i).dosval.at(condbandedge);
        condlower = condbandedge;
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction band region edge = " << dcondlower << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        j = condbandedge - 1;
        while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) < dcondlower) && (j > 0)) {
          condlower = j;
          dcondlower = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
          j--;
        }
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction band region edge = " << dcondlower << endl;
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction band region edge energy = " << AGL_data.AGL_edos_properties.at(i).energy.at(condlower) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 	
        dcondupper = AGL_data.AGL_edos_properties.at(i).dosval.at(condbandedge);
        condupper = condbandedge;
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction band region edge = " << dcondupper << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        j = condbandedge + 1;
        while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) > dcondupper) && (j < AGL_data.AGL_edos_properties.at(i).dosval.size())) {
          condupper = j;
          dcondupper = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
          j++;
        }
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction band region edge = " << dcondupper << endl;
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction band region edge energy = " << AGL_data.AGL_edos_properties.at(i).energy.at(condupper) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 	
        // Linear fits to points on either side of transition to zero
        // [OBSOLETE] vallower = valbandedge - 1;
        // [OBSOLETE] valupper = valbandedge;
        // [OBSOLETE] condlower = condbandedge;
        // [OBSOLETE] condupper = condbandedge + 1;
        valregionedge = vallower;
        condregionedge = condupper;
      }
      // Have upper and lower limits for valence and conduction band edges
      // Check if valence and conduction band edge regions touch: if so, then perform a single fit
      if ((condlower - valupper) <= 1) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Performing single fit" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // Find lowest DOS value in fitting region
        valdosmin = AGL_data.AGL_edos_properties.at(i).dosval.at(vallower);
        for (uint j = vallower; j <= condupper; j++) {
          if (AGL_data.AGL_edos_properties.at(i).dosval.at(j) < valdosmin) {
            valdosmin = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
          }
        }
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "valdosmin = " << valdosmin << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // If all DOS points in fitting region are greater than the band gap tolerance, then no restriction on polynomial
        if (valdosmin >= egap_tol) {
          energtofit.clear();
          dostofit.clear();
          // Range of values to fit
          for (uint j = vallower; j <= condupper; j++) {
            // Only fit values within two DOS points or within factor 10 of minimum
            // This prioritizes fitting the bottom of the DOS
            j_min_EF = min_EF - j;
            if ((aurostd::abs(AGL_data.AGL_edos_properties.at(i).dosval.at(j) / AGL_data.AGL_edos_properties.at(i).dosval.at(min_EF)) < 10.0) || (std::abs(j_min_EF) <= 2)) {
              energtofit.push_back(AGL_data.AGL_edos_properties.at(i).energy.at(j));
              dostofit.push_back(AGL_data.AGL_edos_properties.at(i).dosval.at(j));
            }
          }
          if(LVERBOSE) {
            for (uint j = 0; j < energtofit.size(); j++) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Energy to fit = " << energtofit.at(j) << ", dostofit = " << dostofit.at(j) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
          }
          nminlimit = 1;
          nmaxlimit = 0;
          if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "Calling edosfit" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          aglerror = AGL_functions::edosfit(energtofit, dostofit, energtoeval, dospolyeval, nminlimit, nmaxlimit, AGL_data.gaussxm_debug, FileMESSAGE);
          if (aglerror != 0) {
            if (aglerror == 2) {
              polyfitfail = true;
            } else {
              return aglerror;
            }
          }
          if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "edosfit successful" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
        } else {
          // If lowest point in DOS in fitted region is less than gap tolerance, then shrink fitted region until lowest point in fitted polynomial is under the tolerance
          valdosmin = AGL_data.AGL_edos_properties.at(i).dosval.at(vallower);
          for (uint j = vallower; j <= valupper; j++) {
            if (AGL_data.AGL_edos_properties.at(i).dosval.at(j) < valdosmin) {
              valdosmin = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
            }
          }
          if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: vallower = " << vallower << endl;
            aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valupper = " << valupper << endl;
            aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: condlower = " << condlower << endl;
            aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: condupper = " << condupper << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          energtofit.clear();
          dostofit.clear();
          // Range of values to fit
          for (uint j = vallower; j <= condupper; j++) {
            // Only fit values within one DOS point or within factor 10 of minimum
            // This prioritizes fitting the bottom of the DOS
            // [OBSOLETE] j_min_EF = min_EF - j;
            // [OBSOLETE] if ((aurostd::abs(AGL_data.AGL_edos_properties.at(i).dosval.at(j) / AGL_data.AGL_edos_properties.at(i).dosval.at(min_EF)) < 10.0) || (aurostd::abs(j_min_EF) <= 2)) {
            energtofit.push_back(AGL_data.AGL_edos_properties.at(i).energy.at(j));
            dostofit.push_back(AGL_data.AGL_edos_properties.at(i).dosval.at(j));
            // [OBSOLETE] }
          }
          nminlimit = 1;
          nmaxlimit = 0;
          if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "Calling edosfit" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if(LVERBOSE) {
            for (uint j = 0; j < energtofit.size(); j++) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Energy to fit = " << energtofit.at(j) << ", dostofit = " << dostofit.at(j) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
          }
          aglerror = AGL_functions::edosfit(energtofit, dostofit, energtoeval, dospolyeval, nminlimit, nmaxlimit, AGL_data.gaussxm_debug, FileMESSAGE);
          if (aglerror != 0) {
            if (aglerror == 2) {
              polyfitfail = true;
            } else {
              return aglerror;
            }
          } else if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "edosfit successful" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if (!polyfitfail) {
            valpolymin = dospolyeval.at(0);
            for (uint j = 0; j < dospolyeval.size(); j++) {
              if (dospolyeval.at(j) < valpolymin) {
                valpolymin = dospolyeval.at(j);
              }
            }
          }
          // If lowest point in DOS in fitted region is less than gap tolerance and lowest point in fitted polynomial is not,
          // the shrink fitted region until lowest point in fitted polynomial is under the tolerance
          itermax = condupper - vallower + 5;
          iter = 0;
          while ((valpolymin > egap_tol) && ((condupper - vallower) > 1) && (!polyfitfail)) {
            // [OBSOLETE] if (AGL_data.AGL_edos_properties.at(i).dosval.at(vallower) > egap_tol) {
            if (vallower < (AGL_data.AGL_edos_properties.at(i).dosval.size()-1)) {
              if ((AGL_data.AGL_edos_properties.at(i).dosval.at(vallower) > egap_tol) && (AGL_data.AGL_edos_properties.at(i).dosval.at(vallower+1) > egap_tol)) {
                vallower = vallower + 1;
              }
            }
            // [OBSOLETE] } else {
            // [OBSOLETE]   if (AGL_data.AGL_edos_properties.at(i).dosval.at(vallower) > egap_tol) {
            // [OBSOLETE] vallower = vallower + 1;
            // [OBSOLETE]   }	      
            // [OBSOLETE] }
            // [OBSOLETE] if (AGL_data.AGL_edos_properties.at(i).dosval.at(condupper) > egap_tol) {
            if (condupper > 0) {
              if ((AGL_data.AGL_edos_properties.at(i).dosval.at(condupper) > egap_tol) && (AGL_data.AGL_edos_properties.at(i).dosval.at(condupper-1) > egap_tol)) {
                condupper = condupper - 1;
              }
            }
            // [OBSOLETE] } else {
            // [OBSOLETE]   if (AGL_data.AGL_edos_properties.at(i).dosval.at(condupper) > egap_tol) {
            // [OBSOLETE] condupper = condupper - 1;
            // [OBSOLETE]   }	    
            // [OBSOLETE] }
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: vallower = " << vallower << endl;
            aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: condupper = " << condupper << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            energtofit.clear();
            dostofit.clear();
            // Range of values to fit
            for (uint j = vallower; j <= condupper; j++) {
              energtofit.push_back(AGL_data.AGL_edos_properties.at(i).energy.at(j));
              dostofit.push_back(AGL_data.AGL_edos_properties.at(i).dosval.at(j));
            }
            nminlimit = 1;
            nmaxlimit = 0;
            if (LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Calling edosfit" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            aglerror = AGL_functions::edosfit(energtofit, dostofit, energtoeval, dospolyeval, nminlimit, nmaxlimit, AGL_data.gaussxm_debug, FileMESSAGE);
            if (aglerror != 0) {
              if (aglerror == 2) {
                polyfitfail = true;
                aurostd::StringstreamClean(aus);
                aus << _AGLSTR_MESSAGE_ + "edosfit unsuccessful" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              } else {
                return aglerror;
              }
            }
            if (!polyfitfail) {
              if (LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AGLSTR_MESSAGE_ + "edosfit successful" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              for (uint j = 0; j < dostofit.size(); j++) {
                if (dostofit.at(j) < valpolymin) {
                  valpolymin = dostofit.at(j);
                }
              }
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valpolymin = " << valpolymin << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            iter++;
            if ((iter > itermax) && (valpolymin > egap_tol)) {
              polyfitfail = true;
            }
          }
        }
        // Initialize ingap to false for this DOS
        ingap = false;
        if (!polyfitfail) {
          // Test for band gap threshold to find band gap
          for (uint j = 0; j < dospolyeval.size(); j++) {
            if(LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Energy = " << energtoeval.at(j) << ", dospolyeval = " << dospolyeval.at(j) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if (ingap) {
              condbandmin = energtoeval.at(j);
              if (dospolyeval.at(j) > egap_tol) {
                aurostd::StringstreamClean(aus);
                aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction band minimum (polynomial) = " << condbandmin << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                break;
              }
            } else if (dospolyeval.at(j) < egap_tol) {
              valbandmax = energtoeval.at(j);
              ingap = true;
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence band maximum (polynomial) = " << valbandmax << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
          }
        }
        if (ingap && (!polyfitfail)) {
          // Bandgap found, is equal to difference between conduction band minimum and valence band maximum
          AGL_data.AGL_edos_properties.at(i).edos_band_gap = aurostd::abs(condbandmin - valbandmax);
          if (AGL_data.AGL_edos_properties.at(i).edos_band_gap < 0.0) {
            AGL_data.AGL_edos_properties.at(i).edos_band_gap = 0.0;
          }
          AGL_data.AGL_edos_properties.at(i).edos_gap_set = true;
        } else {
          // Bandgap not found in polynomial, check in original DOS data
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Searching DOS for band gap" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          for (uint j = valregionedge; j <= condregionedge; j++) {
            if(LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "j = " << j << endl;
              aus << _AGLSTR_MESSAGE_ + "DOS = " << AGL_data.AGL_edos_properties.at(i).dosval.at(j) << endl;
              aus << _AGLSTR_MESSAGE_ + "Energy = " << AGL_data.AGL_edos_properties.at(i).energy.at(j) << endl;
              aus << _AGLSTR_MESSAGE_ + "ingap = " << ingap << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if (ingap) {
              condbandmin = AGL_data.AGL_edos_properties.at(i).energy.at(j);
              if (AGL_data.AGL_edos_properties.at(i).dosval.at(j) > egap_tol) {
                aurostd::StringstreamClean(aus);
                aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction band minimum (DOS data) = " << condbandmin << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                break;
              }
            } else if (AGL_data.AGL_edos_properties.at(i).dosval.at(j) < egap_tol) {
              valbandmax = AGL_data.AGL_edos_properties.at(i).energy.at(j);
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence band maximum (DOS data) = " << valbandmax << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              ingap = true;
              if (j == condregionedge) {
                condbandmin = AGL_data.AGL_edos_properties.at(i).energy.at(j);
              }
            }
          }
          if (ingap) {
            // Bandgap found, is equal to difference between conduction band minimum and valence band maximum
            AGL_data.AGL_edos_properties.at(i).edos_band_gap = aurostd::abs(condbandmin - valbandmax);
            if (AGL_data.AGL_edos_properties.at(i).edos_band_gap < 0.0) {
              AGL_data.AGL_edos_properties.at(i).edos_band_gap = 0.0;
            }
            AGL_data.AGL_edos_properties.at(i).edos_gap_set = true;
          } else {
            // Bandgap not found, material is metallic at this pressure
            AGL_data.AGL_edos_properties.at(i).edos_band_gap = 0.0;
            AGL_data.AGL_edos_properties.at(i).edos_gap_set = true;
          }
        }
        if (ingap) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence band maximum = " << valbandmax << endl;
          aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction band minimum = " << condbandmin << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Material is metallic at this pressure" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      } else {
        // Valence and conduction band edge regions are separate
        // Fit separate monotonic polynomials to each region
        // Range of values to fit to get valence band maximum
        // Find lowest DOS value in fitting region
        valdosmin = AGL_data.AGL_edos_properties.at(i).dosval.at(vallower);
        for (uint j = vallower; j <= valupper; j++) {
          if (AGL_data.AGL_edos_properties.at(i).dosval.at(j) < valdosmin) {
            valdosmin = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
          }
        }
        if (LVERBOSE) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: vallower = " << vallower << endl;
          aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valupper = " << valupper << endl;
          aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: condlower = " << condlower << endl;
          aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: condupper = " << condupper << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // If all DOS points in fitting region are greater than the band gap tolerance, then no restriction on polynomial
        if (valdosmin >= egap_tol) {
          energtofit.clear();
          dostofit.clear();
          for (uint j = vallower; j <= valupper; j++) {
            energtofit.push_back(AGL_data.AGL_edos_properties.at(i).energy.at(j));
            dostofit.push_back(AGL_data.AGL_edos_properties.at(i).dosval.at(j));
          }
          nminlimit = 0;
          nmaxlimit = 0;
          if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "Calling edosfit" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          aglerror = AGL_functions::edosfit(energtofit, dostofit, energtoeval, dospolyeval, nminlimit, nmaxlimit, AGL_data.gaussxm_debug, FileMESSAGE);
          if (aglerror != 0) {
            if (aglerror == 2) {
              polyfitfail = true;
            } else {
              return aglerror;
            }
          } else if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "edosfit successful" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
        } else {
          if ((valupper - vallower) < 2) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: vallower = " << vallower << endl;
            aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valupper = " << valupper << endl;		    
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          energtofit.clear();
          dostofit.clear();
          // Range of values to fit
          for (uint j = vallower; j <= valupper; j++) {
            energtofit.push_back(AGL_data.AGL_edos_properties.at(i).energy.at(j));
            dostofit.push_back(AGL_data.AGL_edos_properties.at(i).dosval.at(j));
          }
          nminlimit = 0;
          nmaxlimit = 0;
          if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "Calling edosfit" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          aglerror = AGL_functions::edosfit(energtofit, dostofit, energtoeval, dospolyeval, nminlimit, nmaxlimit, AGL_data.gaussxm_debug, FileMESSAGE);
          if (aglerror != 0) {
            if (aglerror == 2) {
              polyfitfail = true;
            } else {
              return aglerror;
            }
          } else if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "edosfit successful" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if (!polyfitfail) {
            valpolymin = dospolyeval.at(0);
            for (uint j = 0; j < dospolyeval.size(); j++) {
              if (dospolyeval.at(j) < valpolymin) {
                valpolymin = dospolyeval.at(j);
              }
            }
          }
          // If lowest point in DOS in fitted region is less than gap tolerance and lowest point in fitted polynomial is not,
          // then shrink fitted region until lowest point in fitted polynomial is under the tolerance
          itermax = valupper - vallower + 5;
          iter = 0;
          while ((valpolymin > egap_tol) && ((valupper - vallower) > 1) && (!polyfitfail)) {
            // [OBSOLETE] if (AGL_data.AGL_edos_properties.at(i).dosval.at(vallower) > egap_tol) {
            if (vallower < (AGL_data.AGL_edos_properties.at(i).dosval.size()-1)) {
              if ((AGL_data.AGL_edos_properties.at(i).dosval.at(vallower) > egap_tol) && (AGL_data.AGL_edos_properties.at(i).dosval.at(vallower+1) > egap_tol)) {
                vallower = vallower + 1;
              }
            }
            // [OBSOLETE]  } else {
            // [OBSOLETE] if (AGL_data.AGL_edos_properties.at(i).dosval.at(vallower) > egap_tol) {
            // [OBSOLETE] vallower = vallower + 1;				
            // [OBSOLETE]   }
            // [OBSOLETE] }
            // [OBSOLETE] if (AGL_data.AGL_edos_properties.at(i).dosval.at(valupper) > egap_tol) {
            if (valupper > 0) {
              if ((AGL_data.AGL_edos_properties.at(i).dosval.at(valupper) > egap_tol) && (AGL_data.AGL_edos_properties.at(i).dosval.at(valupper-1) > egap_tol)) {
                valupper = valupper - 1;
              }
            }
            // [OBSOLETE] } else {
            // [OBSOLETE]   if (AGL_data.AGL_edos_properties.at(i).dosval.at(valupper) > egap_tol) {
            // [OBSOLETE] valupper = valupper - 1;
            // [OBSOLETE]   }
            // [OBSOLETE] }
            energtofit.clear();
            dostofit.clear();
            // Range of values to fit
            for (uint j = vallower; j <= valupper; j++) {
              energtofit.push_back(AGL_data.AGL_edos_properties.at(i).energy.at(j));
              dostofit.push_back(AGL_data.AGL_edos_properties.at(i).dosval.at(j));
            }
            nminlimit = 0;
            nmaxlimit = 0;
            if (LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Calling edosfit" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            aglerror = AGL_functions::edosfit(energtofit, dostofit, energtoeval, dospolyeval, nminlimit, nmaxlimit, AGL_data.gaussxm_debug, FileMESSAGE);
            if (aglerror != 0) {
              if (aglerror == 2) {
                polyfitfail = true;
              } else {
                return aglerror;
              }
            } else if (LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "edosfit successful" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if (!polyfitfail) {
              for (uint j = 0; j < dospolyeval.size(); j++) {
                if (dospolyeval.at(j) < valpolymin) {
                  valpolymin = dospolyeval.at(j);
                }
              }
            }
            iter++;
            if ((iter > itermax) && (valpolymin > egap_tol)) {
              polyfitfail = true;
            }
          }
        }
        if (polyfitfail) {
          // Polynomial fit failed: test original DOS data to find band edge
          for (uint j = 0; j < dostofit.size(); j++) {
            if(LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Energy = " << energtofit.at(j) << ", dostofit = " << dostofit.at(j) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if (dostofit.at(j) < egap_tol) {
              valbandmax = energtofit.at(j);
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence band maximum (DOS) = " << valbandmax << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              break;
            } else if (j == (dostofit.size() - 1)) {
              for (uint k = valregionedge; k <= condregionedge; k++) {
                valbandmax = AGL_data.AGL_edos_properties.at(i).energy.at(k);
                if (AGL_data.AGL_edos_properties.at(i).dosval.at(k) < egap_tol) {
                  aurostd::StringstreamClean(aus);
                  aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence band maximum (DOS data) = " << valbandmax << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                  break;
                } 
              }
            }
          }
        } else {
          // Polynomial fit successful: test for band gap threshold to find valence band maximum
          for (uint j = 0; j < dospolyeval.size(); j++) {
            if(LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Energy = " << energtoeval.at(j) << ", dospolyeval = " << dospolyeval.at(j) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if (dospolyeval.at(j) < egap_tol) {
              valbandmax = energtoeval.at(j);
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence band maximum (polynomial) = " << valbandmax << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              break;
            } else if (j == (dospolyeval.size() - 1)) {
              for (uint k = valregionedge; k <= condregionedge; k++) {
                valbandmax = AGL_data.AGL_edos_properties.at(i).energy.at(k);
                if (AGL_data.AGL_edos_properties.at(i).dosval.at(k) < egap_tol) {
                  aurostd::StringstreamClean(aus);
                  aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence band maximum (DOS data) = " << valbandmax << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                  break;
                } 
              }
            }
          }
        }
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valence band maximum = " << valbandmax << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	  
        // Range of values to fit to get conduction band minimum
        // Find lowest DOS value in fitting region
        condosmin = AGL_data.AGL_edos_properties.at(i).dosval.at(condlower);
        for (uint j = condlower; j <= condupper; j++) {
          if (AGL_data.AGL_edos_properties.at(i).dosval.at(j) < condosmin) {
            condosmin = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
          }
        }
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "condosmin= " << condosmin << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	  
        // If all DOS points in fitting region are greater than the band gap tolerance, then no restriction on polynomial
        if (condosmin >= egap_tol) {
          energtofit.clear();
          dostofit.clear();
          for (uint j = condlower; j <= condupper; j++) {
            energtofit.push_back(AGL_data.AGL_edos_properties.at(i).energy.at(j));
            dostofit.push_back(AGL_data.AGL_edos_properties.at(i).dosval.at(j));
          }
          nminlimit = 0;
          nmaxlimit = 0;
          if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "Calling edosfit" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          aglerror = AGL_functions::edosfit(energtofit, dostofit, energtoeval, dospolyeval, nminlimit, nmaxlimit, AGL_data.gaussxm_debug, FileMESSAGE);
          if (aglerror != 0) {
            if (aglerror == 2) {
              polyfitfail = true;
            } else {
              return aglerror;
            }
          } else if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "edosfit successful" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
        } else {
          energtofit.clear();
          dostofit.clear();
          // Range of values to fit
          for (uint j = condlower; j <= condupper; j++) {
            energtofit.push_back(AGL_data.AGL_edos_properties.at(i).energy.at(j));
            dostofit.push_back(AGL_data.AGL_edos_properties.at(i).dosval.at(j));
          }
          nminlimit = 0;
          nmaxlimit = 0;
          if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "Calling edosfit" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          aglerror = AGL_functions::edosfit(energtofit, dostofit, energtoeval, dospolyeval, nminlimit, nmaxlimit, AGL_data.gaussxm_debug, FileMESSAGE);
          if (aglerror != 0) {
            if (aglerror == 2) {
              polyfitfail = true;
            } else {
              return aglerror;
            }	    
          } else if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "edosfit successful" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if (!polyfitfail) {
            condpolymin = dospolyeval.at(0);
            for (uint j = 0; j < dospolyeval.size(); j++) {
              if (dospolyeval.at(j) < condpolymin) {
                condpolymin = dospolyeval.at(j);
              }
            }
          }
          // If lowest point in DOS in fitted region is less than gap tolerance and lowest point in fitted polynomial is not,
          // then shrink fitted region until lowest point in fitted polynomial is under the tolerance
          itermax = condupper - condlower + 5;
          iter = 0;
          while ((condpolymin > egap_tol) && ((condupper - condlower) > 1) && (!polyfitfail)) {
            // [OBSOLETE] if (AGL_data.AGL_edos_properties.at(i).dosval.at(condlower) > egap_tol) { 
            if (condlower < (AGL_data.AGL_edos_properties.at(i).dosval.size()-1)) {
              if ((AGL_data.AGL_edos_properties.at(i).dosval.at(condlower) > egap_tol) && (AGL_data.AGL_edos_properties.at(i).dosval.at(condlower+1) > egap_tol)) {
                condlower = condlower + 1;
              }
            }
            // [OBSOLETE] } else {
            // [OBSOLETE]   if (AGL_data.AGL_edos_properties.at(i).dosval.at(condlower) > egap_tol) {
            // [OBSOLETE] condlower = condlower + 1;
            // [OBSOLETE]  }
            // [OBSOLETE] }
            // [OBSOLETE] if (AGL_data.AGL_edos_properties.at(i).dosval.at(condupper) > egap_tol) {
            if (condupper > 0) {
              if ((AGL_data.AGL_edos_properties.at(i).dosval.at(condupper) > egap_tol) && (AGL_data.AGL_edos_properties.at(i).dosval.at(condupper-1) > egap_tol)) {
                condupper = condupper - 1;
              }
            }
            // [OBSOLETE] } else {
            // [OBSOLETE]   if (AGL_data.AGL_edos_properties.at(i).dosval.at(condupper) > egap_tol) {
            // [OBSOLETE]   condupper = condupper - 1;
            // [OBSOLETE]   }
            // [OBSOLETE] }	    
            energtofit.clear();
            dostofit.clear();
            // Range of values to fit
            for (uint j = condlower; j <= condupper; j++) {
              energtofit.push_back(AGL_data.AGL_edos_properties.at(i).energy.at(j));
              dostofit.push_back(AGL_data.AGL_edos_properties.at(i).dosval.at(j));
            }
            if(LVERBOSE) {
              for (uint j = 0; j < energtofit.size(); j++) {
                aurostd::StringstreamClean(aus);
                aus << _AGLSTR_MESSAGE_ + "Energy to fit = " << energtofit.at(j) << ", dostofit = " << dostofit.at(j) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
            }
            nminlimit = 0;
            nmaxlimit = 0;
            if (LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Calling edosfit" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            aglerror = AGL_functions::edosfit(energtofit, dostofit, energtoeval, dospolyeval, nminlimit, nmaxlimit, AGL_data.gaussxm_debug, FileMESSAGE);
            if (aglerror != 0) {
              if (aglerror == 2) {
                polyfitfail = true;
              } else {
                return aglerror;
              }
            } else if (LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "edosfit successful" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if (!polyfitfail) {
              for (uint j = 0; j < dospolyeval.size(); j++) {
                if (dospolyeval.at(j) < condpolymin) {
                  condpolymin = dospolyeval.at(j);
                }
              }
            }
            if ((iter > itermax) && (condpolymin > egap_tol)) {
              polyfitfail = true;
            }
          }
        }
        if (polyfitfail) {
          // Polynomial fit failed: test original DOS data to find band edge
          for (uint j = 0; j < dostofit.size(); j++) {
            if(LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Energy = " << energtofit.at(j) << ", dostofit = " << dostofit.at(j) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if (dostofit.at(j) > egap_tol) {
              condbandmin = energtofit.at(j);
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction band minimum (DOS) = " << condbandmin << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	    
              break;
            } else if (j == (dostofit.size() - 1)) {
              for (uint k = valregionedge + 1; k <= condregionedge; k++) {
                condbandmin = AGL_data.AGL_edos_properties.at(i).energy.at(k);
                if (AGL_data.AGL_edos_properties.at(i).dosval.at(k) > egap_tol) {
                  aurostd::StringstreamClean(aus);
                  aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction band minimum (DOS data) = " << condbandmin << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	    
                  break;
                }
              }
            }
          }
        } else {
          // Polynomial fit successful: test for band gap threshold to find conduction band minimum
          for (uint j = 0; j < dospolyeval.size(); j++) {
            if(LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Energy = " << energtoeval.at(j) << ", dospolyeval = " << dospolyeval.at(j) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if (dospolyeval.at(j) > egap_tol) {
              condbandmin = energtoeval.at(j);
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction band minimum (polynomial) = " << condbandmin << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	    
              break;
            } else if (j == (dospolyeval.size() - 1)) {
              for (uint k = valregionedge + 1; k <= condregionedge; k++) {
                condbandmin = AGL_data.AGL_edos_properties.at(i).energy.at(k);
                if (AGL_data.AGL_edos_properties.at(i).dosval.at(k) > egap_tol) {
                  aurostd::StringstreamClean(aus);
                  aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction band minimum (DOS data) = " << condbandmin << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	    
                  break;
                }
              }
            }
          }     
        }
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: conduction band minimum = " << condbandmin << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: vallower = " << vallower << endl;
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: valupper = " << valupper << endl;
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: condlower = " << condlower << endl;
        aus << _AGLSTR_MESSAGE_ + "Electronic DOS gap: condupper = " << condupper << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
        // Calculate band gap
        AGL_data.AGL_edos_properties.at(i).edos_band_gap = aurostd::abs(condbandmin - valbandmax);
        if (AGL_data.AGL_edos_properties.at(i).edos_band_gap < 0.0) {
          AGL_data.AGL_edos_properties.at(i).edos_band_gap = 0.0;
        }
        AGL_data.AGL_edos_properties.at(i).edos_gap_set = true;
      }
      if (AGL_data.AGL_edos_properties.at(i).edos_band_gap > 0.0) {
        // There is a finite band gap, so DOS at Fermi energy should be zero
        AGL_data.AGL_edos_properties.at(i).dosval_EF = 0.0;
        AGL_data.AGL_edos_properties.at(i).dosval_EF_set = true;
      } else {
        // Otherwise, material is metallic at this pressure so DOS should have finite value at Fermi energy
        // Interpolate DOS to get value at Fermi energy
        // Find region around Fermi energy where DOS is monotonic
        if (AGL_data.AGL_edos_properties.at(i).dosval.at(below_EF) < AGL_data.AGL_edos_properties.at(i).dosval.at(above_EF)) {
          dvallower = AGL_data.AGL_edos_properties.at(i).dosval.at(below_EF);
          vallower = below_EF;
          uint j = below_EF - 1;
          while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) < dvallower) && (j > 0)) {
            vallower = j;
            dvallower = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
            j--;
          }
          dvalupper = AGL_data.AGL_edos_properties.at(i).dosval.at(above_EF);
          valupper = above_EF;
          j = above_EF - 1;
          while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) > dvalupper) && (j > 0)) {
            valupper = j;
            dvalupper = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
            j--;
          }
        } else {
          dvallower = AGL_data.AGL_edos_properties.at(i).dosval.at(below_EF);
          vallower = below_EF;
          uint j = below_EF - 1;
          while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) > dvallower) && (j > 0)) {
            vallower = j;
            dvallower = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
            j--;
          }
          dvalupper = AGL_data.AGL_edos_properties.at(i).dosval.at(above_EF);
          valupper = above_EF;
          j = above_EF - 1;
          while ((AGL_data.AGL_edos_properties.at(i).dosval.at(j) < dvalupper) && (j > 0)) {
            valupper = j;
            dvalupper = AGL_data.AGL_edos_properties.at(i).dosval.at(j);
            j--;
          }
        }
        energtofit.clear();
        dostofit.clear();
        // Range of values to fit
        for (uint j = vallower; j <= valupper; j++) {
          // Only fit values within one DOS point or within factor 10 of points on each side of Fermi level
          // This prioritizes fitting the bottom of the DOS
          j_below_EF = below_EF - j;
          j_above_EF = above_EF - j;
          if ((aurostd::abs(AGL_data.AGL_edos_properties.at(i).dosval.at(j) / AGL_data.AGL_edos_properties.at(i).dosval.at(below_EF)) < 10.0) || (aurostd::abs(AGL_data.AGL_edos_properties.at(i).dosval.at(j) / AGL_data.AGL_edos_properties.at(i).dosval.at(above_EF)) < 10.0) || (std::abs(j_below_EF) < 2) || (std::abs(j_above_EF) < 2)) {
            energtofit.push_back(AGL_data.AGL_edos_properties.at(i).energy.at(j));
            dostofit.push_back(AGL_data.AGL_edos_properties.at(i).dosval.at(j));
          }
        }
        if(LVERBOSE) {
          for (uint j = 0; j < energtofit.size(); j++) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "Energy to fit = " << energtofit.at(j) << ", dostofit = " << dostofit.at(j) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
        }
        // Fit by polynomial
        nminlimit = 1;
        nmaxlimit = 0;
        if (LVERBOSE) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Calling edosfit" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        aglerror = AGL_functions::edosfit(energtofit, dostofit, energtoeval, dospolyeval, nminlimit, nmaxlimit, AGL_data.gaussxm_debug, FileMESSAGE);
        if (aglerror != 0) {
          if (aglerror == 2) {
            polyfitfail = true;
          } else {
            return aglerror;
          }
        } else {
          if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "edosfit successful" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          // Evaluate DOS at Fermi energy
          energdiffmin = aurostd::abs(energtoeval.at(0) - AGL_data.AGL_edos_properties.at(i).Fermi_energy);
          ienergdiffmin = 0;
          for (uint j = 0; j < energtoeval.size(); j++) {
            if (aurostd::abs(energtoeval.at(j) - AGL_data.AGL_edos_properties.at(i).Fermi_energy) < energdiffmin) {
              energdiffmin = aurostd::abs(energtoeval.at(j) - AGL_data.AGL_edos_properties.at(i).Fermi_energy);
              ienergdiffmin = j;
            }
          }
          AGL_data.AGL_edos_properties.at(i).dosval_EF = dospolyeval.at(ienergdiffmin);
          if (AGL_data.AGL_edos_properties.at(i).dosval_EF < 0.0) {
            AGL_data.AGL_edos_properties.at(i).dosval_EF = 0.0;
          }
          AGL_data.AGL_edos_properties.at(i).dosval_EF_set = true;
        }
      }
      if (polyfitfail) {
        // Polynomial fit failed: find DOS value at Fermi energy in calculated DOS
        energdiffmin = aurostd::abs(AGL_data.AGL_edos_properties.at(i).energy.at(0) - AGL_data.AGL_edos_properties.at(i).Fermi_energy);
        ienergdiffmin = 0;
        for (uint j = 0; j < AGL_data.AGL_edos_properties.at(i).energy.size(); j++) {
          if (aurostd::abs(AGL_data.AGL_edos_properties.at(i).energy.at(j) - AGL_data.AGL_edos_properties.at(i).Fermi_energy) < energdiffmin) {
            energdiffmin = aurostd::abs(AGL_data.AGL_edos_properties.at(i).energy.at(j) - AGL_data.AGL_edos_properties.at(i).Fermi_energy);
            ienergdiffmin = j;
          }
        }
        AGL_data.AGL_edos_properties.at(i).dosval_EF = AGL_data.AGL_edos_properties.at(i).dosval.at(ienergdiffmin);
        if (AGL_data.AGL_edos_properties.at(i).dosval_EF < 0.0) {
          AGL_data.AGL_edos_properties.at(i).dosval_EF = 0.0;
        }
        AGL_data.AGL_edos_properties.at(i).dosval_EF_set = true;	
      }
      pressure_list.push_back(AGL_data.AGL_edos_properties.at(i).pressure_external);
      edosgap_list.push_back(AGL_data.AGL_edos_properties.at(i).edos_band_gap);
      dosvalEF_list.push_back(AGL_data.AGL_edos_properties.at(i).dosval_EF);
    }
    // Get derivatives of band gap and DOS(EF) as a function of pressure at zero pressure
    // Find also minimum and maximum values of band gap and corresponding pressures
    // First, find electronic band gap at pressure extrema and at value closest to zero
    double presabsmin = aurostd::abs(pressure_list.at(0));
    double presmax = pressure_list.at(0);
    double presmin = pressure_list.at(0);
    uint ipresabsmin = 0;
    uint ipresmax = 0;
    uint ipresmin = 0;
    double egapmin = edosgap_list.at(0);
    double egapmax = edosgap_list.at(0);
    double egapminpres = pressure_list.at(0);
    double egapmaxpres = pressure_list.at(0);
    double edosmin = dosvalEF_list.at(0);
    double edosmax = dosvalEF_list.at(0);
    double edosminpres = pressure_list.at(0);
    double edosmaxpres = pressure_list.at(0);
    for (uint i = 0; i < pressure_list.size(); i++) {
      if (pressure_list.at(i) > presmax) {
        presmax = pressure_list.at(i);
        ipresmax = i;
      }
      if (pressure_list.at(i) < presmin) {
        presmin = pressure_list.at(i);
        ipresmin = i;
      }
      if (aurostd::abs(pressure_list.at(i)) < presabsmin) {
        presabsmin = aurostd::abs(pressure_list.at(i));
        ipresabsmin = i;
      }
      if (edosgap_list.at(i) < egapmin) {
        egapmin = edosgap_list.at(i);
        egapminpres = pressure_list.at(i);
      }
      if (edosgap_list.at(i) > egapmax) {
        egapmax = edosgap_list.at(i);
        egapmaxpres = pressure_list.at(i);
      }
      if (dosvalEF_list.at(i) < edosmin) {
        edosmin = dosvalEF_list.at(i);
        edosminpres = pressure_list.at(i);
      }
      if (dosvalEF_list.at(i) > edosmax) {
        edosmax = dosvalEF_list.at(i);
        edosmaxpres = pressure_list.at(i);
      }
    }
    AGL_data.egap_min = egapmin;
    AGL_data.egap_max = egapmax;
    AGL_data.egap_min_pressure = egapminpres;
    AGL_data.egap_max_pressure = egapmaxpres;    
    AGL_data.dosvalEF_min = edosmin;
    AGL_data.dosvalEF_max = edosmax;
    AGL_data.dosvalEF_min_pressure = edosminpres;
    AGL_data.dosvalEF_max_pressure = edosmaxpres;    
    double presabsmingap = edosgap_list.at(ipresabsmin);
    double presmingap = edosgap_list.at(ipresmin);
    double presmaxgap = edosgap_list.at(ipresmax);
    double egapmaxdiff = 0.0;
    bool egapmaxdiffpoint = false;
    if (LVERBOSE) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Pressure minimum = " << pressure_list.at(ipresmin) << endl;
      aus << _AGLSTR_MESSAGE_ + "Pressure maximum = " << pressure_list.at(ipresmax) << endl;
      aus << _AGLSTR_MESSAGE_ + "Pressure absolute minimum = " << pressure_list.at(ipresabsmin) << endl;
      aus << _AGLSTR_MESSAGE_ + "Electronic band gap at pressure minimum = " << edosgap_list.at(ipresmin) << endl;
      aus << _AGLSTR_MESSAGE_ + "Electronic band gap at pressure maximum = " << edosgap_list.at(ipresmax) << endl;
      aus << _AGLSTR_MESSAGE_ + "Electronic band gap at absolute pressure minimum = " << edosgap_list.at(ipresabsmin) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
    }
    // Check if gap at zero pressure is greater or less than gap at both pressure extrema, or if it is in between
    // This will determine maximum number of permitted turning points
    if ((presabsmingap > presmingap) && (presabsmingap > presmaxgap)) {
      // Local maximum at zero pressure
      nminlimit = 0;
      nmaxlimit = 1;
    } else if ((presabsmingap < presmingap) && (presabsmingap < presmaxgap)) {
      // Local minimum at zero pressure
      nminlimit = 1;
      nmaxlimit = 0;
    } else {
      // Gap monotonically increases or decreases with pressure
      nminlimit = 0;
      nmaxlimit = 0;   
    }
    // Check if there is a band gap point above or below both endpoints that is nearer to zero pressure than the endpoint pressures
    for (uint i = 0; i < pressure_list.size(); i++) {
      egapmaxdiffpoint = false;
      if (aurostd::abs(edosgap_list.at(i) - edosgap_list.at(ipresmin)) > egapmaxdiff) {
        egapmaxdiff = aurostd::abs(edosgap_list.at(i) - edosgap_list.at(ipresmin));
        egapmaxdiffpoint = true;
      } 
      if (aurostd::abs(edosgap_list.at(i) - edosgap_list.at(ipresmax)) > egapmaxdiff) {
        egapmaxdiff = aurostd::abs(edosgap_list.at(i) - edosgap_list.at(ipresmax));
        egapmaxdiffpoint = true;
      }
      if (LVERBOSE) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Pressure = " << pressure_list.at(i) << endl;
        aus << _AGLSTR_MESSAGE_ + "Electronic band gap = " << edosgap_list.at(i) << endl;
        aus << _AGLSTR_MESSAGE_ + "Pressure min distance = " << aurostd::abs(pressure_list.at(i) - pressure_list.at(ipresmin)) << endl;
        aus << _AGLSTR_MESSAGE_ + "Pressure max distance = " << aurostd::abs(pressure_list.at(i) - pressure_list.at(ipresmax)) << endl;
        aus << _AGLSTR_MESSAGE_ + "Pressure absolute min distance = " << aurostd::abs(pressure_list.at(i) - pressure_list.at(ipresmin)) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }	
      if (pressure_list.at(i) < pressure_list.at(ipresabsmin)) {
        if (aurostd::abs(pressure_list.at(i) - pressure_list.at(ipresabsmin)) < aurostd::abs(pressure_list.at(i) - pressure_list.at(ipresmin))) {
          if ((edosgap_list.at(i) > edosgap_list.at(ipresmin)) && (edosgap_list.at(i) > edosgap_list.at(ipresmax)) && egapmaxdiffpoint) {
            nminlimit = 0;
            nmaxlimit = 1;
          } else if ((edosgap_list.at(i) < edosgap_list.at(ipresmin)) && (edosgap_list.at(i) < edosgap_list.at(ipresmax)) && egapmaxdiffpoint) {
            nminlimit = 1;
            nmaxlimit = 0;
          }
        }
      } else if (pressure_list.at(i) > pressure_list.at(ipresabsmin)) {
        if (aurostd::abs(pressure_list.at(i) - pressure_list.at(ipresabsmin)) < aurostd::abs(pressure_list.at(i) - pressure_list.at(ipresmax))) {
          if ((edosgap_list.at(i) > edosgap_list.at(ipresmin)) && (edosgap_list.at(i) > edosgap_list.at(ipresmax)) && egapmaxdiffpoint) {
            nminlimit = 0;
            nmaxlimit = 1;
          } else if ((edosgap_list.at(i) < edosgap_list.at(ipresmin)) && (edosgap_list.at(i) < edosgap_list.at(ipresmax))  && egapmaxdiffpoint) {
            nminlimit = 1;
            nmaxlimit = 0;
          }
        }
      }
    }
    if(LVERBOSE) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Number of minima limit = " << nminlimit << endl;
      aus << _AGLSTR_MESSAGE_ + "Number of maxima limit = " << nmaxlimit << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
    }
    // Get derivative of band gap as a function of pressure at zero pressure
    aglerror = AGL_functions::edosgap_pressure_fit(pressure_list, edosgap_list, AGL_data.edosgap_pressure, edosgapfit_list, nminlimit, nmaxlimit, AGL_data.gaussxm_debug, FileMESSAGE);
    if (aglerror != 0) {
      return aglerror;
    }
    // Store band gap values evaluated from polynomial
    double presmindiff;
    double prestol = 0.1;
    uint ipresval;
    for (uint j = 0; j < AGL_data.AGL_edos_properties.size(); j++) {
      // First find matching pressure
      presmindiff = aurostd::abs(AGL_data.AGL_edos_properties.at(j).pressure_external - pressure_list.at(0));
      ipresval = 0;
      for (uint k = 0; k < pressure_list.size(); k++) {
        if (presmindiff > aurostd::abs(AGL_data.AGL_edos_properties.at(j).pressure_external - pressure_list.at(k))) {
          presmindiff = aurostd::abs(AGL_data.AGL_edos_properties.at(j).pressure_external - pressure_list.at(k));
          ipresval = k;
        }
      }
      // If the difference is less than the pressure tolerance, then save this value of the band gap
      if (presmindiff < prestol) {
        AGL_data.AGL_edos_properties.at(j).edos_band_gap_polyval = edosgapfit_list.at(ipresval);
        AGL_data.AGL_edos_properties.at(j).edos_gap_poly_set = true;
      }
    }
    if (LVERBOSE) {
      string ofilepgapfitname = AGL_data.dirpathname + "/AGL_pressure_egap_fit.out";
      stringstream ofilepgapfitss;
      for (uint j = 0; j < pressure_list.size(); j++) {
        ofilepgapfitss << pressure_list.at(j) << "\t" << edosgapfit_list.at(j) << endl;
      }
      if(!aurostd::stringstream2file(ofilepgapfitss, ofilepgapfitname, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Unable to open file " << ofilepgapfitname.c_str() <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }
    }
    // Get derivative of DOS(EF) as a function of pressure at zero pressure
    // First find number of turning points
    double presabsmindos = dosvalEF_list.at(ipresabsmin);
    double presmindos = dosvalEF_list.at(ipresmin);
    double presmaxdos = dosvalEF_list.at(ipresmax);
    double edosmaxdiff = 0.0;
    bool edosmaxdiffpoint = false;
    // Check if DOS at zero pressure is greater or less than DOS at both pressure extrema, or if it is in between
    // This will determine maximum number of permitted turning points
    if ((presabsmindos > presmindos) && (presabsmindos > presmaxdos)) {
      // Local maximum at zero pressure
      nminlimit = 0;
      nmaxlimit = 1;
    } else if ((presabsmindos < presmindos) && (presabsmindos < presmaxdos)) {
      // Local minimum at zero pressure
      nminlimit = 1;
      nmaxlimit = 0;
    } else {
      // DOS monotonically increases or decreases with pressure
      nminlimit = 0;
      nmaxlimit = 0;   
    }
    // Check if there is a DOS value point above or below both endpoints that is nearer to zero pressure than the endpoint pressures
    for (uint i = 0; i < pressure_list.size(); i++) {
      edosmaxdiffpoint = false;
      if (aurostd::abs(dosvalEF_list.at(i) - dosvalEF_list.at(ipresmin)) > edosmaxdiff) {
        edosmaxdiff = aurostd::abs(dosvalEF_list.at(i) - dosvalEF_list.at(ipresmin));
        edosmaxdiffpoint = true;
      } 
      if (aurostd::abs(dosvalEF_list.at(i) - dosvalEF_list.at(ipresmax)) > edosmaxdiff) {
        edosmaxdiff = aurostd::abs(dosvalEF_list.at(i) - dosvalEF_list.at(ipresmax));
        edosmaxdiffpoint = true;
      }
      if (pressure_list.at(i) < pressure_list.at(ipresabsmin)) {
        if (aurostd::abs(pressure_list.at(i) - pressure_list.at(ipresabsmin)) < aurostd::abs(pressure_list.at(i) - pressure_list.at(ipresmin))) {
          if ((dosvalEF_list.at(i) > dosvalEF_list.at(ipresmin)) && (dosvalEF_list.at(i) > dosvalEF_list.at(ipresmax)) && edosmaxdiffpoint) {
            nminlimit = 0;
            nmaxlimit = 1;
          } else if ((dosvalEF_list.at(i) < dosvalEF_list.at(ipresmin)) && (dosvalEF_list.at(i) < dosvalEF_list.at(ipresmax)) && edosmaxdiffpoint) {
            nminlimit = 1;
            nmaxlimit = 0;
          }
        }
      } else if (pressure_list.at(i) > pressure_list.at(ipresabsmin)) {
        if (aurostd::abs(pressure_list.at(i) - pressure_list.at(ipresabsmin)) < aurostd::abs(pressure_list.at(i) - pressure_list.at(ipresmax))) {
          if ((dosvalEF_list.at(i) > dosvalEF_list.at(ipresmin)) && (dosvalEF_list.at(i) > dosvalEF_list.at(ipresmax)) && edosmaxdiffpoint) {
            nminlimit = 0;
            nmaxlimit = 1;
          } else if ((dosvalEF_list.at(i) < dosvalEF_list.at(ipresmin)) && (dosvalEF_list.at(i) < dosvalEF_list.at(ipresmax)) && edosmaxdiffpoint) {
            nminlimit = 1;
            nmaxlimit = 0;
          }
        }
      }
    }	

    aglerror = AGL_functions::edosgap_pressure_fit(pressure_list, dosvalEF_list, AGL_data.dosvalEF_pressure, dosvalEFfit_list, nminlimit, nmaxlimit, AGL_data.gaussxm_debug, FileMESSAGE);
    if (aglerror != 0) {
      return aglerror;
    }
    for (uint j = 0; j < AGL_data.AGL_edos_properties.size(); j++) {
      // First find matching pressure
      presmindiff = aurostd::abs(AGL_data.AGL_edos_properties.at(j).pressure_external - pressure_list.at(0));
      ipresval = 0;
      if (LVERBOSE) {
        aurostd::StringstreamClean(aus);
        aus << "Pressure point j = " << j << endl;    
        aus << "Pressure value j = " << AGL_data.AGL_edos_properties.at(j).pressure_external << endl;    
        aus << "Pressure value ipresval = " << pressure_list.at(0) << endl;    
        aus << "presmindiff = " << presmindiff << endl;      
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      for (uint k = 0; k < pressure_list.size(); k++) {
        if (presmindiff > aurostd::abs(AGL_data.AGL_edos_properties.at(j).pressure_external - pressure_list.at(k))) {
          presmindiff = aurostd::abs(AGL_data.AGL_edos_properties.at(j).pressure_external - pressure_list.at(k));
          ipresval = k;
          if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << "Pressure point ipresval = " << ipresval << endl;    
            aus << "Pressure value ipresval = " << pressure_list.at(ipresval) << endl;    
            aus << "presmindiff = " << presmindiff << endl;      
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
        }
      }
      // If the difference is less than the pressure tolerance, then save this value of the band gap
      if (presmindiff < prestol) {
        AGL_data.AGL_edos_properties.at(j).dosval_EF_polyval = dosvalEFfit_list.at(ipresval);
        AGL_data.AGL_edos_properties.at(j).dosval_EF_poly_set = true;      
      }
      if (LVERBOSE) {
        aurostd::StringstreamClean(aus);
        aus << "Pressure point j = " << j << endl;    
        aus << "Pressure value j = " << AGL_data.AGL_edos_properties.at(j).pressure_external << endl;    
        aus << "Pressure point ipresval = " << j << endl;    
        aus << "Pressure value ipresval = " << pressure_list.at(ipresval) << endl;    
        aus << "DOS value j = " << AGL_data.AGL_edos_properties.at(j).dosval_EF_polyval << endl;      
        aus << "DOS value ipresval = " << dosvalEFfit_list.at(ipresval) << endl;      
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }
    if (LVERBOSE) {
      string ofilepdosfitname = AGL_data.dirpathname + "/AGL_pressure_edos_fit.out";
      stringstream ofilepdosfitss;
      aurostd::StringstreamClean(aus);
      for (uint j = 0; j < pressure_list.size(); j++) {
        ofilepdosfitss << pressure_list.at(j) << "\t" << dosvalEFfit_list.at(j) << endl;
        aus << pressure_list.at(j) << "\t" << dosvalEFfit_list.at(j) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	  
      if(!aurostd::stringstream2file(ofilepdosfitss, ofilepdosfitname, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Unable to open file " << ofilepdosfitname.c_str() <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }
    }
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::edosfit
// *************************************************************************** 
namespace AGL_functions {
  //
  // edosfit: Interpolate electronic density of states with polynomial function
  //
  uint edosfit(vector<double>& energtofit, vector<double>& dostofit, vector<double>& energtoeval, vector<double>& dospolyeval, uint& nminlimit, uint& nmaxlimit, bool gxmdebug, ofstream& FileMESSAGE) {
    bool LVERBOSE=(FALSE || XHOST.DEBUG);
    uint pferr = 0;
    // [OBSOLETE] uint nmaxlimit = 0;
    uint nmin, nmax;
    ostringstream aus;
    uint first_entry_tofit = 0;
    uint last_entry_tofit = dostofit.size() - 1;
    double rms = 0.0;
    double energval, dosval, energrange, denergpoints;
    int npolycoeffmax;
    uint nenergpoints;
    // [OBSOLETE] bool ingap = false;
    vector<double> polycoeffwork, weight;
    if (dostofit.size() > 1) {
      npolycoeffmax = dostofit.size() - 1;
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Not enough points to fit" << endl;
      for (uint i = 0; i < dostofit.size(); i++) {
        aus << energtofit.at(i) << "\t" << dostofit.at(i) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 2;
    } 
    polycoeffwork.resize(npolycoeffmax);
    for (uint i = 0; i < polycoeffwork.size(); i++) {
      polycoeffwork.at(i) = 0.0;
    }
    if(LVERBOSE) {
      aurostd::StringstreamClean(aus);
      for (uint i = 0; i < dostofit.size(); i++) {
        aus << energtofit.at(i) << "\t" << dostofit.at(i) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 
    // Initialize the weights
    weight.resize(dostofit.size());
    for (uint i = 0; i < weight.size(); i++) {
      weight.at(i) = 1.0;
    }
    // Calculate energy range
    energrange = energtofit.at(last_entry_tofit) - energtofit.at(first_entry_tofit);
    if (energrange < 0.0) {
      energrange = energrange * -1.0;
    }
    // Calculate number of 1meV intervals
    denergpoints = energrange / 0.001;
    // Convert to unsigned int
    nenergpoints = static_cast<unsigned int>(denergpoints + 1.5);
    // Calculate energy grid at 1meV intervals
    energval = energtofit.at(0);
    energtoeval.clear();
    dospolyeval.clear();
    for (uint j = 0; j < nenergpoints; j++) {
      energtoeval.push_back(energval);
      energval = energval + 0.001;
    }
    // Fit polynomials to energy-edos data
    // Loop over decreasing polynomial orders until polynomial with appropriate number of extrema in relevant interval is obtaining
    for (int i = npolycoeffmax; i > 0; i--) {
      // Reinitialize polynomial coeffecients to zero
      for (uint j = 0; j < polycoeffwork.size(); j++) {
        polycoeffwork.at(j) = 0.0;
      }
      // Reinitialize number of maxima and minima to zero
      nmin = 0;
      nmax = 0;
      // Fit energy-edos data to obtain polynomial
      uint pferr = AGL_functions::polynom_fit (first_entry_tofit, last_entry_tofit, energtofit, dostofit, weight, rms, i, polycoeffwork, gxmdebug, FileMESSAGE);
      if(pferr != 0) {
        return 2;
      }
      // Evaluate polynomial at 1meV intervals
      dospolyeval.clear();
      for (uint j = 0; j < nenergpoints; j++) {
        uint pferr = AGL_functions::polynom_eval (energtoeval.at(j), polycoeffwork, dosval, 0);
        if(pferr != 0) {
          return 2;
        }
        dospolyeval.push_back(dosval);
      }
      // Count the number of local minima and maxima
      for (uint j = 1; j < (dospolyeval.size() - 1); j++) {
        if ((dospolyeval.at(j-1) < dospolyeval.at(j)) && (dospolyeval.at(j+1) < dospolyeval.at(j))) {
          nmax++;
        }
        if ((dospolyeval.at(j-1) > dospolyeval.at(j)) && (dospolyeval.at(j+1) > dospolyeval.at(j))) {
          nmin++;
        }
      }
      if ((nmin <= nminlimit) && (nmax <= nmaxlimit)) {
        break;
      }
    }
    if(LVERBOSE) {
      for (uint j = 0; j < polycoeffwork.size(); j++) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Polynomial coefficient = " << polycoeffwork.at(j) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }
    return pferr;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::edosgap_pressure_fit
// ***************************************************************************
namespace AGL_functions {
  //
  // edosgap_pressure_fit: Fit pressure versus electronic density of states data with polynomial
  // Derivative of polynomial at zero pressure is a screening descriptor for effect of pressure on electronic properties
  //
  uint edosgap_pressure_fit(vector<double>& pressure, vector<double>& edosgap, double& edosgap_pressure, vector<double>& edosgapfit, uint& nminlimit, uint& nmaxlimit, bool gxmdebug, ofstream& FileMESSAGE) {
    bool LVERBOSE=(FALSE || XHOST.DEBUG);
    uint pferr = 0;
    // [OBSOLETE] uint nmaxlimit = 0;
    uint nmin, nmax;
    ostringstream aus;
    uint first_entry_tofit = 0;
    uint last_entry_tofit = pressure.size() - 1;
    int nminpolyorder;
    double rms = 0.0;
    double zero_pressure = 0.0;
    double edosgapval;
    int npolycoeffmax;
    vector<double> polycoeffwork, weight, edosgappolyeval;
    double egap_tol = 1.0e-6;
    if (edosgap.size() > 1) {
      npolycoeffmax = edosgap.size() - 1;
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Not enough points to fit" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 2;
    } 
    polycoeffwork.resize(npolycoeffmax);
    for (uint i = 0; i < polycoeffwork.size(); i++) {
      polycoeffwork.at(i) = 0.0;
    }
    if(LVERBOSE) {
      aurostd::StringstreamClean(aus);
      for (uint i = 0; i < edosgap.size(); i++) {
        aus << pressure.at(i) << "\t" << edosgap.at(i) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Remove zero regions at beginning and end of pressure range from region to fit
    uint k = 0;
    while ((edosgap.at(k) < egap_tol) && (k < (edosgap.size() - 1))) {
      k++;
      first_entry_tofit = k;
    }
    k = edosgap.size() - 1;
    while ((edosgap.at(k) < egap_tol) && (k > 0)) {
      k--;
      last_entry_tofit = k;
    }
    if(LVERBOSE) {
      aurostd::StringstreamClean(aus);
      aus << "First entry to fit = " << first_entry_tofit << endl;
      aus << "Last entry to fit = " << last_entry_tofit << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Initialize the weights
    weight.resize(edosgap.size());
    for (uint i = 0; i < weight.size(); i++) {
      weight.at(i) = 1.0;
    }
    // Fit polynomials to pressure-edosgap data
    if (first_entry_tofit < last_entry_tofit) {
      nminpolyorder = last_entry_tofit - first_entry_tofit;
      if (nminpolyorder > 2) {
        nminpolyorder = 2;
      }
      // Loop over decreasing polynomial orders until polynomial with appropriate number of extrema in relevant interval is obtained
      // [OBSOLETE] for (int i = npolycoeffmax; i > 0; i--)
      for (int i = npolycoeffmax; i > nminpolyorder; i--) {
        // Reinitialize polynomial coeffecients to zero
        for (uint j = 0; j < polycoeffwork.size(); j++) {
          polycoeffwork.at(j) = 0.0;
        }
        // Reinitialize number of maxima and minima to zero
        nmin = 0;
        nmax = 0;
        // Fit energy-edos data to obtain polynomial
        pferr = AGL_functions::polynom_fit (first_entry_tofit, last_entry_tofit, pressure, edosgap, weight, rms, i, polycoeffwork, gxmdebug, FileMESSAGE);
        if(pferr != 0) {
          return 2;
        }
        // Evaluate polynomial
        edosgappolyeval.clear();
        for (uint j = first_entry_tofit; j <= last_entry_tofit; j++) {
          pferr = AGL_functions::polynom_eval (pressure.at(j), polycoeffwork, edosgapval, 0);
          if(pferr != 0) {
            return 2;
          }
          edosgappolyeval.push_back(edosgapval);
        }
        // Count the number of local minima and maxima
        for (uint j = 1; j < (edosgappolyeval.size() - 1); j++) {
          if ((edosgappolyeval.at(j-1) < edosgappolyeval.at(j)) && (edosgappolyeval.at(j+1) < edosgappolyeval.at(j))) {
            nmax++;
          }
          if ((edosgappolyeval.at(j-1) > edosgappolyeval.at(j)) && (edosgappolyeval.at(j+1) > edosgappolyeval.at(j))) {
            nmin++;
          }
        }
        if ((nmin <= nminlimit) && (nmax <= nmaxlimit)) {
          break;
        }
      }
      edosgapfit.clear();
      for (uint j = 0; j < pressure.size(); j++) {
        if ((j < first_entry_tofit) || (j > last_entry_tofit)) {
          edosgapval = 0.0;
        } else if (edosgap.at(j) < egap_tol) {
          if ((j >= 1) && (j < (edosgap.size() - 1))) {
            if ((edosgap.at(j-1) < egap_tol) && (edosgap.at(j+1) < egap_tol)) {
              edosgapval = 0.0;
            } else {
              pferr = AGL_functions::polynom_eval (pressure.at(j), polycoeffwork, edosgapval, 0);
              if(pferr != 0) {
                return 2;
              }	  
            }
          } else if ((j < 1) && (j < (edosgap.size() - 1))) {
            if ((edosgap.at(j+1) < egap_tol)) {
              edosgapval = 0.0;
            } else {
              pferr = AGL_functions::polynom_eval (pressure.at(j), polycoeffwork, edosgapval, 0);
              if(pferr != 0) {
                return 2;
              }	  
            }
          } else if ((j >= 1) && (j >= (edosgap.size() - 1))) {
            if ((edosgap.at(j-1) < egap_tol)) {
              edosgapval = 0.0;
            } else {
              pferr = AGL_functions::polynom_eval (pressure.at(j), polycoeffwork, edosgapval, 0);
              if(pferr != 0) {
                return 2;
              }	  
            }
          } else {
            edosgapval = 0.0;	
          }
        } else {
          uint pferr = AGL_functions::polynom_eval (pressure.at(j), polycoeffwork, edosgapval, 0);
          if(pferr != 0) {
            return 2;
          }
        }
        if (edosgapval < 0.0) {
          edosgapval = 0.0;
        }
        if (LVERBOSE) {
          aurostd::StringstreamClean(aus);
          aus << "First entry to fit = " << first_entry_tofit << endl;
          aus << "Last entry to fit = " << last_entry_tofit << endl;
          aus << "Pressure point = " << j << endl;
          aus << "Gap/DOS value = " << edosgapval << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        edosgapfit.push_back(edosgapval);
      }
    } else {
      // There are no non-zero points to fit
      edosgapfit.clear();   
      for (uint j = 0; j < pressure.size(); j++) {
        edosgapfit.push_back(0.0);
      }
      for (uint j = 0; j < polycoeffwork.size(); j++) {
        polycoeffwork.at(j) = 0.0;
      }
    }
    if (LVERBOSE) {
      aurostd::StringstreamClean(aus);
      for (uint i = 0; i < edosgapfit.size(); i++) {
        aus << pressure.at(i) << "\t" << edosgapfit.at(i) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Evaluate first derivative at zero pressure
    uint ipresabsmin = 0;
    double presabsmin = pressure.at(0);
    for (uint i = 0; i < pressure.size(); i++) {
      if (aurostd::abs(pressure.at(i)) < presabsmin) {
        presabsmin = aurostd::abs(pressure.at(i));
        ipresabsmin = i;
      }
    }
    if (LVERBOSE) {
      aurostd::StringstreamClean(aus);
      aus << "Zero pressure value = " << presabsmin << endl;
      aus << "Zero pressure point = " << ipresabsmin << endl;      
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    if (edosgap.at(ipresabsmin) < egap_tol) {
      if ((ipresabsmin >= 1) && (ipresabsmin < (edosgap.size() - 1))) {
        if ((edosgap.at(ipresabsmin-1) < egap_tol) && (edosgap.at(ipresabsmin+1) < egap_tol)) {
          edosgap_pressure = 0.0;
        } else {
          pferr = AGL_functions::polynom_eval (zero_pressure, polycoeffwork, edosgap_pressure, 1);
          if(pferr != 0) {
            return 2;
          }	  
        }
      } else if ((ipresabsmin < 1) && (ipresabsmin < (edosgap.size() - 1))) {
        if ((edosgap.at(ipresabsmin+1) < egap_tol)) {
          edosgap_pressure = 0.0;
        } else {
          pferr = AGL_functions::polynom_eval (zero_pressure, polycoeffwork, edosgap_pressure, 1);
          if(pferr != 0) {
            return 2;
          }	  
        }
      } else if ((ipresabsmin >= 1) && (ipresabsmin >= (edosgap.size() - 1))) {
        if ((edosgap.at(ipresabsmin-1) < egap_tol)) {
          edosgap_pressure = 0.0;
        } else {
          pferr = AGL_functions::polynom_eval (zero_pressure, polycoeffwork, edosgap_pressure, 1);
          if(pferr != 0) {
            return 2;
          }	  
        }
      } else {
        edosgap_pressure = 0.0;	
      }
    } else {
      pferr = AGL_functions::polynom_eval (zero_pressure, polycoeffwork, edosgap_pressure, 1);
      if(pferr != 0) {
        return 2;
      }	  
    }
    return pferr;
  }
} // namespace AGL_functions

#endif  // _AFLOW_AGL_ELECTRONIC_CPP
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
