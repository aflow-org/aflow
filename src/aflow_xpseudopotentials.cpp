// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2007-2021
#ifndef _AFLOW_XPSEUDOPOTENTIAL_CPP
#define _AFLOW_XPSEUDOPOTENTIAL_CPP
#include "aflow.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// _XPSEUDOPOTENTIAL_PROTOTYPES_
/*
   ./aflow --scrub=POTCAR --FILE /common/VASP/potpaw_PBE/current/Mo_pv/POTCAR
   ./aflow --scrub=OUTCAR --FILE /common/LIB3/LIB/AgCdCo/TFCC001.ABC/OUTCAR.relax2.xz
   ./aflow --pseudopotentials_check=/tmp/POTCAR1
   ./aflow --pseudopotentials_check=/tmp/POTCAR2
   ./aflow --pseudopotentials_check=/tmp/OUTCAR.relax2
   ./aflow --pseudopotentials_check=/common/LIB3/LIB/AgCdCo/TFCC001.ABC/OUTCAR.relax2.xz

   ./aflow --use_aflow.in=aflow.in --beep --force --showPID --lib2raw="/common/LIB3/LIB/TeW_pvY_sv/TFCC001.ABC"

#!/bin/sh
#echo "$1"
#echo "$2"
STR1=`cat "/common/LIB1/LIB/$1/A1/aflow.in" | grep AUID | head -1 | sed "s/\[VASP_POTCAR_AUID\]/if(AUID==\""/g`
STR2="\") {"$2"}   //   "$1
echo $STR1$STR2

*/

std::vector<xPOTCAR> vxpseudopotential;        // store starting from ONE

#define PSEUDOPOTENTIAL_GENERATOR_pad 70

bool xPOTCAR_FixBoot(xPOTCAR& xPOT) {
  bool fix=FALSE;
  if(!xPOT.vENMAX.size()) { fix=TRUE; xPOT.vENMAX.push_back(NNN);}                                    // if not identified for BOOT STRAP
  if(!xPOT.vENMIN.size()) { fix=TRUE; xPOT.vENMIN.push_back(NNN);}                                    // if not identified for BOOT STRAP
  if(!xPOT.vPOMASS.size()) { fix=TRUE; xPOT.vPOMASS.push_back(NNN);}                                  // if not identified for BOOT STRAP
  if(!xPOT.vZVAL.size()) { fix=TRUE; xPOT.vZVAL.push_back(NNN);}                                      // if not identified for BOOT STRAP
  if(!xPOT.vEATOM.size()) { fix=TRUE; xPOT.vEATOM.push_back(NNN);}                                    // if not identified for BOOT STRAP
  if(!xPOT.vRCORE.size()) { fix=TRUE; xPOT.vRCORE.push_back(NNN);}                                    // if not identified for BOOT STRAP
  if(!xPOT.vRWIGS.size()) { fix=TRUE; xPOT.vRWIGS.push_back(NNN);}                                    // if not identified for BOOT STRAP
  if(!xPOT.vEAUG.size()) { fix=TRUE; xPOT.vEAUG.push_back(NNN);}                                      // if not identified for BOOT STRAP
  if(!xPOT.vRAUG.size()) { fix=TRUE; xPOT.vRAUG.push_back(NNN);}                                      // if not identified for BOOT STRAP
  if(!xPOT.vRMAX.size()) { fix=TRUE; xPOT.vRMAX.push_back(NNN);}                                      // if not identified for BOOT STRAP
  if(!xPOT.vTITEL.size()) { fix=TRUE; xPOT.vTITEL.push_back("N/A");}                                  // if not identified for BOOT STRAP
  if(!xPOT.vLEXCH.size()) { fix=TRUE; xPOT.vLEXCH.push_back("N/A");}                                  // if not identified for BOOT STRAP
  if(!xPOT.species.size()) { fix=TRUE; xPOT.species.push_back("N/A");}                                // if not identified for BOOT STRAP
  if(!xPOT.species_Z.size()) { fix=TRUE; xPOT.species_Z.push_back(0);}                                // if not identified for BOOT STRAP
  if(!xPOT.species_pp.size()) { fix=TRUE; xPOT.species_pp.push_back("N/A");}                          // if not identified for BOOT STRAP
  if(!xPOT.species_pp_type.size()) { fix=TRUE; xPOT.species_pp_type.push_back("N/A");}                // if not identified for BOOT STRAP
  if(!xPOT.species_pp_version.size()) { fix=TRUE; xPOT.species_pp_version.push_back("N/A");}          // if not identified for BOOT STRAP
  if(!xPOT.species_pp_AUID.size()) { fix=TRUE; xPOT.species_pp_AUID.push_back("N/A");}                // if not identified for BOOT STRAP
  if(!xPOT.species_pp_groundstate_energy.size()) { fix=TRUE; xPOT.species_pp_groundstate_energy.push_back(NNN);}          // if not identified for BOOT STRAP
  if(!xPOT.species_pp_groundstate_structure.size()) { fix=TRUE; xPOT.species_pp_groundstate_structure.push_back("N/A");}  // if not identified for BOOT STRAP
  return fix;
}

xPOTCAR xPOTCAR_Finder(vector<string>& species_pp_AUID,vector<string>& species_pp_AUID_collisions,const string& TITEL,const string& LEXCH,const double& EATOM,const double& RMAX,bool LVERBOSE) {
  xPOTCAR xPOT;
  bool found=FALSE;
  for(uint ipp=0;ipp<vxpseudopotential.size();ipp++) {
    bool test=
      (TITEL==vxpseudopotential.at(ipp).vTITEL.at(0)) &&
      (LEXCH==vxpseudopotential.at(ipp).vLEXCH.at(0)) &&
      (aurostd::abs(EATOM-vxpseudopotential.at(ipp).vEATOM.at(0))<0.00001) &&
      (aurostd::abs(RMAX-vxpseudopotential.at(ipp).vRMAX.at(0))<0.0001);
    if(test && found) {
      if(!aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_GGA/potcar.Apr00/Xe") &&
          !aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_GGA/potcar.Apr00/Kr") &&
          !aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_GGA/potcar.Apr00/Li_pv") &&
          !aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_LDA/potcar.Apr00/Li_pv")) {
        //if(LVERBOSE)
        cerr << XPID << "xPOTCAR::xPOTCAR_Finder: COLLISION: POTCAR " << vxpseudopotential.at(ipp).filename << " " << xPOT.filename << endl;
        species_pp_AUID_collisions.push_back(vxpseudopotential.at(ipp).AUID);
      }
    }
    if(test && !found) {
      found=TRUE;
      if(LVERBOSE) cerr << XPID << "xPOTCAR::xPOTCAR_Finder: FOUND: POTCAR=" << vxpseudopotential.at(ipp).filename << endl;
      species_pp_AUID.push_back(vxpseudopotential.at(ipp).AUID);
      xPOT=vxpseudopotential.at(ipp);
    }
  }
  if(!found) {
    if(vxpseudopotential.size()) {
      // if(LVERBOSE)
      cerr << XPID << "xPOTCAR::xPOTCAR_Finder: NOT FOUND: TITEL=" << TITEL << " LEXCH=" <<  LEXCH << " EATOM=" <<  EATOM << " RMAX =" <<  RMAX << endl;
      //    cout << XPID << "xPOTCAR::xPOTCAR_Finder: NOT FOUND: TITEL=" << TITEL << " LEXCH=" <<  LEXCH << " EATOM=" <<  EATOM << " RMAX =" <<  RMAX << endl;
      // throw aurostd::xerror(__AFLOW_FILE__,XPID+"xPOTCAR::xPOTCAR_Finder():","Throw for debugging purposes.",_GENERIC_ERROR_);
    }
    species_pp_AUID.push_back("N/A");
  }

  xPOTCAR_FixBoot(xPOT);
  return xPOT;
}

xPOTCAR xPOTCAR_Finder(const string& AUID,bool LVERBOSE) {
  xPOTCAR xPOT;
  bool found=FALSE;
  for(uint ipp=0;ipp<vxpseudopotential.size();ipp++) {
    bool test= (AUID==vxpseudopotential.at(ipp).AUID);
    if(test && found) {
      if(!aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_GGA/potcar.Apr00/Xe") &&
          !aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_GGA/potcar.Apr00/Kr") &&
          !aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_GGA/potcar.Apr00/Li_pv") &&
          !aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_LDA/potcar.Apr00/Li_pv")) {
        //if(LVERBOSE)
        cerr << XPID << "xPOTCAR::xPOTCAR_Finder: COLLISION: POTCAR " << vxpseudopotential.at(ipp).filename << " " << xPOT.filename << endl;
      }
    }
    if(test && !found) {
      found=TRUE;
      if(LVERBOSE) cerr << XPID << "xPOTCAR::xPOTCAR_Finder: FOUND: POTCAR=" << vxpseudopotential.at(ipp).filename << endl;
      xPOT=vxpseudopotential.at(ipp);
    }
  }
  if(!found) {
    if(vxpseudopotential.size()) {
      // if(LVERBOSE)
      cerr << XPID << "xPOTCAR::xPOTCAR_Finder: NOT FOUND: AUID=" << AUID << endl;
      //   cout << XPID << "xPOTCAR::xPOTCAR_Finder: NOT FOUND: AUID=" << AUID << endl;
      //   throw aurostd::xerror(__AFLOW_FILE__,XPID+"xPOTCAR::xPOTCAR_Finder():","Throw for debugging purposes.",_GENERIC_ERROR_);
    }
  }

  xPOTCAR_FixBoot(xPOT);
  return xPOT;
}

bool xPOTCAR_PURE_Printer(xPOTCAR& xPOT,ostream& oss,bool LVERBOSE) {
  if(XHOST.PSEUDOPOTENTIAL_GENERATOR && xPOT.species.size()==1) {  //SC20200326
    string comment="";
    // objects/functions for references energies  //SC20200326
    //   vector<string> tokens;
    //   aurostd::string2tokens(xPOT.species_pp_version.at(0),tokens,":");
    // [OBSOLETE] if(tokens.size()>0) {vdate.clear();vdate.push_back(tokens.at(tokens.size()-1));};

    double groundstate_energy=NNN;
    string groundstate_structure="N/A_"+xPOT.AUID;

    double volume_atom;
    double spin_atom;

    bool found=xPOTCAR_EnthalpyReference_AUID(xPOT.AUID,"",groundstate_structure,groundstate_energy,volume_atom,spin_atom);

    if(found)  cerr << XPID << "xPOTCAR_PURE_Printer:     FOUND AUID=" << xPOT.AUID << endl;
    if(!found) cerr << XPID << "xPOTCAR_PURE_Printer: NOT_FOUND AUID=" << xPOT.AUID << endl;

    comment=xPOT.species_pp_version.at(0);
    xPOT.vTITEL.at(0)=aurostd::RemoveWhiteSpaces(xPOT.vTITEL.at(0));
    xPOT.vLEXCH.at(0)=aurostd::RemoveWhiteSpaces(xPOT.vLEXCH.at(0));
    if(LVERBOSE) {cerr << XPID << "xPOTCAR_PURE_Printer: vTITEL.at(0)=" << xPOT.vTITEL.at(0) << endl;}   //SC20200326
    if(LVERBOSE) {cerr << XPID << "xPOTCAR_PURE_Printer: species.at(0)=" << xPOT.species.at(0) << endl;}   //SC20200326
    if(LVERBOSE) {cerr << XPID << "xPOTCAR_PURE_Printer: species_Z.at(0)=" << xPOT.species_Z.at(0) << endl;}    //SC20200326
    oss << "  " << endl;
    oss << "  // ******************************************************************************************************************************************************** " << endl;
    //   oss << "  // [AFLOW]START=" << comment << " " << endl;
    oss << "  // " << comment << " " << comment << " " << comment << " " << comment << " " << endl;
    oss << "  // " << xPOT.filename << endl;    //SC20200326
    oss << "  " << aurostd::PaddedPOST("{",PSEUDOPOTENTIAL_GENERATOR_pad) << "      // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("xPOTCAR x;",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.filename=\""+xPOT.filename+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.AUID=\""+xPOT.AUID+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vTITEL.push_back(\""+xPOT.vTITEL.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.pp_type=\""+xPOT.pp_type+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species.push_back(\""+xPOT.species.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_Z.push_back("+aurostd::utype2string<int>(xPOT.species_Z.at(0))+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp.push_back(\""+xPOT.species_pp.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_type.push_back(\""+xPOT.species_pp_type.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_version.push_back(\""+xPOT.species_pp_version.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_AUID.push_back(\""+xPOT.AUID+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    if(!found) oss << "    " << aurostd::PaddedPOST("x.species_pp_groundstate_energy.push_back("+aurostd::utype2string<double>(groundstate_energy,10)+");//"+xPOT.AUID,PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    if(found)  oss << "    " << aurostd::PaddedPOST("x.species_pp_groundstate_energy.push_back("+aurostd::utype2string<double>(groundstate_energy,10)+");//",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_groundstate_structure.push_back(\""+groundstate_structure+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.POTCAR_PAW="+aurostd::bool2string(xPOT.POTCAR_PAW)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.POTCAR_TYPE=\""+xPOT.POTCAR_TYPE+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.POTCAR_KINETIC="+aurostd::bool2string(xPOT.POTCAR_KINETIC)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.POTCAR_GW="+aurostd::bool2string(xPOT.POTCAR_GW)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.POTCAR_AE="+aurostd::bool2string(xPOT.POTCAR_AE)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // EATOM RCORE RWIGS EAUG RAUG ENMAX ENMIN POMASS ZVAL RMAX LEXCH
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vENMAX.push_back("+aurostd::utype2string<double>(xPOT.vENMAX.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vENMIN.push_back("+aurostd::utype2string<double>(xPOT.vENMIN.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vPOMASS.push_back("+aurostd::utype2string<double>(xPOT.vPOMASS.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vZVAL.push_back("+aurostd::utype2string<double>(xPOT.vZVAL.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vEATOM.push_back("+aurostd::utype2string<double>(xPOT.vEATOM.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vRCORE.push_back("+aurostd::utype2string<double>(xPOT.vRCORE.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vRWIGS.push_back("+aurostd::utype2string<double>(xPOT.vRWIGS.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vEAUG.push_back("+aurostd::utype2string<double>(xPOT.vEAUG.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vRAUG.push_back("+aurostd::utype2string<double>(xPOT.vRAUG.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vRMAX.push_back("+aurostd::utype2string<double>(xPOT.vRMAX.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vLEXCH.push_back(\""+xPOT.vLEXCH.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vxpseudopotential.push_back(x);",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "  " << aurostd::PaddedPOST("}",PSEUDOPOTENTIAL_GENERATOR_pad) << "      // " << comment << endl;    //SC20200326
    oss << "  // " << xPOT.filename << endl;    //SC20200326
    //   oss << "  // [AFLOW]STOP=" << comment << " " << endl;
    oss << "  // ******************************************************************************************************************************************************** " << endl;
    oss << endl;    //SC20200326
    return TRUE;
  }
  return FALSE;
}

ostream& operator<<(ostream& oss,const xPOTCAR& xPOT) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  for(uint i=0;i<xPOT.species.size();i++) {
    string comment=xPOT.species_pp_version.at(i);
    oss << "  // ******************************************************************************************************************************************************** " << endl;
    oss << "  // [AFLOW]START=" << comment << " " << endl;
    oss << "  // " << comment << " " << comment << " " << comment << " " << comment << " " << endl;
    oss << "  // " << xPOT.filename << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("filename=\""+xPOT.filename+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("AUID=\""+xPOT.AUID+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vTITEL.push_back(\""+xPOT.vTITEL.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("pp_type=\""+xPOT.pp_type+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species.push_back(\""+xPOT.species.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species_Z.push_back("+aurostd::utype2string<int>(xPOT.species_Z.at(i))+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp.push_back(\""+xPOT.species_pp.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_type.push_back(\""+xPOT.species_pp_type.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_version.push_back(\""+xPOT.species_pp_version.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_AUID.push_back(\""+xPOT.species_pp_AUID.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_groundstate_energy.push_back("+aurostd::utype2string<double>(xPOT.species_pp_groundstate_energy.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_groundstate_structure.push_back(\""+xPOT.species_pp_groundstate_structure.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_PAW="+aurostd::bool2string(xPOT.POTCAR_PAW)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_TYPE=\""+xPOT.POTCAR_TYPE+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_KINETIC="+aurostd::bool2string(xPOT.POTCAR_KINETIC)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_GW="+aurostd::bool2string(xPOT.POTCAR_GW)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_AE="+aurostd::bool2string(xPOT.POTCAR_AE)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vENMAX.push_back("+aurostd::utype2string<double>(xPOT.vENMAX.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vENMIN.push_back("+aurostd::utype2string<double>(xPOT.vENMIN.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vPOMASS.push_back("+aurostd::utype2string<double>(xPOT.vPOMASS.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vZVAL.push_back("+aurostd::utype2string<double>(xPOT.vZVAL.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vEATOM.push_back("+aurostd::utype2string<double>(xPOT.vEATOM.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vRCORE.push_back("+aurostd::utype2string<double>(xPOT.vRCORE.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vRWIGS.push_back("+aurostd::utype2string<double>(xPOT.vRWIGS.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vEAUG.push_back("+aurostd::utype2string<double>(xPOT.vEAUG.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;   //SC20200326
    oss << "    " << aurostd::PaddedPOST("vRAUG.push_back("+aurostd::utype2string<double>(xPOT.vRAUG.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;   //SC20200326
    oss << "    " << aurostd::PaddedPOST("vRMAX.push_back("+aurostd::utype2string<double>(xPOT.vRMAX.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vLEXCH.push_back(\""+xPOT.vLEXCH.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "  // " << xPOT.filename << endl;    //SC20200326
    oss << "  // [AFLOW]STOP=" << comment << " " << endl;
    oss << "  // ******************************************************************************************************************************************************** " << endl;
  }
  return oss;
}

uint xPOTCAR_Initialize(void) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  uint COLLISIONS=0;
  if(LDEBUG) cerr << "xpseudopotentials_initialize [BEGIN]" << endl;
  long double start=aurostd::get_useconds();
  //  if(LDEBUG) cerr << "xpseudopotentials_initialize: " << start << endl;
  if(LDEBUG) cerr << "vxpseudopotential.size()=" << vxpseudopotential.size() << endl;
#include "aflow_xpseudopotentials_data.cpp"
  long double stop=aurostd::get_useconds();
  if(LDEBUG) cerr << "xpseudopotentials_initialize: " << (stop-start)/1000  << endl;
  if(LDEBUG) cerr << "vxpseudopotential.size()=" << vxpseudopotential.size() << endl;
  if(LDEBUG) cerr << "xpseudopotentials_initialize [END]" << endl;

  for(uint i=0;i<vxpseudopotential.size();i++) {
    for(uint j=i+1;j<vxpseudopotential.size();j++) {
      if(vxpseudopotential.at(i).AUID!=vxpseudopotential.at(j).AUID) {
        if(vxpseudopotential.at(i).vTITEL.at(0)==vxpseudopotential.at(j).vTITEL.at(0)) {
          if(vxpseudopotential.at(i).vLEXCH.at(0)==vxpseudopotential.at(j).vLEXCH.at(0)) {
            if(abs(vxpseudopotential.at(i).vEATOM.at(0)-vxpseudopotential.at(j).vEATOM.at(0))<0.00001) {
              if(abs(vxpseudopotential.at(i).vRMAX.at(0)-vxpseudopotential.at(j).vRMAX.at(0))<0.0001) {
                cerr << "COLLISION" << " "
                  << "FILENAME " << vxpseudopotential.at(i).filename << " " << vxpseudopotential.at(j).filename << " "
                  << "TITEL " << vxpseudopotential.at(i).vTITEL.at(0) << " "
                  << "LEXCH " << vxpseudopotential.at(i).vLEXCH.at(0) << " "
                  << "EATOM " << vxpseudopotential.at(i).vEATOM.at(0) << " "
                  << "RMAX " << vxpseudopotential.at(i).vRMAX.at(0) << " "
                  << endl;
                COLLISIONS++;
              }
            }
          }
        }
      }
    }
  }
  //  cerr << "COLLISIONS=" << COLLISIONS << endl;
  return COLLISIONS;
}

bool xPOTCAR_EnthalpyReference_AUID(string AUID,string METAGGA) {
  string groundstate_structure="";
  double groundstate_energy=0.0,volume_atom=0.0,spin_atom=0.0;
  return  xPOTCAR_EnthalpyReference_AUID(AUID,METAGGA,groundstate_structure,groundstate_energy,volume_atom,spin_atom);
}

bool xPOTCAR_EnthalpyReference_AUID(string AUID,string METAGGA,string& groundstate_structure,double& groundstate_energy,double& volume_atom,double& spin_atom) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  bool VERBOSE=0;//TRUE;
  if(LDEBUG) cerr << XPID << "xPOTCAR_EnthalpyReference_AUID: [BEGIN]" << endl;
  bool found=FALSE;
  if(LDEBUG) cout <<"ERROR (xPOTCAR_EnthalpyReference_AUID): AUID=[" << AUID << "]" << endl; // throw aurostd::xerror(__AFLOW_FILE__,XPID+"xPOTCAR_EnthalpyReference_AUID():","Throw for debugging purposes.",_GENERIC_ERROR_);
  if(LDEBUG) cout <<"ERROR (xPOTCAR_EnthalpyReference_AUID): METAGGA=[" << METAGGA << "]" << endl; // throw aurostd::xerror(__AFLOW_FILE__,XPID+"xPOTCAR_EnthalpyReference_AUID():","Throw for debugging purposes.",_GENERIC_ERROR_);

  bool nKIN=FALSE,SCAN=FALSE;
  if(METAGGA.empty() || METAGGA=="none" || METAGGA=="NONE") {nKIN=TRUE;SCAN=FALSE;}
  if(METAGGA=="SCAN" || METAGGA=="scan") {nKIN=FALSE;SCAN=TRUE;}

  if(!XHOST.PSEUDOPOTENTIAL_GENERATOR) if(VERBOSE) cout <<"xPOTCAR_EnthalpyReference_AUID: AUID=[" << AUID << "]  METAGGA=[" << METAGGA << "]" << endl;

  // Ac
  // find /common/LIB1/LIB/Ac* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="06e8eff4a6b826a0" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.031959;volume_atom=44.76122;spin_atom=0.0;} // Ac_s:PAW_PBE:06Sep2000
  if(AUID=="335f20e70f06b78b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.093931;volume_atom=45.40978;spin_atom=0.0;} // Ac:PAW_PBE:06Sep2000
  if(AUID=="5ece90a899819730" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.662485;volume_atom=40.48842;spin_atom=0.0;} // Ac:PAW_LDA_KIN:12Apr2000
  if(AUID=="863169f4af8ed89f" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.637417;volume_atom=40.48951;spin_atom=0.0;} // Ac_s:PAW_LDA:12Apr2000
  if(AUID=="b35f9d074448e332" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.04764;volume_atom=45.53578;spin_atom=0.0;} // Ac:PAW_PBE_KIN:06Sep2000
  if(AUID=="b35f9d074448e332" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-90.49116;volume_atom=43.98308;spin_atom=0.0;} // Ac:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="b5bf833f40cc220c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.045405;volume_atom=44.99825;spin_atom=0.0;} // Ac_s:PAW_GGA:11Apr2000
  if(AUID=="f04c8b0982efae26" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.064376;volume_atom=45.1122;spin_atom=0.0;} // Ac:PAW_GGA:11Apr2000
  if(AUID=="f4d16e5d594ed3ef" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.648738;volume_atom=40.45141;spin_atom=0.0;} // Ac:PAW_LDA:12Apr2000

  // Ag
  // find /common/LIB1/LIB/Ag* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="42cafd70767f07e7" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.727438;volume_atom=17.63908;spin_atom=0.0;} // Ag_GW:PAW_PBE_KIN:06Mar2008
  if(AUID=="4c76b3914a101101" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.711702;volume_atom=15.6433;spin_atom=0.0;} // Ag_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="639bb9917452d998" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.744496;volume_atom=16.20044;spin_atom=0.0;} // Ag:LDA:01Apr2000
  if(AUID=="6f991220b199a503" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.763214;volume_atom=15.83572;spin_atom=0.0;} // Ag_pv:PAW_LDA_KIN:09Dec2005
  if(AUID=="9144029176631616" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.745994;volume_atom=16.06125;spin_atom=0.0;} // Ag:PAW_LDA:17Apr2000
  if(AUID=="b38744d8573af1ff" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.796279;volume_atom=15.89569;spin_atom=0.0;} // Ag_GW:PAW_LDA_KIN:06Mar2008
  if(AUID=="c2b59d1e03816af0" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.6531552;volume_atom=17.21056;spin_atom=0.0;} // Ag_sv_GW:PAW_PBE_KIN:05Dec2013
  if(AUID=="cc0097ef8847f6c4" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.718367;volume_atom=17.65138;spin_atom=0.0;} // Ag:PAW_PBE_KIN:02Apr2005
  if(AUID=="cc0097ef8847f6c4" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-31.38398;volume_atom=16.9292;spin_atom=0.0;} // Ag:PAW_PBE_KIN:SCAN:02Apr2005
  if(AUID=="dac6684dedca50c0" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.720421;volume_atom=17.98532;spin_atom=0.0;} // Ag:GGA:01Apr2000
  if(AUID=="dc0f72d3f1bd620d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.726841;volume_atom=17.78507;spin_atom=0.0;} // Ag:PAW_GGA:18Jul2000
  if(AUID=="e7b3b0854117d479" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.784394;volume_atom=15.91371;spin_atom=0.0;} // Ag:PAW_LDA_KIN:06Sep2000
  if(AUID=="eca65d7d992efb48" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.82746;volume_atom=17.85122;spin_atom=0.0;} // Ag:PAW_PBE:06Sep2000
  if(AUID=="f4452364bfcfa5d4" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.701376;volume_atom=17.48231;spin_atom=0.0;} // Ag_pv:PAW_PBE_KIN:09Dec2005
  if(AUID=="f4452364bfcfa5d4" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-31.3515;volume_atom=16.84626;spin_atom=0.0;} // Ag_pv:PAW_PBE_KIN:SCAN:09Dec2005

  // Al
  // find /common/LIB1/LIB/Al* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="2510ee85c48b48c2" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.187269;volume_atom=15.73444;spin_atom=0.0;} // Al_sv_GW:PAW_LDA:2Feb2008
  if(AUID=="2b253cb3349b4835" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.191843;volume_atom=15.80178;spin_atom=0.0;} // Al:PAW_LDA:17Apr2000
  if(AUID=="35865a9b2b516833" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.744696;volume_atom=16.47165;spin_atom=0.0;} // Al:PAW_PBE_KIN:04Jan2001
  if(AUID=="35865a9b2b516833" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.748459;volume_atom=16.1309;spin_atom=0.0;} // Al:PAW_PBE_KIN:SCAN:04Jan2001
  if(AUID=="423324aa45d0242f" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.743564;volume_atom=16.47192;spin_atom=0.0;} // Al:PAW_PBE:04Jan2001
  if(AUID=="624cbf44af81780c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.194643;volume_atom=15.74631;spin_atom=0.0;} // Al:LDA:01Apr2000
  if(AUID=="81d77d1303d46a24" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.190916;volume_atom=15.80575;spin_atom=0.0;} // Al:PAW_LDA_KIN:17Apr2000
  if(AUID=="9ed43f446f226080" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.751211;volume_atom=16.45294;spin_atom=0.0;} // Al_GW:PAW_PBE_KIN:19Mar2012
  if(AUID=="a0a3e933f8f5dcc7" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.720925;volume_atom=16.36401;spin_atom=0.0;} // Al:GGA:01Apr2000
  if(AUID=="b9cdcfa63ddf1ea9" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.212536;volume_atom=15.12083;spin_atom=0.0;} // Al_h:PAW_GGA:08Apr2002
  if(AUID=="de9492ed3ab234ca" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.800213;volume_atom=14.87718;spin_atom=0.0;} // Al_h:PAW_PBE:08Apr2002
  if(AUID=="e7afb3f96db37614" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.691939;volume_atom=16.52734;spin_atom=0.0;} // Al:GGA:01Apr2000
  if(AUID=="ec638285f61e73c3" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.231;volume_atom=15.60734;spin_atom=0.0;} // Al:LDA:01Apr2000
  if(AUID=="f22666edc17781c8" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.198066;volume_atom=15.7813;spin_atom=0.0;} // Al_GW:PAW_LDA_KIN:19Mar2012
  if(AUID=="f2de0641bce99445" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.695804;volume_atom=16.56004;spin_atom=0.0;} // Al:PAW_GGA:05Jan2001
  if(AUID=="ff45f346e944433a" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.764606;volume_atom=16.23023;spin_atom=0.0;} // Al_sv_GW:PAW_PBE:2Feb2008

  // Ar
  // find /common/LIB1/LIB/Ar* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="33d4ba6d5d91a625" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.1416475;volume_atom=29.31217;spin_atom=0.0;} // Ar_GW:PAW_LDA_KIN:21Aug2013
  if(AUID=="4dc09f7b5e65769e" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.1445373;volume_atom=29.62463;spin_atom=0.0;} // Ar:PAW_LDA_KIN:07Sep2000
  if(AUID=="5adcf9fd506a6f8f" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.07957;volume_atom=50.79861;spin_atom=0.0;} // Ar:PAW_GGA:06Sep2000
  if(AUID=="5ea4f89868f4ce2e" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.145017;volume_atom=29.6839;spin_atom=0.0;} // Ar:PAW_LDA:07Sep2000
  if(AUID=="76ed84ac7c922a74" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.065575;volume_atom=44.08257;spin_atom=0.0;} // Ar:PAW_PBE:07Sep2000
  if(AUID=="c50dbc48fec4befc" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.04020414;volume_atom=44.08101;spin_atom=0.0;} // Ar:PAW_PBE_KIN:07Sep2000
  if(AUID=="c50dbc48fec4befc" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.917285;volume_atom=36.21949;spin_atom=0.0;} // Ar:PAW_PBE_KIN:SCAN:07Sep2000
  if(AUID=="efa3f1a7f900b645" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.0431514;volume_atom=48.27225;spin_atom=0.0;} // Ar_GW:PAW_PBE_KIN:21Aug2013

  // As
  if(AUID=="f1de4433a638eae7" && nKIN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-4.65243;volume_atom=22.6706;spin_atom=0.0;} // As:PAW_PBE:06Sep2000
  if(AUID=="1170cc199ae6ad79" && SCAN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-20.2444;volume_atom=21.5214;spin_atom=0.0;} // As_d:PAW_PBE_KIN:SCAN:11Apr2003
  if(AUID=="b7b512c9cb6a6957" && SCAN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-20.1938;volume_atom=21.8855;spin_atom=0.0;} // As:PAW_PBE_KIN:SCAN:22Sep2009
  if(AUID=="1170cc199ae6ad79" && nKIN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-4.67878;volume_atom=22.5925;spin_atom=0.0;} // As_d:PAW_PBE_KIN:11Apr2003

  // Au
  // find /common/LIB1/LIB/Au* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="101ad6a2dfb494d2" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.771285;volume_atom=16.27679;spin_atom=0.0;} // Au_sv_GW:PAW_LDA_KIN:13Sep2013
  if(AUID=="23da4a7f5773027a" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.463428;volume_atom=16.50809;spin_atom=0.0;} // Au_GW:PAW_LDA_KIN:23Mar2010
  if(AUID=="3477c9645a808524" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.272119;volume_atom=18.03681;spin_atom=0.0;} // Au:PAW_PBE:06Sep2000
  if(AUID=="57972e71335a410c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.203148;volume_atom=18.04886;spin_atom=0.0;} // Au:PAW_GGA:18Jul2000
  if(AUID=="6b9f89c79a20072f" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.448969;volume_atom=16.53801;spin_atom=0.0;} // Au:PAW_LDA_KIN:04Oct2007
  if(AUID=="6db06caa45317191" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.593645;volume_atom=17.46777;spin_atom=0.0;} // Au_sv_GW:PAW_PBE_KIN:13Sep2013
  if(AUID=="834bf85c5ba24d77" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.396411;volume_atom=16.66323;spin_atom=0.0;} // Au:PAW_LDA:04Feb1998
  if(AUID=="9f4c25db32e99a1b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.233364;volume_atom=17.79449;spin_atom=0.0;} // Au_GW:PAW_PBE_KIN:23Mar2010
  if(AUID=="c42e6bdfe7e7024e" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.38948;volume_atom=16.76765;spin_atom=0.0;} // Au:LDA:01Apr2000
  if(AUID=="e2936a3049f7d935" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.223078;volume_atom=17.83091;spin_atom=0.0;} // Au:PAW_PBE_KIN:04Oct2007
  if(AUID=="e2936a3049f7d935" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-69.20261;volume_atom=17.00782;spin_atom=0.0;} // Au:PAW_PBE_KIN:SCAN:04Oct2007
  if(AUID=="fdbc4d4856dea2ca" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.193035;volume_atom=18.19974;spin_atom=0.0;} // Au:GGA:01Apr2000

  // B
  if(AUID=="00b4dfcc28b5887b" && nKIN) {found=TRUE;groundstate_structure="ICSD_56992";groundstate_energy=-6.691806;volume_atom=7.248621;spin_atom=0.0;} // B_h:PAW_GGA:18Jul2000
  if(AUID=="76781ebe8489383f" && nKIN) {found=TRUE;groundstate_structure="ICSD_56992";groundstate_energy=-6.693257;volume_atom=7.241668;spin_atom=0.0;} // B_h:PAW_PBE:07Sep2000
  if(AUID=="70110ee6c6cbaf90" && nKIN) {found=TRUE;groundstate_structure="ICSD_56992";groundstate_energy=-7.469139;volume_atom=6.998381;spin_atom=0.0;} // B_h:PAW_LDA:17Apr2000
  if(AUID=="a7226ec232e2ee27" && nKIN) {found=TRUE;groundstate_structure="ICSD_56992";groundstate_energy=-6.69039;volume_atom=7.248991;spin_atom=0.0;} // B_s:PAW_PBE:22Jan2003

  // Ba
  if(AUID=="062068333d4e7e80" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.92400;volume_atom=63.1759;spin_atom=0.0;} // Ba_sv:PAW_PBE:06Sep2000
  if(AUID=="0ea886f29a8f18e6" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.92401;volume_atom=62.1307;spin_atom=0.0;} // Ba_sv:PAW_GGA:14Apr2000
  if(AUID=="b5d1dafcf2f8d798" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-2.19536;volume_atom=55.5803;spin_atom=0.0;} // Ba:LDA:01Apr2000
  if(AUID=="45efe879c954c89e" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.89600;volume_atom=62.4601;spin_atom=0.0;} // Ba:GGA:01Apr2000
  if(AUID=="7cbe610c465db194" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-2.23839;volume_atom=53.9574;spin_atom=0.0;} // Ba_sv:PAW_LDA:17Apr2000
  if(AUID=="63fd0fe091a069b0" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-39.5316;volume_atom=63.8147;spin_atom=0.0;} // Ba_sv:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="63fd0fe091a069b0" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.90889;volume_atom=63.2532;spin_atom=0.0;} // Ba_sv:PAW_PBE_KIN:06Sep2000
  if(AUID=="5b4abc32c2c70aaf" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-4.74296;volume_atom=62.53586;spin_atom=0.0;} // Ba_sv_GW:PAW_PBE_KIN:23Mar2010
  if(AUID=="91a1497765539c20" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-2.250323;volume_atom=53.98184;spin_atom=0.0;} // Ba_sv:PAW_LDA_KIN:17Apr2000
  if(AUID=="63fd0fe091a069b0" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.908961;volume_atom=63.22138;spin_atom=0.0;} // Ba_sv:PAW_PBE_KIN:06Sep2000
  if(AUID=="8c5d3fb64ef288e5" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-5.12931;volume_atom=52.92512;spin_atom=0.0;} // Ba_sv_GW:PAW_LDA_KIN:23Mar2010

  // Be
  if(AUID=="abd9038ab359fe8e" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.19956;volume_atom=7.60317;spin_atom=0.0;} // Be_sv:PAW_LDA:23Feb1998
  if(AUID=="80136f0733d30e11" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-3.74102;volume_atom=7.92019;spin_atom=0.0;} // Be_sv:PAW_PBE:06Sep2000
  if(AUID=="250ccb49102b1ab7" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-3.75366;volume_atom=7.92077;spin_atom=0.0;} // Be:PAW_PBE:06Sep2000
  if(AUID=="362b70d9e55ede03" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-3.72800;volume_atom=7.91229;spin_atom=0.0;} // Be:PAW_GGA:11Feb1998
  if(AUID=="0134da61389a9f60" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-3.73743;volume_atom=7.84295;spin_atom=0.0;} // Be:GGA:01Apr2000
  if(AUID=="653b8d1063d6b4e0" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.19996;volume_atom=7.57850;spin_atom=0.0;} // Be:LDA:01Apr2000
  if(AUID=="d89972c63683236a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.21025;volume_atom=7.57396;spin_atom=0.0;} // Be:PAW_LDA:02Feb1998
  if(AUID=="b72167bbb69596d1" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.36112;volume_atom=7.92170;spin_atom=0.0;} // Be:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="b72167bbb69596d1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-3.7651;volume_atom=7.91695;spin_atom=0.0;}  // Be:PAW_PBE_KIN:06Sep2000
  if(AUID=="ca7aaa39015dc4f6" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-3.76207;volume_atom=7.93251;spin_atom=0.0;} // Be_sv:PAW_PBE:06Sep2000
  if(AUID=="3c1c54261ede6034" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.206423;volume_atom=7.566058;spin_atom=0.0;} // Be_GW:PAW_LDA_KIN:04Mar2008
  if(AUID=="9c7dbe99a9ffabaf" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-3.70811;volume_atom=7.922699;spin_atom=0.0;} // Be_GW:PAW_PBE_KIN:04Mar2008
  if(AUID=="54e866d928c8d2d6" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.208449;volume_atom=7.572612;spin_atom=0.0;} // Be:PAW_LDA_KIN:02Feb1998
  if(AUID=="22e3a620c52a6e6e" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.29211;volume_atom=7.507654;spin_atom=0.0;} // Be_sv_GW:PAW_LDA:31Mar2010
  if(AUID=="b68d3dd58972282d" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.99249;volume_atom=7.864782;spin_atom=0.0;} // Be_sv_GW:PAW_PBE:31Mar2010
  if(AUID=="90caaee04c9d7fc2" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.203062;volume_atom=7.599556;spin_atom=0.0;} // Be_sv:PAW_LDA:26Mar2012

  // Bi
  if(AUID=="8b1158914bff5430" && nKIN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-3.87274;volume_atom=36.8656;spin_atom=0.0;} // Bi:PAW_PBE:08Apr2002
  if(AUID=="bba01b714d812459" && nKIN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-4.03716;volume_atom=36.2786;spin_atom=0.0;} // Bi_d:PAW_PBE:06Sep2000
  if(AUID=="2efd860034c70e39" && SCAN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-78.5112;volume_atom=35.6425;spin_atom=0.0;} // Bi:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="46860e6586f9539f" && SCAN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-78.5664;volume_atom=34.5043;spin_atom=0.0;} // Bi_d:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="46860e6586f9539f" && nKIN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-3.88235;volume_atom=36.2564;spin_atom=0.0;} // Bi_d:PAW_PBE_KIN:06Sep2000
  if(AUID=="2efd860034c70e39" && nKIN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-3.87117;volume_atom=36.8949;spin_atom=0.0;} // Bi:PAW_PBE_KIN:08Apr2002

  // Br
  if(AUID=="c53210bda32c7827" && nKIN) {found=TRUE;groundstate_structure="A11";groundstate_energy=-1.5898;volume_atom=45.5907;spin_atom=0.0;} // Br:PAW_PBE:06Sep2000
  if(AUID=="42d23b031e31632a" && nKIN) {found=TRUE;groundstate_structure="A11";groundstate_energy=-1.569501;volume_atom=39.0509;spin_atom=0.0;} // Br:PAW_PBE_KIN:06Sep2000

  // C
  if(AUID=="f3800ea7e3e20d86" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-9.19568;volume_atom=10.6949;spin_atom=0.0;} // C_h:PAW_PBE:20Dec2001
  if(AUID=="a720f5418b5ad14f" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-9.22034;volume_atom=10.4453;spin_atom=0.0;} // C:PAW_PBE:08Apr2002
  if(AUID=="8d30df5759a1e14d" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-9.19508;volume_atom=10.5847;spin_atom=0.0;} // C_s:PAW_PBE:06Sep2000
  if(AUID=="8397f0cab3a0348d" && SCAN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-10.1261;volume_atom=8.99501;spin_atom=0.0;} // C_h:PAW_PBE_KIN:SCAN:06Feb2004
  if(AUID=="f9bce748fcbf37b3" && SCAN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-10.0941;volume_atom=9.02062;spin_atom=0.0;} // C:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="f9bce748fcbf37b3" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-9.225506;volume_atom=11.01637;spin_atom=0.0;} // C:PAW_PBE_KIN:08Apr2002
  if(AUID=="00d46f7c41fd745f" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-10.11818;volume_atom=8.628322;spin_atom=0.0;} // C_h:PAW_LDA:17Apr2000
  if(AUID=="31c540de1efdd804" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-9.268609;volume_atom=11.13834;spin_atom=0.0;} // C:GGA:01Apr2000
  if(AUID=="3a89016e2e639e47" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-9.216902;volume_atom=11.17504;spin_atom=0.0;} // C_h:PAW_GGA:18Jul2000
  if(AUID=="3ff7c326f3a2d786" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-9.141556;volume_atom=10.85588;spin_atom=0.0;} // C_s:GGA:01Apr2000
  if(AUID=="459fb1e14999a1a8" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-10.09531;volume_atom=8.513899;spin_atom=0.0;} // C_s:PAW_LDA:04May1998
  if(AUID=="58ee909441d7673b" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-9.240065;volume_atom=11.17413;spin_atom=0.0;} // C:PAW_GGA:05Jan2001
  if(AUID=="73c35c1982ecc74b" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-10.00952;volume_atom=8.650105;spin_atom=0.0;} // C_s:LDA:01Apr2000
  if(AUID=="97cfc186c185e847" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-9.221835;volume_atom=11.20526;spin_atom=0.0;} // C_s:PAW_GGA:08Oct1999
  if(AUID=="9e8ca933cd7d349f" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-10.11561;volume_atom=8.569782;spin_atom=0.0;} // C:PAW_LDA:31May2000
  if(AUID=="d951d3c4d963900a" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-10.15636;volume_atom=8.521433;spin_atom=0.0;} // C:LDA:01Apr2000

  // Ca (19)
  // find /common/LIB1/LIB/Ca* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="1bd33f6e3f1fd4b0" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.08587;volume_atom=37.34595;spin_atom=0.0;} // Ca_sv_GW:PAW_LDA_KIN:31Mar2010
  if(AUID=="26e7ac00e1d3cd41" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.210568;volume_atom=37.71927;spin_atom=0.0;} // Ca_pv:PAW_LDA:05May1998
  if(AUID=="3bdd926d04dbca1e" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.183257;volume_atom=37.8791;spin_atom=0.0;} // Ca_sv:PAW_LDA:17Apr2000
  if(AUID=="3e3d86ebf8f9e84a" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.931989;volume_atom=41.58018;spin_atom=0.0;} // Ca:GGA:01Apr2000
  if(AUID=="529aebe19b9289ea" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.914223;volume_atom=41.75753;spin_atom=0.0;} // Ca_pv:PAW_PBE_KIN:06Sep2000
  if(AUID=="529aebe19b9289ea" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-9.663724;volume_atom=43.00712;spin_atom=0.0;} // Ca_pv:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="560b47e7f6be795e" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.917229;volume_atom=41.82487;spin_atom=0.0;} // Ca_sv:PAW_GGA:04May1998
  if(AUID=="57183fbccf86bacc" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.216359;volume_atom=37.6043;spin_atom=0.0;} // Ca_pv:PAW_LDA_KIN:05May1998
  if(AUID=="572c23d8cdb07c91" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.947546;volume_atom=41.33207;spin_atom=0.0;} // Ca:GGA:01Apr2000
  if(AUID=="5ecdefc3dd3c9a74" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.848079;volume_atom=41.49652;spin_atom=0.0;} // Ca_sv_GW:PAW_PBE_KIN:31Mar2010
  if(AUID=="7f030fc2723afab5" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.907499;volume_atom=42.55465;spin_atom=0.0;} // Ca:PAW_GGA:10Feb1998
  if(AUID=="818406d11e64bddc" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.930055;volume_atom=42.18691;spin_atom=0.0;} // Ca_sv:PAW_PBE_KIN:06Sep2000
  if(AUID=="818406d11e64bddc" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-9.662664;volume_atom=42.94814;spin_atom=0.0;} // Ca_sv:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="a5df1335919961ad" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.001074;volume_atom=42.15212;spin_atom=0.0;} // Ca_sv:PAW_PBE:06Sep2000
  if(AUID=="ae8076243474d29f" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.20423;volume_atom=37.85516;spin_atom=0.0;} // Ca_sv:PAW_LDA_KIN:17Apr2000
  if(AUID=="b216760dca178acb" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.976382;volume_atom=41.79727;spin_atom=0.0;} // Ca_pv:PAW_PBE:06Sep2000
  if(AUID=="d0c3b753999f51e0" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.911633;volume_atom=41.56355;spin_atom=0.0;} // Ca_pv:PAW_GGA:05May1998
  if(AUID=="d85b55fde567529b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.242847;volume_atom=37.32481;spin_atom=0.0;} // Ca:LDA:01Apr2000
  if(AUID=="e65bce3da30c5d46" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.182359;volume_atom=38.29302;spin_atom=0.0;} // Ca:LDA:01Apr2000 (this is _pv)

  // Cd
  if(AUID=="0b0758f9de1eb816" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-0.761895;volume_atom=22.4448;spin_atom=0.0;} // Cd:PAW_GGA:04May1998
  if(AUID=="3388a058d6fd344d" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-0.906233;volume_atom=22.4532;spin_atom=0.0;} // Cd:PAW_PBE:06Sep2000
  if(AUID=="c5c179cf112468e6" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.503460;volume_atom=19.7986;spin_atom=0.0;} // Cd:PAW_LDA:03Mar1998
  if(AUID=="07875e5a0789cd01" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.530730;volume_atom=20.0760;spin_atom=0.0;} // Cd:LDA:01Apr2000
  if(AUID=="f0b4c3233e3b19d8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-0.778952;volume_atom=22.8753;spin_atom=0.0;} // Cd:GGA:01Apr2000
  if(AUID=="fdf368e153b383d1" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-30.40360;volume_atom=21.1426;spin_atom=0.0;} // Cd:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="fdf368e153b383d1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-0.743299;volume_atom=22.4275;spin_atom=0.0;} // Cd:PAW_PBE_KIN:06Sep2000
  if(AUID=="63ded7f1e98bd74e" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.53873;volume_atom=19.7295;spin_atom=0.0;} // Cd_f_GW:PAW_LDA_KIN:18May2010
  if(AUID=="8e504fb676bf06a4" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-0.7555161;volume_atom=22.34894;spin_atom=0.0;} // Cd_f_GW:PAW_PBE_KIN:18May2010
  if(AUID=="ab3a0195ab049ec8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.518227;volume_atom=19.8174;spin_atom=0.0;} // Cd:PAW_LDA_KIN:03Mar1998
  if(AUID=="74bf44fed86520ba" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.523063;volume_atom=19.13637;spin_atom=0.0;} // Cd_sv_GW:PAW_LDA_KIN:16Apr2014
  if(AUID=="43daba3719e5e3b4" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-0.7795479;volume_atom=21.42407;spin_atom=0.0;} // Cd_sv_GW:PAW_PBE_KIN:16Apr2014

  // Ce
  if(AUID=="510a8ebc02942ebd" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.256439;volume_atom=27.3549;spin_atom=0.0;} // Ce_s:PAW_GGA:11May2000
  if(AUID=="27449258606dd200" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.645601;volume_atom=25.52952;spin_atom=0.0;} // Ce_s:PAW_LDA:17Apr2000
  if(AUID=="6e846dc10a52b11b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.930997;volume_atom=26.05794;spin_atom=0.0;} // Ce:PAW_PBE:28Sep2000
  if(AUID=="43bef7da77f3563c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.899931;volume_atom=22.79678;spin_atom=0.0;} // Ce:PAW_LDA:28Sep2000
  if(AUID=="af507d55827ec15c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.964473;volume_atom=26.03088;spin_atom=0.0;} // Ce:PAW_GGA:29Sep2000
  if(AUID=="abfa73bca84b87e0" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.77123;volume_atom=37.644;spin_atom=0.0;} // Ce_3:PAW_PBE:06Sep2000
  if(AUID=="c30b547637edabc3" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.74004;volume_atom=37.32093;spin_atom=0.0;} // Ce_3:PAW_GGA:11May2000
  if(AUID=="24ba1850e358faff" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-44.56046;volume_atom=24.57593;spin_atom=0.0;} // Ce_h:PAW_PBE_KIN:SCAN:03Mar2005
  if(AUID=="a896fb0aac886619" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-43.13392;volume_atom=37.35713;spin_atom=0.0;} // Ce_3:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="bc711278cc7e92ff" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-44.61209;volume_atom=24.63577;spin_atom=0.0;} // Ce:PAW_PBE_KIN:SCAN:23Dec2003

  // Cl
  if(AUID=="6bf17162620b7ce3" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-1.787206;volume_atom=500;spin_atom=0.0;} // Cl:PAW_PBE:17Jan2003
  if(AUID=="6dca7e86aecb68cc" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-1.776937;volume_atom=500;spin_atom=0.0;} // Cl_h:PAW_PBE:08Apr2002
  if(AUID=="9d680135d96b4029" && SCAN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-7.906078;volume_atom=500;spin_atom=0.0;} // Cl:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="9d680135d96b4029" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-1.787335;volume_atom=500;spin_atom=0.0;} // Cl:PAW_PBE_KIN:06Sep2000

  // Co
  if(AUID=="41f6635e51e063f9" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.13760;volume_atom=10.0906;spin_atom=1.51284;} // Co:LDA:01Apr2000
  if(AUID=="4fb418aa1d0f1607" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.05079;volume_atom=10.9889;spin_atom=1.62675;} // Co:GGA:01Apr2000
  if(AUID=="f3c42ce518194ac4" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.99076;volume_atom=10.8013;spin_atom=1.57160;} // Co:PAW_GGA:03Mar1998
  if(AUID=="9c4653e42c9f870b" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.10876;volume_atom=10.8414;spin_atom=1.60252;} // Co:PAW_PBE:06Sep2000
  if(AUID=="23d97019de916d58" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.11138;volume_atom=9.98379;spin_atom=1.49449;} // Co:PAW_LDA:03Mar1998
  if(AUID=="4cc3ed3cb4d06e66" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-17.2118;volume_atom=10.3970;spin_atom=1.72628;} // Co:PAW_PBE_KIN:SCAN:02Aug2007
  if(AUID=="6d93f7d0ea628829" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-17.2429;volume_atom=10.4421;spin_atom=1.72940;} // Co_pv:PAW_PBE_KIN:SCAN:23Apr2009
  if(AUID=="12d6a4385e456ff2" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-17.3028;volume_atom=10.2840;spin_atom=1.72626;} // Co_sv:PAW_PBE_KIN:SCAN:23Jul2007
  if(AUID=="12d6a4385e456ff2" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.11042;volume_atom=10.6401;spin_atom=1.58914;} // Co_sv:PAW_PBE_KIN:23Jul2007
  if(AUID=="57bd04e1209ca57e" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.053049;volume_atom=10.77001;spin_atom=1.594416;} // Co_GW:PAW_PBE_KIN:31Mar2010
  if(AUID=="55caba6c967c2367" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.144426;volume_atom=9.940874;spin_atom=1.490024;} // Co_GW:PAW_LDA_KIN:31Mar2010
  if(AUID=="6637b9e694868663" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.129802;volume_atom=9.946954;spin_atom=1.489165;} // Co:PAW_LDA_KIN:26Mar2009
  if(AUID=="4cc3ed3cb4d06e66" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.037893;volume_atom=10.7764;spin_atom=1.590394;} // Co:PAW_PBE_KIN:02Aug2007
  if(AUID=="6d93f7d0ea628829" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.030628;volume_atom=10.78008;spin_atom=1.593403;} // Co_pv:PAW_PBE_KIN:23Apr2009
  if(AUID=="9273ad55118191e1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.116928;volume_atom=9.953821;spin_atom=1.489631;} // Co_pv:PAW_LDA_KIN:15Mar2012
  if(AUID=="552096566047fff9" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.199133;volume_atom=9.840014;spin_atom=1.486248;} // Co_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="86d7a01ba3c1b75e" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.184125;volume_atom=10.60544;spin_atom=1.588372;} // Co_sv_GW:PAW_PBE_KIN:05Dec2013
  if(AUID=="ac83e5a2e9a83028" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.1571;volume_atom=9.866361;spin_atom=1.488013;} // Co_sv:PAW_LDA_KIN:23Jul2007

  // Cr
  if(AUID=="e56aa7dff1571851" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.6156;volume_atom=10.8356;spin_atom=0.0;} // Cr_pv:PAW_LDA:07Sep2000
  if(AUID=="411ed623d1b5054d" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.51285;volume_atom=11.3963;spin_atom=0.0;} // Cr:PAW_PBE:06Sep2000
  if(AUID=="1365c4634d5c2ea1" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.46687;volume_atom=11.3926;spin_atom=0.0;} // Cr:PAW_GGA:03Mar1998
  if(AUID=="f020a26862c0a532" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.62940;volume_atom=11.5175;spin_atom=0.0;} // Cr_pv:PAW_PBE:07Sep2000
  if(AUID=="7ab8ddc788e0ece9" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.6036;volume_atom=10.7989;spin_atom=0.0;} // Cr:LDA:01Apr2000
  if(AUID=="235696d86e852bef" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.43370;volume_atom=11.5367;spin_atom=0.0;} // Cr:GGA:01Apr2000
  if(AUID=="caa6b9695b37b7c7" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.46956;volume_atom=11.5096;spin_atom=0.0;} // Cr_pv:PAW_GGA:07Sep2000
  if(AUID=="46e25663256eec1e" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.6345;volume_atom=10.7192;spin_atom=0.0;} // Cr:PAW_LDA:03Mar1998
  if(AUID=="46684a319ca854a2" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-18.2521;volume_atom=11.1136;spin_atom=0.0;} // Cr:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="c349e20347a8f61f" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-18.3863;volume_atom=11.0185;spin_atom=0.0;} // Cr_sv:PAW_PBE_KIN:SCAN:23Jul2007
  if(AUID=="d84af42941e8d51d" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-18.3160;volume_atom=11.2542;spin_atom=0.0;} // Cr_pv:PAW_PBE_KIN:SCAN:02Aug2007
  if(AUID=="46684a319ca854a2" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.48959;volume_atom=11.3875;spin_atom=0.0;} // Cr:PAW_PBE_KIN:06Sep2000
  if(AUID=="d84af42941e8d51d" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.51042;volume_atom=11.4657;spin_atom=0.0;} // Cr_pv:PAW_PBE_KIN:02Aug2007
  if(AUID=="86c325658e1c2bb2" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.63447;volume_atom=10.72272;spin_atom=0.0;} // Cr:PAW_LDA_KIN:03Mar1998
  if(AUID=="f5a28ee5ce72c4f2" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.63269;volume_atom=10.79429;spin_atom=0.0;} // Cr_pv:PAW_LDA_KIN:02Aug2007
  if(AUID=="81df1668e8d9304c" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.35907;volume_atom=10.76206;spin_atom=0.0;} // Cr_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="145bfb04c18c9c72" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.29901;volume_atom=11.42431;spin_atom=0.0;} // Cr_sv_GW:PAW_PBE_KIN:05Dec2013
  if(AUID=="c349e20347a8f61f" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.601694;volume_atom=11.38384;spin_atom=0.0;} // Cr_sv:PAW_PBE_KIN:23Jul2007
  if(AUID=="79d9b9059234f4d3" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.173256;volume_atom=11.96508;spin_atom=0.0;} // Cr_sv:PAW_LDA_KIN:23Jul2007

  // Cs
  if(AUID=="42d5fdd268b1a941" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-2.707642;volume_atom=94.24478;spin_atom=0.0;} // Cs_sv_GW:PAW_LDA_KIN:23Mar2010
  if(AUID=="47d2b2b79b3bb561" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-0.9922865;volume_atom=94.77585;spin_atom=0.0;} // Cs_sv:PAW_LDA:07Sep2000
  if(AUID=="485f7a26842607ea" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-0.8521117;volume_atom=115.6232;spin_atom=0.0;} // Cs_sv:PAW_PBE_KIN:08Apr2002
  if(AUID=="a7b37bbd60c921b0" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-0.9924374;volume_atom=94.68309;spin_atom=0.0;} // Cs_sv:PAW_LDA_KIN:07Sep2000
  if(AUID=="e33c8a4092f5c27c" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-2.526557;volume_atom=114.9539;spin_atom=0.0;} // Cs_sv_GW:PAW_PBE_KIN:23Mar2010
  if(AUID=="e054d996abee22c5" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-0.852451;volume_atom=115.5724;spin_atom=0.0;} // Cs_sv:PAW_PBE:08Apr2002

  // Cu
  // find /common/LIB1/LIB/Cu* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="1794dd824ade7bfb" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.61717;volume_atom=10.7855;spin_atom=0.0;} // Cu_pv:PAW_LDA:19Apr2000
  if(AUID=="cd77b2d6e432ab5d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.72886;volume_atom=11.9621;spin_atom=0.0;} // Cu:PAW_GGA:05Jan2001
  if(AUID=="82741f3fd74a85b4" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.76087;volume_atom=12.0520;spin_atom=0.0;} // Cu:GGA:01Apr2000
  if(AUID=="62263cf2c885649c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.09701;volume_atom=11.8086;spin_atom=0.0;} // Cu_pv:PAW_PBE:06Sep2000
  if(AUID=="c7a1b91d209a648b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.71935;volume_atom=11.9654;spin_atom=0.0;} // Cu:PAW_PBE:05Jan2001
  if(AUID=="920126d2172870b5" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.70071;volume_atom=10.8985;spin_atom=0.0;} // Cu:PAW_LDA:03Mar1998
  if(AUID=="88d5d10238695c60" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.75059;volume_atom=10.9687;spin_atom=0.0;} // Cu:LDA:01Apr2000
  if(AUID=="a81d4798925d371d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.63758;volume_atom=11.7904;spin_atom=0.0;} // Cu_pv:PAW_GGA:19Apr2000
  if(AUID=="2ce0157f813571ea" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-15.1713;volume_atom=11.1286;spin_atom=0.0;} // Cu_pv:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="34032e61f91e67dc" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-15.1334;volume_atom=11.2002;spin_atom=0.0;} // Cu:PAW_PBE_KIN:SCAN:22Jun2005
  if(AUID=="34032e61f91e67dc" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.72818;volume_atom=11.9295;spin_atom=0.0;} // Cu:PAW_PBE_KIN:22Jun2005
  if(AUID=="2ce0157f813571ea" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.7502;volume_atom=11.7714;spin_atom=0.0;} // Cu_pv:PAW_PBE_KIN:06Sep2000
  if(AUID=="425c1e69fb1d3d12" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.750682;volume_atom=11.83027;spin_atom=0.0;} // Cu_GW:PAW_PBE_KIN:19May2006
  if(AUID=="d15bb62d47c7b9be" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.740633;volume_atom=10.812;spin_atom=0.0;} // Cu_GW:PAW_LDA_KIN:19May2006
  if(AUID=="611c1c4066ac8e2c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.709788;volume_atom=10.88644;spin_atom=0.0;} // Cu:PAW_LDA_KIN:22Jun2005
  if(AUID=="aed87d6c59ccf64b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.740238;volume_atom=10.78904;spin_atom=0.0;} // Cu_pv:PAW_LDA_KIN:19Apr2000
  if(AUID=="9180a6efde161dbb" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.83218;volume_atom=10.71468;spin_atom=0.0;} // Cu_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="f99eb576b351d71c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.765172;volume_atom=11.6499;spin_atom=0.0;} // Cu_sv_GW:PAW_PBE_KIN:05Dec2013

  // Dy
  if(AUID=="0373199190a8abb7" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.58813;volume_atom=31.8570;spin_atom=0.0;} // Dy_3:PAW_PBE:06Sep2000
  if(AUID=="57f7009ab6dbe961" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.52172;volume_atom=31.5074;spin_atom=0.0;} // Dy_3:PAW_GGA:10May2000
  if(AUID=="b7037854f4315aa3" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-47.8191;volume_atom=30.5976;spin_atom=0.0;} // Dy_3:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="06df93ed854cd141" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.233;volume_atom=31.1247;spin_atom=0.0;} // Dy:PAW_PBE_KIN:23Dec2003
  if(AUID=="b7037854f4315aa3" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.527491;volume_atom=31.7888;spin_atom=0.0;} // Dy_3:PAW_PBE_KIN:06Sep2000

  // Eu
  if(AUID=="9b801fd8fd52e116" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.23775;volume_atom=43.37685;spin_atom=7.26175;} // Eu:PAW_PBE:08Apr2002
  if(AUID=="2453fc28e3c2ebe1" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.28213;volume_atom=42.56765;spin_atom=7.243493;} // Eu:PAW_PBE_KIN:23Dec2003

  // F
  if(AUID=="3defbdd6f1dad3f2" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-1.790868;volume_atom=500;spin_atom=0.0;} // F_s:PAW_PBE:06Sep2000
  if(AUID=="4c13232a0a03cdc4" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-1.780263;volume_atom=500;spin_atom=0.0;} // F:PAW_PBE_KIN:08Apr2002
  if(AUID=="4c13232a0a03cdc4" && SCAN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-3.16947;volume_atom=500;spin_atom=0.0;} // F:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="675e493e51b0f2f2" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-1.859526;volume_atom=500;spin_atom=0.0;} // F:PAW_PBE:08Apr2002
  if(AUID=="306f88e141704ce5" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-1.908987;volume_atom=500;spin_atom=0.0;} // F_h:PAW_PBE:07Sep2000

  // Fe
  if(AUID=="53451df4122ca435" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.27793;volume_atom=10.5121;spin_atom=2.03940;} // Fe:LDA:01Apr2000
  if(AUID=="8a9cdad8b2fd617c" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.16576;volume_atom=11.2687;spin_atom=2.15865;} // Fe:PAW_GGA:03Mar1998
  if(AUID=="ca7f0085c55a7aa1" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.30490;volume_atom=11.6137;spin_atom=2.31462;} // Fe:GGA:01Apr2000
  if(AUID=="f5bec9ed423a0ff3" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.31138;volume_atom=11.3320;spin_atom=2.19316;} // Fe:PAW_PBE:06Sep2000
  if(AUID=="24dd99c81a522c7f" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.23151;volume_atom=10.3385;spin_atom=1.95457;} // Fe:PAW_LDA:03Mar1998
  if(AUID=="8c4dfe9b38461b3c" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.11742;volume_atom=11.2608;spin_atom=2.15340;} // Fe_pv:PAW_GGA:06May1998
  if(AUID=="15e65f5e50047695" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.45502;volume_atom=11.3277;spin_atom=2.19491;} // Fe_pv:PAW_PBE:06Sep2000
  if(AUID=="ea29acbd64d22bd0" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.17964;volume_atom=10.3583;spin_atom=1.95509;} // Fe_pv:PAW_LDA:03Mar1998
  if(AUID=="46955a9bf291f348" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.27136;volume_atom=11.2061;spin_atom=2.14058;} // Fe_sv:PAW_GGA:14Sep2000
  if(AUID=="4d78633077b48d89" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-18.0562;volume_atom=11.5112;spin_atom=2.62602;} // Fe:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="c3ba95aef6c439d4" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-18.1013;volume_atom=11.6239;spin_atom=2.66996;} // Fe_pv:PAW_PBE_KIN:SCAN:02Aug2007
  if(AUID=="dbda82e42dec4a56" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-18.1874;volume_atom=11.3325;spin_atom=2.62763;} // Fe_sv:PAW_PBE_KIN:SCAN:23Jul2007
  if(AUID=="4d78633077b48d89" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.23903;volume_atom=11.3121;spin_atom=2.18884;} // Fe:PAW_PBE_KIN:06Sep2000
  if(AUID=="c3ba95aef6c439d4" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.256;volume_atom=11.2646;spin_atom=2.18544;} // Fe_pv:PAW_PBE_KIN:02Aug2007
  if(AUID=="d702aa64f30158e8" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.260077;volume_atom=10.28232;spin_atom=1.943251;} // Fe_GW:PAW_LDA_KIN:31Mar2010
  if(AUID=="e976c0bff7b57001" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.255884;volume_atom=11.2649;spin_atom=2.197398;} // Fe_GW:PAW_PBE_KIN:31Mar2010
  if(AUID=="1dc7f4fd5097d6bd" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.227033;volume_atom=10.34655;spin_atom=1.961236;} // Fe:PAW_LDA_KIN:03Mar1998
  if(AUID=="5053a33b2b8de4ce" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.250277;volume_atom=10.31224;spin_atom=1.948558;} // Fe_pv:PAW_LDA_KIN:02Aug2007
  if(AUID=="e46d2489183a0133" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.51587;volume_atom=10.17803;spin_atom=1.908769;} // Fe_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="e1d42714d0665eef" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.516331;volume_atom=11.05797;spin_atom=2.146284;} // Fe_sv_GW:PAW_PBE_KIN:05Dec2013
  if(AUID=="dbda82e42dec4a56" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.338977;volume_atom=11.13428;spin_atom=2.159969;} // Fe_sv:PAW_PBE_KIN:23Jul2007
  if(AUID=="b6b806587202e51a" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.302782;volume_atom=10.21304;spin_atom=1.922156;} // Fe_sv:PAW_LDA_KIN:23Jul2007

  // Ga
  if(AUID=="5a3627680584037d" && nKIN) {found=TRUE;groundstate_structure="A11";groundstate_energy=-2.901519;volume_atom=19.96062;spin_atom=0.0;} // Ga_h:PAW_PBE_KIN:09Apr2002
  if(AUID=="567d70c3d85e2c7f" && nKIN) {found=TRUE;groundstate_structure="A11";groundstate_energy=-2.901526;volume_atom=19.96071;spin_atom=0.0;} // Ga_h:PAW_PBE:09Apr2002
  if(AUID=="47e95a8d720b4f46" && nKIN) {found=TRUE;groundstate_structure="A11";groundstate_energy=-2.926373;volume_atom=19.96145;spin_atom=0.0;} // Ga_h:PAW_RPBE:09Apr2002

  // Ge
  if(AUID=="e58c5d01a1b37232" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-4.49261;volume_atom=24.1781;spin_atom=0.0;} // Ge:PAW_PBE:05Jan2001
  if(AUID=="e8adcacdfca69c78" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.17465;volume_atom=22.4994;spin_atom=0.0;} // Ge:PAW_LDA:03Mar1998
  if(AUID=="fd03d8c50dd891a8" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-4.52888;volume_atom=23.8391;spin_atom=0.0;} // Ge:GGA:01Apr2000
  if(AUID=="23df10d879b3b58b" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-4.50372;volume_atom=23.6953;spin_atom=0.0;} // Ge_h:PAW_PBE:09Apr2002
  if(AUID=="1a9537cf1728245b" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-4.62213;volume_atom=23.7855;spin_atom=0.0;} // Ge_d:PAW_PBE:06Sep2000
  if(AUID=="60aa6eae51da5ce4" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.19244;volume_atom=22.2190;spin_atom=0.0;} // Ge:LDA:01Apr2000
  if(AUID=="8f36ae9587237414" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-4.52998;volume_atom=23.6635;spin_atom=0.0;} // Ge_h:PAW_RPBE:09Apr2002
  if(AUID=="4538ff07d45d4093" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.13954;volume_atom=22.1566;spin_atom=0.0;} // Ge_d:PAW_LDA:03Mar1998
  if(AUID=="584c68479b3b04b1" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-4.48878;volume_atom=23.7431;spin_atom=0.0;} // Ge_d:PAW_GGA:03Mar1998
  if(AUID=="efd53a76a6f5032e" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.18648;volume_atom=22.0633;spin_atom=0.0;} // Ge_h:PAW_LDA:21Jan2003
  if(AUID=="18fa5cbdea403800" && SCAN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-19.0570;volume_atom=22.4931;spin_atom=0.0;} // Ge_d:PAW_PBE_KIN:SCAN:03Jul2007
  if(AUID=="571d96beea462000" && SCAN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-19.0597;volume_atom=22.5845;spin_atom=0.0;} // Ge_h:PAW_PBE_KIN:SCAN:09Apr2002
  if(AUID=="39a1162598519002" && SCAN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-19.0406;volume_atom=22.9647;spin_atom=0.0;} // Ge:PAW_PBE_KIN:SCAN:05Jan2001
  if(AUID=="18fa5cbdea403800" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-4.51764;volume_atom=23.7152;spin_atom=0.0;} // Ge_d:PAW_PBE_KIN:03Jul2007
  if(AUID=="571d96beea462000" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-4.50383;volume_atom=23.6937;spin_atom=0.0;} // Ge_h:PAW_PBE_KIN:09Apr2002
  if(AUID=="39a1162598519002" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-4.49199;volume_atom=24.1778;spin_atom=0.0;} // Ge:PAW_PBE_KIN:05Jan2001

  // H
  if(AUID=="48b40c1ec0fc696e" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-3.369226;volume_atom=500;spin_atom=0.0;} // H_h:PAW_PBE:07Sep2000
  if(AUID=="772cd4ac80b01737" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-3.386066;volume_atom=500;spin_atom=0.0;} // H:PAW_PBE:15Jun2001
  if(AUID=="772cd4ac80b01737" && SCAN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-3.448563;volume_atom=500;spin_atom=0.0;} // H:PAW_PBE:SCAN:15Jun2001
  if(AUID=="b4cc088636d527f0" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-3.386348;volume_atom=500;spin_atom=0.0;} // H:PAW_PBE:15Jun2001

  // He
  if(AUID=="c135275f4a0c77c9" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=0.481074;volume_atom=5.48651;spin_atom=0.0;} // He:PAW_GGA:05Jan2001
  if(AUID=="60575574eefd1ae1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-0.04401;volume_atom=9.55685;spin_atom=0.0;} // He:PAW_LDA:07Sep2000
  if(AUID=="d3ff93e340e01503" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=0.236864;volume_atom=6.96149;spin_atom=0.0;} // He:PAW_PBE:05Jan2001

  // Hf
  if(AUID=="bb7552942568c756" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.88402;volume_atom=21.9665;spin_atom=0.0;} // Hf:GGA:01Apr2000
  if(AUID=="1af9281ed7be5121" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.7693;volume_atom=20.3825;spin_atom=0.0;} // Hf:LDA:01Apr2000
  if(AUID=="1703a3c3d5bdc240" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.87902;volume_atom=22.3346;spin_atom=0.0;} // Hf:PAW_GGA:6May2002
  if(AUID=="2a020b71ae180392" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.95711;volume_atom=22.2563;spin_atom=0.0;} // Hf:PAW_PBE:20Jan2003
  if(AUID=="5e77afe83e49db7f" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.7507;volume_atom=20.7174;spin_atom=0.0;} // Hf:PAW_LDA:21Jan2003
  if(AUID=="d56380d31571155d" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.83075;volume_atom=22.3351;spin_atom=0.0;} // Hf_pv:PAW_GGA:17Apr2000
  if(AUID=="95ae304e242d499a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.6877;volume_atom=20.7631;spin_atom=0.0;} // Hf_pv:PAW_LDA:17Apr2000
  if(AUID=="15832a2c336c16e0" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.95294;volume_atom=22.4349;spin_atom=0.0;} // Hf_pv:PAW_PBE:06Sep2000
  if(AUID=="b5f0d256277acd2b" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-62.7979;volume_atom=21.3757;spin_atom=0.0;} // Hf:PAW_PBE_KIN:SCAN:20Jan2003
  if(AUID=="d78f974a4747b7ff" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-65.4899;volume_atom=21.3256;spin_atom=0.0;} // Hf_sv:PAW_PBE_KIN:SCAN:10Jan2008
  if(AUID=="8be591141a9fad1a" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-62.7249;volume_atom=21.4698;spin_atom=0.0;} // Hf_pv:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="8be591141a9fad1a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.92235;volume_atom=22.3816;spin_atom=0.0;} // Hf_pv:PAW_PBE_KIN:06Sep2000
  if(AUID=="9f5b6eb8e91b7bb6" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.74498;volume_atom=20.72519;spin_atom=0.0;} // Hf:PAW_LDA_KIN:21Jan2003
  if(AUID=="b5f0d256277acd2b" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.958298;volume_atom=22.25371;spin_atom=0.0;} // Hf:PAW_PBE_KIN:20Jan2003
  if(AUID=="117e9d5b1780b640" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.69392;volume_atom=20.79667;spin_atom=0.0;} // Hf_pv:PAW_LDA_KIN:17Apr2000
  if(AUID=="db6f396ed2ce01c6" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.76232;volume_atom=22.33873;spin_atom=0.0;} // Hf_sv_GW:PAW_PBE_KIN:16Jan2015
  if(AUID=="a43077ac3bc180ef" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-13.543;volume_atom=20.73522;spin_atom=0.0;} // Hf_sv_GW:PAW_LDA_KIN:16Jan2015

  // Hg - requires rechecking
  if(AUID=="85cf414c5473e942" && nKIN) {found=TRUE;groundstate_structure="A10(A3)";groundstate_energy=-0.917393;volume_atom=22.3518;spin_atom=0.0;} // Hg:PAW_LDA:04Feb1998
  if(AUID=="bc158f9698c492b9" && nKIN) {found=TRUE;groundstate_structure="A10(A3)";groundstate_energy=-0.232554;volume_atom=29.8458;spin_atom=0.0;} // Hg:GGA:01Apr2000
  if(AUID=="96247dcf8cda9590" && nKIN) {found=TRUE;groundstate_structure="A10(A3)";groundstate_energy=-0.300828;volume_atom=29.3560;spin_atom=0.0;} // Hg:PAW_PBE:06Sep2000
  if(AUID=="df9134b9ac23884a" && nKIN) {found=TRUE;groundstate_structure="A10(A3)";groundstate_energy=-0.197681;volume_atom=28.7964;spin_atom=0.0;} // Hg:PAW_GGA:04May1998
  if(AUID=="03e3289690e1f860" && nKIN) {found=TRUE;groundstate_structure="A10(A3)";groundstate_energy=-0.958689;volume_atom=22.7055;spin_atom=0.0;} // Hg:LDA:01Apr2000
  if(AUID=="56eba2b878545a76" && SCAN) {found=TRUE;groundstate_structure="A10(A3)";groundstate_energy=-68.11140;volume_atom=26.8710;spin_atom=0.0;} // Hg:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="56eba2b878545a76" && nKIN) {found=TRUE;groundstate_structure="A10(A3)";groundstate_energy=-0.158807;volume_atom=27.701;spin_atom=0.0;} // Hg:PAW_PBE_KIN:06Sep2000

  // Ho
  if(AUID=="ac7fc7278038dec4" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.56866;volume_atom=31.3589;spin_atom=0.0;} // Ho_3:PAW_PBE:06Sep2000
  if(AUID=="20bea66adf91e4a2" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.49473;volume_atom=31.0514;spin_atom=0.0;} // Ho_3:PAW_GGA:10May2000
  if(AUID=="141af062a506aee7" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-49.0018;volume_atom=30.1155;spin_atom=0.0;} // Ho_3:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="141af062a506aee7" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.507265;volume_atom=31.29005;spin_atom=0.0;} // Ho_3:PAW_PBE_KIN:06Sep2000

  // K
  if(AUID=="1faff02cfb784b29" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.02995;volume_atom=73.8776;spin_atom=0.0;} // K:GGA:01Apr2000
  if(AUID=="bc9d6fb83de4a3e6" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.17040;volume_atom=64.3539;spin_atom=0.0;} // K:LDA:01Apr2000
  if(AUID=="e9c22b07bc0ec1ed" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.02890;volume_atom=73.4124;spin_atom=0.0;} // K:GGA:01Apr2000
  if(AUID=="90f35bb229eef165" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.02436;volume_atom=73.3337;spin_atom=0.0;} // K_pv:PAW_GGA:11Feb1998
  if(AUID=="828b9d8df4c3212b" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.16884;volume_atom=63.8721;spin_atom=0.0;} // K:LDA:01Apr2000
  if(AUID=="57bec8a3e025c3e4" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.16059;volume_atom=63.6556;spin_atom=0.0;} // K_pv:PAW_LDA:02Feb1998
  if(AUID=="941fbf66b21efd05" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.02692;volume_atom=72.4893;spin_atom=0.0;} // K_pv:PAW_PBE:17Jan2003
  if(AUID=="5c717bdf563a1964" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.03858;volume_atom=72.8681;spin_atom=0.0;} // K_sv:PAW_GGA:04May1998
  if(AUID=="cc26596520b9b525" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.09647;volume_atom=73.1151;spin_atom=0.0;} // K_sv:PAW_PBE:06Sep2000
  if(AUID=="85c370c6dd5f769e" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.15723;volume_atom=63.4202;spin_atom=0.0;} // K_sv:PAW_LDA:24Mar1998
  if(AUID=="ef04627355867a05" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.29912;volume_atom=71.3848;spin_atom=0.0;} // K_pv:PAW_PBE_KIN:SCAN:17Jan2003
  if(AUID=="3ac2263c1999fd5c" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.33540;volume_atom=74.0621;spin_atom=0.0;} // K_sv:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="ef04627355867a05" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.02694;volume_atom=72.4597;spin_atom=0.0;} // K_pv:PAW_PBE_KIN:17Jan2003
  if(AUID=="3ac2263c1999fd5c" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.04665;volume_atom=72.9004;spin_atom=0.0;} // K_sv:PAW_PBE_KIN:06Sep2000
  if(AUID=="cd771e7ca0b75f0c" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.161109;volume_atom=63.0932;spin_atom=0.0;} // K_pv:PAW_LDA_KIN:22Nov2003
  if(AUID=="ffbcf5d09569ab63" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-3.508309;volume_atom=72.28369;spin_atom=0.0;} // K_sv_GW:PAW_PBE_KIN:31Mar2010
  if(AUID=="c231ba698148ee95" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.164799;volume_atom=63.68187;spin_atom=0.0;} // K_sv:PAW_LDA_KIN:22Nov2003
  if(AUID=="f21b701c5a994872" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-3.679968;volume_atom=63.10559;spin_atom=0.0;} // K_sv_GW:PAW_LDA_KIN:31Mar2010

  // I
  if(AUID=="89deb5ffd67a7166" && SCAN) {found=TRUE;groundstate_structure="A14";groundstate_energy=-36.69694;volume_atom=42.95589;spin_atom=0.0;} // I:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="7e60ff94db40760b" && nKIN) {found=TRUE;groundstate_structure="A14";groundstate_energy=-1.51735;volume_atom=49.44386;spin_atom=0.0;} // I:PAW_PBE:08Apr2002
  if(AUID=="89deb5ffd67a7166" && nKIN) {found=TRUE;groundstate_structure="A14";groundstate_energy=-1.516921;volume_atom=49.59526;spin_atom=0.0;} // I:PAW_PBE_KIN:08Apr2002

  // In
  if(AUID=="5ce1eee07a5df3a7" && nKIN) {found=TRUE;groundstate_structure="A6";groundstate_energy=-2.72115;volume_atom=27.1064;spin_atom=0.0;} // In_d:PAW_PBE:06Sep2000
  if(AUID=="b26165280c5a6d7c" && SCAN) {found=TRUE;groundstate_structure="A6";groundstate_energy=-33.3068;volume_atom=26.1419;spin_atom=0.0;} // In:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="ff28085cf0ed0c70" && SCAN) {found=TRUE;groundstate_structure="A6";groundstate_energy=-33.3094;volume_atom=25.5954;spin_atom=0.0;} // In_d:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="ff28085cf0ed0c70" && nKIN) {found=TRUE;groundstate_structure="A6";groundstate_energy=-2.55842;volume_atom=27.0463;spin_atom=0.0;} // In_d:PAW_PBE_KIN:06Sep2000

  // Ir
  // find /common/LIB1/LIB/Ir* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="78171971d3315e28" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-10.3309;volume_atom=13.8873;spin_atom=0.0;} // Ir:PAW_LDA:10Feb1998
  if(AUID=="8714595dd67f5cc6" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.72200;volume_atom=14.6380;spin_atom=0.0;} // Ir:GGA:01Apr2000
  if(AUID=="0d6f43382108d13d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-10.2667;volume_atom=13.9078;spin_atom=0.0;} // Ir:LDA:01Apr2000
  if(AUID=="b33a03214533a06a" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.79336;volume_atom=14.5807;spin_atom=0.0;} // Ir:PAW_GGA:04May1998
  if(AUID=="9c41e654d31db16a" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.85711;volume_atom=14.5260;spin_atom=0.0;} // Ir:PAW_PBE:06Sep2000
  if(AUID=="be8b84b10c2543f2" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-71.2587;volume_atom=13.5393;spin_atom=0.0;} // Ir:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="be8b84b10c2543f2" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.84886;volume_atom=14.4792;spin_atom=0.0;} // Ir:PAW_PBE_KIN:06Sep2000
  if(AUID=="f143f80475599f27" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-10.3273;volume_atom=13.89488;spin_atom=0.0;} // Ir:PAW_LDA_KIN:10Feb1998
  if(AUID=="c40fbe0470edb1a3" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-10.89579;volume_atom=13.78158;spin_atom=0.0;} // Ir_sv_GW:PAW_LDA_KIN:23Mar2010
  if(AUID=="016aa529c50df352" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-9.406343;volume_atom=14.36741;spin_atom=0.0;} // Ir_sv_GW:PAW_PBE_KIN:23Mar2010

  // La
  // find /common/LIB1/LIB/La* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="d0888fff3f31a114" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.56862;volume_atom=31.6652;spin_atom=0.0;} // La:PAW_LDA:17Apr2000
  if(AUID=="ec88cfb1842bf27d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.90610;volume_atom=36.2172;spin_atom=0.0;} // La:PAW_GGA:14Apr2000
  if(AUID=="8353806ee117d7fd" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.91782;volume_atom=36.5244;spin_atom=0.0;} // La:PAW_PBE:06Sep2000
  if(AUID=="d04dad60260af779" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.86223;volume_atom=37.0366;spin_atom=0.0;} // La_s:PAW_PBE:06Sep2000
  if(AUID=="b2abe18f66fcdafe" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.86314;volume_atom=36.6893;spin_atom=0.0;} // La_s:PAW_GGA:17Apr2000
  if(AUID=="0549ee555a0b18c8" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.53774;volume_atom=32.3324;spin_atom=0.0;} // La_s:PAW_LDA:17Apr2000
  if(AUID=="450ae928112e37cc" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-43.0405;volume_atom=37.4721;spin_atom=0.0;} // La_s:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="4098c299a18eef08" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-43.1348;volume_atom=36.4040;spin_atom=0.0;} // La:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="4098c299a18eef08" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.87431;volume_atom=36.72;spin_atom=0.0;} // La:PAW_PBE_KIN:06Sep2000
  if(AUID=="47faddba936a625c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.376577;volume_atom=31.48441;spin_atom=0.0;} // La_GW:PAW_LDA_KIN:16May2012
  if(AUID=="ff51171bdf527343" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.699423;volume_atom=36.54313;spin_atom=0.0;} // La_GW:PAW_PBE_KIN:16May2012
  if(AUID=="1bafc710bd53cb57" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.583786;volume_atom=31.67771;spin_atom=0.0;} // La:PAW_LDA_KIN:17Apr2000
  if(AUID=="54eb962a4ef4ea32" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.543443;volume_atom=32.34424;spin_atom=0.0;} // La_s:PAW_LDA_KIN:17Apr2000
  if(AUID=="450ae928112e37cc" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.820638;volume_atom=37.1566;spin_atom=0.0;} // La_s:PAW_PBE_KIN:06Sep2000

  // Li
  // find /common/LIB1/LIB/Li* -name "A2" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="d51fe7240453d7c4" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-1.897390;volume_atom=20.3189;spin_atom=0.0;} // Li:PAW_PBE:17Jan2003
  if(AUID=="67021bad74ddb4c5" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-1.897720;volume_atom=20.3318;spin_atom=0.0;} // Li_sv:PAW_GGA:23Jan2001
  if(AUID=="99a2f30dac0e41d6" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-1.905030;volume_atom=20.2771;spin_atom=0.0;} // Li_sv:PAW_PBE:23Jan2001
  if(AUID=="e3ad72d8dbbd284c" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-2.04188;volume_atom=19.0006;spin_atom=0.0;} // Li:PAW_LDA:21Jan2003
  if(AUID=="1c81ce7af52fff4e" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-2.04908;volume_atom=18.8562;spin_atom=0.0;} // Li:LDA:01Apr2000
  if(AUID=="6076711c64467c4f" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-1.90718;volume_atom=19.7515;spin_atom=0.0;} // Li:GGA:01Apr2000
  if(AUID=="c680d72aeb02882e" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-2.04290;volume_atom=19.0034;spin_atom=0.0;} // Li_sv:PAW_LDA:19Jan2001
  if(AUID=="1c81ce7af52fff4e" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-2.04908;volume_atom=18.8562;spin_atom=0.0;} // Li:LDA:01Apr2000
  if(AUID=="0666c5aca67cb28d" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-1.89268;volume_atom=20.3860;spin_atom=0.0;} // Li:PAW_GGA:21Jan2003
  if(AUID=="6076711c64467c4f" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-1.90718;volume_atom=19.7515;spin_atom=0.0;} // Li:GGA:01Apr2000
  if(AUID=="84026814dcd62d7a" && SCAN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-2.33938;volume_atom=20.8216;spin_atom=0.0;} // Li:PAW_PBE_KIN:SCAN:17Jan2003
  if(AUID=="84026814dcd62d7a" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.89756;volume_atom=20.3191;spin_atom=0.0;} // Li:PAW_PBE_KIN:17Jan2003
  if(AUID=="fbd8e4887832ff54" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.90441;volume_atom=19.9894;spin_atom=0.0;} // Li_sv:PAW_PBE:10Sep2004
  if(AUID=="e3a6007a49bfc154" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-3.795153;volume_atom=18.96645;spin_atom=0.0;} // Li_AE_GW:PAW_LDA:25Mar2010
  if(AUID=="66ecb6ec79d94c89" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-3.689522;volume_atom=20.19741;spin_atom=0.0;} // Li_AE_GW:PAW_PBE:25Mar2010
  if(AUID=="78d965ff372a7e11" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-2.041035;volume_atom=18.95202;spin_atom=0.0;} // Li_GW:PAW_LDA_KIN:11May2007
  if(AUID=="98109ece1c9fcc5d" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.894638;volume_atom=20.32433;spin_atom=0.0;} // Li_GW:PAW_PBE_KIN:11May2007
  if(AUID=="157d4a0ea0712143" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-3.795337;volume_atom=18.9764;spin_atom=0.0;} // Li_sv_GW:PAW_LDA:25Mar2010
  if(AUID=="15ba50edd26f214b" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-3.689587;volume_atom=20.20269;spin_atom=0.0;} // Li_sv_GW:PAW_PBE:25Mar2010
  if(AUID=="2d76bbafaaad041b" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-2.044447;volume_atom=18.75312;spin_atom=0.0;} // Li_sv:PAW_LDA:10Sep2004
  if(AUID=="faf632cf08a133d7" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-2.04007;volume_atom=19.03729;spin_atom=0.0;} // Li:PAW_LDA_KIN:15Mar2012

  // Mg
  if(AUID=="285c78eec2b06a80" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.51618;volume_atom=22.8934;spin_atom=0.0;} // Mg:GGA:01Apr2000
  if(AUID=="d8bb30571ef22203" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.52082;volume_atom=22.7555;spin_atom=0.0;} // Mg:GGA:01Apr2000
  if(AUID=="78efd19d72d98cb8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.78695;volume_atom=21.4123;spin_atom=0.0;} // Mg:LDA:01Apr2000
  if(AUID=="61fec4224ab4cab8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.78165;volume_atom=21.5278;spin_atom=0.0;} // Mg:LDA:01Apr2000
  if(AUID=="de98b20ad5c678cf" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.52338;volume_atom=22.8721;spin_atom=0.0;} // Mg:PAW_GGA:05Jan2001
  if(AUID=="7d88d84366cb7141" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.54148;volume_atom=22.8463;spin_atom=0.0;} // Mg:PAW_PBE:05Jan2001
  if(AUID=="b3d0cc0187ea8112" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.78256;volume_atom=21.5714;spin_atom=0.0;} // Mg:PAW_LDA:02Mar1998
  if(AUID=="7a5a79e9587d1be6" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.48173;volume_atom=22.8793;spin_atom=0.0;} // Mg:GGA:01Apr2000
  if(AUID=="1308594682e50447" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.78486;volume_atom=21.5561;spin_atom=0.0;} // Mg:LDA:01Apr2000
  if(AUID=="3e32a29daafd48c5" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.47720;volume_atom=22.7969;spin_atom=0.0;} // Mg_pv:PAW_GGA:10Feb1998
  if(AUID=="c2b2c35b8213a287" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.75631;volume_atom=21.4740;spin_atom=0.0;} // Mg_pv:PAW_LDA:02Mar1998
  if(AUID=="393ee52ce57367b8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.59348;volume_atom=22.7983;spin_atom=0.0;} // Mg_pv:PAW_PBE:06Sep2000
  if(AUID=="2c5dfb19cf1554e5" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.96828;volume_atom=22.6560;spin_atom=0.0;} // Mg:PAW_PBE_KIN:SCAN:13Apr2007
  if(AUID=="26735ac174693885" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.98465;volume_atom=22.7620;spin_atom=0.0;} // Mg_pv:PAW_PBE_KIN:SCAN:13Apr2007
  if(AUID=="2c5dfb19cf1554e5" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.5056;volume_atom=22.8575;spin_atom=0.0;} // Mg:PAW_PBE_KIN:13Apr2007
  if(AUID=="26735ac174693885" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.504;volume_atom=22.672;spin_atom=0.0;} // Mg_pv:PAW_PBE_KIN:13Apr2007
  if(AUID=="10cfb435085a0f1b" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.777176;volume_atom=21.57801;spin_atom=0.0;} // Mg_GW:PAW_LDA_KIN:13Apr2007
  if(AUID=="de8b5248bd86b871" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.502108;volume_atom=22.88404;spin_atom=0.0;} // Mg_GW:PAW_PBE_KIN:13Apr2007
  if(AUID=="590af3dee2a4fe48" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.809508;volume_atom=21.42695;spin_atom=0.0;} // Mg:PAW_LDA_KIN:13Apr2007
  if(AUID=="013b0f18810c46a4" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.767867;volume_atom=21.39455;spin_atom=0.0;} // Mg_pv:PAW_LDA_KIN:13Apr2007
  if(AUID=="149749856a0edabd" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.499509;volume_atom=22.63839;spin_atom=0.0;} // Mg_pv_GW:PAW_PBE_KIN:20Apr2010
  if(AUID=="069806e2f2e95d29" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.763302;volume_atom=21.36101;spin_atom=0.0;} // Mg_pv_GW:PAW_LDA_KIN:20Apr2010
  if(AUID=="e17bc5af0b84ca63" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-16.05968;volume_atom=21.27808;spin_atom=0.0;} // Mg_sv_GW:PAW_LDA:20Apr2010
  if(AUID=="9d41ca07323da062" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.519768;volume_atom=22.50192;spin_atom=0.0;} // Mg_sv:PAW_PBE:12Apr2007
  if(AUID=="ece8e1f37bb686ef" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.776121;volume_atom=21.29943;spin_atom=0.0;} // Mg_sv:PAW_LDA:15Mar2012
  if(AUID=="d42c93affd861a85" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-15.65218;volume_atom=22.59257;spin_atom=0.0;} // Mg_sv_GW:PAW_PBE_KIN:20Apr2010

  // Mn
  if(AUID=="8f2e213b3592f62e" && nKIN) {found=TRUE;groundstate_structure="A12";groundstate_energy=-9.02786;volume_atom=10.7275;spin_atom=0.0;} // Mn:PAW_PBE:06Sep2000
  if(AUID=="99be850476e2dfb3" && nKIN) {found=TRUE;groundstate_structure="A12";groundstate_energy=-9.15350;volume_atom=10.6973;spin_atom=0.0;} // Mn_pv:PAW_PBE:07Sep2000
  if(AUID=="47732958260bfe13" && nKIN) {found=TRUE;groundstate_structure="A12";groundstate_energy=-8.98011;volume_atom=10.7161;spin_atom=0.0;} // Mn:PAW_PBE_KIN:06Sep2000
  if(AUID=="392051866cb4fe60" && nKIN) {found=TRUE;groundstate_structure="A12";groundstate_energy=-8.99108;volume_atom=10.7527;spin_atom=0.0;} // Mn_pv:PAW_PBE_KIN:02Aug2007

  // Mo
  if(AUID=="8f74857b7edbc55b" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.8344;volume_atom=15.6219;spin_atom=0.0;} // Mo:GGA:01Apr2000
  if(AUID=="4690ce5c08300f84" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.9104;volume_atom=15.6337;spin_atom=0.0;} // Mo:PAW_GGA:08Jan2002
  if(AUID=="c56388dcfe8555c0" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.1515;volume_atom=14.8749;spin_atom=0.0;} // Mo:LDA:01Apr2000
  if(AUID=="c90275c91471b466" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.9465;volume_atom=15.5844;spin_atom=0.0;} // Mo:PAW_PBE:08Apr2002
  if(AUID=="3e6c7aec9ade11fe" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.5202;volume_atom=16.1255;spin_atom=0.0;} // Mo:GGA:01Apr2000
  if(AUID=="097174f08d663402" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.7663;volume_atom=15.3630;spin_atom=0.0;} // Mo:LDA:01Apr2000
  if(AUID=="173298095bee51b0" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.8042;volume_atom=15.9061;spin_atom=0.0;} // Mo_pv:PAW_GGA:08Jan2002
  if(AUID=="517b9d282b85b101" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.0734;volume_atom=15.1997;spin_atom=0.0;} // Mo_pv:PAW_LDA:08Jan2002
  if(AUID=="444a5f723c4c0576" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.8439;volume_atom=15.8555;spin_atom=0.0;} // Mo_pv:PAW_PBE:08Apr2002
  if(AUID=="027c5a9d41c6523f" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.1895;volume_atom=15.0334;spin_atom=0.0;} // Mo_sv:PAW_LDA:15Nov2001
  if(AUID=="5001e3854b3663c0" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-34.8382;volume_atom=15.5390;spin_atom=0.0;} // Mo_pv:PAW_PBE_KIN:SCAN:04Feb2005
  if(AUID=="0808a3790788f802" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-34.7737;volume_atom=15.3662;spin_atom=0.0;} // Mo:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="eb2ad6b5ae2c88f9" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-34.8748;volume_atom=15.5342;spin_atom=0.0;} // Mo_sv:PAW_PBE_KIN:SCAN:02Feb2006
  if(AUID=="0808a3790788f802" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.9454;volume_atom=15.5866;spin_atom=0.0;} // Mo:PAW_PBE_KIN:08Apr2002
  if(AUID=="5001e3854b3663c0" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.9234;volume_atom=15.7371;spin_atom=0.0;} // Mo_pv:PAW_PBE_KIN:04Feb2005
  if(AUID=="eb2ad6b5ae2c88f9" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.9308;volume_atom=15.7292;spin_atom=0.0;} // Mo_sv:PAW_PBE_KIN:02Feb2006
  if(AUID=="25f5fda9fc8c9cff" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.21278;volume_atom=14.94573;spin_atom=0.0;} // Mo:PAW_LDA_KIN:08Apr2002
  if(AUID=="bd15acc21b11b902" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.16272;volume_atom=15.06315;spin_atom=0.0;} // Mo_pv:PAW_LDA_KIN:04Feb2005
  if(AUID=="9e3e9147c43cf059" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.71705;volume_atom=15.69242;spin_atom=0.0;} // Mo_sv_GW:PAW_PBE_KIN:05Dec2013
  if(AUID=="fc1f92a69942b43d" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.95249;volume_atom=15.01251;spin_atom=0.0;} // Mo_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="0f8c2e9ac364f9e3" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.13006;volume_atom=15.06347;spin_atom=0.0;} // Mo_sv:PAW_LDA_KIN:29Jan2005

  // N
  if(AUID=="3cd4003f27559a49" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-8.318698;volume_atom=500;spin_atom=0.0;} // N:PAW_PBE:08Apr2002
  if(AUID=="be869d8257082ad4" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-8.324198;volume_atom=500;spin_atom=0.0;} // N_h:PAW_PBE:13Feb2001
  if(AUID=="4ce4637079c90d89" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-8.319179;volume_atom=500;spin_atom=0.0;} // N:PAW_PBE_KIN:08Apr2002
  if(AUID=="4ce4637079c90d89" && SCAN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-9.252821;volume_atom=500;spin_atom=0.0;} // N:PAW_PBE_KIN:SCAN:08Apr2002

  // Na
  if(AUID=="bd75e1ab544c0ed0" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.31096;volume_atom=36.2541;spin_atom=0.0;} // Na_pv:PAW_PBE:05Jan2001
  if(AUID=="192d4f6f863806ad" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.31443;volume_atom=35.0106;spin_atom=0.0;} // Na_sv:PAW_GGA:28Sep2000
  if(AUID=="637322e42c89bd5e" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.31537;volume_atom=35.6060;spin_atom=0.0;} // Na_sv:PAW_PBE:28Sep2000
  if(AUID=="2215677fd2f98aea" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.45287;volume_atom=33.2473;spin_atom=0.0;} // Na:PAW_LDA:24Mar1998
  if(AUID=="e7d71dba59886669" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.30398;volume_atom=37.1851;spin_atom=0.0;} // Na:PAW_GGA:05Jan2001
  if(AUID=="79f12c5d16cf98b8" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.30640;volume_atom=36.7352;spin_atom=0.0;} // Na:PAW_PBE:08Apr2002
  if(AUID=="f44fe5066a3bc573" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.30382;volume_atom=36.4114;spin_atom=0.0;} // Na_pv:PAW_GGA:05Jan2001
  if(AUID=="189d870b4537263b" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.45445;volume_atom=31.4092;spin_atom=0.0;} // Na_sv:PAW_LDA:28Sep2000
  if(AUID=="b994b2842f3b80ab" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.45431;volume_atom=32.6970;spin_atom=0.0;} // Na_pv:PAW_LDA:28Sep2000
  if(AUID=="86f3cf1cb9c67ee6" && SCAN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-4.20653;volume_atom=37.2429;spin_atom=0.0;} // Na_pv:PAW_PBE_KIN:SCAN:19Sep2006
  if(AUID=="63e83f9cf84ac969" && SCAN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-4.18698;volume_atom=36.5426;spin_atom=0.0;} // Na:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="86f3cf1cb9c67ee6" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.311;volume_atom=36.1448;spin_atom=0.0;} // Na_pv:PAW_PBE_KIN:19Sep2006
  if(AUID=="e720028a65d0fbc7" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.452684;volume_atom=33.25262;spin_atom=0.0;} // Na:PAW_LDA_KIN:24Mar1998
  if(AUID=="63e83f9cf84ac969" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.305819;volume_atom=36.73294;spin_atom=0.0;} // Na:PAW_PBE_KIN:08Apr2002
  if(AUID=="0e2d9945823de573" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.452836;volume_atom=32.52883;spin_atom=0.0;} // Na_pv:PAW_LDA_KIN:19Dec2003
  if(AUID=="6e7e25a6de2a42c0" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-3.564566;volume_atom=31.81405;spin_atom=0.0;} // Na_sv_GW:PAW_LDA_KIN:11May2015
  if(AUID=="60a75938718c0f9d" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.453047;volume_atom=31.60811;spin_atom=0.0;} // Na_sv:PAW_LDA:28Sep2000
  if(AUID=="3455418d35853cd4" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-3.390494;volume_atom=35.54009;spin_atom=0.0;} // Na_sv_GW:PAW_PBE_KIN:11May2015
  if(AUID=="1c8fc7a7d87f9aa3" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.313068;volume_atom=35.11102;spin_atom=0.0;} // Na_sv:PAW_PBE:28Sep2000

  // Nb
  if(AUID=="31fbbfcce7b1cf49" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.3250;volume_atom=16.8059;spin_atom=0.0;} // Nb:LDA:01Apr2000
  if(AUID=="724f3dcc974f7f25" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.1730;volume_atom=17.8196;spin_atom=0.0;} // Nb:GGA:01Apr2000
  if(AUID=="e1ae346e0f041fdc" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.99099;volume_atom=18.3604;spin_atom=0.0;} // Nb:GGA:01Apr2000
  if(AUID=="31ac6525497e7fb5" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.0804;volume_atom=17.3804;spin_atom=0.0;} // Nb:LDA:01Apr2000
  if(AUID=="3df5731180246390" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.0550;volume_atom=18.2510;spin_atom=0.0;} // Nb_pv:PAW_GGA:09Jan2002
  if(AUID=="73b4ad5a036338fa" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.1463;volume_atom=17.3220;spin_atom=0.0;} // Nb_pv:PAW_LDA:09Jan2002
  if(AUID=="d9711dcad08e9fe9" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.0830;volume_atom=18.2314;spin_atom=0.0;} // Nb_pv:PAW_PBE:08Apr2002
  if(AUID=="79b1b1f832ef1be0" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.1943;volume_atom=18.0754;spin_atom=0.0;} // Nb_sv:PAW_GGA:14Nov2001
  if(AUID=="141a7bea5919e0b3" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.2253;volume_atom=18.0525;spin_atom=0.0;} // Nb_sv:PAW_PBE:17Jan2003
  if(AUID=="9fdcd7998786bd34" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.2626;volume_atom=17.0887;spin_atom=0.0;} // Nb_sv:PAW_LDA:15Nov2001
  if(AUID=="2d2f42f890bd3354" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-33.0916;volume_atom=18.0829;spin_atom=0.0;} // Nb_pv:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="b383f09c3645c1a5" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-33.2389;volume_atom=17.8122;spin_atom=0.0;} // Nb_sv:PAW_PBE_KIN:SCAN:25May2007
  if(AUID=="b383f09c3645c1a5" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.2147;volume_atom=18.0366;spin_atom=0.0;} // Nb_sv:PAW_PBE_KIN:25May2007
  if(AUID=="2d2f42f890bd3354" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.081;volume_atom=18.23947;spin_atom=0.0;} // Nb_pv:PAW_PBE_KIN:08Apr2002
  if(AUID=="125b6a029adf24ea" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.14547;volume_atom=17.3083;spin_atom=0.0;} // Nb_pv:PAW_LDA_KIN:09Jan2002
  if(AUID=="9b8caffbad647df6" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.54389;volume_atom=17.00623;spin_atom=0.0;} // Nb_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="7a3805e4c05d950e" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.47372;volume_atom=17.9751;spin_atom=0.0;} // Nb_sv_GW:PAW_PBE_KIN:05Dec2013
  if(AUID=="7b5ec33ec65ced5a" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.25399;volume_atom=17.05046;spin_atom=0.0;} // Nb_sv:PAW_LDA_KIN:25May2007

  // Ne
  // find /common/LIB1/LIB/Ne* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="bf15cfac49144396" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.048593;volume_atom=13.7731;spin_atom=0.0;} // Ne:PAW_LDA:07Sep2000
  if(AUID=="78bf25649465a255" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.103775;volume_atom=13.5480;spin_atom=0.0;} // Ne:LDA:01Apr2000
  if(AUID=="4bda2c3e931df0db" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.102230;volume_atom=13.5480;spin_atom=0.0;} // Ne:LDA:01Apr2000
  if(AUID=="445c00adf152283d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.066215;volume_atom=20.2857;spin_atom=0.0;} // Ne:PAW_GGA:05Jan2001
  if(AUID=="91b32db4af57705e" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.032434;volume_atom=20.1045;spin_atom=0.0;} // Ne:PAW_PBE:05Jan2001
  if(AUID=="326c9a1ee646c679" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.394710;volume_atom=15.4322;spin_atom=0.0;} // Ne:PAW_PBE_KIN:SCAN:05Jan2001
  if(AUID=="1db85f57d622abef" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.04866122;volume_atom=13.79651;spin_atom=0.0;} // Ne:PAW_LDA_KIN:07Sep2000
  if(AUID=="326c9a1ee646c679" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.03233834;volume_atom=20.8374;spin_atom=0.0344651;} // Ne:PAW_PBE_KIN:05Jan2001
  if(AUID=="3444109214f68eab" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.04640846;volume_atom=13.06904;spin_atom=0.0;} // Ne_GW:PAW_LDA_KIN:22Sep2014
  if(AUID=="f46c79223af4d753" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.03232092;volume_atom=20.4637;spin_atom=0.0330588;} // Ne_GW:PAW_PBE_KIN:02Oct2006
  if(AUID=="42320e2cde3e888a" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.02750396;volume_atom=22.17472;spin_atom=0.0389689;} // Ne_GW:PAW_PBE_KIN:22Sep2014
  if(AUID=="687e59556bdcdb00" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.04245677;volume_atom=14.02031;spin_atom=0.0;} // Ne_s_GW:PAW_LDA_KIN:02Oct2006

  // Ni
  // find /common/LIB1/LIB/Ni* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="1141cd9c7118591b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.48921;volume_atom=10.76129;spin_atom=0.6315579;} // Ni_pv:PAW_PBE_KIN:06Sep2000
  if(AUID=="1141cd9c7118591b" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-16.13729;volume_atom=10.27715;spin_atom=0.7384191;} // Ni_pv:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="134ab84006e2f008" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.586101;volume_atom=10.11318;spin_atom=0.5863435;} // Ni:LDA:01Apr2000
  if(AUID=="1dc196d6de4f6b37" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.571077;volume_atom=10.9016;spin_atom=0.6273349;} // Ni:PAW_PBE:06Sep2000
  if(AUID=="1f50d23b4fbbe9d3" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.46999;volume_atom=10.8419;spin_atom=0.6269131;} // Ni:PAW_PBE_KIN:02Aug2007
  if(AUID=="1f50d23b4fbbe9d3" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-16.06428;volume_atom=10.35633;spin_atom=0.6963174;} // Ni:PAW_PBE_KIN:SCAN:02Aug2007
  if(AUID=="3d3b495397078a92" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.59795;volume_atom=9.959794;spin_atom=0.5959346;} // Ni_pv:PAW_LDA_KIN:19Apr2000
  if(AUID=="510ec32ec7300763" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.406212;volume_atom=10.64451;spin_atom=0.6311128;} // Ni_sv_GW:PAW_PBE_KIN:05Dec2013
  if(AUID=="55fd1c17c8c96b19" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.391172;volume_atom=10.76669;spin_atom=0.6146924;} // Ni_pv:PAW_GGA:19Apr2000
  if(AUID=="5737f6f3ba793dbe" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.567688;volume_atom=10.02489;spin_atom=0.5781327;} // Ni:PAW_LDA:03Mar1998
  if(AUID=="5fe58cf81946c37d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.582183;volume_atom=9.995466;spin_atom=0.5849204;} // Ni:PAW_LDA_KIN:02Aug2007
  if(AUID=="bee457b2befa4114" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.501635;volume_atom=9.950662;spin_atom=0.5836186;} // Ni_pv:PAW_LDA:19Apr2000
  if(AUID=="bf5d7aa5efb4320a" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.444359;volume_atom=9.878276;spin_atom=0.5938513;} // Ni_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="ead949d4e6fa93d1" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.777832;volume_atom=10.78229;spin_atom=0.6338126;} // Ni_pv:PAW_PBE:06Sep2000
  if(AUID=="eafb70a88aeb4012" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.476548;volume_atom=10.99426;spin_atom=0.6368326;} // Ni:GGA:01Apr2000
  if(AUID=="f316c1e20cc629f5" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.461334;volume_atom=10.87603;spin_atom=0.6086769;} // Ni:PAW_GGA:03Mar1998
  if(AUID=="f57c7aee2773cf5c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.599106;volume_atom=9.979222;spin_atom=0.5910477;} // Ni_GW:PAW_LDA_KIN:31Mar2010
  if(AUID=="fd7da3695dd95380" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.482908;volume_atom=10.80456;spin_atom=0.6293194;} // Ni_GW:PAW_PBE_KIN:31Mar2010

  // O
  if(AUID=="85ded82734544fa9" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-4.931145;volume_atom=500;spin_atom=1.0;} // O:PAW_PBE:08Apr2002
  if(AUID=="5c530b0eaed90a00" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-4.931071;volume_atom=500;spin_atom=1.0;} // O:PAW_PBE_KIN:08Apr2002
  if(AUID=="5c530b0eaed90a00" && SCAN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-6.022492;volume_atom=500;spin_atom=1.0;} // O:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="847e34315251cef5" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-4.701027;volume_atom=500;spin_atom=1.0;} // O_s:PAW_PBE:07Sep2000
  if(AUID=="60ee78e0c6ceac3d" && nKIN) {found=TRUE;groundstate_structure="diatom";groundstate_energy=-5.017967;volume_atom=500;spin_atom=1.0;} // O_h:PAW_PBE:20Dec2001

  // Os
  if(AUID=="4cdd79622cbed39c" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.1466;volume_atom=14.3178;spin_atom=0.0;} // Os:PAW_GGA:06Feb2003
  if(AUID=="458230748f9ee8c4" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.1170;volume_atom=14.3175;spin_atom=0.0;} // Os:GGA:01Apr2000
  if(AUID=="8c0f432aa073c44a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.6831;volume_atom=13.6759;spin_atom=0.0;} // Os:LDA:01Apr2000
  if(AUID=="b1954be3210262f1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.2440;volume_atom=14.2137;spin_atom=0.0;} // Os:PAW_PBE:17Jan2003
  if(AUID=="d38b52b3db36dc96" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.7031;volume_atom=13.7273;spin_atom=0.0;} // Os:PAW_LDA:21Jan2003
  if(AUID=="216ecda62b9b32bf" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.2193;volume_atom=14.2979;spin_atom=0.0;} // Os_pv:PAW_PBE:20Jan2003
  if(AUID=="d2665c78be0bcda4" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.1086;volume_atom=14.3251;spin_atom=0.0;} // Os_pv:PAW_GGA:10Feb1998
  if(AUID=="6948d5b374665a4d" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.6613;volume_atom=13.7168;spin_atom=0.0;} // Os_pv:PAW_LDA:10Feb1998
  if(AUID=="0b4b416ff63ebb16" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-71.8359;volume_atom=13.6516;spin_atom=0.0;} // Os:PAW_PBE_KIN:SCAN:17Jan2003
  if(AUID=="6fba18d27313ce5e" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-71.7960;volume_atom=13.7403;spin_atom=0.0;} // Os_pv:PAW_PBE_KIN:SCAN:20Jan2003
  if(AUID=="6fba18d27313ce5e" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.2224;volume_atom=14.2969;spin_atom=0.0;} // Os_pv:PAW_PBE_KIN:20Jan2003
  if(AUID=="3d2f9f5957efb0e2" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.63657;volume_atom=13.65506;spin_atom=0.0;} // Os:PAW_LDA_KIN:21Jan2003
  if(AUID=="0b4b416ff63ebb16" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.24261;volume_atom=14.21246;spin_atom=0.0;} // Os:PAW_PBE_KIN:17Jan2003
  if(AUID=="9ce762cd7cd9b03f" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.66875;volume_atom=13.79887;spin_atom=0.0;} // Os_pv:PAW_LDA_KIN:28Mar2012
  if(AUID=="e444b49bdb1f5f2b" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-13.80236;volume_atom=13.68508;spin_atom=0.0;} // Os_sv_GW:PAW_LDA_KIN:23Mar2010
  if(AUID=="63300b11641f5a42" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.35418;volume_atom=14.19354;spin_atom=0.0;} // Os_sv_GW:PAW_PBE_KIN:23Mar2010

  // P
  if(AUID=="0df9932192e3546b" && nKIN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-5.32414;volume_atom=16.0927;spin_atom=0.0;} // P:PAW_PBE:17Jan2003
  if(AUID=="d64d49a8028ca2bc" && SCAN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-10.3005;volume_atom=16.0865;spin_atom=0.0;} // P_h:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="3cafdbeb02f9e951" && SCAN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-10.3006;volume_atom=16.0927;spin_atom=0.0;} // P:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="3cafdbeb02f9e951" && nKIN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-5.32303;volume_atom=16.0811;spin_atom=0.0;} // P:PAW_PBE_KIN:06Sep2000

  // Pb
  // find /common/LIB1/LIB/Pb* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="1c02af03a0cb18b2" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.59467;volume_atom=31.8712;spin_atom=0.0;} // Pb_d:GGA:01Apr2000
  if(AUID=="58341640c6c638ed" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.26624;volume_atom=29.0840;spin_atom=0.0;} // Pb_d:LDA:01Apr2000
  if(AUID=="7d26482cab914656" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.54563;volume_atom=31.4146;spin_atom=0.0;} // Pb_d:PAW_GGA:04May1998
  if(AUID=="0e47cf73a95cbfc2" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.70470;volume_atom=31.5298;spin_atom=0.0;} // Pb_d:PAW_PBE:06Sep2000
  if(AUID=="1010cb02b8b71033" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.20301;volume_atom=28.7142;spin_atom=0.0;} // Pb_d:PAW_LDA:30Apr1998
  if(AUID=="2f64d0d7f1c21e2c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.57056;volume_atom=31.7950;spin_atom=0.0;} // Pb:PAW_PBE:08Apr2002
  if(AUID=="899670f2ebab4069" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.58825;volume_atom=31.8825;spin_atom=0.0;} // Pb:GGA:01Apr2000
  if(AUID=="9568f688bb71b2b0" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.60201;volume_atom=31.7526;spin_atom=0.0;} // Pb:PAW_GGA:25Jul2001
  if(AUID=="626de6bdde387007" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.25933;volume_atom=29.0453;spin_atom=0.0;} // Pb:LDA:01Apr2000
  if(AUID=="381b125c2959cfe6" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.28081;volume_atom=28.9939;spin_atom=0.0;} // Pb:PAW_LDA:02Sep2001
  if(AUID=="727d5248f9fb2e47" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-75.7142;volume_atom=30.6165;spin_atom=0.0;} // Pb:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="3143a75eeaa29ea5" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-75.7300;volume_atom=30.1590;spin_atom=0.0;} // Pb_d:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="3143a75eeaa29ea5" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.56172;volume_atom=31.4742;spin_atom=0.0;} // Pb_d:PAW_PBE_KIN:06Sep2000
  if(AUID=="fff14192b2b4f176" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.554114;volume_atom=31.3173;spin_atom=0.0;} // Pb_d_GW:PAW_PBE_KIN:14Apr2014
  if(AUID=="33aab1fdb40d3950" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.28143;volume_atom=28.96858;spin_atom=0.0;} // Pb:PAW_LDA_KIN:02Sep2001
  if(AUID=="a50c38205409299f" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.252156;volume_atom=28.72551;spin_atom=0.0;} // Pb_d:PAW_LDA_KIN:30Apr1998
  if(AUID=="e7d3c828acba821d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.250507;volume_atom=28.62015;spin_atom=0.0;} // Pb_d_GW:PAW_LDA_KIN:14Apr2014
  if(AUID=="727d5248f9fb2e47" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.570961;volume_atom=31.81152;spin_atom=0.0;} // Pb:PAW_PBE_KIN:08Apr2002
  if(AUID=="3bda0fd1d9069479" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.578003;volume_atom=30.50902;spin_atom=0.0;} // Pb_sv_GW:PAW_PBE_KIN:04Apr2014
  if(AUID=="10185f68827f501a" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.261472;volume_atom=27.99589;spin_atom=0.0;} // Pb_sv_GW:PAW_LDA_KIN:04Apr2014

  // Pd
  // find /common/LIB1/LIB/Pd* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="2ef0a2bc50307dc4" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.216941;volume_atom=15.22565;spin_atom=0.2877056;} // Pd:PAW_PBE_KIN:04Jan2005
  if(AUID=="2ef0a2bc50307dc4" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-32.81275;volume_atom=14.75361;spin_atom=0.4598601;} // Pd:PAW_PBE_KIN:SCAN:04Jan2005
  if(AUID=="310e1fd8b9174c88" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.495885;volume_atom=14.0938;spin_atom=0.0;} // Pd:PAW_LDA_KIN:04Jan2005
  if(AUID=="35af3e709779fcbb" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.477963;volume_atom=14.05396;spin_atom=0.0;} // Pd_pv:PAW_LDA_KIN:28Jan2005
  if(AUID=="35cb7e1d80725a17" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.204;volume_atom=15.15592;spin_atom=0.2780651;} // Pd_pv:PAW_PBE_KIN:28Jan2005
  if(AUID=="35cb7e1d80725a17" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-32.77616;volume_atom=14.64568;spin_atom=0.4412058;} // Pd_pv:PAW_PBE_KIN:SCAN:28Jan2005
  if(AUID=="3b6d075d41b9645f" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.440518;volume_atom=14.23686;spin_atom=0.0;} // Pd:PAW_LDA:09Oct1998
  if(AUID=="4158852ee0320965" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.359429;volume_atom=14.25476;spin_atom=0.0;} // Pd_pv:PAW_LDA:17Apr2000
  if(AUID=="46106d1a706186be" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.138558;volume_atom=15.49124;spin_atom=0.0;} // Pd_pv:PAW_GGA:04Mar1998
  if(AUID=="4af3f66c14a34cee" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.381966;volume_atom=15.39383;spin_atom=0.0;} // Pd_pv:PAW_PBE:06Sep2000
  if(AUID=="552d7d192268351d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.428931;volume_atom=14.30396;spin_atom=0.0;} // Pd:LDA:01Apr2000
  if(AUID=="6ade515555bfea29" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.178475;volume_atom=15.34327;spin_atom=0.0;} // Pd:PAW_PBE:05Jan2001
  if(AUID=="78fd29f76da976e8" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.759958;volume_atom=15.025;spin_atom=0.2621222;} // Pd_sv_GW:PAW_PBE_KIN:05Dec2013
  if(AUID=="7e8d64cc173f9f0c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.212492;volume_atom=15.38239;spin_atom=0.0;} // Pd:PAW_GGA:05Jan2001
  if(AUID=="854a0d2f9191bbb2" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.198015;volume_atom=15.4798;spin_atom=0.0;} // Pd:GGA:01Apr2000
  if(AUID=="863baeb4a12a2226" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.238536;volume_atom=15.15262;spin_atom=0.2751329;} // Pd_GW:PAW_PBE_KIN:06Mar2008
  if(AUID=="d97ddcc6e311c95c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.524692;volume_atom=14.021;spin_atom=0.0;} // Pd_GW:PAW_LDA_KIN:06Mar2008
  if(AUID=="f8d2ab794f3a7628" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.030837;volume_atom=13.95344;spin_atom=0.0;} // Pd_sv_GW:PAW_LDA_KIN:05Dec2013

  // Pt
  // find /common/LIB1/LIB/Pt* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="4934ee186867e17b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.01467;volume_atom=15.8490;spin_atom=0.0;} // Pt:GGA:01Apr2000
  if(AUID=="59081d68fabd053c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.05482;volume_atom=15.6527;spin_atom=0.0;} // Pt:PAW_PBE:05Jan2001
  if(AUID=="129e32c460592b3d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.42489;volume_atom=14.9080;spin_atom=0.0;} // Pt:LDA:01Apr2000
  if(AUID=="2f9c34be5eee2bc7" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.45079;volume_atom=14.8389;spin_atom=0.0;} // Pt:PAW_LDA:17Apr2000
  if(AUID=="05e3fa0ef96c01d7" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.04404;volume_atom=15.7566;spin_atom=0.0;} // Pt:PAW_GGA:05Jan2001
  if(AUID=="c67ecad4e9259793" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-70.2902;volume_atom=14.7627;spin_atom=0.0;} // Pt:PAW_PBE_KIN:SCAN:04Feb2005
  if(AUID=="35f9da8fbc175cc5" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-70.1284;volume_atom=14.8589;spin_atom=0.0;} // Pt_pv:PAW_PBE_KIN:SCAN:12Dec2005
  if(AUID=="c67ecad4e9259793" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.09746;volume_atom=15.5516;spin_atom=0.0;} // Pt:PAW_PBE_KIN:04Feb2005
  if(AUID=="35f9da8fbc175cc5" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.08134;volume_atom=15.4851;spin_atom=0.0;} // Pt_pv:PAW_PBE_KIN:12Dec2005
  if(AUID=="901eab88fea883a0" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.529267;volume_atom=14.71255;spin_atom=0.0;} // Pt_f_GW:PAW_LDA_KIN:10Mar2009
  if(AUID=="fc52fca180d8c423" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.11373;volume_atom=15.52079;spin_atom=0.0;} // Pt_GW:PAW_PBE_KIN:10Mar2009
  if(AUID=="693825ba67350fc5" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.495628;volume_atom=14.6916;spin_atom=0.0;} // Pt_pv:PAW_LDA_KIN:12Dec2005
  if(AUID=="15e78a872bbd5702" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.302459;volume_atom=14.59476;spin_atom=0.0;} // Pt_sv_GW:PAW_LDA_KIN:23Mar2010
  if(AUID=="cffa53570eebb464" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.054165;volume_atom=15.43248;spin_atom=0.0;} // Pt_sv_GW:PAW_PBE_KIN:23Mar2010
  if(AUID=="d5cf0a3860dcc78c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.507509;volume_atom=14.74523;spin_atom=0.0;} // Pt_new:PAW_LDA_KIN:15Nov2007

  // Rb
  if(AUID=="4a3726ab01a44892" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-0.962354;volume_atom=90.2922;spin_atom=0.0;} // Rb_sv:PAW_PBE:06Sep2000

  // Re
  // find /common/LIB1/LIB/Re* -name "A3" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="1e59564cdfe6dd0a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.34746;volume_atom=14.96078;spin_atom=0.0;} // Re:PAW_GGA:05Jan2001
  if(AUID=="28be1918e906c6ea" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-14.04642;volume_atom=14.86271;spin_atom=0.0;} // Re_sv_GW:PAW_PBE_KIN:23Mar2010
  if(AUID=="2e6c1204bb1f83d9" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-13.75744;volume_atom=14.33156;spin_atom=0.0;} // Re:PAW_LDA_KIN:21Jan2003
  if(AUID=="330065ce03f40a48" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-13.71319;volume_atom=14.41776;spin_atom=0.0;} // Re_pv:PAW_LDA_KIN:11Feb1998
  if(AUID=="3f9c2349a9115560" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.38064;volume_atom=14.96051;spin_atom=0.0;} // Re_pv:PAW_PBE_KIN:06Sep2000
  if(AUID=="3f9c2349a9115560" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-71.03393;volume_atom=14.34273;spin_atom=0.0;} // Re_pv:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="598a3a038e3785b3" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-13.65369;volume_atom=14.41401;spin_atom=0.0;} // Re_pv:PAW_LDA:11Feb1998
  if(AUID=="6c3f2ec3842638b8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.41343;volume_atom=14.8627;spin_atom=0.0;} // Re:PAW_PBE_KIN:17Jan2003
  if(AUID=="6c3f2ec3842638b8" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-71.09379;volume_atom=14.28247;spin_atom=0.0;} // Re:PAW_PBE_KIN:SCAN:17Jan2003
  if(AUID=="7289c4501ec2c23a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-13.7534;volume_atom=14.33325;spin_atom=0.0;} // Re:PAW_LDA:21Jan2003
  if(AUID=="754454a117ea525b" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.2208;volume_atom=15.03947;spin_atom=0.0;} // Re_pv:PAW_GGA:11Feb1998
  if(AUID=="7c2099bc993d23a3" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.30544;volume_atom=14.92092;spin_atom=0.0;} // Re:GGA:01Apr2000
  if(AUID=="98510e8bfa14709e" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-13.77736;volume_atom=14.25279;spin_atom=0.0;} // Re:LDA:01Apr2000
  if(AUID=="d72276b27490b853" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.43246;volume_atom=14.99852;spin_atom=0.0;} // Re_pv:PAW_PBE:06Sep2000
  if(AUID=="df4c77c476ac2743" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-15.37206;volume_atom=14.30951;spin_atom=0.0;} // Re_sv_GW:PAW_LDA_KIN:23Mar2010
  if(AUID=="e5efc178b4648996" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.41127;volume_atom=14.86526;spin_atom=0.0;} // Re:PAW_PBE:17Jan2003

  // Rh
  // find /common/LIB1/LIB/Rh* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="e9aa6eb55d4587f1" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.13631;volume_atom=14.1993;spin_atom=0.0;} // Rh_pv:PAW_GGA:17Apr2000
  if(AUID=="c6d7ff86b44b5d4c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.52654;volume_atom=13.3693;spin_atom=0.0;} // Rh_pv:PAW_LDA:17Apr2000
  if(AUID=="b1547f3221c0c660" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.34058;volume_atom=14.1854;spin_atom=0.0;} // Rh_pv:PAW_PBE:06Sep2000
  if(AUID=="883dfe9dc76b6f31" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.58863;volume_atom=13.3502;spin_atom=0.0;} // Rh:LDA:01Apr2000
  if(AUID=="b02edafce0da47bc" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.17866;volume_atom=14.2151;spin_atom=0.0;} // Rh:GGA:01Apr2000
  if(AUID=="e06760ebe9e0e518" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.27044;volume_atom=14.1233;spin_atom=0.0;} // Rh:PAW_PBE:06Sep2000
  if(AUID=="2297d0177b7db39d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.22290;volume_atom=14.1495;spin_atom=0.0;} // Rh:PAW_GGA:04May1998
  if(AUID=="b03cdc726a115d0b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.62793;volume_atom=13.3185;spin_atom=0.0;} // Rh:PAW_LDA:03Mar1998
  if(AUID=="dc3fbc9801231714" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-33.8609;volume_atom=13.4868;spin_atom=0.0;} // Rh:PAW_PBE_KIN:SCAN:04Feb2005
  if(AUID=="98605c42de12449b" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-33.7489;volume_atom=13.4857;spin_atom=0.0;} // Rh_pv:PAW_PBE_KIN:SCAN:25Jan2005
  if(AUID=="98605c42de12449b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.25357;volume_atom=13.9311;spin_atom=0.0;} // Rh_pv:PAW_PBE_KIN:25Jan2005
  if(AUID=="10b1a1345dd201d9" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.270549;volume_atom=13.93035;spin_atom=0.0;} // Rh_GW:PAW_PBE_KIN:06Mar2008
  if(AUID=="f67899f1f39d8de1" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.707898;volume_atom=13.16461;spin_atom=0.0;} // Rh_GW:PAW_LDA_KIN:06Mar2008
  if(AUID=="8dec99af5c2e1b15" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.711322;volume_atom=13.16041;spin_atom=0.0;} // Rh:PAW_LDA_KIN:04Feb2005
  if(AUID=="dc3fbc9801231714" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.275773;volume_atom=13.92373;spin_atom=0.0;} // Rh:PAW_PBE_KIN:04Feb2005
  if(AUID=="c4a1fff764b83ce4" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.682509;volume_atom=13.17032;spin_atom=0.0;} // Rh_pv:PAW_LDA_KIN:25Jan2005
  if(AUID=="786e475ef59bd566" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.824001;volume_atom=13.11256;spin_atom=0.0;} // Rh_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="84e12fc493054c95" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.405578;volume_atom=13.85661;spin_atom=0.0;} // Rh_sv_GW:PAW_PBE_KIN:05Dec2013

  // Ru
  if(AUID=="e2e4aaaa46da8131" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.20414;volume_atom=13.8077;spin_atom=0.0;} // Ru:PAW_PBE:06Sep2000
  if(AUID=="1fa7a7266ef21653" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.09581;volume_atom=13.8865;spin_atom=0.0;} // Ru:GGA:01Apr2000
  if(AUID=="24117781e96ec111" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.08258;volume_atom=13.9181;spin_atom=0.0;} // Ru_pv:PAW_GGA:10Feb1998
  if(AUID=="d45ec3c5462bef36" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.5801;volume_atom=13.1589;spin_atom=0.0;} // Ru:LDA:01Apr2000
  if(AUID=="48750268bde7c56d" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.16053;volume_atom=13.8456;spin_atom=0.0;} // Ru:PAW_GGA:03Mar1998
  if(AUID=="64d0c8f1c723d0dc" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.6345;volume_atom=13.1584;spin_atom=0.0;} // Ru:PAW_LDA:17Apr2000
  if(AUID=="04267bcea324606c" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.5430;volume_atom=13.2301;spin_atom=0.0;} // Ru_pv:PAW_LDA:03Mar1998
  if(AUID=="a3d039d8d88a4a86" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.27135;volume_atom=13.8821;spin_atom=0.0;} // Ru_pv:PAW_PBE:06Sep2000
  if(AUID=="c6fbedbf10339a1d" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.16320;volume_atom=13.8588;spin_atom=0.0;} // Ru_sv:PAW_GGA:02Oct2001
  if(AUID=="6a824616a7311021" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-34.8327;volume_atom=13.3205;spin_atom=0.0;} // Ru_pv:PAW_PBE_KIN:SCAN:28Jan2005
  if(AUID=="ee46e79f573bbba1" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-34.9124;volume_atom=13.0842;spin_atom=0.0;} // Ru:PAW_PBE_KIN:SCAN:04Feb2005
  if(AUID=="1550a87ae48844e0" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-34.8652;volume_atom=13.3100;spin_atom=0.0;} // Ru_sv:PAW_GGA_KIN:SCAN:28Jan2005
  if(AUID=="6a824616a7311021" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.23339;volume_atom=13.6814;spin_atom=0.0;} // Ru_pv:PAW_PBE_KIN:28Jan2005
  if(AUID=="75e129c7fca73918" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.72705;volume_atom=13.01646;spin_atom=0.0;} // Ru:PAW_LDA_KIN:04Feb2005
  if(AUID=="ee46e79f573bbba1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.256405;volume_atom=13.62536;spin_atom=0.0;} // Ru:PAW_PBE_KIN:04Feb2005
  if(AUID=="d8a5660388ccadd7" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.69806;volume_atom=13.05715;spin_atom=0.0;} // Ru_pv:PAW_LDA_KIN:28Jan2005
  if(AUID=="be15e79ed6e7404c" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.40363;volume_atom=13.01084;spin_atom=0.0;} // Ru_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="c476a1b5dba68ade" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.946162;volume_atom=13.62664;spin_atom=0.0;} // Ru_sv_GW:PAW_PBE_KIN:05Dec2013
  if(AUID=="1be9ee7614ff5941" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.7046;volume_atom=13.06641;spin_atom=0.0;} // Ru_sv:PAW_LDA_KIN:28Jan2005
  if(AUID=="1550a87ae48844e0" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.283335;volume_atom=13.69816;spin_atom=0.0;} // Ru_sv:PAW_GGA_KIN:28Jan2005

  // S
  //  if(AUID=="a18b03f7564daaaa" && SCAN) {found=TRUE;groundstate_structure="A16";groundstate_energy=-9.649852;volume_atom=26.64947;spin_atom=0.0;} // S:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="a18b03f7564daaaa" && nKIN) {found=TRUE;groundstate_structure="A16";groundstate_energy=-4.12508;volume_atom=35.81378;spin_atom=0.0;} // S:PAW_PBE_KIN:06Sep2000
  if(AUID=="e75b388c58799a1c" && nKIN) {found=TRUE;groundstate_structure="A16";groundstate_energy=-4.126363;volume_atom=36.36556;spin_atom=0.0;} // S:PAW_PBE:17Jan2003
  if(AUID=="a18b03f7564daaaa" && SCAN) {found=TRUE;groundstate_structure="A16";groundstate_energy=-9.64992;volume_atom=26.63135;spin_atom=0.0;} // S:PAW_PBE_KIN:SCAN:06Sep2000

  // Sb
  if(AUID=="6ec70613ee3528b2" && nKIN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-3.88932;volume_atom=27.1685;spin_atom=0.0;} // Sb:PAW_PBE:06Sep2000
  if(AUID=="bbf81761207d2817" && SCAN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-37.1553;volume_atom=28.6389;spin_atom=0.0;} // Sb:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="bbf81761207d2817" && nKIN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-4.13485;volume_atom=31.7896;spin_atom=0.0;} // Sb:PAW_PBE_KIN:06Sep2000

  // Sc
  if(AUID=="4e27e0eebb4ca9ea" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.69169;volume_atom=22.2639;spin_atom=0.0;} // Sc:LDA:01Apr2000
  if(AUID=="11f9d8ee4231562f" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.13734;volume_atom=24.1795;spin_atom=0.0;} // Sc:GGA:01Apr2000
  if(AUID=="1d27da013e36ce9b" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.18297;volume_atom=24.0170;spin_atom=0.0;} // Sc:PAW_GGA:08Aug2001
  if(AUID=="90cbeb0ff7fd69f7" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.15503;volume_atom=24.8607;spin_atom=0.0;} // Sc:GGA:01Apr2000
  if(AUID=="5b645c3a2198a968" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.67171;volume_atom=22.9075;spin_atom=0.0;} // Sc:LDA:01Apr2000
  if(AUID=="a563cb0b31f21d93" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.22496;volume_atom=24.2228;spin_atom=0.0;} // Sc_sv:PAW_GGA:07Sep2000
  if(AUID=="4b48722b5ae8f9e9" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.33212;volume_atom=24.4214;spin_atom=0.0;} // Sc_sv:PAW_PBE:07Sep2000
  if(AUID=="dd64f21bfc2e4cab" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.72459;volume_atom=22.3699;spin_atom=0.0;} // Sc_sv:PAW_LDA:07Sep2000
  if(AUID=="69bd66903db1e199" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-14.1593;volume_atom=24.2183;spin_atom=0.619075;} // Sc:PAW_PBE_KIN:SCAN:04Feb2005  MAGNETIC ??
  if(AUID=="83d4f13df3294a60" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-14.2457;volume_atom=24.7117;spin_atom=0.597575;} // Sc_sv:PAW_PBE_KIN:SCAN:07Sep2000  MAGNETIC ??
  if(AUID=="83d4f13df3294a60" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.24799;volume_atom=24.4193;spin_atom=0.0;} // Sc_sv:PAW_PBE_KIN:07Sep2000
  if(AUID=="2bfc025403470bd0" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.729086;volume_atom=22.22153;spin_atom=0.0;} // Sc:PAW_LDA_KIN:26Mar2012
  if(AUID=="69bd66903db1e199" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.202602;volume_atom=24.1799;spin_atom=0.0;} // Sc:PAW_PBE_KIN:04Feb2005
  if(AUID=="94509e69a97fd6ec" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.171505;volume_atom=22.30684;spin_atom=0.0;} // Sc_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="109c2eef02dea969" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.725005;volume_atom=22.36042;spin_atom=0.0;} // Sc_sv:PAW_LDA_KIN:07Sep2000
  if(AUID=="bbde08d55886a708" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.695936;volume_atom=24.35525;spin_atom=0.0;} // Sc_sv_GW:PAW_PBE_KIN:05Dec2013

  // Se
  if(AUID=="40060a482797271e" && nKIN) {found=TRUE;groundstate_structure="A8";groundstate_energy=-3.48266;volume_atom=29.6441;spin_atom=0.0;} // Se:PAW_PBE:06Sep2000
  if(AUID=="00143077dee7333f" && SCAN) {found=TRUE;groundstate_structure="A8";groundstate_energy=-20.0853;volume_atom=28.3006;spin_atom=0.0;} // Se:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="00143077dee7333f" && nKIN) {found=TRUE;groundstate_structure="A8";groundstate_energy=-3.49791;volume_atom=29.7337;spin_atom=0.0;} // Se:PAW_PBE_KIN:06Sep2000

  // Si
  if(AUID=="096bd93da38a71ba" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.46219;volume_atom=20.3474;spin_atom=0.0;} // Si_h:PAW_GGA:08Apr2002
  if(AUID=="b19b07eb7ec794be" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.43180;volume_atom=20.2984;spin_atom=0.0;} // Si:GGA:01Apr2000
  if(AUID=="5d164d220fd7ceee" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.44184;volume_atom=20.3640;spin_atom=0.0;} // Si_h:PAW_PBE:08Apr2002
  if(AUID=="69684f2b4007eeda" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.44713;volume_atom=20.2871;spin_atom=0.0;} // Si:GGA:01Apr2000
  if(AUID=="530c405087cb5e2b" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.99306;volume_atom=19.5385;spin_atom=0.0;} // Si:LDA:01Apr2000
  if(AUID=="940bc6e41f3ea1ba" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.96148;volume_atom=19.6993;spin_atom=0.0;} // Si_h:PAW_LDA:21Jan2003
  if(AUID=="6b54cf2aa9fd187c" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.97452;volume_atom=19.5596;spin_atom=0.0;} // Si:LDA:01Apr2000
  if(AUID=="2e636f3ea3dd411d" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.95931;volume_atom=19.7009;spin_atom=0.0;} // Si:PAW_LDA:02Apr1999
  if(AUID=="3ea77bd2a2af2ad4" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.42373;volume_atom=20.4310;spin_atom=0.0;} // Si:PAW_PBE:05Jan2001
  if(AUID=="d47b3f10236456be" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.43030;volume_atom=20.4209;spin_atom=0.0;} // Si:PAW_GGA:05Jan2001
  if(AUID=="afc7580b135e98ce" && SCAN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-10.0057;volume_atom=19.9704;spin_atom=0.0;} // Si:PAW_PBE_KIN:SCAN:05Jan2001
  if(AUID=="afc7580b135e98ce" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.42376;volume_atom=20.4329;spin_atom=0.0;} // Si:PAW_PBE_KIN:05Jan2001

  // Sm
  if(AUID=="d3011b597251aacb" && nKIN) {found=TRUE;groundstate_structure="C19";groundstate_energy=-4.712267;volume_atom=33.91257;spin_atom=0.0;} // Sm_3:PAW_PBE:07Sep2000
  if(AUID=="8d38efd1df55cc4b" && nKIN) {found=TRUE;groundstate_structure="C19";groundstate_energy=-4.63869;volume_atom=33.59555;spin_atom=0.0;} // Sm_3:PAW_GGA:11May2000

  // Sn
  if(AUID=="ec3f84d587dd2653" && nKIN) {found=TRUE;groundstate_structure="A5";groundstate_energy=-3.96372;volume_atom=28.1306;spin_atom=0.0;} // Sn_d:PAW_PBE:06Sep2000
  if(AUID=="37c5d78c5699ee73" && nKIN) {found=TRUE;groundstate_structure="A5";groundstate_energy=-3.79514;volume_atom=28.3410;spin_atom=0.0;} // Sn:PAW_PBE:08Apr2002
  if(AUID=="4ddcdf83c7589054" && SCAN) {found=TRUE;groundstate_structure="A5";groundstate_energy=-35.8121;volume_atom=27.4595;spin_atom=0.0;} // Sn:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="6431a791389bb4f7" && SCAN) {found=TRUE;groundstate_structure="A5";groundstate_energy=-35.7984;volume_atom=27.2303;spin_atom=0.0;} // Sn_d:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="6431a791389bb4f7" && nKIN) {found=TRUE;groundstate_structure="A5";groundstate_energy=-3.80374;volume_atom=28.1095;spin_atom=0.0;} // Sn_d:PAW_PBE_KIN:06Sep2000
  if(AUID=="4ddcdf83c7589054" && nKIN) {found=TRUE;groundstate_structure="A5";groundstate_energy=-3.79431;volume_atom=28.3399;spin_atom=0.0;} // Sn:PAW_PBE_KIN:08Apr2002
  if(AUID=="70e23373b95459e4" && nKIN) {found=TRUE;groundstate_structure="A5";groundstate_energy=-3.827754;volume_atom=28.28456;spin_atom=0.0;} // Sn:GGA:01Apr2000

  // Sr
  // find /common/LIB1/LIB/Sr* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="82ec420862ed9432" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.73980;volume_atom=40.5798;spin_atom=0.0;} // Sr:GGA:01Apr2000
  if(AUID=="6fb0d22509f8b232" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.92701;volume_atom=41.8890;spin_atom=0.0;} // Sr:LDA:01Apr2000
  if(AUID=="24554a4dcb7ec11b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.61166;volume_atom=54.1148;spin_atom=0.0;} // Sr_pv:GGA:01Apr2000
  if(AUID=="b3f00f0671094a1c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.90232;volume_atom=48.0657;spin_atom=0.0;} // Sr_pv:LDA:01Apr2000
  if(AUID=="bbfa6aec08364643" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.62474;volume_atom=53.7508;spin_atom=0.0;} // Sr_sv:PAW_GGA:10Feb1998
  if(AUID=="3951da71de7f3413" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.88167;volume_atom=48.2867;spin_atom=0.0;} // Sr_sv:PAW_LDA:10Feb1998
  if(AUID=="68b7a60c49c97739" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.68354;volume_atom=54.6082;spin_atom=0.0631431;} /// Sr_sv:PAW_PBE:07Sep2000
  if(AUID=="403b3e0afc9132cc" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-22.3021;volume_atom=54.0609;spin_atom=0.156893;} // Sr_sv:PAW_PBE_KIN:SCAN:07Sep2000
  if(AUID=="403b3e0afc9132cc" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.63703;volume_atom=54.4729;spin_atom=0.0527548;} // Sr_sv:PAW_PBE_KIN:07Sep2000
  if(AUID=="413055480ee6327d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.65192;volume_atom=48.08185;spin_atom=0.0;} // Sr_sv_GW:PAW_LDA_KIN:23Mar2010
  if(AUID=="a447636c05ea63c7" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.354383;volume_atom=54.30419;spin_atom=0.0617999;} // Sr_sv_GW:PAW_PBE_KIN:23Mar2010
  if(AUID=="0ef5f931a158c7cf" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.898481;volume_atom=48.29546;spin_atom=0.0;} // Sr_sv:PAW_LDA_KIN:10Feb1998

  // Ta
  if(AUID=="ecb6cb1d273c3dcf" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.8043;volume_atom=17.8706;spin_atom=0.0;} // Ta:GGA:01Apr2000
  if(AUID=="c88d4f9a05e45d41" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.7802;volume_atom=18.0957;spin_atom=0.0;} // Ta:PAW_GGA:06Feb2003
  if(AUID=="d1561ec02d179810" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.8929;volume_atom=17.1116;spin_atom=0.0;} // Ta:PAW_LDA:21Jan2003
  if(AUID=="d83777bdbc002be9" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.8600;volume_atom=18.0798;spin_atom=0.0;} // Ta:PAW_PBE:17Jan2003
  if(AUID=="de2a3eb2f27e9633" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.7319;volume_atom=18.2562;spin_atom=0.0;} // Ta_pv:PAW_GGA:07Sep2000
  if(AUID=="926ce9924f7a5fd2" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.8489;volume_atom=18.2778;spin_atom=0.0;} // Ta_pv:PAW_PBE:07Sep2000
  if(AUID=="6d000fcf1ccde6da" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.8285;volume_atom=17.2465;spin_atom=0.0;} // Ta_pv:PAW_LDA:07Sep2000
  if(AUID=="096d94d0306813d6" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-66.6819;volume_atom=17.5029;spin_atom=0.0;} // Ta:PAW_PBE_KIN:SCAN:17Jan2003
  if(AUID=="2287798c5502b426" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-66.6142;volume_atom=17.5162;spin_atom=0.0;} // Ta_pv:PAW_PBE_KIN:SCAN:07Sep2000
  if(AUID=="2287798c5502b426" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.8085;volume_atom=18.2463;spin_atom=0.0;} // Ta_pv:PAW_PBE_KIN:07Sep2000
  if(AUID=="5fb21f2d50782c0a" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.89255;volume_atom=17.10657;spin_atom=0.0;} // Ta:PAW_LDA_KIN:21Jan2003
  if(AUID=="096d94d0306813d6" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.85995;volume_atom=18.08551;spin_atom=0.0;} // Ta:PAW_PBE_KIN:17Jan2003
  if(AUID=="c09ef008a3f0e9d6" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.82661;volume_atom=17.24153;spin_atom=0.0;} // Ta_pv:PAW_LDA_KIN:07Sep2000
  if(AUID=="2d534c8fcbaedf16" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-15.3079;volume_atom=17.13619;spin_atom=0.0;} // Ta_sv_GW:PAW_LDA_KIN:23Mar2010
  if(AUID=="bcd3b7b7689d66d7" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-14.29729;volume_atom=18.14686;spin_atom=0.0;} // Ta_sv_GW:PAW_PBE_KIN:23Mar2010

  // Tc
  if(AUID=="fd329b53fcbb5dbd" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.2088;volume_atom=14.4754;spin_atom=0.0;} // Tc:GGA:01Apr2000
  if(AUID=="aabb0d26abfc5599" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.3047;volume_atom=14.4300;spin_atom=0.0;} // Tc:PAW_PBE:17Jan2003
  if(AUID=="34c9772cc010d2cb" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.6391;volume_atom=13.7601;spin_atom=0.0;} // Tc:LDA:01Apr2000
  if(AUID=="80004c8d11aa64f2" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.2801;volume_atom=14.4848;spin_atom=0.0;} // Tc:PAW_GGA:21Dec2000
  if(AUID=="0eb01e3427b59ddc" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.2038;volume_atom=14.5814;spin_atom=0.0;} // Tc_pv:PAW_GGA:20Feb1998
  if(AUID=="43f93be60f147933" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.6972;volume_atom=13.8130;spin_atom=0.0;} // Tc:PAW_LDA:21Jan2003
  if(AUID=="c33964edb4f74ebf" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.3597;volume_atom=14.5515;spin_atom=0.0;} // Tc_pv:PAW_PBE:06Sep2000
  if(AUID=="aaf18f3b4fa31411" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.6002;volume_atom=13.9073;spin_atom=0.0;} // Tc_pv:PAW_LDA:03Mar1998
  if(AUID=="5751220a43dddd9b" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-35.0200;volume_atom=13.9455;spin_atom=0.0;} // Tc:PAW_PBE_KIN:SCAN:04Feb2005
  if(AUID=="5c4bf0cceb477531" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-35.0765;volume_atom=14.0558;spin_atom=0.0;} // Tc_pv:PAW_PBE_KIN:SCAN:04Feb2005
  if(AUID=="3816a915d663bd59" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-35.0857;volume_atom=14.0351;spin_atom=0.0;} // Tc_sv:PAW_PBE_KIN:SCAN:23Mar2010
  if(AUID=="5c4bf0cceb477531" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.352;volume_atom=14.3541;spin_atom=0.0;} // Tc_pv:PAW_PBE_KIN:04Feb2005
  if(AUID=="7f49cada4f0d0649" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.77448;volume_atom=13.65925;spin_atom=0.0;} // Tc:PAW_LDA_KIN:04Feb2005
  if(AUID=="5751220a43dddd9b" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.37483;volume_atom=14.27513;spin_atom=0.0;} // Tc:PAW_PBE_KIN:04Feb2005
  if(AUID=="09d8862331cbb74e" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.73743;volume_atom=13.73272;spin_atom=0.0;} // Tc_pv:PAW_LDA_KIN:04Feb2005
  if(AUID=="69d305a31cd8c1fa" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.61325;volume_atom=14.32109;spin_atom=0.0;} // Tc_sv_GW:PAW_PBE_KIN:05Dec2013
  if(AUID=="1c51ade7c2153ec3" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.99266;volume_atom=13.70234;spin_atom=0.0;} // Tc_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="89fc48286bd60f36" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.70304;volume_atom=13.72423;spin_atom=0.0;} // Tc_sv:PAW_LDA_KIN:23Mar2010
  if(AUID=="3816a915d663bd59" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.35915;volume_atom=14.35312;spin_atom=0.0;} // Tc_sv:PAW_PBE_KIN:23Mar2010

  // Te
  if(AUID=="96c99e8488c592ae" && nKIN) {found=TRUE;groundstate_structure="A8";groundstate_energy=-3.14141;volume_atom=34.8509;spin_atom=0.0;} // Te:PAW_PBE:08Apr2002
  if(AUID=="8d2f6568d1d06421" && SCAN) {found=TRUE;groundstate_structure="A8";groundstate_energy=-37.2600;volume_atom=33.9467;spin_atom=0.0;} // Te:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="8d2f6568d1d06421" && nKIN) {found=TRUE;groundstate_structure="A8";groundstate_energy=-3.14191;volume_atom=34.9412;spin_atom=0.0;} // Te:PAW_PBE_KIN:08Apr2002

  // Ti
  if(AUID=="e97d34c6b4cd62df" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.53337;volume_atom=15.9367;spin_atom=0.0;} // Ti:LDA:01Apr2000
  if(AUID=="8498e164009c60f7" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.73332;volume_atom=17.2006;spin_atom=0.0;} // Ti:GGA:01Apr2000
  if(AUID=="264e8b27aa49e864" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.53301;volume_atom=15.9106;spin_atom=0.0;} // Ti:PAW_LDA:03Oct2001
  if(AUID=="47d8726ba10b866f" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.74356;volume_atom=17.0410;spin_atom=0.0;} // Ti:PAW_GGA:08Aug2001
  if(AUID=="e18697ab9588a353" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.76396;volume_atom=17.1060;spin_atom=0.0;} // Ti:PAW_PBE:08Apr2002
  if(AUID=="257a6b2b7e91023a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.75293;volume_atom=17.5203;spin_atom=0.0;} // Ti:GGA:01Apr2000
  if(AUID=="8e962327621a06f8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.51437;volume_atom=16.2957;spin_atom=0.0;} // Ti:LDA:01Apr2000
  if(AUID=="8fecbeeb2db133cb" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.77528;volume_atom=17.1860;spin_atom=0.0;} // Ti_pv:PAW_GGA:07Sep2000
  if(AUID=="c8c8567f19acfb26" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.80066;volume_atom=17.1803;spin_atom=0.0;} // Ti_sv:PAW_GGA:07Sep2000
  if(AUID=="567ef720cface091" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.54336;volume_atom=16.0264;spin_atom=0.0;} // Ti_pv:PAW_LDA:07Sep2000
  if(AUID=="3f1fd5d5748fffbc" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.89056;volume_atom=17.2662;spin_atom=0.0;} // Ti_pv:PAW_PBE:07Sep2000
  if(AUID=="9cdd346558c528c6" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.93863;volume_atom=17.2475;spin_atom=0.0;} // Ti_sv:PAW_PBE:07Sep2000
  if(AUID=="5b9d5db647f6fc85" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.53889;volume_atom=16.0094;spin_atom=0.0;} // Ti_sv:PAW_LDA:07Sep2000
  if(AUID=="e19a9be220c9f36d" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-15.9523;volume_atom=16.9114;spin_atom=0.0;} // Ti:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="571036b918ee987a" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-15.9827;volume_atom=17.0821;spin_atom=0.0;} // Ti_pv:PAW_PBE_KIN:SCAN:07Sep2000
  if(AUID=="9695cfc266c5809c" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-16.0198;volume_atom=17.2317;spin_atom=0.0;} // Ti_sv:PAW_PBE_KIN:SCAN:26Sep2005
  if(AUID=="9695cfc266c5809c" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.82563;volume_atom=17.2087;spin_atom=0.0;} // Ti_sv:PAW_PBE_KIN:26Sep2005
  if(AUID=="6ac0bee4033365b9" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.533261;volume_atom=15.91457;spin_atom=0.0;} // Ti:PAW_LDA_KIN:03Oct2001
  if(AUID=="e19a9be220c9f36d" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.763801;volume_atom=17.10648;spin_atom=0.0;} // Ti:PAW_PBE_KIN:08Apr2002
  if(AUID=="571036b918ee987a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.801842;volume_atom=17.25905;spin_atom=0.0;} // Ti_pv:PAW_PBE_KIN:07Sep2000
  if(AUID=="97a9b01629f53a13" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.544074;volume_atom=16.02674;spin_atom=0.0;} // Ti_pv:PAW_LDA_KIN:07Sep2000
  if(AUID=="7fe2ca0acf938fa1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.75661;volume_atom=15.951;spin_atom=0.0;} // Ti_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="60802bd648654a15" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.05503;volume_atom=17.18451;spin_atom=0.0;} // Ti_sv_GW:PAW_PBE_KIN:05Dec2013
  if(AUID=="24ae6dcc846ed9d1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.540703;volume_atom=15.97375;spin_atom=0.0;} // Ti_sv:PAW_LDA_KIN:26Sep2005

  // Tl
  if(AUID=="9f2d6c60104e71a3" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.25921;volume_atom=30.8473;spin_atom=0.0;} // Tl_d:GGA:01Apr2000
  if(AUID=="06fbb1e6c3452edb" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.93920;volume_atom=27.1009;spin_atom=0.0;} // Tl_d:LDA:01Apr2000
  if(AUID=="958e09d75d34dcd0" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.23399;volume_atom=30.4705;spin_atom=0.0;} // Tl_d:PAW_GGA:11Feb1998
  if(AUID=="1a0488ce9d93df98" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.90255;volume_atom=26.8803;spin_atom=0.0;} // Tl_d:PAW_LDA:11Feb1998
  if(AUID=="b88a3147de55f4a8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.24350;volume_atom=31.1190;spin_atom=0.0;} // Tl:PAW_PBE:08Apr2002
  if(AUID=="f473cbb4831e0bc0" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.36274;volume_atom=30.7601;spin_atom=0.0;} // Tl_d:PAW_PBE:06Sep2000
  if(AUID=="1298edb6f5cd44f1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.28269;volume_atom=30.9588;spin_atom=0.0;} // Tl:PAW_GGA:25Jul2001
  if(AUID=="2ad3380a31e82af9" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.95325;volume_atom=27.3640;spin_atom=0.0;} // Tl:PAW_LDA:03Oct2001
  if(AUID=="0e2337aac65480f2" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-72.3276;volume_atom=29.3500;spin_atom=0.0;} // Tl:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="4f0fac6e26fd86a9" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-72.3189;volume_atom=29.0596;spin_atom=0.0;} // Tl_d:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="4f0fac6e26fd86a9" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.22536;volume_atom=30.6548;spin_atom=0.0;} // Tl_d:PAW_PBE_KIN:06Sep2000
  if(AUID=="6ac0bee4033365b9" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.533261;volume_atom=15.91457;spin_atom=0.0;} // Ti:PAW_LDA_KIN:03Oct2001
  if(AUID=="e19a9be220c9f36d" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.763801;volume_atom=17.10648;spin_atom=0.0;} // Ti:PAW_PBE_KIN:08Apr2002
  if(AUID=="571036b918ee987a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.801842;volume_atom=17.25905;spin_atom=0.0;} // Ti_pv:PAW_PBE_KIN:07Sep2000
  if(AUID=="97a9b01629f53a13" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.544074;volume_atom=16.02674;spin_atom=0.0;} // Ti_pv:PAW_LDA_KIN:07Sep2000
  if(AUID=="7fe2ca0acf938fa1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.75661;volume_atom=15.951;spin_atom=0.0;} // Ti_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="60802bd648654a15" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.05503;volume_atom=17.18451;spin_atom=0.0;} // Ti_sv_GW:PAW_PBE_KIN:05Dec2013
  if(AUID=="24ae6dcc846ed9d1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.540703;volume_atom=15.97375;spin_atom=0.0;} // Ti_sv:PAW_LDA_KIN:26Sep2005

  // V
  if(AUID=="aee485cb636ab31d" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.90719;volume_atom=13.3310;spin_atom=0.0;} // V:GGA:01Apr2000
  if(AUID=="40048762e7327671" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.95900;volume_atom=12.3899;spin_atom=0.0;} // V:LDA:01Apr2000
  if(AUID=="dce68ed9d844d970" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.92373;volume_atom=13.1803;spin_atom=0.0;} // V:PAW_GGA:07Aug2001
  if(AUID=="a741261c8d2045c9" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.95883;volume_atom=12.3350;spin_atom=0.0;} // V:PAW_LDA:07Aug2001
  if(AUID=="183846be5b03c64d" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.94402;volume_atom=13.2022;spin_atom=0.0;} // V:PAW_PBE:08Apr2002
  if(AUID=="ea4e243d3a0747f9" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.99426;volume_atom=12.6635;spin_atom=0.0;} // V:LDA:01Apr2000
  if(AUID=="1dd37802a9ba101a" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.92582;volume_atom=13.3629;spin_atom=0.0;} // V_pv:PAW_GGA:07Sep2000
  if(AUID=="d16570ec1e90da05" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.99489;volume_atom=13.5698;spin_atom=0.0;} // V:GGA:01Apr2000
  if(AUID=="0725afc42e47381e" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.93464;volume_atom=12.5085;spin_atom=0.0;} // V_pv:PAW_LDA:07Sep2000
  if(AUID=="b7d3211b63a6261e" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.07822;volume_atom=13.3934;spin_atom=0.0;} // V_pv:PAW_PBE:07Sep2000
  if(AUID=="6555390e553533e4" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.95881;volume_atom=13.3628;spin_atom=0.0;} // V_sv:PAW_GGA:14Sep2000
  if(AUID=="c3b5d34259b9b5f9" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.93559;volume_atom=12.5083;spin_atom=0.0;} // V_sv:PAW_LDA:07Sep2000
  if(AUID=="94da8dda3e65c762" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.11496;volume_atom=13.3971;spin_atom=0.0;} // V_sv:PAW_PBE:07Sep2000
  if(AUID=="62ba629379463606" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-17.3619;volume_atom=12.8510;spin_atom=0.0;} // V:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="46928c16c2f8a488" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-17.4590;volume_atom=13.0021;spin_atom=0.0;} // V_sv:PAW_PBE_KIN:SCAN:02Aug2007
  if(AUID=="ff7f68300ba3e357" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-17.4221;volume_atom=13.0707;spin_atom=0.0;} // V_pv:PAW_PBE_KIN:SCAN:07Sep2000
  if(AUID=="46928c16c2f8a488" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.98932;volume_atom=13.3687;spin_atom=0.0;} // V_sv:PAW_PBE_KIN:02Aug2007
  if(AUID=="ea63bd55cec65d9b" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.959163;volume_atom=12.33326;spin_atom=0.0;} // V:PAW_LDA_KIN:07Aug2001
  if(AUID=="62ba629379463606" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.943553;volume_atom=13.20028;spin_atom=0.0;} // V:PAW_PBE_KIN:08Apr2002
  if(AUID=="48af57fd67c93591" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.936945;volume_atom=12.5009;spin_atom=0.0;} // V_pv:PAW_LDA_KIN:07Sep2000
  if(AUID=="ff7f68300ba3e357" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.953831;volume_atom=13.38845;spin_atom=0.0;} // V_pv:PAW_PBE_KIN:07Sep2000
  if(AUID=="243544817f53020d" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.90621;volume_atom=12.41552;spin_atom=0.0;} // V_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="72baf87216b94021" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.97486;volume_atom=13.29934;spin_atom=0.0;} // V_sv_GW:PAW_PBE_KIN:05Dec2013
  if(AUID=="af90f4bc7d752a11" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.950845;volume_atom=12.47145;spin_atom=0.0;} // V_sv:PAW_LDA_KIN:02Aug2007

  // W
  if(AUID=="fc6e347dc7ddd7ab" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.8623;volume_atom=15.9358;spin_atom=0.0;} // W:GGA:01Apr2000
  if(AUID=="4a36ee2cfb9fae34" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.9177;volume_atom=15.9658;spin_atom=0.0;} // W:PAW_GGA:21Dec2000
  if(AUID=="6a5db51b820aa7a5" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-14.2099;volume_atom=15.2787;spin_atom=0.0;} // W:PAW_LDA:19Jan2001
  if(AUID=="c895658f3e6606f4" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-14.1883;volume_atom=15.1794;spin_atom=0.0;} // W:LDA:01Apr2000
  if(AUID=="22a83c4c0d26e6eb" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.7816;volume_atom=16.1779;spin_atom=0.0;} // W_pv:PAW_GGA:15Jul1998
  if(AUID=="e05803166f7d0e3e" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-14.0476;volume_atom=15.4745;spin_atom=0.0;} // W_pv:PAW_LDA:22Jul1998
  if(AUID=="cad81c4b63d2d96e" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.9546;volume_atom=16.1649;spin_atom=0.0;} // W_pv:PAW_PBE:06Sep2000
  if(AUID=="1febde78746c14e4" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-13.0124;volume_atom=15.9045;spin_atom=0.0;} // W:PAW_PBE:08Apr2002
  if(AUID=="b99d5a263927c2ea" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-69.9743;volume_atom=15.5134;spin_atom=0.0;} // W:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="1d8c010cef748b9c" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-69.6572;volume_atom=15.6161;spin_atom=0.0;} // W_sv:PAW_PBE_KIN:SCAN:04Sep2015
  if(AUID=="1d8c010cef748b9c" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.9534;volume_atom=16.1168;spin_atom=0.0;} // W_sv:PAW_PBE_KIN:04Sep2015
  if(AUID=="78d442c6ee8e1864" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-14.20311;volume_atom=15.28103;spin_atom=0.0;} // W:PAW_LDA_KIN:19Jan2001
  if(AUID=="b99d5a263927c2ea" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-13.01182;volume_atom=15.91018;spin_atom=0.0;} // W:PAW_PBE_KIN:08Apr2002
  if(AUID=="cdaff66d42eb21cb" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-15.01969;volume_atom=16.04622;spin_atom=0.0;} // W_sv_GW:PAW_PBE_KIN:23Mar2010
  if(AUID=="61a24681ec6520ee" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-16.171;volume_atom=15.38559;spin_atom=0.0;} // W_sv_GW:PAW_LDA_KIN:23Mar2010
  if(AUID=="c5c6cbb458696867" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-14.10342;volume_atom=15.46532;spin_atom=0.0;} // W_sv:PAW_LDA_KIN:04Sep2015

  // Xe
  if(AUID=="cde8ed26c95905b5" && nKIN) {found=TRUE;groundstate_structure="isolated";groundstate_energy=-0.008491;volume_atom=3375;spin_atom=0.0;} // Xe:PAW_PBE:07Sep2000

  // Y
  if(AUID=="377b8f487d291e0e" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.33584;volume_atom=32.5652;spin_atom=0.0;} // Y:GGA:01Apr2000
  if(AUID=="a0835c00b2263c8b" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.40419;volume_atom=31.9468;spin_atom=0.0;} // Y:GGA:01Apr2000
  if(AUID=="1844bc40a30b31c7" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.90740;volume_atom=29.8089;spin_atom=0.0;} // Y_sv:PAW_LDA:10Feb1998
  if(AUID=="20813d90f7aff62f" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.38136;volume_atom=32.4851;spin_atom=0.0;} // Y_sv:PAW_GGA:10Feb1998
  if(AUID=="18ed87089235423b" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.98897;volume_atom=29.2790;spin_atom=0.0;} // Y:LDA:01Apr2000
  if(AUID=="5c8efcbfc77486c0" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.89015;volume_atom=29.8469;spin_atom=0.0;} // Y:LDA:01Apr2000
  if(AUID=="b07bf0a656aa22bc" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.46317;volume_atom=32.7578;spin_atom=0.0;} // Y_sv:PAW_PBE:06Sep2000
  if(AUID=="2fb876577258940d" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-27.8500;volume_atom=33.3097;spin_atom=0.602783;} // Y_sv:PAW_PBE_KIN:SCAN:25May2007   MAGNETIC ??
  if(AUID=="2fb876577258940d" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.42905;volume_atom=32.5612;spin_atom=0.0;} // Y_sv:PAW_PBE_KIN:25May2007
  if(AUID=="2db4b684eb247746" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.111379;volume_atom=29.6223;spin_atom=0.0;} // Y_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="3189f0e60a077125" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.570338;volume_atom=32.51187;spin_atom=0.0;} // Y_sv_GW:PAW_PBE_KIN:05Dec2013
  if(AUID=="1eae9b15c3704e56" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.942676;volume_atom=29.62325;spin_atom=0.0;} // Y_sv:PAW_LDA_KIN:25May2007

  // Yb
  // find /common/LIB1/LIB/Yb* -name "A1" | sed "s/\/common/aflow --force --lib2raw=&/g" > /tmp/xpot && aflow --multi=8 --file=/tmp/xpot | grep XXX | sort | uniq
  if(AUID=="4e35fa207bf4aaa9" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.51978;volume_atom=39.7380;spin_atom=0.0;} // Yb_2:PAW_PBE:06Sep2000
  if(AUID=="b1f8bc4d83e97833" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.44398;volume_atom=39.1468;spin_atom=0.0;} // Yb_2:PAW_GGA:10May2000
  if(AUID=="de1c22a384a8d945" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.67034;volume_atom=36.9163;spin_atom=0.0;} // Yb:PAW_PBE:24Feb2003
  if(AUID=="a2a66bbb705fd114" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-50.5319;volume_atom=37.1321;spin_atom=0.0;} // Yb_2:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="d4d1484e34700879" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-53.5119;volume_atom=28.3498;spin_atom=0.0;} // Yb_3:PAW_PBE_KIN:SCAN:08Jul2013
  if(AUID=="a2a66bbb705fd114" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.45707;volume_atom=39.26736;spin_atom=0.0;} // Yb_2:PAW_PBE_KIN:06Sep2000
  if(AUID=="d4d1484e34700879" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.433541;volume_atom=29.52432;spin_atom=0.0;} // Yb_3:PAW_PBE_KIN:08Jul2013

  // Zn
  if(AUID=="c765a5877b95887a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.12141;volume_atom=15.3617;spin_atom=0.0;} // Zn:GGA:01Apr2000
  if(AUID=="e75d47ef817f48d8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.11328;volume_atom=15.2205;spin_atom=0.0;} // Zn:PAW_GGA:03Mar1998
  if(AUID=="4b8b99f638e1a173" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.86830;volume_atom=13.5768;spin_atom=0.0;} // Zn:PAW_LDA:03Mar1998
  if(AUID=="d321d33f5619ef83" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.89186;volume_atom=13.5675;spin_atom=0.0;} // Zn:LDA:01Apr2000
  if(AUID=="d4ad14ef9da00329" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.26581;volume_atom=15.1693;spin_atom=0.0;} // Zn:PAW_PBE:06Sep2000
  if(AUID=="7ed9d78dd30715fd" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-13.4842;volume_atom=14.1685;spin_atom=0.0;} // Zn:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="7ed9d78dd30715fd" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.10624;volume_atom=15.1313;spin_atom=0.0;} // Zn:PAW_PBE_KIN:06Sep2000
  if(AUID=="264df91d482e7eb1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.876314;volume_atom=13.3412;spin_atom=0.0;} // Zn_GW:PAW_LDA_KIN:09Oct2010
  if(AUID=="710533d0171c7210" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.875921;volume_atom=13.45657;spin_atom=0.0;} // Zn:PAW_LDA_KIN:03Mar1998
  if(AUID=="c9d470aabc061907" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.103282;volume_atom=15.01631;spin_atom=0.0;} // Zn_GW:PAW_PBE_KIN:09Oct2010
  if(AUID=="02e332a61a92f598" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.860839;volume_atom=13.05282;spin_atom=0.0;} // Zn_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="53d23506232b67ac" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.124756;volume_atom=14.5293;spin_atom=0.0;} // Zn_sv_GW:PAW_PBE_KIN:05Dec2013

  // Zr
  if(AUID=="c87b485bdb0c4683" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.36116;volume_atom=21.3443;spin_atom=0.0;} // Zr:LDA:01Apr2000
  if(AUID=="014ba38906f94787" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.47756;volume_atom=23.4023;spin_atom=0.0;} // Zr:PAW_PBE:08Apr2002
  if(AUID=="05c7a473e8b0fed1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.45764;volume_atom=23.3305;spin_atom=0.0;} // Zr:PAW_GGA:08Aug2001
  if(AUID=="a2c98479cda145dd" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.49915;volume_atom=22.8558;spin_atom=0.0;} // Zr:GGA:01Apr2000
  if(AUID=="844977898c624f61" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.22027;volume_atom=21.8813;spin_atom=0.0;} // Zr_sv:PAW_LDA:10Feb1998
  if(AUID=="ae40f9ef1519f7f7" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.40052;volume_atom=23.4147;spin_atom=0.0;} // Zr:GGA:01Apr2000
  if(AUID=="7471a45b48d5cbca" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.43218;volume_atom=23.3391;spin_atom=0.0;} // Zr_sv:PAW_GGA:10Feb1998
  if(AUID=="3f478742b0dd98d0" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.21967;volume_atom=21.9211;spin_atom=0.0;} // Zr:LDA:01Apr2000
  if(AUID=="8f81b69844f3963f" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.54365;volume_atom=23.4268;spin_atom=0.0;} // Zr_sv:PAW_PBE:07Sep2000
  if(AUID=="a3ef4a907c7b4c96" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-30.6972;volume_atom=23.1863;spin_atom=0.0;} // Zr_sv:PAW_PBE_KIN:SCAN:04Jan2005
  if(AUID=="a3ef4a907c7b4c96" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.51242;volume_atom=23.1783;spin_atom=0.0;} // Zr_sv:PAW_PBE_KIN:04Jan2005
  if(AUID=="4b771ea128970dd1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.26381;volume_atom=23.18791;spin_atom=0.0;} // Zr_sv_GW:PAW_PBE_KIN:05Dec2013
  if(AUID=="1d3446a7830fa99e" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.06212;volume_atom=21.62537;spin_atom=0.0;} // Zr_sv_GW:PAW_LDA_KIN:05Dec2013
  if(AUID=="d9f9fc183b6567ac" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.283964;volume_atom=21.62691;spin_atom=0.0;} // Zr_sv:PAW_LDA_KIN:04Jan2005


  // ./xgo Sm_3:PAW_GGA:11May2000 "found=TRUE;groundstate_structure=\"ICSD_246657\";groundstate_energy=-4.621400;volume_atom=33.447633;spin_atom=0.0;"// FIX
  // ./xgo Sm_3:PAW_GGA:11May2000 && 0 "found=TRUE;groundstate_structure=\"ICSD_652637\";groundstate_energy=-4.64136;volume_atom=33.5075;spin_atom=0.0;"// IT HAS LDAU
  // ./xgo Sm_3:PAW_PBE:07Sep2000 && 0 "found=TRUE;groundstate_structure=\"A1\";groundstate_energy=-4.7062;volume_atom=33.8339;spin_atom=0.0;"// IT HAS LDAU
  // ./xgo Ce "found=TRUE;groundstate_structure=\"A1\";groundstate_energy=-5.92998;volume_atom=26.0579;spin_atom=0.0;"
  // ./xgo Ce "found=TRUE;groundstate_structure=\"ICSD_2284-mS4\";groundstate_energy=-5.93013;volume_atom=26.0697;spin_atom=0.0;"
  // ./xgo Cl_h:PAW_PBE:08Apr2002 "found=TRUE;groundstate_structure=\"A11\";groundstate_energy=-1.8156;volume_atom=37.3299;spin_atom=0.0;"WAITING

  if(!found) { volume_atom=999999,spin_atom=999999;} // some defaults
  //  if(!found) cerr <<"ERROR (xPOTCAR_EnthalpyReference_AUID): NOT FOUND: AUID=" << AUID << endl;// throw aurostd::xerror(__AFLOW_FILE__,XPID+"xPOTCAR_EnthalpyReference_AUID():","Throw for debugging purposes.",_GENERIC_ERROR_);
  if(LDEBUG && !found) cout <<"ERROR (xPOTCAR_EnthalpyReference_AUID): NOT FOUND: AUID=" << AUID << endl;// throw aurostd::xerror(__AFLOW_FILE__,XPID+"xPOTCAR_EnthalpyReference_AUID():","Throw for debugging purposes.",_GENERIC_ERROR_);
  if(LDEBUG &&  found) cout <<"ERROR (xPOTCAR_EnthalpyReference_AUID): FOUND: AUID=" << AUID << endl;// throw aurostd::xerror(__AFLOW_FILE__,XPID+"xPOTCAR_EnthalpyReference_AUID():","Throw for debugging purposes.",_GENERIC_ERROR_);
  if(LDEBUG) cerr << XPID << "xPOTCAR_EnthalpyReference_AUID: [END]" << endl;
  return found;
};



#endif // _AFLOW_XPSEUDOPOTENTIAL_CPP

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************
