// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2007-2019
#ifndef _AFLOW_XPROTOTYPE_CPP
#define _AFLOW_XPROTOTYPE_CPP
#include "aflow.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// _XPROTOTYPE_

namespace xprototype {
  // constructors
  xprototype::xprototype() {free();}  
  // destructor
  xprototype::~xprototype() {free();} 
  // free
  
  void xprototype::free() { 
    verbose=FALSE;
    // add all others
    catalog=""; 
    volume=NNN; 
    label="";   
    parameter_list.clear();
    parameter_values.clear();
    parameter_set_id=""; 
    weblink="";      
    stoichiometry.clear();
    Pearson_symbol="";
    space_group_number=NNN;
    space_group_symbol_H_M="";
    space_group_symbol_Hall="";
    space_group_symbol_Schoenflies="";
    for(uint i=0;i<Wyckoff_letters.size();i++) { Wyckoff_letters.at(i).clear(); } Wyckoff_letters.clear();
    for(uint i=0;i<Wyckoff_site_symmetries.size();i++) { Wyckoff_site_symmetries.at(i).clear(); } Wyckoff_site_symmetries.clear();
    for(uint i=0;i<Wyckoff_multiplicities.size();i++) { Wyckoff_multiplicities.at(i).clear(); } Wyckoff_multiplicities.clear();
    prototype_material=""; 
    common_name=""; 
    mineral_name=""; 
    phase=""; 
    strukturbericht="";
    similar_materials.clear();
    comments.clear();
    title="";   
   }
  
  const xprototype& xprototype::operator=(const xprototype& b) {      // operator=
    if(this!=&b) {free();copy(b);}  
    return *this;
  }

  void xprototype::copy(const xprototype& b) {  //copy PRIVATE  
    // will populate
    verbose=b.verbose;
    catalog=b.catalog;
    volume=b.volume;
    label=b.label;
    parameter_list.clear(); for(uint i=0;i<b.parameter_list.size();i++) { parameter_list.push_back(b.parameter_list.at(i)); }
    parameter_values.clear(); for(uint i=0;i<b.parameter_values.size();i++) { parameter_values.push_back(b.parameter_values.at(i)); }
    parameter_set_id=b.parameter_set_id;
    weblink=b.weblink;
    stoichiometry.clear();for(uint i=0;i<b.stoichiometry.size();i++) { stoichiometry.push_back(b.stoichiometry.at(i)); }
    Pearson_symbol=b.Pearson_symbol;
    space_group_number=b.space_group_number;
    space_group_symbol_H_M=b.space_group_symbol_H_M;
    space_group_symbol_Hall=b.space_group_symbol_Hall;
    space_group_symbol_Schoenflies=b.space_group_symbol_Schoenflies;
    for(uint i=0;i<Wyckoff_letters.size();i++) { Wyckoff_letters.at(i).clear(); }
    Wyckoff_letters.clear(); for(uint i=0;i<b.Wyckoff_letters.size();i++) { Wyckoff_letters.push_back(b.Wyckoff_letters.at(i)); }
    for(uint i=0;i<Wyckoff_site_symmetries.size();i++) { Wyckoff_site_symmetries.at(i).clear(); }
    Wyckoff_site_symmetries.clear(); for(uint i=0;i<b.Wyckoff_site_symmetries.size();i++) { Wyckoff_site_symmetries.push_back(b.Wyckoff_site_symmetries.at(i)); }
    for(uint i=0;i<Wyckoff_multiplicities.size();i++) { Wyckoff_multiplicities.at(i).clear(); }
    Wyckoff_multiplicities.clear(); for(uint i=0;i<b.Wyckoff_multiplicities.size();i++) { Wyckoff_multiplicities.push_back(b.Wyckoff_multiplicities.at(i)); }
    prototype_material=b.prototype_material;
    common_name=b.common_name;
    mineral_name=b.mineral_name;
    phase=b.phase;
    strukturbericht=b.strukturbericht;
    similar_materials.clear();for(uint i=0;i<b.similar_materials.size();i++) { similar_materials.push_back(b.similar_materials.at(i)); }
    comments.clear();for(uint i=0;i<b.comments.size();i++) { comments.push_back(b.comments.at(i)); }
    title=b.title;
 }

  void xprototype::clear(){
    xprototype a;(*this)=a;
  }

  ostream& operator<<(ostream& oss,const xprototype& prototype) {
    oss.setf(std::ios::fixed,std::ios::floatfield);
    oss.precision(10);
    oss << "verbose=" << prototype.verbose << endl;
    oss << "catalog=" << prototype.catalog << endl;
    oss << "volume=" << prototype.volume << endl;
    oss << "label=" << prototype.label << endl;
    oss << "parameter_list="; for(uint i=0;i<prototype.parameter_list.size();i++) { oss << prototype.parameter_list.at(i) << " "; } oss << endl;
    oss << "parameter_values="; for(uint i=0;i<prototype.parameter_values.size();i++) { oss << prototype.parameter_values.at(i) << " "; } oss << endl;
    oss << "parameter_set_id=" << prototype.parameter_set_id << endl;
    oss << "weblink=" << prototype.weblink << endl;
    oss << "stoichiometry=";for(uint i=0;i<prototype.stoichiometry.size();i++) { oss << prototype.stoichiometry.at(i) << " "; } oss << endl;
    oss << "Pearson_symbol=" << prototype.Pearson_symbol << endl;
    oss << "space_group_number=" << prototype.space_group_number << endl;
    oss << "space_group_symbol_H_M=" << prototype.space_group_symbol_H_M << endl;
    oss << "space_group_symbol_Hall=" << prototype.space_group_symbol_Hall << endl;
    oss << "space_group_symbol_Schoenflies=" << prototype.space_group_symbol_Schoenflies << endl;
    oss << "Wyckoff_letters="; for(uint i=0;i<prototype.Wyckoff_letters.size();i++) { for(uint j=0;i<prototype.Wyckoff_letters.at(i).size();j++) { oss << prototype.Wyckoff_letters.at(i).at(j) << " "; } oss << endl; }
    oss << "Wyckoff_site_symmetries="; for(uint i=0;i<prototype.Wyckoff_site_symmetries.size();i++) { for(uint j=0;i<prototype.Wyckoff_site_symmetries.at(i).size();j++) { oss << prototype.Wyckoff_site_symmetries.at(i).at(j) << " "; } oss << endl; }
    oss << "Wyckoff_multiplicities="; for(uint i=0;i<prototype.Wyckoff_multiplicities.size();i++) { for(uint j=0;i<prototype.Wyckoff_multiplicities.at(i).size();j++) { oss << prototype.Wyckoff_multiplicities.at(i).at(j) << " "; } oss << endl; }
    oss << "prototype_material=" << prototype.prototype_material << endl;
    oss << "common_name=" << prototype.common_name << endl;
    oss << "mineral_name=" << prototype.mineral_name << endl;
    oss << "phase=" << prototype.phase << endl;
    oss << "strukturbericht=" << prototype.strukturbericht << endl;
    oss << "similar_materials=";for(uint i=0;i<prototype.similar_materials.size();i++) { oss << prototype.similar_materials.at(i) << " "; } oss << endl;
    oss << "comments=";for(uint i=0;i<prototype.comments.size();i++) { oss << prototype.comments.at(i) << " "; } oss << endl;
    oss << "title=" << prototype.title << endl;
    return oss;
  }

  // ********************************************************************************************************************************************************

  void xprototype::populate(const string& prototype) {  
    free();
    // DEFAULT
    verbose=FALSE;

    // OFFSET
    if(prototype=="") {
      xprototype a;
      (*this)=a;
      return; 
    }
    // ********************************************************************************************************************************************************
    // A2BC4D_tI16_121_d_a_i_b A2BC4D_tI16_121_d_a_i_b A2BC4D_tI16_121_d_a_i_b A2BC4D_tI16_121_d_a_i_b A2BC4D_tI16_121_d_a_i_b A2BC4D_tI16_121_d_a_i_b
    else if(prototype=="A2BC4D_tI16_121_d_a_i_b" || aurostd::toupper(prototype)==aurostd::toupper("stannite")) { // A2BC4D_tI16_121_d_a_i_b
      // stannite example
      catalog="anrl";
      volume=1;
      label="A2BC4D_tI16_121_d_a_i_b";
      parameter_list.push_back("a");parameter_list.push_back("c/a");parameter_list.push_back("x4");parameter_list.push_back("z4");
      parameter_values.push_back(5.46);parameter_values.push_back(1.96428571429);parameter_values.push_back(0.245);parameter_values.push_back(0.132);
      parameter_set_id="001";
      weblink=_AFLOW_PROTOTYPE_ENCYCLOPEDIA_+"A2BC4D_tI16_121_d_a_i_b.html";
      stoichiometry.push_back(2);stoichiometry.push_back(1);stoichiometry.push_back(4);stoichiometry.push_back(1);
      // symmetry
      Pearson_symbol="tI16";
      space_group_number=121;
      space_group_symbol_H_M="I-42m"; // or lookup table
      space_group_symbol_Hall="I -4 2"; // or lookup table
      space_group_symbol_Schoenflies="D_{2d}^{11}"; // or lookup table
      vector<string> Wyckoff_letters_tmp;
      // letters
      Wyckoff_letters_tmp.push_back("d");Wyckoff_letters.push_back(Wyckoff_letters_tmp);Wyckoff_letters_tmp.clear();
      Wyckoff_letters_tmp.push_back("a");Wyckoff_letters.push_back(Wyckoff_letters_tmp);Wyckoff_letters_tmp.clear();
      Wyckoff_letters_tmp.push_back("i");Wyckoff_letters.push_back(Wyckoff_letters_tmp);Wyckoff_letters_tmp.clear();
      Wyckoff_letters_tmp.push_back("b");Wyckoff_letters.push_back(Wyckoff_letters_tmp);Wyckoff_letters_tmp.clear();
      // site symmetries
      vector<string> Wyckoff_site_symmetries_tmp;
      Wyckoff_site_symmetries_tmp.push_back("-4..");Wyckoff_site_symmetries.push_back(Wyckoff_site_symmetries_tmp);Wyckoff_site_symmetries_tmp.clear();
      Wyckoff_site_symmetries_tmp.push_back("-42m");Wyckoff_site_symmetries.push_back(Wyckoff_site_symmetries_tmp);Wyckoff_site_symmetries_tmp.clear();
      Wyckoff_site_symmetries_tmp.push_back("..m");Wyckoff_site_symmetries.push_back(Wyckoff_site_symmetries_tmp);Wyckoff_site_symmetries_tmp.clear();
      Wyckoff_site_symmetries_tmp.push_back("-42m");Wyckoff_site_symmetries.push_back(Wyckoff_site_symmetries_tmp);Wyckoff_site_symmetries_tmp.clear();
      vector<uint> Wyckoff_multiplicities_tmp;
      Wyckoff_multiplicities_tmp.push_back(4);Wyckoff_multiplicities.push_back(Wyckoff_multiplicities_tmp);Wyckoff_multiplicities_tmp.clear();
      Wyckoff_multiplicities_tmp.push_back(2);Wyckoff_multiplicities.push_back(Wyckoff_multiplicities_tmp);Wyckoff_multiplicities_tmp.clear();
      Wyckoff_multiplicities_tmp.push_back(8);Wyckoff_multiplicities.push_back(Wyckoff_multiplicities_tmp);Wyckoff_multiplicities_tmp.clear();
      Wyckoff_multiplicities_tmp.push_back(2);Wyckoff_multiplicities.push_back(Wyckoff_multiplicities_tmp);Wyckoff_multiplicities_tmp.clear();
      // designations
      prototype_material="Cu2FeS4Sn";
      common_name="";
      mineral_name="stannite";
      phase="";
      strukturbericht="H2_{6}";
      similar_materials.push_back("Cu2CdSe4Sn");similar_materials.push_back("CoCu2S4Sn");similar_materials.push_back("Cu2GeHgS4");similar_materials.push_back("Cu2HgS4Sn");similar_materials.push_back("Ag2FeS4Sn");
      comments.push_back("If $c=2a$, $x=1/4$, and $z=3/8$, the atoms are on the sites of the diamond ($A4$) structure. If, in addition, the Cu, Fe, and Sn atoms are replaced by a single atom type, the crystal reduces to the zincblende ($B3$) structure.");
      title="Stannite (Cu2FeS4Sn, H2_{6}) Structure";    
      return; 
    }
    // ********************************************************************************************************************************************************
    
    throw aurostd::xerror(__AFLOW_FILE__,"xprototype::xprototype():","Prototype does not exist: "+prototype,_FILE_NOT_FOUND_);
  }
} // namespace xprototype

#endif // _AFLOW_XPROTOTYPE_CPP

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************
