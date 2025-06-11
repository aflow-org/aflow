// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo
// fixed for XZ - SC 2018-2019
// fixed for new AUID language (SC 2019)
// fixed for tree search on the AUID directories (SC 2019) super-speed

#ifndef _AFLOWLIB_WEB_INTERFACE_CPP_
#define _AFLOWLIB_WEB_INTERFACE_CPP_
#include "aflow.h"
#include "aflow_pocc.h" //CO20200624
#include "aflow_symmetry_spacegroup.h" //DX20200929
#include "aflowlib.h" //SD20221207

const string _DEVIL_PROTOTYPES_STRING_ = "64,65,549,550,f8269,f9083,f8819";

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ***************************************************************************
namespace aflowlib {
  //  class _aflowlib_entry 

  _aflowlib_entry::_aflowlib_entry() {  // constructor PUBLIC
    free();
  }
  _aflowlib_entry::~_aflowlib_entry() { // destructor PUBLIC
    free();

    data_api.clear();
    data_source.clear();
    for(uint i=0;i<vLDAU.size();i++){vLDAU[i].clear();}vLDAU.clear(); //CO20210104 clear()
  }

  void _aflowlib_entry::copy(const _aflowlib_entry& b) { // copy PRIVATE
    entry=b.entry;ventry.clear();for(uint i=0;i<b.ventry.size();i++) ventry.push_back(b.ventry.at(i));
    auid=b.auid;
    vauid=b.vauid;vauid.clear();for(uint i=0;i<b.vauid.size();i++) vauid.push_back(b.vauid.at(i));
    aurl=b.aurl;vaurl.clear();for(uint i=0;i<b.vaurl.size();i++) vaurl.push_back(b.vaurl.at(i));
    system_name=b.system_name;
    vkeywords.clear();for(uint i=0;i<b.vkeywords.size();i++) vkeywords.push_back(b.vkeywords.at(i));
    aflowlib_date=b.aflowlib_date;vaflowlib_date.clear();for(uint i=0;i<b.vaflowlib_date.size();i++) vaflowlib_date.push_back(b.vaflowlib_date.at(i)); //CO20200624 - adding LOCK date
    aflowlib_version=b.aflowlib_version;
    aflowlib_entries=b.aflowlib_entries;
    vaflowlib_entries.clear();for(uint i=0;i<b.vaflowlib_entries.size();i++) vaflowlib_entries.push_back(b.vaflowlib_entries.at(i));
    aflowlib_entries_number=b.aflowlib_entries_number;
    aflow_version=b.aflow_version;
    catalog=b.catalog;
    data_api=b.data_api;
    data_source=b.data_source;
    data_language=b.data_language;
    error_status=b.error_status;
    author=b.author;vauthor.clear();for(uint i=0;i<b.vauthor.size();i++) vauthor.push_back(b.vauthor.at(i));
    calculation_cores=b.calculation_cores;calculation_memory=b.calculation_memory;calculation_time=b.calculation_time;
    corresponding=b.corresponding;vcorresponding.clear();for(uint i=0;i<b.vcorresponding.size();i++) vcorresponding.push_back(b.vcorresponding.at(i));
    loop=b.loop;vloop.clear();for(uint i=0;i<b.vloop.size();i++) vloop.push_back(b.vloop.at(i));
    node_CPU_Cores=b.node_CPU_Cores;node_CPU_MHz=b.node_CPU_MHz;node_CPU_Model=b.node_CPU_Model;node_RAM_GB=b.node_RAM_GB;
    Bravais_lattice_orig=b.Bravais_lattice_orig;Bravais_lattice_relax=b.Bravais_lattice_relax;
    code=b.code;
    composition=b.composition;vcomposition.clear();for(uint i=0;i<b.vcomposition.size();i++) vcomposition.push_back(b.vcomposition.at(i));
    compound=b.compound;
    density=b.density;
    density_orig=b.density_orig; //DX20190124 - add original crystal info
    dft_type=b.dft_type;vdft_type.clear();for(uint i=0;i<b.vdft_type.size();i++) vdft_type.push_back(b.vdft_type.at(i));
    eentropy_cell=b.eentropy_cell;eentropy_atom=b.eentropy_atom;
    Egap=b.Egap;Egap_fit=b.Egap_fit;
    energy_cell=b.energy_cell;energy_atom=b.energy_atom;energy_atom_relax1=b.energy_atom_relax1;
    energy_cutoff=b.energy_cutoff;
    delta_electronic_energy_convergence=b.delta_electronic_energy_convergence;
    delta_electronic_energy_threshold=b.delta_electronic_energy_threshold;
    nkpoints=b.nkpoints;
    nkpoints_irreducible=b.nkpoints_irreducible;
    kppra=b.kppra;
    kpoints=b.kpoints;
    kpoints_nnn_relax=b.kpoints_nnn_relax;
    kpoints_nnn_static=b.kpoints_nnn_static;
    kpoints_pairs=b.kpoints_pairs;
    kpoints_bands_path_grid=b.kpoints_bands_path_grid;
    enthalpy_cell=b.enthalpy_cell;enthalpy_atom=b.enthalpy_atom;
    enthalpy_formation_cell=b.enthalpy_formation_cell;enthalpy_formation_atom=b.enthalpy_formation_atom;
    enthalpy_formation_cce_300K_cell=b.enthalpy_formation_cce_300K_cell;enthalpy_formation_cce_300K_atom=b.enthalpy_formation_cce_300K_atom;  //CO20200624
    enthalpy_formation_cce_0K_cell=b.enthalpy_formation_cce_0K_cell;enthalpy_formation_cce_0K_atom=b.enthalpy_formation_cce_0K_atom;  //CO20200624
    entropy_forming_ability=b.entropy_forming_ability;  //CO20200624
    entropic_temperature=b.entropic_temperature;
    files=b.files;vfiles.clear();for(uint i=0;i<b.vfiles.size();i++) vfiles.push_back(b.vfiles.at(i));
    files_LIB=b.files_LIB;vfiles_LIB.clear();for(uint i=0;i<b.vfiles_LIB.size();i++) vfiles_LIB.push_back(b.vfiles_LIB.at(i));
    files_RAW=b.files_RAW;vfiles_RAW.clear();for(uint i=0;i<b.vfiles_RAW.size();i++) vfiles_RAW.push_back(b.vfiles_RAW.at(i));
    files_WEB=b.files_WEB;vfiles_WEB.clear();for(uint i=0;i<b.vfiles_WEB.size();i++) vfiles_WEB.push_back(b.vfiles_WEB.at(i));
    forces=b.forces;vforces.clear();for(uint i=0;i<b.vforces.size();i++) vforces.push_back(b.vforces.at(i));
    Egap_type=b.Egap_type;
    geometry=b.geometry;vgeometry.clear();for(uint i=0;i<b.vgeometry.size();i++) vgeometry.push_back(b.vgeometry.at(i));
    geometry_orig=b.geometry_orig;vgeometry_orig.clear();for(uint i=0;i<b.vgeometry_orig.size();i++) vgeometry_orig.push_back(b.vgeometry_orig.at(i)); //DX20190124 - add original crystal info
    lattice_system_orig=b.lattice_system_orig;lattice_variation_orig=b.lattice_variation_orig;
    lattice_system_relax=b.lattice_system_relax;lattice_variation_relax=b.lattice_variation_relax;
    ldau_TLUJ=b.ldau_TLUJ;
    vLDAU=b.vLDAU;  //ME20190129
    natoms=b.natoms;
    natoms_orig=b.natoms_orig; //DX20190124 - add original crystal info
    nbondxx=b.nbondxx;vnbondxx.clear();for(uint i=0;i<b.vnbondxx.size();i++) vnbondxx.push_back(b.vnbondxx.at(i));
    nspecies=b.nspecies;
    Pearson_symbol_orig=b.Pearson_symbol_orig;Pearson_symbol_relax=b.Pearson_symbol_relax;
    positions_cartesian=b.positions_cartesian;vpositions_cartesian.clear();for(uint i=0;i<b.vpositions_cartesian.size();i++) vpositions_cartesian.push_back(b.vpositions_cartesian.at(i));
    positions_fractional=b.positions_fractional;vpositions_fractional.clear();for(uint i=0;i<b.vpositions_fractional.size();i++) vpositions_fractional.push_back(b.vpositions_fractional.at(i));
    pressure=b.pressure;
    stress_tensor=b.stress_tensor;vstress_tensor.clear();for(uint i=0;i<b.vstress_tensor.size();i++) vstress_tensor.push_back(b.vstress_tensor.at(i));
    pressure_residual=b.pressure_residual;
    Pulay_stress=b.Pulay_stress;
    prototype=b.prototype;
    PV_cell=b.PV_cell;PV_atom=b.PV_atom;
    scintillation_attenuation_length=b.scintillation_attenuation_length;
    sg=b.sg;sg2=b.sg2;vsg.clear();for(uint i=0;i<b.vsg.size();i++){vsg.push_back(b.vsg[i]);} vsg2.clear();for(uint i=0;i<b.vsg2.size();i++){vsg2.push_back(b.vsg2[i]);}  //CO20171202
    spacegroup_orig=b.spacegroup_orig;spacegroup_relax=b.spacegroup_relax;
    species=b.species;vspecies.clear();for(uint i=0;i<b.vspecies.size();i++) vspecies.push_back(b.vspecies.at(i));
    species_pp=b.species_pp;vspecies_pp.clear();for(uint i=0;i<b.vspecies_pp.size();i++) vspecies_pp.push_back(b.vspecies_pp.at(i));
    species_pp_version=b.species_pp_version;vspecies_pp_version.clear();for(uint i=0;i<b.vspecies_pp_version.size();i++) vspecies_pp_version.push_back(b.vspecies_pp_version.at(i));
    species_pp_ZVAL=b.species_pp_ZVAL;vspecies_pp_ZVAL.clear();for(uint i=0;i<b.vspecies_pp_ZVAL.size();i++) vspecies_pp_ZVAL.push_back(b.vspecies_pp_ZVAL.at(i));
    species_pp_AUID=b.species_pp_AUID;vspecies_pp_AUID.clear();for(uint i=0;i<b.vspecies_pp_AUID.size();i++) vspecies_pp_AUID.push_back(b.vspecies_pp_AUID.at(i));
    METAGGA=b.METAGGA;
    spin_cell=b.spin_cell;spin_atom=b.spin_atom;
    spinD=b.spinD;vspinD.clear();for(uint i=0;i<b.vspinD.size();i++) vspinD.push_back(b.vspinD.at(i));
    spinD_magmom_orig=b.spinD_magmom_orig;vspinD_magmom_orig.clear();for(uint i=0;i<b.vspinD_magmom_orig.size();i++) vspinD_magmom_orig.push_back(b.vspinD_magmom_orig.at(i));
    spinF=b.spinF;
    sponsor=b.sponsor;vsponsor.clear();for(uint i=0;i<b.vsponsor.size();i++) vsponsor.push_back(b.vsponsor.at(i));
    stoichiometry=b.stoichiometry;vstoichiometry.clear();for(uint i=0;i<b.vstoichiometry.size();i++) vstoichiometry.push_back(b.vstoichiometry.at(i));
    valence_cell_std=b.valence_cell_std;valence_cell_iupac=b.valence_cell_iupac;
    volume_cell=b.volume_cell;volume_atom=b.volume_atom;
    volume_cell_orig=b.volume_cell_orig;volume_atom_orig=b.volume_atom_orig; //DX20190124 - add original crystal info
    //DX20190124 - added original symmetry info - START
    // SYMMETRY
    crystal_family_orig=b.crystal_family_orig;
    crystal_system_orig=b.crystal_system_orig;
    crystal_class_orig=b.crystal_class_orig;
    point_group_Hermann_Mauguin_orig=b.point_group_Hermann_Mauguin_orig;
    point_group_Schoenflies_orig=b.point_group_Schoenflies_orig;
    point_group_orbifold_orig=b.point_group_orbifold_orig;
    point_group_type_orig=b.point_group_type_orig;
    point_group_order_orig=b.point_group_order_orig;
    point_group_structure_orig=b.point_group_structure_orig;
    Bravais_lattice_lattice_type_orig=b.Bravais_lattice_lattice_type_orig;
    Bravais_lattice_lattice_variation_type_orig=b.Bravais_lattice_lattice_variation_type_orig;
    Bravais_lattice_lattice_system_orig=b.Bravais_lattice_lattice_system_orig;
    Bravais_superlattice_lattice_type_orig=b.Bravais_superlattice_lattice_type_orig;
    Bravais_superlattice_lattice_variation_type_orig=b.Bravais_superlattice_lattice_variation_type_orig;
    Bravais_superlattice_lattice_system_orig=b.Bravais_superlattice_lattice_system_orig;
    Pearson_symbol_superlattice_orig=b.Pearson_symbol_superlattice_orig;
    reciprocal_geometry_orig=b.reciprocal_geometry_orig;vreciprocal_geometry_orig.clear();for(uint i=0;i<b.vreciprocal_geometry_orig.size();i++) vreciprocal_geometry_orig.push_back(b.vreciprocal_geometry_orig.at(i));
    reciprocal_volume_cell_orig=b.reciprocal_volume_cell_orig;
    reciprocal_lattice_type_orig=b.reciprocal_lattice_type_orig;
    reciprocal_lattice_variation_type_orig=b.reciprocal_lattice_variation_type_orig;
    Wyckoff_letters_orig=b.Wyckoff_letters_orig;
    Wyckoff_multiplicities_orig=b.Wyckoff_multiplicities_orig;
    Wyckoff_site_symmetries_orig=b.Wyckoff_site_symmetries_orig;
    //DX20190124 - added original symmetry info - END
    //DX20180823 - added more symmetry info - START
    // SYMMETRY
    crystal_family=b.crystal_family;
    crystal_system=b.crystal_system;
    crystal_class=b.crystal_class;
    point_group_Hermann_Mauguin=b.point_group_Hermann_Mauguin;
    point_group_Schoenflies=b.point_group_Schoenflies;
    point_group_orbifold=b.point_group_orbifold;
    point_group_type=b.point_group_type;
    point_group_order=b.point_group_order;
    point_group_structure=b.point_group_structure;
    Bravais_lattice_lattice_type=b.Bravais_lattice_lattice_type;
    Bravais_lattice_lattice_variation_type=b.Bravais_lattice_lattice_variation_type;
    Bravais_lattice_lattice_system=b.Bravais_lattice_lattice_system;
    Bravais_superlattice_lattice_type=b.Bravais_superlattice_lattice_type;
    Bravais_superlattice_lattice_variation_type=b.Bravais_superlattice_lattice_variation_type;
    Bravais_superlattice_lattice_system=b.Bravais_superlattice_lattice_system;
    Pearson_symbol_superlattice=b.Pearson_symbol_superlattice;
    reciprocal_geometry_relax=b.reciprocal_geometry_relax;vreciprocal_geometry_relax.clear();for(uint i=0;i<b.vreciprocal_geometry_relax.size();i++) vreciprocal_geometry_relax.push_back(b.vreciprocal_geometry_relax.at(i));  //CO20220719 _relax
    reciprocal_volume_cell=b.reciprocal_volume_cell; //DX20190124 - fix typo, add reciprocal
    reciprocal_lattice_type=b.reciprocal_lattice_type;
    reciprocal_lattice_variation_type=b.reciprocal_lattice_variation_type;
    Wyckoff_letters=b.Wyckoff_letters;
    Wyckoff_multiplicities=b.Wyckoff_multiplicities;
    Wyckoff_site_symmetries=b.Wyckoff_site_symmetries;
    //DX20180823 - added more symmetry info - END
    //DX20190209 - added anrl info - START
    aflow_prototype_label_orig=b.aflow_prototype_label_orig;
    aflow_prototype_params_list_orig=b.aflow_prototype_params_list_orig;
    aflow_prototype_params_values_orig=b.aflow_prototype_params_values_orig;
    aflow_prototype_label_relax=b.aflow_prototype_label_relax;
    aflow_prototype_params_list_relax=b.aflow_prototype_params_list_relax;
    aflow_prototype_params_values_relax=b.aflow_prototype_params_values_relax;
    //DX20190209 - added anrl info - END
    pocc_parameters=b.pocc_parameters;  //CO20200731
    // AGL/AEL
    agl_thermal_conductivity_300K=b.agl_thermal_conductivity_300K;
    agl_debye=b.agl_debye;
    agl_acoustic_debye=b.agl_acoustic_debye;
    agl_gruneisen=b.agl_gruneisen;
    agl_heat_capacity_Cv_300K=b.agl_heat_capacity_Cv_300K;
    agl_heat_capacity_Cp_300K=b.agl_heat_capacity_Cp_300K;
    agl_thermal_expansion_300K=b.agl_thermal_expansion_300K;
    agl_bulk_modulus_static_300K=b.agl_bulk_modulus_static_300K;
    agl_bulk_modulus_isothermal_300K=b.agl_bulk_modulus_isothermal_300K;
    agl_poisson_ratio_source=b.agl_poisson_ratio_source; //CT20181212
    agl_vibrational_free_energy_300K_cell=b.agl_vibrational_free_energy_300K_cell; //CT20181212
    agl_vibrational_free_energy_300K_atom=b.agl_vibrational_free_energy_300K_atom; //CT20181212
    agl_vibrational_entropy_300K_cell=b.agl_vibrational_entropy_300K_cell; //CT20181212
    agl_vibrational_entropy_300K_atom=b.agl_vibrational_entropy_300K_atom; //CT20181212
    ael_poisson_ratio=b.ael_poisson_ratio;
    ael_bulk_modulus_voigt=b.ael_bulk_modulus_voigt;
    ael_bulk_modulus_reuss=b.ael_bulk_modulus_reuss;
    ael_shear_modulus_voigt=b.ael_shear_modulus_voigt;
    ael_shear_modulus_reuss=b.ael_shear_modulus_reuss;
    ael_bulk_modulus_vrh=b.ael_bulk_modulus_vrh;
    ael_shear_modulus_vrh=b.ael_shear_modulus_vrh;
    ael_elastic_anisotropy=b.ael_elastic_anisotropy; //CO20181129
    ael_youngs_modulus_vrh=b.ael_youngs_modulus_vrh; //CT20181212
    ael_speed_sound_transverse=b.ael_speed_sound_transverse; //CT20181212
    ael_speed_sound_longitudinal=b.ael_speed_sound_longitudinal; //CT20181212
    ael_speed_sound_average=b.ael_speed_sound_average; //CT20181212
    ael_pughs_modulus_ratio=b.ael_pughs_modulus_ratio; //CT20181212
    ael_debye_temperature=b.ael_debye_temperature; //CT20181212
    ael_applied_pressure=b.ael_applied_pressure; //CT20181212
    ael_average_external_pressure=b.ael_average_external_pressure; //CT20181212
    ael_stiffness_tensor = b.ael_stiffness_tensor;  //ME20191105
    ael_compliance_tensor = b.ael_compliance_tensor;  //ME20191105
    // APL // ME20210927
    energy_free_vibrational_cell_apl_300K = b.energy_free_vibrational_cell_apl_300K;
    energy_free_vibrational_atom_apl_300K = b.energy_free_vibrational_atom_apl_300K;
    entropy_vibrational_cell_apl_300K = b.entropy_vibrational_cell_apl_300K;
    entropy_vibrational_atom_apl_300K = b.entropy_vibrational_atom_apl_300K;
    energy_internal_vibrational_cell_apl_300K = b.energy_internal_vibrational_cell_apl_300K;
    energy_internal_vibrational_atom_apl_300K = b.energy_internal_vibrational_atom_apl_300K;
    energy_zero_point_cell_apl = b.energy_zero_point_cell_apl;
    energy_zero_point_atom_apl = b.energy_zero_point_atom_apl;
    heat_capacity_Cv_cell_apl_300K = b.heat_capacity_Cv_cell_apl_300K;
    heat_capacity_Cv_atom_apl_300K = b.heat_capacity_Cv_atom_apl_300K;
    // QHA
    gruneisen_qha = b.gruneisen_qha; //AS20200901
    gruneisen_qha_300K = b.gruneisen_qha_300K; //AS20200903
    thermal_expansion_qha_300K = b.thermal_expansion_qha_300K; //AS20200901
    modulus_bulk_qha_300K = b.modulus_bulk_qha_300K; //AS20200901
    modulus_bulk_derivative_pressure_qha_300K = b.modulus_bulk_derivative_pressure_qha_300K; //AS20201008
    heat_capacity_Cv_atom_qha_300K = b.heat_capacity_Cv_atom_qha_300K; //AS20201008
    heat_capacity_Cv_cell_qha_300K = b.heat_capacity_Cv_cell_qha_300K; //AS20201207
    heat_capacity_Cp_atom_qha_300K = b.heat_capacity_Cp_atom_qha_300K; //AS20201008
    heat_capacity_Cp_cell_qha_300K = b.heat_capacity_Cp_cell_qha_300K; //AS20201207
    volume_atom_qha_300K = b.volume_atom_qha_300K; //AS20201008
    energy_free_atom_qha_300K = b.energy_free_atom_qha_300K; //AS20201008
    energy_free_cell_qha_300K = b.energy_free_cell_qha_300K; //AS20201207
    // BADER
    bader_net_charges=b.bader_net_charges;vbader_net_charges.clear();for(uint i=0;i<b.vbader_net_charges.size();i++) vbader_net_charges.push_back(b.vbader_net_charges.at(i));
    bader_atomic_volumes=b.bader_atomic_volumes;vbader_atomic_volumes.clear();for(uint i=0;i<b.vbader_atomic_volumes.size();i++) vbader_atomic_volumes.push_back(b.vbader_atomic_volumes.at(i));
    // legacy
    server=b.server;
    vserver.clear();for(uint i=0;i<b.vserver.size();i++) vserver.push_back(b.vserver.at(i));
    vserverdir.clear();for(uint i=0;i<b.vserverdir.size();i++) vserverdir.push_back(b.vserverdir.at(i));
    icsd=b.icsd;
    stoich=b.stoich;vstoich.clear();for(uint i=0;i<b.vstoich.size();i++) vstoich.push_back(b.vstoich.at(i));
    // apennsy
    structure_name=b.structure_name;  // apennsy
    structure_description=b.structure_description;  // apennsy
    distance_gnd=b.distance_gnd;  // apennsy
    distance_tie=b.distance_tie;  // apennsy
    pureA=b.pureA;pureB=b.pureB;  // apennsy
    fcc=b.fcc;bcc=b.bcc;hcp=b.hcp;  // apennsy
    stoich_a=b.stoich_a;stoich_b=b.stoich_b;  // apennsy
    bond_aa=b.bond_aa;bond_ab=b.bond_ab;bond_bb=b.bond_bb;  // apennsy
    vNsgroup.clear();for(uint i=0;i<b.vNsgroup.size();i++) vNsgroup.push_back(b.vNsgroup.at(i));  // apennsy
    vsgroup.clear();for(uint i=0;i<b.vsgroup.size();i++) vsgroup.push_back(b.vsgroup.at(i));  // apennsy
    vstr.clear();for(uint i=0;i<b.vstr.size();i++) vstr.push_back(b.vstr.at(i));  // apennsy
    // details from EntryLoader //HE20220913
    el_source_type=b.el_source_type;
    el_source=b.el_source;
  }


  const _aflowlib_entry& _aflowlib_entry::operator=(const _aflowlib_entry& b) {  // operator= PUBLIC
    if(this!=&b) {free(); copy(b);}
    return *this;
  }

  _aflowlib_entry::_aflowlib_entry(const _aflowlib_entry& b) { // copy PUBLIC
    //  free();*this=b;
    copy(b);
  }

  void _aflowlib_entry::free() { // free PRIVATE
    //[CO20200829 - OBSOLETE]ventry.clear();
    //[CO20200829 - OBSOLETE]vauid.clear();
    //[CO20200829 - OBSOLETE]vaurl.clear();
    //[CO20200829 - OBSOLETE]vaflowlib_entries.clear();
    //[CO20200829 - OBSOLETE]vkeywords.clear();
    //[CO20200829 - OBSOLETE]vauthor.clear();
    //[CO20200829 - OBSOLETE]vcorresponding.clear();
    //[CO20200829 - OBSOLETE]vloop.clear();
    //[CO20200829 - OBSOLETE]vcomposition.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vdft_type.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vfiles.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vfiles_LIB.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vfiles_RAW.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vfiles_WEB.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vforces.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vgeometry.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vgeometry_orig.clear(); // clear all vectors //DX20190124 - add original crystal info
    //[CO20200829 - OBSOLETE]vstress_tensor.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vnbondxx.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vpositions_cartesian.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vpositions_fractional.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vspecies.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vspecies_pp.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vspecies_pp_version.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vspecies_pp_ZVAL.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vspecies_pp_AUID.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vspinD.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vspinD_magmom_orig.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vsponsor.clear();
    //[CO20200829 - OBSOLETE]vstoichiometry.clear(); // clear all vectors
    //[CO20200829 - OBSOLETE]vreciprocal_geometry_relax.clear(); // clear all vectors //DX20180824 - added reciprocal lattice parameters  //CO20220719 _relax
    //[CO20200829 - OBSOLETE]// BADER
    //[CO20200829 - OBSOLETE]vbader_net_charges.clear();
    //[CO20200829 - OBSOLETE]vbader_atomic_volumes.clear();
    //[CO20200829 - OBSOLETE]// legacy
    //[CO20200829 - OBSOLETE]vserver.clear();  // clear all vectors
    //[CO20200829 - OBSOLETE]for(uint i=0;i<vserverdir.size();i++)
    //[CO20200829 - OBSOLETE]  vserverdir.at(i).clear();
    //[CO20200829 - OBSOLETE]vserverdir.clear();  // clear all vectors
    //[CO20200829 - OBSOLETE]vstoich.clear(); // clear all vectors

    entry.clear();ventry.clear();
    auid.clear();
    vauid.clear();vauid.clear();
    aurl.clear();vaurl.clear();
    system_name.clear();
    keywords.clear();vkeywords.clear();
    aflowlib_date.clear();vaflowlib_date.clear(); //CO20200624 - adding LOCK date
    aflowlib_version.clear();
    aflowlib_entries.clear();vaflowlib_entries.clear();
    aflowlib_entries_number=0;
    aflow_version.clear();
    catalog.clear();
    data_api="aapi1.2"; // new version of the API
    data_source="aflowlib";
    data_language="";
    error_status.clear();
    author.clear();vauthor.clear();
    calculation_cores=1;
    calculation_memory=AUROSTD_NAN;
    calculation_time=AUROSTD_NAN;
    corresponding.clear();vcorresponding.clear();
    loop.clear();vloop.clear();
    node_CPU_Cores=AUROSTD_NAN;node_CPU_MHz=AUROSTD_NAN;node_CPU_Model.clear();node_RAM_GB=AUROSTD_NAN;
    Bravais_lattice_orig.clear();Bravais_lattice_relax.clear();
    code.clear();
    composition.clear();vcomposition.clear();
    compound.clear();
    density=AUROSTD_NAN;
    density_orig=AUROSTD_NAN; //DX20190124 - add original crystal info
    dft_type.clear();vdft_type.clear();
    eentropy_cell=AUROSTD_NAN;eentropy_atom=AUROSTD_NAN;
    Egap=AUROSTD_NAN;Egap_fit=AUROSTD_NAN;
    Egap_type.clear();
    energy_cell=AUROSTD_NAN;energy_atom=AUROSTD_NAN;energy_atom_relax1=AUROSTD_NAN;
    energy_cutoff=AUROSTD_NAN;
    delta_electronic_energy_convergence=AUROSTD_NAN;
    delta_electronic_energy_threshold=AUROSTD_NAN;
    nkpoints=0;nkpoints_irreducible=0;kppra=0;
    kpoints.clear();
    kpoints_nnn_relax.clear();kpoints_nnn_static.clear();
    kpoints_pairs.clear();
    kpoints_bands_path_grid=0;
    enthalpy_cell=AUROSTD_NAN;enthalpy_atom=AUROSTD_NAN;
    enthalpy_formation_cell=AUROSTD_NAN;enthalpy_formation_atom=AUROSTD_NAN;
    enthalpy_formation_cce_300K_cell=AUROSTD_NAN;enthalpy_formation_cce_300K_atom=AUROSTD_NAN;  //CO20200624
    enthalpy_formation_cce_0K_cell=AUROSTD_NAN;enthalpy_formation_cce_0K_atom=AUROSTD_NAN;  //CO20200624
    entropy_forming_ability=AUROSTD_NAN;  //CO20200624
    entropic_temperature=AUROSTD_NAN;
    files.clear();vfiles.clear();
    files_LIB.clear();vfiles_LIB.clear();
    files_RAW.clear();vfiles_RAW.clear();
    files_WEB.clear();vfiles_WEB.clear();
    forces.clear();vforces.clear();
    geometry.clear();vgeometry.clear();
    geometry_orig.clear();vgeometry_orig.clear(); //DX20190124 - add original crystal info
    lattice_system_orig.clear();lattice_variation_orig.clear();lattice_system_relax.clear();lattice_variation_relax.clear();
    ldau_TLUJ.clear();
    for(uint i=0;i<vLDAU.size();i++){vLDAU[i].clear();}vLDAU.clear(); vLDAU.resize(4);  //ME20190129  //CO20210104 clear()
    natoms=AUROSTD_NAN;
    natoms_orig=AUROSTD_NAN; //DX20190124 - add original crystal info
    nbondxx.clear();vnbondxx.clear();
    nspecies=AUROSTD_NAN;
    Pearson_symbol_orig.clear();Pearson_symbol_relax.clear();
    positions_cartesian.clear();vpositions_cartesian.clear();
    positions_fractional.clear();vpositions_fractional.clear();
    pressure=AUROSTD_NAN;
    stress_tensor.clear();vstress_tensor.clear();
    pressure_residual=AUROSTD_NAN;
    Pulay_stress=AUROSTD_NAN;
    prototype.clear();
    PV_cell=AUROSTD_NAN;PV_atom=AUROSTD_NAN;
    scintillation_attenuation_length=AUROSTD_NAN;
    sg.clear();sg2.clear();vsg.clear();vsg2.clear();  //CO20171202
    spacegroup_orig=AUROSTD_NAN;spacegroup_relax=AUROSTD_NAN; //CO20201111
    species.clear();vspecies.clear();
    species_pp.clear();vspecies_pp.clear();
    species_pp_version.clear();vspecies_pp_version.clear();
    species_pp_ZVAL.clear();vspecies_pp_ZVAL.clear();
    species_pp_AUID.clear();vspecies_pp_AUID.clear();
    METAGGA.clear();
    spin_cell=AUROSTD_NAN;spin_atom=AUROSTD_NAN;
    spinD.clear();vspinD.clear();
    spinD_magmom_orig.clear();vspinD_magmom_orig.clear();
    spinF=AUROSTD_NAN;
    sponsor.clear();vsponsor.clear();
    stoichiometry.clear();vstoichiometry.clear();
    valence_cell_std=AUROSTD_NAN;valence_cell_iupac=AUROSTD_NAN;
    volume_cell=AUROSTD_NAN;volume_atom=AUROSTD_NAN;
    volume_cell_orig=AUROSTD_NAN;volume_atom_orig=AUROSTD_NAN; //DX20190124 - add original crystal info
    //DX20190124 - added original symmetry info - START
    // SYMMETRY
    crystal_family_orig.clear();
    crystal_system_orig.clear();
    crystal_class_orig.clear();
    point_group_Hermann_Mauguin_orig.clear();
    point_group_Schoenflies_orig.clear();
    point_group_orbifold_orig.clear();
    point_group_type_orig.clear();
    point_group_order_orig=AUROSTD_NAN;
    point_group_structure_orig.clear();
    Bravais_lattice_lattice_type_orig.clear();
    Bravais_lattice_lattice_variation_type_orig.clear();
    Bravais_lattice_lattice_system_orig.clear();
    Bravais_superlattice_lattice_type_orig.clear();
    Bravais_superlattice_lattice_variation_type_orig.clear();
    Bravais_superlattice_lattice_system_orig.clear();
    Pearson_symbol_superlattice_orig.clear();
    reciprocal_geometry_orig.clear();vreciprocal_geometry_orig.clear();
    reciprocal_volume_cell_orig=AUROSTD_NAN;
    reciprocal_lattice_type_orig.clear();
    reciprocal_lattice_variation_type_orig.clear();
    Wyckoff_letters_orig.clear();
    Wyckoff_multiplicities_orig.clear();
    Wyckoff_site_symmetries_orig.clear();
    //DX20190124 - added original symmetry info - END
    //DX20180823 - added more symmetry info - START
    // SYMMETRY
    crystal_family.clear();
    crystal_system.clear();
    crystal_class.clear();
    point_group_Hermann_Mauguin.clear();
    point_group_Schoenflies.clear();
    point_group_orbifold.clear();
    point_group_type.clear();
    point_group_order=AUROSTD_NAN;
    point_group_structure.clear();
    Bravais_lattice_lattice_type.clear();
    Bravais_lattice_lattice_variation_type.clear();
    Bravais_lattice_lattice_system.clear();
    Bravais_superlattice_lattice_type.clear();
    Bravais_superlattice_lattice_variation_type.clear();
    Bravais_superlattice_lattice_system.clear();
    Pearson_symbol_superlattice.clear();
    reciprocal_geometry_relax.clear();vreciprocal_geometry_relax.clear(); //CO20220719 _relax
    reciprocal_volume_cell=AUROSTD_NAN;
    reciprocal_lattice_type.clear();
    reciprocal_lattice_variation_type.clear();
    Wyckoff_letters.clear();
    Wyckoff_multiplicities.clear();
    Wyckoff_site_symmetries.clear();
    //DX20180823 - added more symmetry info - END
    //DX20190209 - added anrl info - START
    aflow_prototype_label_orig.clear();
    aflow_prototype_params_list_orig.clear();
    aflow_prototype_params_values_orig.clear();
    aflow_prototype_label_relax.clear();
    aflow_prototype_params_list_relax.clear();
    aflow_prototype_params_values_relax.clear();
    //DX20190209 - added anrl info - END
    pocc_parameters.clear(); //CO20200731
    // AGL/AEL
    agl_thermal_conductivity_300K=AUROSTD_NAN;
    agl_debye=AUROSTD_NAN;
    agl_acoustic_debye=AUROSTD_NAN;
    agl_gruneisen=AUROSTD_NAN;
    agl_heat_capacity_Cv_300K=AUROSTD_NAN;
    agl_heat_capacity_Cp_300K=AUROSTD_NAN;
    agl_thermal_expansion_300K=AUROSTD_NAN;
    agl_bulk_modulus_static_300K=AUROSTD_NAN;
    agl_bulk_modulus_isothermal_300K=AUROSTD_NAN;
    agl_poisson_ratio_source.clear(); //CT20181212
    agl_vibrational_free_energy_300K_cell=AUROSTD_NAN; //CT20181212
    agl_vibrational_free_energy_300K_atom=AUROSTD_NAN; //CT20181212
    agl_vibrational_entropy_300K_cell=AUROSTD_NAN; //CT20181212
    agl_vibrational_entropy_300K_atom=AUROSTD_NAN; //CT20181212
    ael_poisson_ratio=AUROSTD_NAN;
    ael_bulk_modulus_voigt=AUROSTD_NAN;
    ael_bulk_modulus_reuss=AUROSTD_NAN;
    ael_shear_modulus_voigt=AUROSTD_NAN;
    ael_shear_modulus_reuss=AUROSTD_NAN;
    ael_bulk_modulus_vrh=AUROSTD_NAN;
    ael_shear_modulus_vrh=AUROSTD_NAN;
    ael_elastic_anisotropy=AUROSTD_NAN; //CO20181129
    ael_youngs_modulus_vrh=AUROSTD_NAN; //CT20181212
    ael_speed_sound_transverse=AUROSTD_NAN; //CT20181212
    ael_speed_sound_longitudinal=AUROSTD_NAN; //CT20181212
    ael_speed_sound_average=AUROSTD_NAN; //CT20181212
    ael_pughs_modulus_ratio=AUROSTD_NAN; //CT20181212
    ael_debye_temperature=AUROSTD_NAN; //CT20181212
    ael_applied_pressure=AUROSTD_NAN; //CT20181212
    ael_average_external_pressure=AUROSTD_NAN; //CT20181212
    ael_stiffness_tensor.clear();  //ME20191105
    ael_compliance_tensor.clear();  //ME20191105
    // APL // ME20210927
    energy_free_vibrational_cell_apl_300K = AUROSTD_NAN;
    energy_free_vibrational_atom_apl_300K = AUROSTD_NAN;
    entropy_vibrational_cell_apl_300K = AUROSTD_NAN;
    entropy_vibrational_atom_apl_300K = AUROSTD_NAN;
    energy_internal_vibrational_cell_apl_300K = AUROSTD_NAN;
    energy_internal_vibrational_atom_apl_300K = AUROSTD_NAN;
    energy_zero_point_cell_apl = AUROSTD_NAN;
    energy_zero_point_atom_apl = AUROSTD_NAN;
    heat_capacity_Cv_cell_apl_300K = AUROSTD_NAN;
    heat_capacity_Cv_atom_apl_300K = AUROSTD_NAN;
    // QHA
    gruneisen_qha = AUROSTD_NAN;//AS20200901
    gruneisen_qha_300K = AUROSTD_NAN;//AS20200903
    thermal_expansion_qha_300K = AUROSTD_NAN;//AS20200901
    modulus_bulk_qha_300K = AUROSTD_NAN;//AS20200901
    modulus_bulk_derivative_pressure_qha_300K = AUROSTD_NAN;//AS20201008
    heat_capacity_Cv_atom_qha_300K = AUROSTD_NAN;//AS20201008
    heat_capacity_Cv_cell_qha_300K = AUROSTD_NAN;//AS20201207
    heat_capacity_Cp_atom_qha_300K = AUROSTD_NAN;//AS20201008
    heat_capacity_Cp_cell_qha_300K = AUROSTD_NAN;//AS20201207
    volume_atom_qha_300K = AUROSTD_NAN;//AS20201008
    energy_free_atom_qha_300K = AUROSTD_NAN;//AS20201008
    energy_free_cell_qha_300K = AUROSTD_NAN;//AS20201207
    // BADER
    bader_net_charges.clear();vbader_net_charges.clear();
    bader_atomic_volumes.clear();vbader_atomic_volumes.clear();
    // legacy
    server.clear();vserver.clear();vserverdir.clear();
    icsd.clear();
    stoich.clear();vstoich.clear();
    // apennsy
    structure_name.clear();  // apennsy
    structure_description.clear();  // apennsy
    distance_gnd=AUROSTD_NAN;  // apennsy
    distance_tie=AUROSTD_NAN;  // apennsy
    pureA=FALSE;pureB=FALSE;  // apennsy
    fcc=FALSE;bcc=FALSE;hcp=FALSE;  // apennsy
    stoich_a=AUROSTD_NAN;stoich_b=AUROSTD_NAN;  // apennsy
    bond_aa=AUROSTD_NAN;bond_ab=AUROSTD_NAN;bond_bb=AUROSTD_NAN;  // apennsy
    vNsgroup.clear();  // apennsy
    vsgroup.clear();  // apennsy
    vstr.clear();  // apennsy
    // details from EntryLoader //HE20220913
    el_source_type.clear();
    el_source.clear();
  } 

  void _aflowlib_entry::clear() {  // clear PRIVATE
    //[CO20201220 - too slow]_aflowlib_entry _temp;
    //[CO20201220 - too slow]copy(_temp);
    free();
  }

  _aflowlib_entry::_aflowlib_entry(const string& file) { // constructur from file
    stringstream oss;
    if(!aurostd::FileExist(file)) { //SC20190813
      string message = DEFAULT_FILE_AFLOWLIB_ENTRY_OUT + " not found =" + file;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    string entry;
    aurostd::efile2string(file,entry);
    Load(entry,oss);
  }

  /// @brief Set the value of a member
  /// @param keyword name of the member
  /// @param content value to set
  void _aflowlib_entry::Set(const std::string &keyword, const std::string &content) {
    return SetByHash(aurostd::crc64(keyword), content);
  }

  ///@brief Set the value of members
  ///@param key_hash aurostd::crc64 hash of the member name
  ///@param content used to set the member value
  void _aflowlib_entry::SetByHash(const uint64_t key_hash, const std::string &content) {
    string function = "aflowlib::_aflowlib_entry::SetByHash()";  //ME20191119
    if (content.empty()) return;  //CO20180319
    if (content == "null") return; //CO20180319 - aflux integration!
    vector<string> stokens;
    aurostd::string2tokens(content, stokens, ",");
    //CO20180409 - added the else if's for speed, no need to go through more checks than necessary
    //HE20220404 - replace if else tree with hash switch/case - this replaces 26565 string comparisons with 230 jump table lookups
    switch (key_hash) {
      case (aurostd::ctcrc64("auid")): {
        auid = content;
        vauid.clear();
        aflowlib::auid2vauid(auid, vauid); // create VAUID
        } break;
      case (aurostd::ctcrc64("aurl")): {
        aurl = content;
        aurostd::string2tokens(content, stokens, ":");
        for (uint j = 0; j < stokens.size(); j++) vaurl.push_back(stokens.at(j));
        } break;
      //else if(keyword=="title") {title=content;}  //ME20190129 // OBSOLETE ME20200829 - not used anymore
      case (aurostd::ctcrc64("keywords")): {
        keywords = content;
        aurostd::string2tokens(content, stokens, ",");
        for (uint j = 0; j < stokens.size(); j++) vkeywords.push_back(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("aflowlib_date")): {
        aflowlib_date = content;
        aurostd::string2tokens(content, stokens, ",");
        for (uint j = 0; j < stokens.size(); j++) vaflowlib_date.push_back(stokens.at(j));
        } break; //CO20200624 - adding LOCK date
      case (aurostd::ctcrc64("aflowlib_version")): { aflowlib_version = content; } break;
      case (aurostd::ctcrc64("aflowlib_entries")): {
        aflowlib_entries = content;
        aurostd::string2tokens(content, stokens, ",");
        for (uint j = 0; j < stokens.size(); j++) vaflowlib_entries.push_back(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("aflowlib_entries_number")): { aflowlib_entries_number = aurostd::string2utype<int>(content); } break;
      case (aurostd::ctcrc64("aflow_version")): { aflow_version = content; } break;
      case (aurostd::ctcrc64("catalog")): { catalog = content; } break;
      case (aurostd::ctcrc64("data_api")): { data_api = content; } break;
      case (aurostd::ctcrc64("data_source")): { data_source = content; } break;
      case (aurostd::ctcrc64("data_language")): { data_language = content; } break;
      case (aurostd::ctcrc64("error_status")): { error_status = content; } break;
      case (aurostd::ctcrc64("author")): {
        author = content;
        aurostd::string2tokens(content, stokens, ",");
        for (uint j = 0; j < stokens.size(); j++) vauthor.push_back(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("calculation_cores")): { calculation_cores = aurostd::string2utype<int>(content); } break;
      case (aurostd::ctcrc64("calculation_memory")): { calculation_memory = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("calculation_time")): { calculation_time = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("corresponding")): {
        corresponding = content;
        aurostd::string2tokens(content, stokens, ",");
        for (uint j = 0; j < stokens.size(); j++) vcorresponding.push_back(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("loop")): { vloop.push_back(content); } break;  // CHECK THIS OUT IN THE FITURE
      case (aurostd::ctcrc64("node_CPU_Cores")): { node_CPU_Cores = aurostd::string2utype<int>(content); } break;
      case (aurostd::ctcrc64("node_CPU_MHz")): { node_CPU_MHz = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("node_CPU_Model")): { node_CPU_Model = content; } break;
      case (aurostd::ctcrc64("node_RAM_GB")): { node_RAM_GB = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("Bravais_lattice_orig")): { Bravais_lattice_orig = content; } break;
      case (aurostd::ctcrc64("Bravais_lattice_relax")): { Bravais_lattice_relax = content; } break;
      case (aurostd::ctcrc64("code")): { code = content; } break;
      case (aurostd::ctcrc64("composition")): {
        composition = content;
        for (uint j = 0; j < stokens.size(); j++) vcomposition.push_back(aurostd::string2utype<double>(stokens.at(j)));
        } break;
      case (aurostd::ctcrc64("compound")): { compound = content; } break;
      case (aurostd::ctcrc64("density")): { density = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("density_orig")): { density_orig = aurostd::string2utype<double>(content); } break; //DX20190124 - add original crystal info
      case (aurostd::ctcrc64("dft_type")): {
        dft_type = content;
        aurostd::string2tokens(content, stokens, ",");
        for (uint j = 0; j < stokens.size(); j++) vdft_type.push_back(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("eentropy_cell")): { eentropy_cell = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("eentropy_atom")): { eentropy_atom = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("Egap")): { Egap = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("Egap_fit")): { Egap_fit = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("energy_cell")): { energy_cell = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("energy_atom")): { energy_atom = aurostd::string2utype<double>(content); energy_atom_relax1 = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("energy_cutoff")): { energy_cutoff = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("delta_electronic_energy_convergence")): { delta_electronic_energy_convergence = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("delta_electronic_energy_threshold")): { delta_electronic_energy_threshold = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("nkpoints")): { nkpoints = aurostd::string2utype<uint>(content); } break;
      case (aurostd::ctcrc64("nkpoints_irreducible")): { nkpoints_irreducible = aurostd::string2utype<uint>(content); } break;
      case (aurostd::ctcrc64("kppra")): { kppra = aurostd::string2utype<uint>(content); } break;
      case (aurostd::ctcrc64("kpoints")): { kpoints = content; } break;
      case (aurostd::ctcrc64("kpoints_relax")): {
        vector<int> tokens;
        aurostd::string2tokens(content, tokens, ",");
        kpoints_nnn_relax = aurostd::vector2xvector(tokens);
        } break;  //ME20190129
      case (aurostd::ctcrc64("kpoints_static")): {
        vector<int> tokens;
        aurostd::string2tokens(content, tokens, ",");
        kpoints_nnn_static = aurostd::vector2xvector(tokens);
        } break;  //ME20190129
      case (aurostd::ctcrc64("kpoints_bands_path")): { aurostd::string2tokens(content, kpoints_pairs, ","); } break;  //ME20190129
      case (aurostd::ctcrc64("kpoints_bands_nkpts")): {kpoints_bands_path_grid = aurostd::string2utype<int>(content); } break;  //ME20190129
      case (aurostd::ctcrc64("enthalpy_cell")): { enthalpy_cell = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("enthalpy_atom")): { enthalpy_atom = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("enthalpy_formation_cell")): { enthalpy_formation_cell = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("enthalpy_formation_cce_300K_cell")): { enthalpy_formation_cce_300K_cell = aurostd::string2utype<double>(content); } break; //CO20200624
      case (aurostd::ctcrc64("enthalpy_formation_cce_0K_cell")): { enthalpy_formation_cce_0K_cell = aurostd::string2utype<double>(content); } break; //CO20200624
      case (aurostd::ctcrc64("enthalpy_formation_atom")): { enthalpy_formation_atom = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("enthalpy_formation_cce_300K_atom")): { enthalpy_formation_cce_300K_atom = aurostd::string2utype<double>(content); } break; //CO20200624
      case (aurostd::ctcrc64("enthalpy_formation_cce_0K_atom")): { enthalpy_formation_cce_0K_atom = aurostd::string2utype<double>(content); } break; //CO20200624
      case (aurostd::ctcrc64("entropy_forming_ability")): { entropy_forming_ability = aurostd::string2utype<double>(content); } break; //CO20200624
      case (aurostd::ctcrc64("entropic_temperature")): { entropic_temperature = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("files")): {
        files = content;
        for (uint j = 0; j < stokens.size(); j++) vfiles.push_back(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("files_LIB")): {
        files_LIB = content;
        for (uint j = 0; j < stokens.size(); j++) vfiles_LIB.push_back(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("files_RAW")): {
        files_RAW = content;
        for (uint j = 0; j < stokens.size(); j++) vfiles_RAW.push_back(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("files_WEB")): {
        files_WEB = content;
        for (uint j = 0; j < stokens.size(); j++) vfiles_WEB.push_back(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("forces")): {
        forces = content;
        for (uint j = 0; j < stokens.size(); j++) vforces.push_back(aurostd::string2utype<double>(stokens.at(j)));
        } break;  // FIX
      case (aurostd::ctcrc64("geometry")): {
        geometry = content;
        vgeometry = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; //HE20220404 compacted
        if (stokens.size() == 6)
          for (uint j = 0; j < stokens.size(); j++)
            vgeometry.at(j) = aurostd::string2utype<double>(stokens.at(j));
        } break;
      //DX20190124 - add original crystal info - START
      case (aurostd::ctcrc64("geometry_orig")): {
        geometry_orig = content;
        vgeometry_orig = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; //HE20220404 compacted
        if (stokens.size() == 6)
          for (uint j = 0; j < stokens.size(); j++)
            vgeometry_orig.at(j) = aurostd::string2utype<double>(stokens.at(j));
        } break;
      //DX20190124 - add original crystal info - END
      case (aurostd::ctcrc64("lattice_system_orig")): { lattice_system_orig = content; } break;
      case (aurostd::ctcrc64("lattice_variation_orig")): { lattice_variation_orig = content; } break;
      case (aurostd::ctcrc64("lattice_system_relax")): { lattice_system_relax = content; } break;
      case (aurostd::ctcrc64("lattice_variation_relax")): { lattice_variation_relax = content; } break;
      case (aurostd::ctcrc64("ldau_TLUJ")): { ldau_TLUJ = content; } break;
      case (aurostd::ctcrc64("ldau_type")): { vLDAU[0].push_back(aurostd::string2utype<double>(content)); } break;   //ME20190129
      case (aurostd::ctcrc64("ldau_l")): {
        for (uint j = 0; j < stokens.size(); j++)
          vLDAU[1].push_back(aurostd::string2utype<double>(stokens[j]));
        } break;   //ME20190129
      case (aurostd::ctcrc64("ldau_u")): {
        for (uint j = 0; j < stokens.size(); j++)
          vLDAU[2].push_back(aurostd::string2utype<double>(stokens[j]));
        } break;   //ME20190129
      case (aurostd::ctcrc64("ldau_j")): {
        for (uint j = 0; j < stokens.size(); j++)
          vLDAU[3].push_back(aurostd::string2utype<double>(stokens[j]));
        } break;   //ME20190129
      case (aurostd::ctcrc64("natoms")): { natoms = aurostd::string2utype<int>(content); } break;
      case (aurostd::ctcrc64("natoms_orig")): { natoms_orig = aurostd::string2utype<int>(content); } break;  //DX20190124 - add original crystal info
      case (aurostd::ctcrc64("nbondxx")): {
        nbondxx = content;
        for (uint j = 0; j < stokens.size(); j++) vnbondxx.push_back(aurostd::string2utype<double>(stokens.at(j)));
        } break;
      case (aurostd::ctcrc64("nspecies")): { nspecies = aurostd::string2utype<int>(content); } break;
      case (aurostd::ctcrc64("Pearson_symbol_orig")): { Pearson_symbol_orig = content; } break;
      case (aurostd::ctcrc64("Pearson_symbol_relax")): { Pearson_symbol_relax = content; } break;
      case (aurostd::ctcrc64("positions_cartesian")): {
        positions_cartesian = content;
        for (uint j = 0; j < stokens.size(); j++)
          vpositions_cartesian.push_back(aurostd::string2utype<double>(stokens.at(j)));
        } break;   // FIX
      case (aurostd::ctcrc64("positions_fractional")): {
        positions_fractional = content;
        for (uint j = 0; j < stokens.size(); j++)
          vpositions_fractional.push_back(aurostd::string2utype<double>(stokens.at(j)));
        } break;   // FIX
      case (aurostd::ctcrc64("pressure")): { pressure = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("stress_tensor")): {
        stress_tensor = content;
        vstress_tensor = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        if (stokens.size() == 9)
          for (uint j = 0; j < stokens.size(); j++)
            vstress_tensor.at(j) = aurostd::string2utype<double>(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("pressure_residual")): { pressure_residual = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("Pulay_stress")): { Pulay_stress = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("prototype")): { prototype = content; } break;  // apennsy
      case (aurostd::ctcrc64("PV_cell")): { PV_cell = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("PV_atom")): { PV_atom = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("scintillation_attenuation_length")): { scintillation_attenuation_length = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("sg")): {
        sg = content;
        for (uint j = 0; j < stokens.size(); j++) vsg.push_back(stokens.at(j));
        } break; //CO20180101
      case (aurostd::ctcrc64("sg2")): {
        sg2 = content;
        for (uint j = 0; j < stokens.size(); j++) vsg2.push_back(stokens.at(j));
        } break; //CO20180101
      case (aurostd::ctcrc64("spacegroup_orig")): { spacegroup_orig = aurostd::string2utype<int>(content); } break;  //CO20201111
      case (aurostd::ctcrc64("spacegroup_relax")): { spacegroup_relax = aurostd::string2utype<int>(content); } break;  //CO20201111
      case (aurostd::ctcrc64("species")): {
        species = content;
        for (uint j = 0; j < stokens.size(); j++) vspecies.push_back(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("species_pp")): {
        species_pp = content;
        for (uint j = 0; j < stokens.size(); j++) vspecies_pp.push_back(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("species_pp_version")): {
        species_pp_version = content;
        for (uint j = 0; j < stokens.size(); j++) vspecies_pp_version.push_back(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("species_pp_ZVAL")): {
        species_pp_ZVAL = content;
        for (uint j = 0; j < stokens.size(); j++)
          vspecies_pp_ZVAL.push_back(aurostd::string2utype<double>(stokens.at(j)));
        } break;
      case (aurostd::ctcrc64("species_pp_AUID")): {
        species_pp_AUID = content;
        for (uint j = 0; j < stokens.size(); j++) vspecies_pp_AUID.push_back(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("metagga")): case (aurostd::ctcrc64("METAGGA")): { METAGGA = content; } break;
      case (aurostd::ctcrc64("spin_cell")): { spin_cell = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("spin_atom")): { spin_atom = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("spinD")): {
        spinD = content;
        for (uint j = 0; j < stokens.size(); j++) vspinD.push_back(aurostd::string2utype<double>(stokens.at(j)));
        } break;
      case (aurostd::ctcrc64("spinD_magmom_orig")): {
        spinD_magmom_orig = content;
        for (uint j = 0; j < stokens.size(); j++)
          vspinD_magmom_orig.push_back(aurostd::string2utype<double>(stokens.at(j)));
        } break;
      case (aurostd::ctcrc64("spinF")): { spinF = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("sponsor")): {
        sponsor = content;
        aurostd::string2tokens(content, stokens, ",");
        for (uint j = 0; j < stokens.size(); j++) vsponsor.push_back(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("stoichiometry")): {
        stoichiometry = content;
        for (uint j = 0; j < stokens.size(); j++) vstoichiometry.push_back(aurostd::string2utype<double>(stokens.at(j)));
        } break;
      case (aurostd::ctcrc64("Egap_type")): { Egap_type = content; } break;
      case (aurostd::ctcrc64("valence_cell_std")): { valence_cell_std = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("valence_cell_iupac")): { valence_cell_iupac = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("volume_cell")): { volume_cell = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("volume_atom")): { volume_atom = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("volume_cell_orig")): { volume_cell_orig = aurostd::string2utype<double>(content); } break; //DX20190124 - add original crystal info
      case (aurostd::ctcrc64("volume_atom_orig")): { volume_atom_orig = aurostd::string2utype<double>(content); } break; //DX20190124 - add original crystal info
      // legacy
      case (aurostd::ctcrc64("server")): { vserver.push_back(content); } break;
      case (aurostd::ctcrc64("stoich")): {
        // aurostd::string2tokens(ventry.at(i), tokens, "="); //HE20220404 stoken already defined TODO debug
        // stoich = tokens.at(1);
        // aurostd::string2tokens(stoich, stokens);
        for (uint j = 0; j < stokens.size(); j++) vstoich.push_back(aurostd::string2utype<double>(stokens.at(j)));
        } break;
      //SYMMETRY
      //DX20190124 - added original symmetry info - START
      case (aurostd::ctcrc64("crystal_family_orig")): { crystal_family_orig = content; } break;
      case (aurostd::ctcrc64("crystal_system_orig")): { crystal_system_orig = content; } break;
      case (aurostd::ctcrc64("crystal_class_orig")): { crystal_class_orig = content; } break;
      case (aurostd::ctcrc64("point_group_Hermann_Mauguin_orig")): { point_group_Hermann_Mauguin_orig = content; } break;
      case (aurostd::ctcrc64("point_group_Schoenflies_orig")): { point_group_Schoenflies_orig = content; } break;
      case (aurostd::ctcrc64("point_group_orbifold_orig")): { point_group_orbifold_orig = content; } break;
      case (aurostd::ctcrc64("point_group_type_orig")): { point_group_type_orig = content; } break;
      case (aurostd::ctcrc64("point_group_order_orig")): { point_group_order = aurostd::string2utype<uint>(content); } break;
      case (aurostd::ctcrc64("point_group_structure_orig")): { point_group_structure_orig = content; } break;
      case (aurostd::ctcrc64("Bravais_lattice_lattice_type_orig")): { Bravais_lattice_lattice_type_orig = content; } break;
      case (aurostd::ctcrc64("Bravais_lattice_lattice_variation_type_orig")): { Bravais_lattice_lattice_variation_type_orig = content; } break;
      case (aurostd::ctcrc64("Bravais_lattice_lattice_system_orig")): { Bravais_lattice_lattice_system_orig = content; } break;
      case (aurostd::ctcrc64("Bravais_superlattice_lattice_type_orig")): { Bravais_superlattice_lattice_type_orig = content; } break;
      case (aurostd::ctcrc64("Bravais_superlattice_lattice_variation_type_orig")): { Bravais_superlattice_lattice_variation_type_orig = content; } break;
      case (aurostd::ctcrc64("Bravais_superlattice_lattice_system_orig")): { Bravais_superlattice_lattice_system_orig = content; } break;
      case (aurostd::ctcrc64("Pearson_symbol_superlattice_orig")): { Pearson_symbol_superlattice_orig = content; } break;
      case (aurostd::ctcrc64("reciprocal_geometry_orig")): {
        reciprocal_geometry_orig = content;
        vreciprocal_geometry_orig = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        if (stokens.size() == 6)
          for (uint j = 0; j < stokens.size(); j++)
            vreciprocal_geometry_orig.at(j) = aurostd::string2utype<double>(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("reciprocal_volume_cell_orig")): { reciprocal_volume_cell_orig = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("reciprocal_lattice_type_orig")): { reciprocal_lattice_type_orig = content; } break;
      case (aurostd::ctcrc64("reciprocal_lattice_variation_type_orig")): { reciprocal_lattice_variation_type_orig = content; } break;
      case (aurostd::ctcrc64("Wyckoff_letters_orig")): { Wyckoff_letters_orig = content; } break;
      case (aurostd::ctcrc64("Wyckoff_multiplicities_orig")): { Wyckoff_multiplicities_orig = content; } break;
      case (aurostd::ctcrc64("Wyckoff_site_symmetries_orig")): { Wyckoff_site_symmetries_orig = content; } break;
      //DX20190124 - added original symmetry info - END
      //DX20180823 - added more symmetry info - START
      case (aurostd::ctcrc64("crystal_family")): { crystal_family = content; } break;
      case (aurostd::ctcrc64("crystal_system")): { crystal_system = content; } break;
      case (aurostd::ctcrc64("crystal_class")): { crystal_class = content; } break;
      case (aurostd::ctcrc64("point_group_Hermann_Mauguin")): { point_group_Hermann_Mauguin = content; } break;
      case (aurostd::ctcrc64("point_group_Schoenflies")): { point_group_Schoenflies = content; } break;
      case (aurostd::ctcrc64("point_group_orbifold")): { point_group_orbifold = content; } break;
      case (aurostd::ctcrc64("point_group_type")): { point_group_type = content; } break;
      case (aurostd::ctcrc64("point_group_order")): { point_group_order = aurostd::string2utype<uint>(content); } break;
      case (aurostd::ctcrc64("point_group_structure")): { point_group_structure = content; } break;
      case (aurostd::ctcrc64("Bravais_lattice_lattice_type")): { Bravais_lattice_lattice_type = content; } break;
      case (aurostd::ctcrc64("Bravais_lattice_lattice_variation_type")): { Bravais_lattice_lattice_variation_type = content; } break;
      case (aurostd::ctcrc64("Bravais_lattice_lattice_system")): { Bravais_lattice_lattice_system = content; } break;
      case (aurostd::ctcrc64("Bravais_superlattice_lattice_type")): { Bravais_superlattice_lattice_type = content; } break;
      case (aurostd::ctcrc64("Bravais_superlattice_lattice_variation_type")): { Bravais_superlattice_lattice_variation_type = content; } break;
      case (aurostd::ctcrc64("Bravais_superlattice_lattice_system")): { Bravais_superlattice_lattice_system = content; } break;
      case (aurostd::ctcrc64("Pearson_symbol_superlattice")): { Pearson_symbol_superlattice = content; } break;
      case (aurostd::ctcrc64("reciprocal_geometry_relax")):
      case (aurostd::ctcrc64("reciprocal_geometry")): {
        reciprocal_geometry_relax = content;
        vreciprocal_geometry_relax = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        if (stokens.size() == 6)
          for (uint j = 0; j < stokens.size(); j++)
            vreciprocal_geometry_relax.at(j) = aurostd::string2utype<double>(stokens.at(j));
        } break;
      case (aurostd::ctcrc64("reciprocal_volume_cell")): { reciprocal_volume_cell = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("reciprocal_lattice_type")): { reciprocal_lattice_type = content; } break;
      case (aurostd::ctcrc64("reciprocal_lattice_variation_type")): { reciprocal_lattice_variation_type = content; } break;
      case (aurostd::ctcrc64("Wyckoff_letters")): { Wyckoff_letters = content; } break;
      case (aurostd::ctcrc64("Wyckoff_multiplicities")): { Wyckoff_multiplicities = content; } break;
      case (aurostd::ctcrc64("Wyckoff_site_symmetries")): { Wyckoff_site_symmetries = content; } break;
      //DX20180823 - added more symmetry info - END
      //DX20190209 - added anrl info - START
      case (aurostd::ctcrc64("aflow_prototype_label_orig")): { aflow_prototype_label_orig = content; } break;
      case (aurostd::ctcrc64("aflow_prototype_params_list_orig")): { aflow_prototype_params_list_orig = content; } break;
      case (aurostd::ctcrc64("aflow_prototype_params_values_orig")): { aflow_prototype_params_values_orig = content; } break;
      case (aurostd::ctcrc64("aflow_prototype_label_relax")): { aflow_prototype_label_relax = content; } break;
      case (aurostd::ctcrc64("aflow_prototype_params_list_relax")): { aflow_prototype_params_list_relax = content; } break;
      case (aurostd::ctcrc64("aflow_prototype_params_values_relax")): { aflow_prototype_params_values_relax = content; } break;
      //DX20190209 - added anrl info - END
      case (aurostd::ctcrc64("pocc_parameters")): { pocc_parameters = content; } break;  //CO20200731
      // AGL/AEL
      case (aurostd::ctcrc64("agl_thermal_conductivity_300K")): { agl_thermal_conductivity_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("agl_debye")): { agl_debye = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("agl_acoustic_debye")): { agl_acoustic_debye = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("agl_gruneisen")): { agl_gruneisen = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("agl_heat_capacity_Cv_300K")): { agl_heat_capacity_Cv_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("agl_heat_capacity_Cp_300K")): { agl_heat_capacity_Cp_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("agl_thermal_expansion_300K")): { agl_thermal_expansion_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("agl_bulk_modulus_static_300K")): { agl_bulk_modulus_static_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("agl_bulk_modulus_isothermal_300K")): { agl_bulk_modulus_isothermal_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("agl_poisson_ratio_source")): { agl_poisson_ratio_source = content; } break; //CT20181212
      case (aurostd::ctcrc64("agl_vibrational_free_energy_300K_cell")): { agl_vibrational_free_energy_300K_cell = aurostd::string2utype<double>(content); } break; //CT20181212
      case (aurostd::ctcrc64("agl_vibrational_free_energy_300K_atom")): { agl_vibrational_free_energy_300K_atom = aurostd::string2utype<double>(content); } break; //CT20181212
      case (aurostd::ctcrc64("agl_vibrational_entropy_300K_cell")): { agl_vibrational_entropy_300K_cell = aurostd::string2utype<double>(content); } break; //CT20181212
      case (aurostd::ctcrc64("agl_vibrational_entropy_300K_atom")): {  agl_vibrational_entropy_300K_atom = aurostd::string2utype<double>(content); } break; //CT20181212
      case (aurostd::ctcrc64("ael_poisson_ratio")): { ael_poisson_ratio = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("ael_bulk_modulus_voigt")): { ael_bulk_modulus_voigt = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("ael_bulk_modulus_reuss")): { ael_bulk_modulus_reuss = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("ael_shear_modulus_voigt")): { ael_shear_modulus_voigt = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("ael_shear_modulus_reuss")): { ael_shear_modulus_reuss = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("ael_bulk_modulus_vrh")): { ael_bulk_modulus_vrh = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("ael_shear_modulus_vrh")): { ael_shear_modulus_vrh = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("ael_elastic_anisotropy")): { ael_elastic_anisotropy = aurostd::string2utype<double>(content); } break; //CO20181129
      case (aurostd::ctcrc64("ael_youngs_modulus_vrh")): { ael_youngs_modulus_vrh = aurostd::string2utype<double>(content); } break; //CT20181212
      case (aurostd::ctcrc64("ael_speed_sound_transverse")): { ael_speed_sound_transverse = aurostd::string2utype<double>(content); } break; //CT20181212
      case (aurostd::ctcrc64("ael_speed_sound_longitudinal")): { ael_speed_sound_longitudinal = aurostd::string2utype<double>(content); } break; //CT20181212
      case (aurostd::ctcrc64("ael_speed_sound_average")): { ael_speed_sound_average = aurostd::string2utype<double>(content); } break; //CT20181212
      case (aurostd::ctcrc64("ael_pughs_modulus_ratio")): { ael_pughs_modulus_ratio = aurostd::string2utype<double>(content); } break; //CT20181212
      case (aurostd::ctcrc64("ael_debye_temperature")): { ael_debye_temperature = aurostd::string2utype<double>(content); } break; //CT20181212
      case (aurostd::ctcrc64("ael_applied_pressure")): { ael_applied_pressure = aurostd::string2utype<double>(content); } break; //CT20181212
      case (aurostd::ctcrc64("ael_average_external_pressure")): { ael_average_external_pressure = aurostd::string2utype<double>(content); } break; //CT20181212
      //ME20191105 BEGIN
      case (aurostd::ctcrc64("ael_stiffness_tensor")): {
        xmatrix<double> tensor(6, 6);
        vector<string> rows;
        vector<double> r;
        aurostd::string2tokens(content, rows, ";");
        if (rows.size() != 6) {
          stringstream message;
          message << "Could not read ael_stiffness_tensor: wrong number of rows"
                  << " (found " << rows.size() << ", need 6).";
          throw aurostd::xerror(__AFLOW_FILE__, function, message, _FILE_CORRUPT_);
        } else {
          for (int i = 0; i < 6; i++) {
            aurostd::string2tokens(rows[i], r, ",");
            if (r.size() != 6) {
              stringstream message;
              message << "Could not read ael_stiffness_tensor: wrong number of columns"
                      << " in row " << (i + 1)
                      << " (found " << rows.size() << ", need 6).";
              throw aurostd::xerror(__AFLOW_FILE__, function, message, _FILE_CORRUPT_);
            } else {
              for (int j = 0; j < 6; j++) {
                tensor[i + 1][j + 1] = r[j];
              }
            }
          }
        }
        ael_stiffness_tensor = tensor;
        } break;
      case (aurostd::ctcrc64("ael_compliance_tensor")): {
        xmatrix<double> tensor(6, 6);
        vector<string> rows;
        vector<double> r;
        aurostd::string2tokens(content, rows, ";");
        if (rows.size() != 6) {
          stringstream message;
          message << "Could not read ael_compliance_tensor: wrong number of rows"
                  << " (found " << rows.size() << ", need 6).";
          throw aurostd::xerror(__AFLOW_FILE__, function, message, _FILE_CORRUPT_);
        } else {
          for (int i = 0; i < 6; i++) {
            aurostd::string2tokens(rows[i], r, ",");
            if (r.size() != 6) {
              stringstream message;
              message << "Could not read ael_compliance_tensor: wrong number of columns"
                      << " in row " << (i + 1)
                      << " (found " << rows.size() << ", need 6).";
              throw aurostd::xerror(__AFLOW_FILE__, function, message, _FILE_CORRUPT_);
            } else {
              for (int j = 0; j < 6; j++) {
                tensor[i + 1][j + 1] = r[j];
              }
            }
          }
        }
        ael_compliance_tensor = tensor;
        } break;
      //ME20191105 END
      //ME20210927 BEGIN
      // APL
      case (aurostd::ctcrc64("energy_free_vibrational_cell_apl_300K")): { energy_free_vibrational_cell_apl_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("energy_free_vibrational_atom_apl_300K")): { energy_free_vibrational_atom_apl_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("entropy_vibrational_cell_apl_300K")): { entropy_vibrational_cell_apl_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("entropy_vibrational_atom_apl_300K")): { entropy_vibrational_atom_apl_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("energy_internal_vibrational_cell_apl_300K")): { energy_internal_vibrational_cell_apl_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("energy_internal_vibrational_atom_apl_300K")): { energy_internal_vibrational_atom_apl_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("energy_zero_point_cell_apl")): { energy_zero_point_cell_apl = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("energy_zero_point_atom_apl")): { energy_zero_point_atom_apl = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("heat_capacity_Cv_cell_apl_300K")): { heat_capacity_Cv_cell_apl_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("heat_capacity_Cv_atom_apl_300K")): { heat_capacity_Cv_atom_apl_300K = aurostd::string2utype<double>(content); } break;
      //ME20210927 END
      //AS20200901 BEGIN
      // QHA
      case (aurostd::ctcrc64("gruneisen_qha")): { gruneisen_qha = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("gruneisen_qha_300K")): { gruneisen_qha_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("thermal_expansion_qha_300K")): {
        thermal_expansion_qha_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("modulus_bulk_qha_300K")): { modulus_bulk_qha_300K = aurostd::string2utype<double>(content); } break;
      //AS20200901 END
      //AS20201008 BEGIN
      case (aurostd::ctcrc64("modulus_bulk_derivative_pressure_qha_300K")): { modulus_bulk_derivative_pressure_qha_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("heat_capacity_Cv_atom_qha_300K")): { heat_capacity_Cv_atom_qha_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("heat_capacity_Cv_cell_qha_300K")): { heat_capacity_Cv_cell_qha_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("heat_capacity_Cp_atom_qha_300K")): { heat_capacity_Cp_atom_qha_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("heat_capacity_Cp_cell_qha_300K")): { heat_capacity_Cp_cell_qha_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("volume_atom_qha_300K")): { volume_atom_qha_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("energy_free_atom_qha_300K")): { energy_free_atom_qha_300K = aurostd::string2utype<double>(content); } break;
      case (aurostd::ctcrc64("energy_free_cell_qha_300K")): { energy_free_cell_qha_300K = aurostd::string2utype<double>(content); } break;
      //AS20201008 END
      // BADER
      case (aurostd::ctcrc64("bader_net_charges")): { bader_net_charges = content; aurostd::string2tokens<double>(content, vbader_net_charges, ","); } break;
      case (aurostd::ctcrc64("bader_atomic_volumes")): { bader_atomic_volumes = content; aurostd::string2tokens<double>(content, vbader_atomic_volumes, ","); } break;
    }
  }


  // file2aflowlib
  uint _aflowlib_entry::file2aflowlib(const string& file,ostream& oss) {
    if(!aurostd::FileExist(file)) {cerr << "ERROR - _aflowlib_entry::file2aflowlib: " << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << " not found =" << file << endl;return 0;}
    string entry;
    aurostd::efile2string(file,entry);
    return Load(entry,oss);
  }

  // Load
  uint _aflowlib_entry::Load(const stringstream& stream,ostream& oss) {
    return Load(stream.str(),oss);
  }

  // LoadWeb
  uint _aflowlib_entry::url2aflowlib(const string& _url,ostream& oss,bool verbose) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string url=_url;
    if(url.empty()) {cerr << "ERROR - _aflowlib_entry::url2aflowlib: url.empty()" << endl;return 0;}
    string entry;
    if(aurostd::substring2bool(url,"index") || aurostd::substring2bool(url,"format")) {
      aurostd::StringSubst(url,"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,"");
      if(!aurostd::url2string(url,entry,verbose)){return 0;}   //CO, this is a dud
    } else {
      aurostd::StringSubst(url,"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,"");
      if(!aurostd::url2string(url+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,entry,verbose)){return 0;}  //CO, this is a dud
    }
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " entry=" << entry << endl;} //CO20180528
    return Load(entry,oss);
  }

  // Load variant for EntryLoader //HE20220404
  uint _aflowlib_entry::Load(const vector<uint64_t> & key_hash, const vector<string> & content) {
    clear();
    for (size_t member_index=0; member_index<content.size(); member_index++)  {
      std::string single_content = content[member_index];
      aurostd::StringSubst(single_content, "\"", "");
      aurostd::StringSubst(single_content, "],[", ";");
      aurostd::StringSubst(single_content, "[", "");
      aurostd::StringSubst(single_content, "]", "");
      SetByHash(key_hash[member_index], single_content);
    }
    LoadCleanup();
    return key_hash.size();
  }

  // Load overload
  uint _aflowlib_entry::Load(const string& _entry,ostream& oss) {
    clear(); // start from clean
    entry=_entry; // start from loading it up !
    if(entry.empty()) {cerr << "ERROR - _aflowlib_entry::Load: entry.empty()" << endl;return 0;}
    vector<string> tokens;
    string keyword,content,line;
    aurostd::string2tokens(entry,ventry,"|");
    for(uint i=0;i<ventry.size();i++) {
      line=aurostd::RemoveWhiteSpaces(ventry.at(i));
      aurostd::string2tokens(line,tokens,"=");
      if(tokens.size()>0) {
        keyword=tokens.at(0);
        if(tokens.size()>1) {content=tokens.at(1);} else {continue;} //{content="";}  //CO20180319, content="" screws up string2double(), better to leave as AUROSTD_NAN
        Set(keyword, content);
      }
    }
    LoadCleanup();
    if(0) {
      bool html=FALSE; //TRUE;  //CO20201220
      oss << "Keywords" << endl;
      oss << "auid=" << auid << (html?"<br>":"") << endl;
      oss << "aurl=" << aurl << (html?"<br>":"") << endl;
      //oss << "title=" << title << (html?"<br>":"") << endl;  //ME20190129  // OBSOLETE ME20200829 - not used anymore
      oss << "keywords=" << keywords << (html?"<br>":"") << "  vkeywords= ";for(uint j=0;j<vkeywords.size();j++) oss << vkeywords.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "Optional controls keywords (alphabetic order)" << endl;
      oss << "aflowlib_date=" << aflowlib_date << (html?"<br>":"") << endl; 
      oss << "aflowlib_version=" << aflowlib_version << (html?"<br>":"") << endl; 
      oss << "aflowlib_entries=" << aflowlib_entries << (html?"<br>":"") << endl; 
      oss << "aflowlib_entries_number=" << aflowlib_entries_number << (html?"<br>":"") << endl; 
      oss << "aflow_version=" << aflow_version << (html?"<br>":"") << endl; 
      oss << "catalog=" << catalog << (html?"<br>":"") << endl; 
      oss << "data_api=" << data_api << (html?"<br>":"") << endl; 
      oss << "data_source=" << data_source << (html?"<br>":"") << endl; 
      oss << "data_language=" << data_language << (html?"<br>":"") << endl; 
      oss << "error_status=" << error_status << (html?"<br>":"") << endl; 
      oss << "author=" << author << (html?"<br>":"") << "  vauthor= ";for(uint j=0;j<vauthor.size();j++) oss << vauthor.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "calculation_cores=" << calculation_cores << (html?"<br>":"") << endl; 
      oss << "calculation_memory=" << calculation_memory << (html?"<br>":"") << endl; 
      oss << "calculation_time=" << calculation_time << (html?"<br>":"") << endl; 
      oss << "corresponding=" << corresponding << (html?"<br>":"") << "  vcorresponding= ";for(uint j=0;j<vcorresponding.size();j++) oss << vcorresponding.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "loop=" << loop << (html?"<br>":"") << "  vloop= ";for(uint j=0;j<vloop.size();j++) oss << vloop.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "node_CPU_Cores=" << node_CPU_Cores << (html?"<br>":"") << endl; 
      oss << "node_CPU_MHz=" << node_CPU_MHz << (html?"<br>":"") << endl; 
      oss << "node_CPU_Model=" << node_CPU_Model << (html?"<br>":"") << endl; 
      oss << "node_RAM_GB=" << node_RAM_GB << (html?"<br>":"") << endl; 
      oss << "Optional materials keywords (alphabetic order)" << endl;
      oss << "Bravais_lattice_orig=" << Bravais_lattice_orig << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_relax=" << Bravais_lattice_relax << (html?"<br>":"") << endl;
      oss << "code=" << code << (html?"<br>":"") << endl;
      oss << "composition=" << composition << "  vcomposition= ";for(uint j=0;j<vcomposition.size();j++) oss << vcomposition.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "compound=" << compound << (html?"<br>":"") << endl;
      oss << "density=" << density << (html?"<br>":"") << endl;
      oss << "density_orig=" << density_orig << (html?"<br>":"") << endl; //DX20190124 - add original crystal info
      oss << "dft_type=" << dft_type << (html?"<br>":"") << "  vdft_type= ";for(uint j=0;j<vdft_type.size();j++) oss << vdft_type.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "eentropy_cell=" << eentropy_cell << (html?"<br>":"") << endl; 
      oss << "eentropy_atom=" << eentropy_atom << (html?"<br>":"") << endl; 
      oss << "Egap=" << Egap << (html?"<br>":"") << endl; 
      oss << "Egap_fit=" << Egap_fit << (html?"<br>":"") << endl; 
      oss << "Egap_type=" << Egap_type << (html?"<br>":"") << endl;
      oss << "energy_cell=" << energy_cell << (html?"<br>":"") << endl; 
      oss << "energy_atom=" << energy_atom << (html?"<br>":"") << endl; 
      oss << "energy_cutoff=" << energy_cutoff << (html?"<br>":"") << endl; 
      oss << "delta_electronic_energy_convergence=" << delta_electronic_energy_convergence << (html?"<br>":"") << endl; 
      oss << "delta_electronic_energy_threshold=" << delta_electronic_energy_threshold << (html?"<br>":"") << endl; 
      oss << "nkpoints=" << nkpoints << (html?"<br>":"") << endl; 
      oss << "nkpoints_irreducible=" << nkpoints_irreducible << (html?"<br>":"") << endl; 
      oss << "kppra=" << kppra << (html?"<br>":"") << endl; 
      oss << "kpoints_relax=" << aurostd::joinWDelimiter(kpoints_nnn_relax, ",") << (html?"<br>":"") << endl;  //ME20190129
      oss << "kpoints_static=" << aurostd::joinWDelimiter(kpoints_nnn_static, ",") << (html?"<br>":"") << endl;  //ME20190129
      oss << "kpoints_bands_path=" << aurostd::joinWDelimiter(kpoints_pairs, " | ") << endl;  //ME20190129
      oss << "kpoints_bands_nkpts=" << kpoints_bands_path_grid << (html?"<br>":"") << endl;  //ME20190129
      oss << "kpoints=" << kpoints << (html?"<br>":"") << endl;      
      oss << "enthalpy_cell=" << enthalpy_cell << (html?"<br>":"") << endl; 
      oss << "enthalpy_atom=" << enthalpy_atom << (html?"<br>":"") << endl; 
      oss << "enthalpy_formation_cell=" << enthalpy_formation_cell << (html?"<br>":"") << endl; 
      oss << "enthalpy_formation_cce_300K_cell=" << enthalpy_formation_cce_300K_cell << (html?"<br>":"") << endl;   //CO20200624
      oss << "enthalpy_formation_cce_0K_cell=" << enthalpy_formation_cce_0K_cell << (html?"<br>":"") << endl;   //CO20200624
      oss << "enthalpy_formation_atom=" << enthalpy_formation_atom << (html?"<br>":"") << endl; 
      oss << "enthalpy_formation_cce_300K_atom=" << enthalpy_formation_cce_300K_atom << (html?"<br>":"") << endl;   //CO20200624
      oss << "enthalpy_formation_cce_0K_atom=" << enthalpy_formation_cce_0K_atom << (html?"<br>":"") << endl;   //CO20200624
      oss << "entropy_forming_ability=" << entropy_forming_ability << (html?"<br>":"") << endl;   //CO20200624
      oss << "entropic_temperature=" << entropic_temperature << (html?"<br>":"") << endl; 
      // oss << "files=" << files << "  vfiles= ";for(uint j=0;j<vfiles.size();j++) oss << vfiles.at(j) << " "; oss << (html?"<br>":"") << endl;
      // oss << "files_LIB=" << files_LIB << "  vfiles_LIB= ";for(uint j=0;j<vfiles_LIB.size();j++) oss << vfiles_LIB.at(j) << " "; oss << (html?"<br>":"") << endl;
      // oss << "files_RAW=" << files_RAW << "  vfiles_RAW= ";for(uint j=0;j<vfiles_RAW.size();j++) oss << vfiles_RAW.at(j) << " "; oss << (html?"<br>":"") << endl;
      // oss << "files_WEB=" << files_WEB << "  vfiles_WEB= ";for(uint j=0;j<vfiles_WEB.size();j++) oss << vfiles_WEB.at(j) << " "; oss << (html?"<br>":"") << endl;
      // oss << "forces=" << forces << "  vforces= ";for(uint j=0;j<vforces.size();j++) oss << vforces.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "geometry=" << geometry << "  vgeometry= ";for(uint j=0;j<vgeometry.size();j++) oss << vgeometry.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "geometry_orig=" << geometry_orig << "  vgeometry_orig= ";for(uint j=0;j<vgeometry_orig.size();j++) oss << vgeometry_orig.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "lattice_system_orig=" << lattice_system_orig << (html?"<br>":"") << endl;
      oss << "lattice_variation_orig=" << lattice_variation_orig << (html?"<br>":"") << endl;
      oss << "lattice_system_relax=" << lattice_system_relax << (html?"<br>":"") << endl;
      oss << "lattice_variation_relax=" << lattice_variation_relax << (html?"<br>":"") << endl;
      oss << "ldau_TLUJ=" << ldau_TLUJ << (html?"<br>":"") << endl;      
      if (vLDAU.size()>0 && vLDAU[0].size()) {oss << "ldau_type=" << ((int) vLDAU[0][0]) << (html?"<br>":"") << endl;}  //ME20190129
      if (vLDAU.size()>1 && vLDAU[1].size()) {oss << "ldau_l="; oss << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vLDAU[1], 0), ",") << (html?"<br>":"") << endl;}  //ME20190129
      if (vLDAU.size()>2 && vLDAU[2].size()) {oss << "ldau_u="; oss << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vLDAU[2], 9), ",") << (html?"<br>":"") << endl;}  //ME20190129
      if (vLDAU.size()>3 && vLDAU[3].size()) {oss << "ldau_j="; oss << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vLDAU[3], 9), ",") << (html?"<br>":"") << endl;}  //ME20190129
      oss << "natoms=" << natoms << (html?"<br>":"") << endl;
      oss << "natoms_orig=" << natoms_orig << (html?"<br>":"") << endl; //DX20190124 - add original crystal info
      oss << "nbondxx=" << nbondxx << "  vnbondxx= ";for(uint j=0;j<vnbondxx.size();j++) oss << vnbondxx.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "nspecies=" << nspecies << (html?"<br>":"") << endl;
      oss << "Pearson_symbol_orig=" << Pearson_symbol_orig << (html?"<br>":"") << endl;
      oss << "Pearson_symbol_relax=" << Pearson_symbol_relax << (html?"<br>":"") << endl;
      // oss << "positions_cartesian=" << positions_cartesian << "  vpositions_cartesian= ";for(uint j=0;j<vpositions_cartesian.size();j++) oss << vpositions_cartesian.at(j) << " "; oss << (html?"<br>":"") << endl;
      // oss << "positions_fractional=" << positions_fractional << "  vpositions_fractional= ";for(uint j=0;j<vpositions_fractional.size();j++) oss << vpositions_fractional.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "pressure=" << pressure << (html?"<br>":"") << endl; 
      oss << "stress_tensor=" << stress_tensor << "  vstress_tensor= ";for(uint j=0;j<vstress_tensor.size();j++) oss << vstress_tensor.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "pressure_residual=" << pressure_residual << (html?"<br>":"") << endl; 
      oss << "Pulay_stress=" << Pulay_stress << (html?"<br>":"") << endl; 
      oss << "prototype=" << prototype << (html?"<br>":"") << endl;
      oss << "PV_cell=" << PV_cell << (html?"<br>":"") << endl; 
      oss << "PV_atom=" << PV_atom << (html?"<br>":"") << endl; 
      oss << "scintillation_attenuation_length=" << scintillation_attenuation_length << (html?"<br>":"") << endl;
      oss << "sg=" << sg << (html?"<br>":"") << endl;
      oss << "sg2=" << sg2 << (html?"<br>":"") << endl;
      oss << "spacegroup_orig=" << spacegroup_orig << (html?"<br>":"") << endl;
      oss << "spacegroup_relax=" << spacegroup_relax << (html?"<br>":"") << endl;
      oss << "species=" << species << "  vspecies= ";for(uint j=0;j<vspecies.size();j++) oss << vspecies.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "species_pp=" << species_pp << "  vspecies_pp= ";for(uint j=0;j<vspecies_pp.size();j++) oss << vspecies_pp.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "species_pp_version=" << species_pp_version << "  vspecies_pp_version= ";for(uint j=0;j<vspecies_pp_version.size();j++) oss << vspecies_pp_version.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "species_pp_ZVAL=" << species_pp_ZVAL << "  vspecies_pp_ZVAL= ";for(uint j=0;j<vspecies_pp_ZVAL.size();j++) oss << vspecies_pp_ZVAL.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "species_pp_AUID=" << species_pp_AUID << "  vspecies_pp_AUID= ";for(uint j=0;j<vspecies_pp_AUID.size();j++) oss << vspecies_pp_AUID.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "metagga=" << METAGGA << (html?"<br>":"") << endl;
      oss << "spin_cell=" << spin_cell << (html?"<br>":"") << endl; 
      oss << "spin_atom=" << spin_atom << (html?"<br>":"") << endl; 
      oss << "spinD=" << spinD << "  vspinD= "; for(uint j=0;j<vspinD.size();j++) oss << vspinD.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "spinD_magmom_orig=" << spinD_magmom_orig << "  vspinD_magmom_orig= "; for(uint j=0;j<vspinD_magmom_orig.size();j++) oss << vspinD_magmom_orig.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "spinF=" << spinF << (html?"<br>":"") << endl;
      oss << "sponsor=" << sponsor << (html?"<br>":"") << "  vsponsor= ";for(uint j=0;j<vsponsor.size();j++) oss << vsponsor.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "stoichiometry=" << stoichiometry << "  vstoichiometry= ";for(uint j=0;j<vstoichiometry.size();j++) oss << vstoichiometry.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "valence_cell_std=" << valence_cell_std << (html?"<br>":"") << endl; 
      oss << "valence_cell_iupac=" << valence_cell_iupac << (html?"<br>":"") << endl;      
      oss << "volume_cell=" << volume_cell << (html?"<br>":"") << endl; 
      oss << "volume_atom=" << volume_atom << (html?"<br>":"") << endl; 
      oss << "volume_cell_orig=" << volume_cell_orig << (html?"<br>":"") << endl; //DX20190124 - add original crystal info
      oss << "volume_atom_orig=" << volume_atom_orig << (html?"<br>":"") << endl; //DX20190124 - add original crystal info
      //DX20190124 - added original symmetry info - START
      // SYMMETRY
      oss << "crystal_family_orig=" << crystal_family_orig << (html?"<br>":"") << endl;
      oss << "crystal_system_orig=" << crystal_system_orig << (html?"<br>":"") << endl;
      oss << "crystal_class_orig=" << crystal_class_orig << (html?"<br>":"") << endl;
      oss << "point_group_Hermann_Mauguin_orig=" << point_group_Hermann_Mauguin_orig << (html?"<br>":"") << endl;
      oss << "point_group_Schoenflies_orig=" << point_group_Schoenflies_orig << (html?"<br>":"") << endl;
      oss << "point_group_orbifold_orig=" << point_group_orbifold_orig << (html?"<br>":"") << endl;
      oss << "point_group_type_orig=" << point_group_type_orig << (html?"<br>":"") << endl;
      oss << "point_group_order_orig=" << point_group_order_orig << (html?"<br>":"") << endl;
      oss << "point_group_structure_orig=" << point_group_structure_orig << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_lattice_type_orig=" << Bravais_lattice_lattice_type_orig << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_lattice_variation_type_orig=" << Bravais_lattice_lattice_variation_type_orig << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_lattice_system_orig=" << Bravais_lattice_lattice_system_orig << (html?"<br>":"") << endl;
      oss << "Bravais_superlattice_lattice_type_orig=" << Bravais_superlattice_lattice_type_orig << (html?"<br>":"") << endl;
      oss << "Bravais_superlattice_lattice_variation_type_orig=" << Bravais_superlattice_lattice_variation_type_orig << (html?"<br>":"") << endl;
      oss << "Bravais_superlattice_lattice_system_orig=" << Bravais_superlattice_lattice_system_orig << (html?"<br>":"") << endl;
      oss << "Pearson_symbol_superlattice_orig=" << Pearson_symbol_superlattice_orig << (html?"<br>":"") << endl;
      oss << "reciprocal_geometry_orig=" << reciprocal_geometry_orig << "  vreciprocal_geometry_orig= ";for(uint j=0;j<vreciprocal_geometry_orig.size();j++) oss << vreciprocal_geometry_orig.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "reciprocal_volume_cell_orig=" << reciprocal_volume_cell_orig << (html?"<br>":"") << endl;
      oss << "reciprocal_lattice_type_orig=" << reciprocal_lattice_type_orig << (html?"<br>":"") << endl;
      oss << "reciprocal_lattice_variation_type_orig=" << reciprocal_lattice_variation_type_orig << (html?"<br>":"") << endl;
      oss << "Wyckoff_letters_orig=" << Wyckoff_letters_orig << (html?"<br>":"") << endl;
      oss << "Wyckoff_multiplicities_orig=" << Wyckoff_multiplicities_orig << (html?"<br>":"") << endl;
      oss << "Wyckoff_site_symmetries_orig=" << Wyckoff_site_symmetries_orig << (html?"<br>":"") << endl;
      //DX20190124 - added original symmetry info - END
      //DX20180823 - added more symmetry info - START
      // SYMMETRY
      oss << "crystal_family=" << crystal_family << (html?"<br>":"") << endl;
      oss << "crystal_system=" << crystal_system << (html?"<br>":"") << endl;
      oss << "crystal_class=" << crystal_class << (html?"<br>":"") << endl;
      oss << "point_group_Hermann_Mauguin=" << point_group_Hermann_Mauguin << (html?"<br>":"") << endl;
      oss << "point_group_Schoenflies=" << point_group_Schoenflies << (html?"<br>":"") << endl;
      oss << "point_group_orbifold=" << point_group_orbifold << (html?"<br>":"") << endl;
      oss << "point_group_type=" << point_group_type << (html?"<br>":"") << endl;
      oss << "point_group_order=" << point_group_order << (html?"<br>":"") << endl;
      oss << "point_group_structure=" << point_group_structure << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_lattice_type=" << Bravais_lattice_lattice_type << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_lattice_variation_type=" << Bravais_lattice_lattice_variation_type << (html?"<br>":"") << endl;
      oss << "Bravais_lattice_lattice_system=" << Bravais_lattice_lattice_system << (html?"<br>":"") << endl;
      oss << "Bravais_superlattice_lattice_type=" << Bravais_superlattice_lattice_type << (html?"<br>":"") << endl;
      oss << "Bravais_superlattice_lattice_variation_type=" << Bravais_superlattice_lattice_variation_type << (html?"<br>":"") << endl;
      oss << "Bravais_superlattice_lattice_system=" << Bravais_superlattice_lattice_system << (html?"<br>":"") << endl;
      oss << "Pearson_symbol_superlattice=" << Pearson_symbol_superlattice << (html?"<br>":"") << endl;
      oss << "reciprocal_geometry_relax=" << reciprocal_geometry_relax << "  vreciprocal_geometry_relax= ";for(uint j=0;j<vreciprocal_geometry_relax.size();j++) oss << vreciprocal_geometry_relax.at(j) << " "; oss << (html?"<br>":"") << endl; //CO20220719 _relax
      oss << "reciprocal_volume_cell=" << reciprocal_volume_cell << (html?"<br>":"") << endl;
      oss << "reciprocal_lattice_type=" << reciprocal_lattice_type << (html?"<br>":"") << endl;
      oss << "reciprocal_lattice_variation_type=" << reciprocal_lattice_variation_type << (html?"<br>":"") << endl;
      oss << "Wyckoff_letters=" << Wyckoff_letters << (html?"<br>":"") << endl;
      oss << "Wyckoff_multiplicities=" << Wyckoff_multiplicities << (html?"<br>":"") << endl;
      oss << "Wyckoff_site_symmetries=" << Wyckoff_site_symmetries << (html?"<br>":"") << endl;
      //DX20180823 - added more symmetry info - END
      //DX20190208 - added anrl info - START
      oss << "aflow_prototype_label_orig=" << aflow_prototype_label_orig << (html?"<br>":"") << endl;
      oss << "aflow_prototype_params_list_orig=" << aflow_prototype_params_list_orig << (html?"<br>":"") << endl;
      oss << "aflow_prototype_params_values_orig=" << aflow_prototype_params_values_orig << (html?"<br>":"") << endl;
      oss << "aflow_prototype_label_relax=" << aflow_prototype_label_relax << (html?"<br>":"") << endl;
      oss << "aflow_prototype_params_list_relax=" << aflow_prototype_params_list_relax << (html?"<br>":"") << endl;
      oss << "aflow_prototype_params_values_relax=" << aflow_prototype_params_values_relax << (html?"<br>":"") << endl;
      //DX20190208 - added anrl info - END
      oss << "pocc_parameters=" << pocc_parameters << (html?"<br>":"") << endl;  //CO20200731
      // AGL/AEL
      oss << "agl_thermal_conductivity_300K=" << agl_thermal_conductivity_300K << (html?"<br>":"") << endl; 
      oss << "agl_debye=" << agl_debye << (html?"<br>":"") << endl; 
      oss << "agl_acoustic_debye=" << agl_acoustic_debye << (html?"<br>":"") << endl; 
      oss << "agl_gruneisen=" << agl_gruneisen << (html?"<br>":"") << endl; 
      oss << "agl_heat_capacity_Cv_300K=" << agl_heat_capacity_Cv_300K << (html?"<br>":"") << endl; 
      oss << "agl_heat_capacity_Cp_300K=" << agl_heat_capacity_Cp_300K << (html?"<br>":"") << endl; 
      oss << "agl_thermal_expansion_300K=" << agl_thermal_expansion_300K << (html?"<br>":"") << endl; 
      oss << "agl_bulk_modulus_static_300K=" << agl_bulk_modulus_static_300K << (html?"<br>":"") << endl; 
      oss << "agl_bulk_modulus_isothermal_300K=" << agl_bulk_modulus_isothermal_300K << (html?"<br>":"") << endl; 
      oss << "agl_poisson_ratio_source=" << agl_poisson_ratio_source << (html?"<br>":"") << endl; //CT20181212
      oss << "agl_vibrational_free_energy_300K_cell=" << agl_vibrational_free_energy_300K_cell << (html?"<br>":"") << endl; //CT20181212
      oss << "agl_vibrational_free_energy_300K_atom=" << agl_vibrational_free_energy_300K_atom << (html?"<br>":"") << endl; //CT20181212 
      oss << "agl_vibrational_entropy_300K_cell=" << agl_vibrational_entropy_300K_cell << (html?"<br>":"") << endl; //CT20181212 
      oss << "agl_vibrational_entropy_300K_atom=" << agl_vibrational_entropy_300K_atom << (html?"<br>":"") << endl; //CT20181212 
      oss << "ael_poisson_ratio=" << ael_poisson_ratio << (html?"<br>":"") << endl; 
      oss << "ael_bulk_modulus_voigt=" << ael_bulk_modulus_voigt << (html?"<br>":"") << endl; 
      oss << "ael_bulk_modulus_reuss=" << ael_bulk_modulus_reuss << (html?"<br>":"") << endl; 
      oss << "ael_shear_modulus_voigt=" << ael_shear_modulus_voigt << (html?"<br>":"") << endl; 
      oss << "ael_shear_modulus_reuss=" << ael_shear_modulus_reuss << (html?"<br>":"") << endl; 
      oss << "ael_bulk_modulus_vrh=" << ael_bulk_modulus_vrh << (html?"<br>":"") << endl; 
      oss << "ael_shear_modulus_vrh=" << ael_shear_modulus_vrh << (html?"<br>":"") << endl; 
      oss << "ael_elastic_anisotropy=" << ael_elastic_anisotropy << (html?"<br>":"") << endl; //CO20181129
      oss << "ael_youngs_modulus_vrh=" << ael_youngs_modulus_vrh << (html?"<br>":"") << endl; //CT20181212 
      oss << "ael_speed_sound_transverse=" << ael_speed_sound_transverse << (html?"<br>":"") << endl; //CT20181212 
      oss << "ael_speed_sound_longitudinal=" << ael_speed_sound_longitudinal << (html?"<br>":"") << endl; //CT20181212 
      oss << "ael_speed_sound_average=" << ael_speed_sound_average << (html?"<br>":"") << endl; //CT20181212
      oss << "ael_pughs_modulus_ratio=" << ael_pughs_modulus_ratio << (html?"<br>":"") << endl; //CT20181212 
      oss << "ael_debye_temperature=" << ael_debye_temperature << (html?"<br>":"") << endl; //CT20181212 
      oss << "ael_applied_pressure=" << ael_applied_pressure << (html?"<br>":"") << endl; //CT20181212 
      oss << "ael_average_external_pressure=" << ael_average_external_pressure << (html?"<br>":"") << endl; //CT20181212 
      //ME20191105 BEGIN
      oss << "ael_stiffness_tensor="; for (int i = ael_stiffness_tensor.lrows; i <= ael_stiffness_tensor.urows; i++) {for (int j = ael_stiffness_tensor.lcols; j <= ael_stiffness_tensor.ucols; j++) {oss << ael_stiffness_tensor[i][j] << " ";} oss << (html?"<br>":"") << endl;} //ME20191105  //CO20201220
      oss << "ael_compliance_tensor="; for (int i = ael_compliance_tensor.lrows; i <= ael_compliance_tensor.urows; i++) {for (int j = ael_compliance_tensor.lcols; j <= ael_compliance_tensor.ucols; j++) {oss << ael_compliance_tensor[i][j] << " ";} oss << (html?"<br>":"") << endl;} //ME20191105  //CO20201220
      //ME20191105 END
      //ME20210927 BEGIN
      oss << "energy_free_vibrational_cell_apl_300K=" << energy_free_vibrational_cell_apl_300K << (html?"<br>":"") << endl;
      oss << "energy_free_vibrational_atom_apl_300K=" << energy_free_vibrational_atom_apl_300K << (html?"<br>":"") << endl;
      oss << "entropy_vibrational_cell_apl_300K=" << entropy_vibrational_cell_apl_300K << (html?"<br>":"") << endl;
      oss << "entropy_vibrational_atom_apl_300K=" << entropy_vibrational_atom_apl_300K << (html?"<br>":"") << endl;
      oss << "energy_internal_vibrational_cell_apl_300K=" << energy_internal_vibrational_cell_apl_300K << (html?"<br>":"") << endl;
      oss << "energy_internal_vibrational_atom_apl_300K=" << energy_internal_vibrational_atom_apl_300K << (html?"<br>":"") << endl;
      oss << "energy_zero_point_cell_apl=" << energy_zero_point_cell_apl << (html?"<br>":"") << endl;
      oss << "energy_zero_point_atom_apl=" << energy_zero_point_atom_apl << (html?"<br>":"") << endl;
      oss << "heat_capacity_Cv_cell_apl_300K=" << heat_capacity_Cv_cell_apl_300K << (html?"<br>":"") << endl;
      oss << "heat_capacity_Cv_atom_apl_300K=" << heat_capacity_Cv_atom_apl_300K << (html?"<br>":"") << endl;
      //ME20210927 END
      //AS20200901 BEGIN
      //QHA
      oss << "gruneisen_qha=" << gruneisen_qha << (html?"<br>":"") << endl;
      oss << "gruneisen_qha_300K=" << gruneisen_qha_300K << (html?"<br>":"") << endl;
      oss << "thermal_expansion_qha_300K=" << thermal_expansion_qha_300K << (html?"<br>":"") << endl;
      oss << "modulus_bulk_qha_300K=" << modulus_bulk_qha_300K << (html?"<br>":"") << endl;
      //AS20200901 END
      //AS20201008 BEGIN
      oss << "modulus_bulk_derivative_pressure_qha_300K=" << modulus_bulk_derivative_pressure_qha_300K << (html?"<br>":"") << endl;
      oss << "heat_capacity_Cv_atom_qha_300K=" << heat_capacity_Cv_atom_qha_300K << (html?"<br>":"") << endl;
      oss << "heat_capacity_Cv_cell_qha_300K=" << heat_capacity_Cv_cell_qha_300K << (html?"<br>":"") << endl;
      oss << "heat_capacity_Cp_atom_qha_300K=" << heat_capacity_Cp_atom_qha_300K << (html?"<br>":"") << endl;
      oss << "heat_capacity_Cp_cell_qha_300K=" << heat_capacity_Cp_cell_qha_300K << (html?"<br>":"") << endl;
      oss << "volume_atom_qha_300K=" << volume_atom_qha_300K << (html?"<br>":"") << endl;
      oss << "energy_free_atom_qha_300K=" << energy_free_atom_qha_300K << (html?"<br>":"") << endl;
      oss << "energy_free_cell_qha_300K=" << energy_free_cell_qha_300K << (html?"<br>":"") << endl;
      //AS20201008 END
      // BADER
      oss << "bader_net_charges=" << bader_net_charges << "  vbader_net_charges= ";for(uint j=0;j<vbader_net_charges.size();j++) oss << vbader_net_charges.at(j) << " "; oss << (html?"<br>":"") << endl; 
      oss << "bader_atomic_volumes=" << bader_atomic_volumes << "  vbader_atomic_volumes= ";for(uint j=0;j<vbader_atomic_volumes.size();j++) oss << vbader_atomic_volumes.at(j) << " "; oss << (html?"<br>":"") << endl; 
      // legacy
      oss << "server=" << server << (html?"<br>":"") << "  vserver= ";for(uint j=0;j<vserver.size();j++) oss << vserver.at(j) << " "; oss << (html?"<br>":"") << endl;
      oss << "icsd=" << icsd << (html?"<br>":"") << endl;
      oss << "stoich=" << stoich << "  vstoich= ";for(uint j=0;j<vstoich.size();j++) oss << vstoich.at(j) << " "; oss << (html?"<br>":"") << endl;
    }
    return ventry.size();
  }

  // seperated to allow two different load strategies without repeting the cleanup //HE20220404
  void _aflowlib_entry::LoadCleanup() {
    vector<string> tokens;
    //ME20190129 - FIX vLDAU //CO20210713 - there's only 1 type, but the number is repeated for the number of species
    if (vLDAU[0].size()) vLDAU[0].assign(vLDAU[1].size(), vLDAU[0][0]);
    // FIX LOOP
    loop="";
    vloop.push_back("aflowlib");
    //    for(uint j=0;j<vloop.size();j++) loop+=vloop.at(j)+(j<vloop.size()-1?", ":"");
    for(uint j=0;j<vloop.size();j++) loop+=vloop.at(j)+(j<vloop.size()-1?",":""); // no space
    // FIX SERVER
    server="";
    for(uint j=0;j<vserver.size();j++) {
      server+=vserver.at(j)+(j<vserver.size()-1?", ":"");
      vserverdir.push_back(*(new vector<string>(0))); // space
    }
    // FIX ICSD
    if(aurostd::substring2bool(prototype,"_ICSD_")) {
      aurostd::string2tokens(prototype,tokens,"_");
      icsd=tokens.at(tokens.size()-1);
    }
    // FIX APENNSY
    structure_name=prototype;
    structure_description=prototype;
    distance_gnd=AUROSTD_MAX_DOUBLE; // gotta calculate it
    distance_tie=AUROSTD_MAX_DOUBLE; // gotta calculate it
    pureA=FALSE;pureB=FALSE;
    fcc=FALSE; bcc=FALSE;hcp=FALSE;
    stoich_a=AUROSTD_MAX_DOUBLE;stoich_b=AUROSTD_MAX_DOUBLE;
    bond_aa=AUROSTD_MAX_DOUBLE;bond_ab=AUROSTD_MAX_DOUBLE;bond_bb=AUROSTD_MAX_DOUBLE;
    vNsgroup.clear();vsgroup.clear();vstr.clear();  // apennsy
    //fixing species
    if((species.empty()||vspecies.empty())&&!vspecies_pp.empty()){  //CO+HE20221110 - sometimes species is empty
      species.clear();vspecies.clear();
      for(uint i=0;i<vspecies_pp.size();i++){vspecies.push_back(vspecies_pp[i]);aurostd::VASP_PseudoPotential_CleanName_InPlace(vspecies.back());}
      species=aurostd::joinWDelimiter(vspecies,",");
    }
    // DONE
  }

  // aflowlib2string
  string _aflowlib_entry::aflowlib2string(string mode, bool PRINT_NULL) {
    stringstream sss("");
    //  string eendl="\n";

    // this is the normal aflowlib.out mode
    if(mode=="" || mode=="out" || mode=="OUT") {
      string eendl="";

      if(auid.size()) sss << "" << "aurl=" << aurl << eendl;
      if(auid.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "auid=" << auid << eendl;
      //if(!title.empty()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "title=" << title << eendl;  //ME20190125 // OBSOLETE ME20200829 - not used anymore
      if(data_api.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "data_api=" << data_api << eendl;
      if(data_source.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "data_source=" << data_source << eendl;
      if(data_language.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "data_language=" << data_language << eendl;
      if(error_status.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "error_status=" << error_status << eendl;
      // LOOP
      if(vloop.size()) {
        aurostd::sort(vloop);
        sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "loop=";
        for(uint i=0;i<vloop.size();i++) sss << vloop.at(i) << (i<vloop.size()-1?",":"");
        sss << eendl;
      }
      // MATERIALS
      if(code.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "code=" << code << eendl;
      if(compound.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "compound=" << compound << eendl;
      if(prototype.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "prototype=" << prototype << eendl;
      if(nspecies!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "nspecies=" << nspecies << eendl;
      if(natoms!=AUROSTD_NAN)sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "natoms=" << natoms << eendl;
      if(natoms_orig!=AUROSTD_NAN)sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "natoms_orig=" << natoms_orig << eendl; //DX20190124 - add original crystal info
      if(composition.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "composition=" << composition << eendl;
      if(density!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "density=" << density << eendl;
      if(density_orig!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "density_orig=" << density_orig << eendl; //DX20190124 - add original crystal info
      if(scintillation_attenuation_length!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "scintillation_attenuation_length=" << scintillation_attenuation_length << eendl;
      if(stoichiometry.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "stoichiometry=" << stoichiometry << eendl;
      if(species.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species=" << species << eendl;
      if(species_pp.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species_pp=" << species_pp << eendl;
      if(dft_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "dft_type=" << dft_type << eendl;
      // if(species_pp_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species_pp_type=" << species_pp_type << eendl;
      if(species_pp_version.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species_pp_version=" << species_pp_version << eendl;
      if(species_pp_ZVAL.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species_pp_ZVAL=" << species_pp_ZVAL << eendl;
      if(species_pp_AUID.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "species_pp_AUID=" << species_pp_AUID << eendl;
      if(METAGGA.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "metagga=" << METAGGA << eendl;
      //ME20190124 - add more detailed LDAU information
      if(vLDAU.size()>0 && vLDAU[0].size()>0){sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ldau_type=" << vLDAU[0][0] << eendl;} //ME20190124  //CO20210713  //CO+DX20210726 - note precision 0, from wiki: "The first field indicates the type (T) of the DFT+U corrections: type=1, the rotationally invariant version introduced by Liechtenstein et al.29); type=2, the simplified rotationally invariant version introduced by Dudarev et al.30)."
      if(vLDAU.size()>1 && vLDAU[1].size()>0){sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ldau_l=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vLDAU[1], 0), ",") << eendl;}  //ME20190124  //CO20210713  //CO+DX20210726 - note precision 0, from wiki: "The second field indicates the l-quantum number ({L}, one number for each species separated by ",") for which the on-site interaction is added (-1=neglected, 0=$s$, 1=$p$, 2=$d$, 3=$f$)."
      if(vLDAU.size()>2 && vLDAU[2].size()>0){sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ldau_u=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vLDAU[2], 9), ",") << eendl;}  //ME20190124  //CO20210713  //CO+DX20210726 - note precision 9, from wiki: "The third field lists the effective on-site Coulomb interaction parameters ({U}, one number for each species separated by ",")."
      if(vLDAU.size()>3 && vLDAU[3].size()>0){sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ldau_j=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vLDAU[3], 9), ",") << eendl;}  //ME20190124  //CO20210713  //CO+DX20210726 - note precision 9, from wiki: "The fourth field species the effective on-site exchange interaction parameters ({J}, one number for each species separated by ",")."
      if(ldau_TLUJ.size()) {sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ldau_TLUJ=" << ldau_TLUJ << eendl;} //ME20190124  //CO20210713
      if(valence_cell_iupac!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "valence_cell_iupac=" << valence_cell_iupac << eendl;
      if(valence_cell_std!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "valence_cell_std=" << valence_cell_std << eendl;
      if(volume_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "volume_cell=" << volume_cell << eendl;
      if(volume_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "volume_atom=" << volume_atom << eendl;
      if(volume_cell_orig!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "volume_cell_orig=" << volume_cell_orig << eendl; //DX20190124 - add original crystal info
      if(volume_atom_orig!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "volume_atom_orig=" << volume_atom_orig << eendl; //DX20190124 - add original crystal info
      if(pressure!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "pressure=" << pressure << eendl;
      if(stress_tensor.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "stress_tensor=" << stress_tensor << eendl;
      if(pressure_residual!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "pressure_residual=" << pressure_residual << eendl;
      if(Pulay_stress!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Pulay_stress=" << Pulay_stress << eendl;
      if(geometry.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "geometry=" << geometry << eendl;
      if(geometry_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "geometry_orig=" << geometry_orig << eendl; //DX20190124 - add original crystal info
      if(Egap!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Egap=" << Egap << eendl;
      if(Egap_fit!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Egap_fit=" << Egap_fit << eendl;
      if(Egap_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Egap_type=" << Egap_type << eendl;
      if(energy_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_cell=" << energy_cell << eendl;
      if(energy_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_atom=" << energy_atom << eendl;
      if(energy_cutoff!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_cutoff=" << energy_cutoff << eendl;
      if(delta_electronic_energy_convergence!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "delta_electronic_energy_convergence=" << delta_electronic_energy_convergence << eendl;
      if(delta_electronic_energy_threshold!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "delta_electronic_energy_threshold=" << delta_electronic_energy_threshold << eendl;
      // [NOT_PRINTED]     if(nkpoints!=0) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "nkpoints=" << nkpoints << eendl;
      // [NOT_PRINTED]     if(nkpoints_irreducible!=0) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "nkpoints_irreducible=" << nkpoints_irreducible << eendl;
      // [NOT_PRINTED]     if(kppra!=0) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "kppra=" << kppra << eendl;
      //ME20190124 BEGIN - Add the individual pieces of "kpoints" to the out file
      if ((kpoints_nnn_relax.rows == 3) && (sum(kpoints_nnn_relax) > 0)) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "kpoints_relax=" << aurostd::joinWDelimiter(kpoints_nnn_relax, ",") << eendl;
      if ((kpoints_nnn_static.rows == 3) && (sum(kpoints_nnn_static) > 0)) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "kpoints_static=" << aurostd::joinWDelimiter(kpoints_nnn_static, ",") << eendl;
      if (kpoints_pairs.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "kpoints_bands_path=" << aurostd::joinWDelimiter(kpoints_pairs, ",") << eendl;
      if (kpoints_bands_path_grid > 0) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "kpoints_bands_nkpts=" << ((int) kpoints_bands_path_grid) << eendl;
      //ME20190124 END
      if(kpoints.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "kpoints=" << kpoints << eendl;
      if(enthalpy_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_cell=" << enthalpy_cell << eendl;
      if(enthalpy_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_atom=" << enthalpy_atom << eendl;
      if(eentropy_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "eentropy_cell=" << eentropy_cell << eendl;
      if(eentropy_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "eentropy_atom=" << eentropy_atom << eendl;
      if(enthalpy_formation_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_formation_cell=" << enthalpy_formation_cell << eendl;
      if(enthalpy_formation_cce_300K_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_formation_cce_300K_cell=" << enthalpy_formation_cce_300K_cell << eendl; //CO20200624
      if(enthalpy_formation_cce_0K_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_formation_cce_0K_cell=" << enthalpy_formation_cce_0K_cell << eendl; //CO20200624
      if(enthalpy_formation_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_formation_atom=" << enthalpy_formation_atom << eendl;
      if(enthalpy_formation_cce_300K_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_formation_cce_300K_atom=" << enthalpy_formation_cce_300K_atom << eendl; //CO20200624
      if(enthalpy_formation_cce_0K_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "enthalpy_formation_cce_0K_atom=" << enthalpy_formation_cce_0K_atom << eendl; //CO20200624
      if(entropy_forming_ability!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "entropy_forming_ability=" << entropy_forming_ability << eendl; //CO20200624
      if(entropic_temperature!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "entropic_temperature=" << entropic_temperature << eendl;
      if(PV_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "PV_cell=" << PV_cell << eendl;
      if(PV_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "PV_atom=" << PV_atom << eendl;
      if(spin_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spin_cell=" << spin_cell << eendl;
      if(spin_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spin_atom=" << spin_atom << eendl;
      if(spinD.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spinD=" << spinD << eendl;
      if(spinD_magmom_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spinD_magmom_orig=" << spinD_magmom_orig << eendl;
      if(spinF!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spinF=" << spinF << eendl;
      if(stoich.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "stoich=" << stoich << eendl;
      if(calculation_time!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "calculation_time=" << calculation_time << eendl;
      if(calculation_memory!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "calculation_memory=" << calculation_memory << eendl;
      if(calculation_cores!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "calculation_cores=" << calculation_cores << eendl;
      if(nbondxx.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "nbondxx=" << nbondxx << eendl;
      if(sg.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "sg=" << sg << eendl;
      if(sg2.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "sg2=" << sg2 << eendl;
      if(spacegroup_orig!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spacegroup_orig=" << spacegroup_orig << eendl; //CO20201111
      if(spacegroup_relax!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "spacegroup_relax=" << spacegroup_relax << eendl;  //CO20201111
      if(forces.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "forces=" << forces << eendl;
      if(positions_cartesian.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "positions_cartesian=" << positions_cartesian << eendl;
      if(positions_fractional.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "positions_fractional=" << positions_fractional << eendl;
      if(Bravais_lattice_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_orig=" << Bravais_lattice_orig << eendl;
      if(lattice_variation_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "lattice_variation_orig=" << lattice_variation_orig << eendl;
      if(lattice_system_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "lattice_system_orig=" << lattice_system_orig << eendl;
      if(Pearson_symbol_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Pearson_symbol_orig=" << Pearson_symbol_orig << eendl;
      if(Bravais_lattice_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_relax=" << Bravais_lattice_relax << eendl;
      if(lattice_variation_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "lattice_variation_relax=" << lattice_variation_relax << eendl;
      if(lattice_system_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "lattice_system_relax=" << lattice_system_relax << eendl;
      if(Pearson_symbol_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Pearson_symbol_relax=" << Pearson_symbol_relax << eendl;
      //DX20190124 - added original symmetry info - START
      // SYMMETRY
      if(crystal_family_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "crystal_family_orig=" << crystal_family_orig << eendl;
      if(crystal_system_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "crystal_system_orig=" << crystal_system_orig << eendl;
      if(crystal_class_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "crystal_class_orig=" << crystal_class_orig << eendl;
      if(point_group_Hermann_Mauguin_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_Hermann_Mauguin_orig=" << point_group_Hermann_Mauguin_orig << eendl;
      if(point_group_Schoenflies_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_Schoenflies_orig=" << point_group_Schoenflies_orig << eendl;
      if(point_group_orbifold_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_orbifold_orig=" << point_group_orbifold_orig << eendl;
      if(point_group_type_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_type_orig=" << point_group_type_orig << eendl;
      if(point_group_order_orig!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_order_orig=" << point_group_order_orig << eendl;
      if(point_group_structure_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_structure_orig=" << point_group_structure_orig << eendl;
      if(Bravais_lattice_lattice_type_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_lattice_type_orig=" << Bravais_lattice_lattice_type_orig << eendl;
      if(Bravais_lattice_lattice_variation_type_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_lattice_variation_type_orig=" << Bravais_lattice_lattice_variation_type_orig << eendl;
      if(Bravais_lattice_lattice_system_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_lattice_system_orig=" << Bravais_lattice_lattice_system_orig << eendl;
      if(Bravais_superlattice_lattice_type_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_superlattice_lattice_type_orig=" << Bravais_superlattice_lattice_type_orig << eendl;
      if(Bravais_superlattice_lattice_variation_type_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_superlattice_lattice_variation_type_orig=" << Bravais_superlattice_lattice_variation_type_orig << eendl;
      if(Bravais_superlattice_lattice_system_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_superlattice_lattice_system_orig=" << Bravais_superlattice_lattice_system_orig << eendl;
      if(Pearson_symbol_superlattice_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Pearson_symbol_superlattice_orig=" << Pearson_symbol_superlattice_orig << eendl;
      if(reciprocal_geometry_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_geometry_orig=" << reciprocal_geometry_orig << eendl;
      if(reciprocal_volume_cell_orig!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_volume_cell_orig=" << reciprocal_volume_cell_orig << eendl;
      if(reciprocal_lattice_type_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_lattice_type_orig=" << reciprocal_lattice_type_orig << eendl;
      if(reciprocal_lattice_variation_type_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_lattice_variation_type_orig=" << reciprocal_lattice_variation_type_orig << eendl;
      if(Wyckoff_letters_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Wyckoff_letters_orig=" << Wyckoff_letters_orig << eendl;
      if(Wyckoff_multiplicities_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Wyckoff_multiplicities_orig=" << Wyckoff_multiplicities_orig << eendl;
      if(Wyckoff_site_symmetries_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Wyckoff_site_symmetries_orig=" << Wyckoff_site_symmetries_orig << eendl;
      //DX20190124 - added original symmetry info - END
      //DX20180823 - added more symmetry info - START
      // SYMMETRY
      if(crystal_family.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "crystal_family=" << crystal_family << eendl;
      if(crystal_system.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "crystal_system=" << crystal_system << eendl;
      if(crystal_class.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "crystal_class=" << crystal_class << eendl;
      if(point_group_Hermann_Mauguin.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_Hermann_Mauguin=" << point_group_Hermann_Mauguin << eendl;
      if(point_group_Schoenflies.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_Schoenflies=" << point_group_Schoenflies << eendl;
      if(point_group_orbifold.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_orbifold=" << point_group_orbifold << eendl;
      if(point_group_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_type=" << point_group_type << eendl;
      if(point_group_order!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_order=" << point_group_order << eendl;
      if(point_group_structure.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "point_group_structure=" << point_group_structure << eendl;
      if(Bravais_lattice_lattice_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_lattice_type=" << Bravais_lattice_lattice_type << eendl;
      if(Bravais_lattice_lattice_variation_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_lattice_variation_type=" << Bravais_lattice_lattice_variation_type << eendl;
      if(Bravais_lattice_lattice_system.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_lattice_lattice_system=" << Bravais_lattice_lattice_system << eendl;
      if(Bravais_superlattice_lattice_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_superlattice_lattice_type=" << Bravais_superlattice_lattice_type << eendl;
      if(Bravais_superlattice_lattice_variation_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_superlattice_lattice_variation_type=" << Bravais_superlattice_lattice_variation_type << eendl;
      if(Bravais_superlattice_lattice_system.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Bravais_superlattice_lattice_system=" << Bravais_superlattice_lattice_system << eendl;
      if(Pearson_symbol_superlattice.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Pearson_symbol_superlattice=" << Pearson_symbol_superlattice << eendl;
      if(reciprocal_geometry_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_geometry_relax=" << reciprocal_geometry_relax << eendl; //CO20220719 _relax
      if(reciprocal_volume_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_volume_cell=" << reciprocal_volume_cell << eendl;
      if(reciprocal_lattice_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_lattice_type=" << reciprocal_lattice_type << eendl;
      if(reciprocal_lattice_variation_type.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "reciprocal_lattice_variation_type=" << reciprocal_lattice_variation_type << eendl;
      if(Wyckoff_letters.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Wyckoff_letters=" << Wyckoff_letters << eendl;
      if(Wyckoff_multiplicities.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Wyckoff_multiplicities=" << Wyckoff_multiplicities << eendl;
      if(Wyckoff_site_symmetries.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "Wyckoff_site_symmetries=" << Wyckoff_site_symmetries << eendl;
      //DX20180823 - added more symmetry info - END
      //DX20190208 - added anrl info - START
      if(aflow_prototype_label_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "aflow_prototype_label_orig=" << aflow_prototype_label_orig << eendl;
      if(aflow_prototype_params_list_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "aflow_prototype_params_list_orig=" << aflow_prototype_params_list_orig << eendl;
      if(aflow_prototype_params_values_orig.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "aflow_prototype_params_values_orig=" << aflow_prototype_params_values_orig << eendl;
      if(aflow_prototype_label_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "aflow_prototype_label_relax=" << aflow_prototype_label_relax << eendl;
      if(aflow_prototype_params_list_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "aflow_prototype_params_list_relax=" << aflow_prototype_params_list_relax << eendl;
      if(aflow_prototype_params_values_relax.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "aflow_prototype_params_values_relax=" << aflow_prototype_params_values_relax << eendl;
      //DX20190208 - added anrl info - END
      if(pocc_parameters.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "pocc_parameters=" << pocc_parameters << eendl; //CO20200731
      // AGL/AEL
      if(agl_thermal_conductivity_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_thermal_conductivity_300K=" << agl_thermal_conductivity_300K << eendl;
      if(agl_debye!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_debye=" << agl_debye << eendl;
      if(agl_acoustic_debye!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_acoustic_debye=" << agl_acoustic_debye << eendl;
      if(agl_gruneisen!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_gruneisen=" << agl_gruneisen << eendl;
      if(agl_heat_capacity_Cv_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_heat_capacity_Cv_300K=" << agl_heat_capacity_Cv_300K << eendl;
      if(agl_heat_capacity_Cp_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_heat_capacity_Cp_300K=" << agl_heat_capacity_Cp_300K << eendl;
      if(agl_thermal_expansion_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_thermal_expansion_300K=" << agl_thermal_expansion_300K << eendl;
      if(agl_bulk_modulus_static_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_bulk_modulus_static_300K=" << agl_bulk_modulus_static_300K << eendl;
      if(agl_bulk_modulus_isothermal_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_bulk_modulus_isothermal_300K=" << agl_bulk_modulus_isothermal_300K << eendl;
      if(agl_poisson_ratio_source.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_poisson_ratio_source=" << agl_poisson_ratio_source << eendl; //CT20181212
      if(agl_vibrational_free_energy_300K_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_vibrational_free_energy_300K_cell=" << agl_vibrational_free_energy_300K_cell << eendl; //CT20181212
      if(agl_vibrational_free_energy_300K_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_vibrational_free_energy_300K_atom=" << agl_vibrational_free_energy_300K_atom << eendl; //CT20181212
      if(agl_vibrational_entropy_300K_cell!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_vibrational_entropy_300K_cell=" << agl_vibrational_entropy_300K_cell << eendl; //CT20181212
      if(agl_vibrational_entropy_300K_atom!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "agl_vibrational_entropy_300K_atom=" << agl_vibrational_entropy_300K_atom << eendl; //CT20181212
      if(ael_poisson_ratio!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_poisson_ratio=" << ael_poisson_ratio << eendl;
      if(ael_bulk_modulus_voigt!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_bulk_modulus_voigt=" << ael_bulk_modulus_voigt << eendl;
      if(ael_bulk_modulus_reuss!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_bulk_modulus_reuss=" << ael_bulk_modulus_reuss << eendl;
      if(ael_shear_modulus_voigt!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_shear_modulus_voigt=" << ael_shear_modulus_voigt << eendl;
      if(ael_shear_modulus_reuss!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_shear_modulus_reuss=" << ael_shear_modulus_reuss << eendl;
      if(ael_bulk_modulus_vrh!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_bulk_modulus_vrh=" << ael_bulk_modulus_vrh << eendl;
      if(ael_shear_modulus_vrh!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_shear_modulus_vrh=" << ael_shear_modulus_vrh << eendl;
      if(ael_elastic_anisotropy!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_elastic_anisotropy=" << ael_elastic_anisotropy << eendl; //CO20181129
      if(ael_youngs_modulus_vrh!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_youngs_modulus_vrh=" << ael_youngs_modulus_vrh << eendl; //CT20181212
      if(ael_speed_sound_transverse!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_speed_sound_transverse=" << ael_speed_sound_transverse << eendl; //CT20181212
      if(ael_speed_sound_longitudinal!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_speed_sound_longitudinal=" << ael_speed_sound_longitudinal << eendl; //CT20181212
      if(ael_speed_sound_average!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_speed_sound_average=" << ael_speed_sound_average << eendl; //CT20181212
      if(ael_pughs_modulus_ratio!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_pughs_modulus_ratio=" << ael_pughs_modulus_ratio << eendl; //CT20181212
      if(ael_debye_temperature!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_debye_temperature=" << ael_debye_temperature << eendl; //CT20181212
      if(ael_applied_pressure!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_applied_pressure=" << ael_applied_pressure << eendl; //CT20181212
      if(ael_average_external_pressure!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_average_external_pressure=" << ael_average_external_pressure << eendl; //CT20181212
      //ME20191105 BEGIN
      if ((ael_stiffness_tensor.rows == 6) && (ael_stiffness_tensor.cols == 6)) {
        sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_stiffness_tensor="; //CO20200624 - a lot of empty ael_*_tensor appearing in aflowlib.out
        for (int i = 1; i <= 6; i++) {
          for (int j = 1; j <= 6; j++) sss << ael_stiffness_tensor[i][j] << ((j < 6)?",":"");
          sss << ((i < 6)?";":"") << eendl;
        }
      }
      if ((ael_compliance_tensor.rows == 6) && (ael_compliance_tensor.cols == 6)) {
        sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "ael_compliance_tensor=";  //CO20200624 - a lot of empty ael_*_tensor appearing in aflowlib.out
        for (int i = 1; i <= 6; i++) {
          for (int j = 1; j <= 6; j++) {
            sss << ael_compliance_tensor[i][j] << ((j < 6)?",":"");
          }
          sss << ((i < 6)?";":"") << eendl;
        }
      }
      //ME20191105 END
      //ME20210927 BEGIN
      //APL
      if (energy_free_vibrational_cell_apl_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_free_vibrational_cell_apl_300K=" << energy_free_vibrational_cell_apl_300K;
      if (energy_free_vibrational_atom_apl_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_free_vibrational_atom_apl_300K=" << energy_free_vibrational_atom_apl_300K;
      if (entropy_vibrational_cell_apl_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "entropy_vibrational_cell_apl_300K=" << entropy_vibrational_cell_apl_300K;
      if (entropy_vibrational_atom_apl_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "entropy_vibrational_atom_apl_300K=" << entropy_vibrational_atom_apl_300K;
      if (energy_internal_vibrational_cell_apl_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_internal_vibrational_cell_apl_300K=" << energy_internal_vibrational_cell_apl_300K;
      if (energy_internal_vibrational_atom_apl_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_internal_vibrational_atom_apl_300K=" << energy_internal_vibrational_atom_apl_300K;
      if (energy_zero_point_cell_apl!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_zero_point_cell_apl=" << energy_zero_point_cell_apl;
      if (energy_zero_point_atom_apl!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_zero_point_atom_apl=" << energy_zero_point_atom_apl;
      if (heat_capacity_Cv_cell_apl_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "heat_capacity_Cv_cell_apl_300K=" << heat_capacity_Cv_cell_apl_300K;
      if (heat_capacity_Cv_atom_apl_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "heat_capacity_Cv_atom_apl_300K=" << heat_capacity_Cv_atom_apl_300K;
      //ME20210927 END
      //AS20200901 BEGIN
      //QHA
      if(gruneisen_qha!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "gruneisen_qha=" << gruneisen_qha << eendl; //AS20200901
      if(gruneisen_qha_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "gruneisen_qha_300K=" << gruneisen_qha_300K << eendl; //AS20200901
      if(thermal_expansion_qha_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "thermal_expansion_qha_300K=" << thermal_expansion_qha_300K << eendl; //AS20200901
      if(modulus_bulk_qha_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "modulus_bulk_qha_300K=" << modulus_bulk_qha_300K << eendl; //AS20200901
      //AS20200901 END
      //AS20201008 BEGIN
      if(modulus_bulk_derivative_pressure_qha_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "modulus_bulk_derivative_pressure_qha_300K=" << modulus_bulk_derivative_pressure_qha_300K << eendl; //AS20201008
      if(heat_capacity_Cv_atom_qha_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "heat_capacity_Cv_atom_qha_300K=" << heat_capacity_Cv_atom_qha_300K << eendl; //AS20201008
      if(heat_capacity_Cv_cell_qha_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "heat_capacity_Cv_cell_qha_300K=" << heat_capacity_Cv_cell_qha_300K << eendl; //AS20201207
      if(heat_capacity_Cp_atom_qha_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "heat_capacity_Cp_atom_qha_300K=" << heat_capacity_Cp_atom_qha_300K << eendl; //AS20201008
      if(heat_capacity_Cp_cell_qha_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "heat_capacity_Cp_cell_qha_300K=" << heat_capacity_Cp_cell_qha_300K << eendl; //AS20201207
      if(volume_atom_qha_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "volume_atom_qha_300K=" << volume_atom_qha_300K << eendl; //AS20201008
      if(energy_free_atom_qha_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_free_atom_qha_300K=" << energy_free_atom_qha_300K << eendl; //AS20201008
      if(energy_free_cell_qha_300K!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy_free_cell_qha_300K=" << energy_free_cell_qha_300K << eendl; //AS20201207
      //AS20201008 END
      // BADER
      if(bader_net_charges.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "bader_net_charges=" << bader_net_charges << eendl;
      if(bader_atomic_volumes.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "bader_atomic_volumes=" << bader_atomic_volumes << eendl;
      // FILES
      if(vfiles.size()) {
        aurostd::sort(vfiles);
        sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "files=";
        for(uint i=0;i<vfiles.size();i++) sss << vfiles.at(i) << (i<vfiles.size()-1?",":"");
        sss << eendl;
      }
      // CPUS
      if(node_CPU_Model.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "node_CPU_Model=" << node_CPU_Model << eendl;
      if(node_CPU_Cores!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "node_CPU_Cores=" << node_CPU_Cores << eendl;
      if(node_CPU_MHz!=AUROSTD_NAN) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "node_CPU_MHz=" << node_CPU_MHz << eendl;
      if(node_RAM_GB!=INF) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "node_RAM_GB=" << node_RAM_GB << eendl;
      // VERSION/DATE
      if(aflow_version.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "aflow_version=" << aflow_version << eendl;
      if(catalog.size()) sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "catalog=" << catalog << eendl;
      sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "aflowlib_version=" << string(AFLOW_VERSION) << eendl;
      sss << _AFLOWLIB_ENTRY_SEPARATOR_ << "aflowlib_date=" << aurostd::joinWDelimiter(vaflowlib_date,",") << eendl;  //CO20200624 - adding LOCK date
      sss << endl;

    } // out

    // this is the aflowlib.json mode
    if(mode=="json" || mode=="JSON") {  //CO OPERATE HERE ALL THE STRINGS AS BEFORE
      string eendl=",";
      stringstream sscontent_json;
      vector<string> vcontent_json;
      vector<string> sg_tokens;
      stringstream ss_helper;
      vector<vector<string> > vvs;
      vector<string> vs, vs2, vcontent_tmp; //DX20210129 - added vs2, content_tmp
      bool odd_xvec_count;

      //////////////////////////////////////////////////////////////////////////
      if(auid.size()) {
        sscontent_json << "\"aurl\":\"" << aurl << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"aurl\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(auid.size()) {
        sscontent_json << "\"auid\":\"" << auid << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"auid\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //ME20190125 BEGIN
      //////////////////////////////////////////////////////////////////////////
      // OBSOLETE ME20200829 - not used anymore
      //if(!title.empty()) {
      //  sscontent_json << "\"title\":\"" << title << "\"" << eendl;
      //} else {
      //  if(PRINT_NULL) sscontent_json << "\"title\":null" << eendl;
      //}
      //vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //////////////////////////////////////////////////////////////////////////
      //ME20190125 END

      //////////////////////////////////////////////////////////////////////////
      if(data_api.size()) {
        sscontent_json << "\"data_api\":\"" << data_api << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"data_api\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(data_source.size()) {
        sscontent_json << "\"data_source\":\"" << data_source << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"data_source\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(data_language.size()) {
        sscontent_json << "\"data_language\":\"" << data_language << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"data_language\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(error_status.size()) {
        sscontent_json << "\"error_status\":\"" << error_status << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"error_status\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      // LOOP
      //////////////////////////////////////////////////////////////////////////
      if(vloop.size()) {
        aurostd::sort(vloop);
        sscontent_json << "\"loop\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vloop,"\""),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"loop\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      // MATERIALS
      //////////////////////////////////////////////////////////////////////////
      if(code.size()) {
        sscontent_json << "\"code\":\"" << code << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"code\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(compound.size()) {
        sscontent_json << "\"compound\":\"" << compound << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"compound\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(prototype.size()) {
        sscontent_json << "\"prototype\":\"" << prototype << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"prototype\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(nspecies!=AUROSTD_NAN) {
        sscontent_json << "\"nspecies\":" << nspecies;
      } else {
        if(PRINT_NULL) sscontent_json << "\"nspecies\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(natoms!=AUROSTD_NAN) {
        sscontent_json << "\"natoms\":" << natoms;
      } else {
        if(PRINT_NULL) sscontent_json << "\"natoms\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //DX20190124 - add original crystal info - START
      //////////////////////////////////////////////////////////////////////////
      if(natoms_orig!=AUROSTD_NAN) {
        sscontent_json << "\"natoms_orig\":" << natoms_orig;
      } else {
        if(PRINT_NULL) sscontent_json << "\"natoms_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////
      //DX20190124 - add original crystal info - END

      //////////////////////////////////////////////////////////////////////////
      if(vcomposition.size()) {
        //aflowlib_libraries does not specify precision
        sscontent_json << "\"composition\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vcomposition,9),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"composition\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(density!=AUROSTD_NAN) {
        sscontent_json << "\"density\":" << density;
      } else {
        if(PRINT_NULL) sscontent_json << "\"density\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //DX20190124 - add original crystal info - START
      //////////////////////////////////////////////////////////////////////////
      if(density_orig!=AUROSTD_NAN) {
        sscontent_json << "\"density_orig\":" << density_orig;
      } else {
        if(PRINT_NULL) sscontent_json << "\"density_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////
      //DX20190124 - add original crystal info - END

      //////////////////////////////////////////////////////////////////////////
      if(scintillation_attenuation_length!=AUROSTD_NAN) {
        sscontent_json << "\"scintillation_attenuation_length\":" << scintillation_attenuation_length;
      } else {
        if(PRINT_NULL) sscontent_json << "\"scintillation_attenuation_length\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vstoichiometry.size()) {
        //aflowlib_libraries specifies precision of 9
        sscontent_json << "\"stoichiometry\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vstoichiometry,9),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"stoichiometry\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspecies.size()) {
        sscontent_json << "\"species\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies,"\""),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"species\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspecies_pp.size()) {
        sscontent_json << "\"species_pp\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies_pp,"\""),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"species_pp\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(dft_type.size()) {
        //DX+CO START
        //sscontent_json << "\"dft_type\":\"" << dft_type << "\""; //CO, this is technically a vector (RESTAPI paper)
        sscontent_json << "\"dft_type\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vdft_type,"\""),",") << "]";
        //DX+CO END
      } else {
        if(PRINT_NULL) sscontent_json << "\"dft_type\":null" << dft_type;
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      // if(species_pp_type.size()) sscontent_json << "species_pp_type=" << species_pp_type;

      //////////////////////////////////////////////////////////////////////////
      if(vspecies_pp_version.size()) {
        sscontent_json << "\"species_pp_version\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies_pp_version,"\""),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"species_pp_version\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspecies_pp_ZVAL.size()) {
        //aflowlib_libraries does not specify precision
        sscontent_json << "\"species_pp_ZVAL\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vspecies_pp_ZVAL,9),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"species_pp_ZVAL\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspecies_pp_AUID.size()) {
        sscontent_json << "\"species_pp_AUID\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies_pp_AUID,"\""),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"species_pp_AUID\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(METAGGA.size()) {
        sscontent_json << "\"metagga\":\"" << METAGGA << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"metagga\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      //ME20190124 - Modified to include more detailed LDAU information
      if(vLDAU.size()>0 && vLDAU[0].size()) {  //ldau_TLUJ.size()
        ss_helper.str("");
        vs.clear();
        //only string is available, so we have to parse really fast
        //would be nice if we have vldau_TLUJ already
        //ME20190124 - vLDAU is now available
        //[OBSOLETE ME20190124 - use vLDAU] vector<string> ldau_TLUJ_tokens;
        //[OBSOLETE ME20190124 - use vLDAU] aurostd::string2tokens(ldau_TLUJ,ldau_TLUJ_tokens,";");
        //[OBSOLETE ME20190124 - use vLDAU] if(ldau_TLUJ_tokens.size()==4){
        //[OBSOLETE ME20190124 - use vLDAU] conversion to double ENSURES that these are numbers
        //[OBSOLETE ME20190124 - use vLDAU] non-numbers without "" will break json
        //[OBSOLETE ME20190124 - use vLDAU] int T=aurostd::string2utype<int>(ldau_TLUJ_tokens.at(0));
        int T = (int) vLDAU[0][0]; //ME20190124
        vector<int> L(vLDAU[1].begin(), vLDAU[1].end());  //ME20190124
        vector<double> U = vLDAU[2],J = vLDAU[3];  //ME20190124
        //[OBSOLETE ME20190124 - NOT USED] vector<string> ldau_TLUJ_tokens2;
        //[OBSOLETE ME20190124 - use vLDAU] breaking up not necessary, but a nice check that we don't have hanging commas
        //[OBSOLETE ME20190124 - use vLDAU] the extra space at the end will be removed by joinWDelimiter()
        //[OBSOLETE ME20190124 - use vLDAU] aurostd::string2tokens(ldau_TLUJ_tokens.at(1),L,",");
        //[OBSOLETE ME20190124 - use vLDAU] aurostd::string2tokens(ldau_TLUJ_tokens.at(2),U,",");
        //[OBSOLETE ME20190124 - use vLDAU] aurostd::string2tokens(ldau_TLUJ_tokens.at(3),J,",");
        vs.push_back(aurostd::utype2string(T));
        if(L.size()&&U.size()&&J.size()){
          //no precision needed
          vs.push_back("["+aurostd::joinWDelimiter(L,",")+"]");
          vs.push_back("["+aurostd::joinWDelimiter(aurostd::vecDouble2vecString(U,9),",")+"]");
          vs.push_back("["+aurostd::joinWDelimiter(aurostd::vecDouble2vecString(J,9),",")+"]");
        }
        if(T!=0){ss_helper << aurostd::joinWDelimiter(vs,",");} //CO20210713
        vector<string> ldau_keys;  //ME20190124
        aurostd::string2tokens("ldau_type,ldau_l,ldau_u,ldau_j", ldau_keys, ",");  //ME20190124
        for (uint i = 0; i < ldau_keys.size() && i < vs.size(); i++) {  //ME20190124  //CO20210713
          sscontent_json << "\"" << ldau_keys[i] << "\":" << vs[i];  //ME20190124
          vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);  //ME20190124
        }
        //} ME20190124
        vs.clear();
      }
      if(!ss_helper.str().empty()){ //CO20180216 - !empty() is better for strings than !size()
        sscontent_json << "\"ldau_TLUJ\":[" << ss_helper.str() << "]"; ss_helper.str("");
      } else {
        if(PRINT_NULL) sscontent_json << "\"ldau_TLUJ\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(valence_cell_iupac!=AUROSTD_NAN) {
        sscontent_json << "\"valence_cell_iupac\":" << valence_cell_iupac;
      } else {
        if(PRINT_NULL) sscontent_json << "\"valence_cell_iupac\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(valence_cell_std!=AUROSTD_NAN) {
        sscontent_json << "\"valence_cell_std\":" << valence_cell_std;
      } else {
        if(PRINT_NULL) sscontent_json << "\"valence_cell_std\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(volume_cell!=AUROSTD_NAN) {
        sscontent_json << "\"volume_cell\":" << volume_cell;
      } else {
        if(PRINT_NULL) sscontent_json << "\"volume_cell\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(volume_atom!=AUROSTD_NAN) {
        sscontent_json << "\"volume_atom\":" << volume_atom;
      } else {
        if(PRINT_NULL) sscontent_json << "\"volume_atom\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //DX20190124 - add original crystal info - START
      //////////////////////////////////////////////////////////////////////////
      if(volume_cell_orig!=AUROSTD_NAN) {
        sscontent_json << "\"volume_cell_orig\":" << volume_cell_orig;
      } else {
        if(PRINT_NULL) sscontent_json << "\"volume_cell_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(volume_atom_orig!=AUROSTD_NAN) {
        sscontent_json << "\"volume_atom_orig\":" << volume_atom_orig;
      } else {
        if(PRINT_NULL) sscontent_json << "\"volume_atom_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////
      //DX20190124 - add original crystal info - END

      //////////////////////////////////////////////////////////////////////////
      if(pressure!=AUROSTD_NAN) {
        sscontent_json << "\"pressure\":" << pressure;
      } else {
        if(PRINT_NULL) sscontent_json << "\"pressure\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vstress_tensor.size()) {
        //aflowlib_libraries specifies precision of 7
        sscontent_json << "\"stress_tensor\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vstress_tensor,7),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"stress_tensor\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(pressure_residual!=AUROSTD_NAN) {
        sscontent_json << "\"pressure_residual\":" << pressure_residual;
      } else {
        if(PRINT_NULL) sscontent_json << "\"pressure_residual\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Pulay_stress!=AUROSTD_NAN) {
        sscontent_json << "\"Pulay_stress\":" << Pulay_stress;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Pulay_stress\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vgeometry.size()) {
        //aflowlib_libraries specifies precision of 7
        sscontent_json << "\"geometry\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vgeometry,7),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"geometry\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //DX20190124 - add original crystal info - START
      //////////////////////////////////////////////////////////////////////////
      if(vgeometry_orig.size()) {
        //aflowlib_libraries specifies precision of 7
        sscontent_json << "\"geometry_orig\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vgeometry_orig,7),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"geometry_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////
      //DX20190124 - add original crystal info - END

      //////////////////////////////////////////////////////////////////////////
      if(Egap!=AUROSTD_NAN) {
        sscontent_json << "\"Egap\":" << Egap;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Egap\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Egap_fit!=AUROSTD_NAN) {
        sscontent_json << "\"Egap_fit\":" << Egap_fit;
      } else {
        if(PRINT_NULL) sscontent_json << "\"Egap_fit\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Egap_type.size()) {
        sscontent_json << "\"Egap_type\":\"" << Egap_type << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Egap_type\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(energy_cell!=AUROSTD_NAN) {
        sscontent_json << "\"energy_cell\":" << energy_cell;
      } else {
        if(PRINT_NULL) sscontent_json << "\"energy_cell\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(energy_atom!=AUROSTD_NAN) {
        sscontent_json << "\"energy_atom\":" << energy_atom;
      } else {
        if(PRINT_NULL) sscontent_json << "\"energy_atom\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(energy_cutoff!=AUROSTD_NAN) {
        sscontent_json << "\"energy_cutoff\":" << energy_cutoff;
      } else {
        if(PRINT_NULL) sscontent_json << "\"energy_cutoff\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(delta_electronic_energy_convergence!=AUROSTD_NAN) {
        sscontent_json << "\"delta_electronic_energy_convergence\":" << delta_electronic_energy_convergence;
      } else {
        if(PRINT_NULL) sscontent_json << "\"delta_electronic_energy_convergence\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(delta_electronic_energy_threshold!=AUROSTD_NAN) {
        sscontent_json << "\"delta_electronic_energy_threshold\":" << delta_electronic_energy_threshold;
      } else {
        if(PRINT_NULL) sscontent_json << "\"delta_electronic_energy_threshold\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////
      // [NOT_PRINTED]     if(nkpoints!=0) {
      // [NOT_PRINTED]       sscontent_json << "\"nkpoints\":" << nkpoints;
      // [NOT_PRINTED]     } else {
      // [NOT_PRINTED]       if(PRINT_NULL) sscontent_json << "\"nkpoints\":null";
      // [NOT_PRINTED]     }
      // [NOT_PRINTED]     vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////

      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////
      // [NOT_PRINTED]     if(nkpoints_irreducible!=0) {
      // [NOT_PRINTED]       sscontent_json << "\"nkpoints_irreducible\":" << nkpoints_irreducible;
      // [NOT_PRINTED]     } else {
      // [NOT_PRINTED]       if(PRINT_NULL) sscontent_json << "\"nkpoints_irreducible\":null";
      // [NOT_PRINTED]     }
      // [NOT_PRINTED]     vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////

      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////
      // [NOT_PRINTED]     if(kppra!=0) {
      // [NOT_PRINTED]       sscontent_json << "\"kppra\":" << kppra;
      // [NOT_PRINTED]     } else {
      // [NOT_PRINTED]       if(PRINT_NULL) sscontent_json << "\"kppra\":null";
      // [NOT_PRINTED]     }
      // [NOT_PRINTED]     vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      // [NOT_PRINTED]     //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(kpoints.size()) {
        //this one is a bit complicated, so we will test if the string was created, and then recreate the json array on the spot
        ss_helper.str("");
        vs.clear();
        //ME20190124 - Add the individual pieces of "kpoints" to the json file
        if ((kpoints_nnn_relax.rows==3) && (sum(kpoints_nnn_relax) > 0)) {  //ME20190128
          vs.push_back("["+aurostd::joinWDelimiter(kpoints_nnn_relax,",")+"]");
          sscontent_json << "\"kpoints_relax\":" << vs.back();  //ME20190124
        } else if (PRINT_NULL) {  //ME20190124
          sscontent_json << "\"kpoints_relax\":null";  //ME20190124
        }
        vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);  //ME20190124
        if ((kpoints_nnn_static.rows==3) && (sum(kpoints_nnn_static) > 0)) {  //ME20190128
          vs.push_back("["+aurostd::joinWDelimiter(kpoints_nnn_static,",")+"]");
          if (sum(kpoints_nnn_static) > 0) sscontent_json << "\"kpoints_static\":" << vs.back();  //ME20190124
        } else if (PRINT_NULL) {  //ME20190124
          sscontent_json << "\"kpoints_static\":null";  //ME20190124
        }
        vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);  //ME20190124
        if(kpoints_pairs.size()){
          //first for escape characters in \Gamma or \Sigma
          vector<string> kpoints_pairs_new;
          char issue_c='\\';
          stringstream issue_ss; issue_ss << issue_c;
          string fix_s="\\\\";
          for(uint i=0;i<kpoints_pairs.size();i++){
            kpoints_pairs_new.push_back(aurostd::StringSubst(kpoints_pairs.at(i),issue_ss.str(),fix_s));
          }
          vs.push_back("["+aurostd::joinWDelimiter(aurostd::wrapVecEntries(kpoints_pairs_new,"\""),",")+"]");
          sscontent_json << "\"kpoints_bands_path\":" << vs.back();  //ME20190124
        } else if (PRINT_NULL) {  //ME20190124
          sscontent_json << "\"kpoints_bands_path\":null";  //ME20190124
        }
        vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);  //ME20190124
        ss_helper << aurostd::joinWDelimiter(vs, ",");  //ME20190128
        if(kpoints_bands_path_grid!=0){
          //ME20190128 - This causes kpoints to only be written when the band structure
          // was calculated. This is inconsistent with the aflowlib.out file
          // [OBSOLETE ME20190128] ss_helper << aurostd::joinWDelimiter(vs,",") << "," << aurostd::utype2string(kpoints_bands_path_grid);
          ss_helper << "," << aurostd::utype2string(kpoints_bands_path_grid);
          sscontent_json << "\"kpoints_bands_nkpts\":" << ((int) kpoints_bands_path_grid);  //ME20190124
        } else if (PRINT_NULL) {  //ME20190124
          sscontent_json << "\"kpoints_bands_nkpts\":null";  //ME20190124
        }
        vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);  //ME20190124
      }
      if(!ss_helper.str().empty()){ //CO20180216 - !empty() is better for strings than !size()
        sscontent_json << "\"kpoints\":[" << ss_helper.str() << "]"; ss_helper.str("");
      } else {
        if(PRINT_NULL) sscontent_json << "\"kpoints\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_cell!=AUROSTD_NAN) {
        sscontent_json << "\"enthalpy_cell\":" << enthalpy_cell;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_cell\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_atom!=AUROSTD_NAN) {
        sscontent_json << "\"enthalpy_atom\":" << enthalpy_atom;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_atom\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(eentropy_cell!=AUROSTD_NAN) {
        sscontent_json << "\"eentropy_cell\":" << eentropy_cell;
      } else {
        if(PRINT_NULL) sscontent_json << "\"eentropy_cell\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(eentropy_atom!=AUROSTD_NAN) {
        sscontent_json << "\"eentropy_atom\":" << eentropy_atom;
      } else {
        if(PRINT_NULL) sscontent_json << "\"eentropy_atom\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_formation_cell!=AUROSTD_NAN) {
        sscontent_json << "\"enthalpy_formation_cell\":" << enthalpy_formation_cell;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_formation_cell\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_formation_cce_300K_cell!=AUROSTD_NAN) {  //CO20200624
        sscontent_json << "\"enthalpy_formation_cce_300K_cell\":" << enthalpy_formation_cce_300K_cell;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_formation_cce_300K_cell\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_formation_cce_0K_cell!=AUROSTD_NAN) {  //CO20200624
        sscontent_json << "\"enthalpy_formation_cce_0K_cell\":" << enthalpy_formation_cce_0K_cell;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_formation_cce_0K_cell\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_formation_atom!=AUROSTD_NAN) {
        sscontent_json << "\"enthalpy_formation_atom\":" << enthalpy_formation_atom;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_formation_atom\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_formation_cce_300K_atom!=AUROSTD_NAN) {  //CO20200624
        sscontent_json << "\"enthalpy_formation_cce_300K_atom\":" << enthalpy_formation_cce_300K_atom;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_formation_cce_300K_atom\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(enthalpy_formation_cce_0K_atom!=AUROSTD_NAN) {  //CO20200624
        sscontent_json << "\"enthalpy_formation_cce_0K_atom\":" << enthalpy_formation_cce_0K_atom;
      } else {
        if(PRINT_NULL) sscontent_json << "\"enthalpy_formation_cce_0K_atom\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(entropy_forming_ability!=AUROSTD_NAN) {  //CO20200624
        sscontent_json << "\"entropy_forming_ability\":" << entropy_forming_ability;
      } else {
        if(PRINT_NULL) sscontent_json << "\"entropy_forming_ability\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(entropic_temperature!=AUROSTD_NAN) {
        sscontent_json << "\"entropic_temperature\":" << entropic_temperature;
      } else {
        if(PRINT_NULL) sscontent_json << "\"entropic_temperature\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(PV_cell!=AUROSTD_NAN) {
        sscontent_json << "\"PV_cell\":" << PV_cell;
      } else {
        if(PRINT_NULL) sscontent_json << "\"PV_cell\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(PV_atom!=AUROSTD_NAN) {
        sscontent_json << "\"PV_atom\":" << PV_atom;
      } else {
        if(PRINT_NULL) sscontent_json << "\"PV_atom\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(spin_cell!=AUROSTD_NAN) {
        sscontent_json << "\"spin_cell\":" << spin_cell;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spin_cell\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(spin_atom!=AUROSTD_NAN) {
        sscontent_json << "\"spin_atom\":" << spin_atom;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spin_atom\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspinD.size()) {
        //aflowlib_libraries specifies precision of 5
        sscontent_json << "\"spinD\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vspinD,5),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"spinD\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vspinD_magmom_orig.size()) {
        //aflowlib_libraries specifies precision of 5
        sscontent_json << "\"spinD_magmom_orig\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vspinD_magmom_orig,5),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"spinD_magmom_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(spinF!=AUROSTD_NAN) {
        sscontent_json << "\"spinF\":" << spinF;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spinF\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      // [OBSOLETE] 
      // [OBSOLETE] if(stoich.size()) {
      // [OBSOLETE]   //just use the string SC made
      // [OBSOLETE]   sscontent_json << "\"stoich\":\"" << stoich << "\"";
      // [OBSOLETE]   ////aflowlib_libraries does not specify precision
      // [OBSOLETE]   //sscontent_json << "\"stoich\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vstoich,9),",") << "]";
      // [OBSOLETE] } else {
      // [OBSOLETE]   if(PRINT_NULL) sscontent_json << "\"stoich\":null";
      // [OBSOLETE] }
      // [OBSOLETE] vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(calculation_time!=AUROSTD_NAN) {
        sscontent_json << "\"calculation_time\":" << calculation_time;
      } else {
        if(PRINT_NULL) sscontent_json << "\"calculation_time\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(calculation_memory!=AUROSTD_NAN) {
        sscontent_json << "\"calculation_memory\":" << calculation_memory;
      } else {
        if(PRINT_NULL) sscontent_json << "\"calculation_memory\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(calculation_cores!=AUROSTD_NAN) {
        sscontent_json << "\"calculation_cores\":" << calculation_cores;
      } else {
        if(PRINT_NULL) sscontent_json << "\"calculation_cores\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vnbondxx.size()) {
        //aflowlib_libraries does not specify precision
        sscontent_json << "\"nbondxx\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vnbondxx,9),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"nbondxx\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(sg.size()) {
        aurostd::string2tokens(sg,sg_tokens,",");
        sscontent_json << "\"sg\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(sg_tokens,"\""),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"sg\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(sg2.size()) {
        aurostd::string2tokens(sg2,sg_tokens,",");
        sscontent_json << "\"sg2\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(sg_tokens,"\""),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"sg2\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(spacegroup_orig!=AUROSTD_NAN) {  //CO20201111
        sscontent_json << "\"spacegroup_orig\":" << spacegroup_orig;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spacegroup_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(spacegroup_relax!=AUROSTD_NAN) { //CO20201111
        sscontent_json << "\"spacegroup_relax\":" << spacegroup_relax;
      } else {
        if(PRINT_NULL) sscontent_json << "\"spacegroup_relax\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vforces.size()) {
        ss_helper.str("");
        vs.clear();
        vvs.clear();
        odd_xvec_count=false;
        for(uint i=0;i<vforces.size();i++){
          if(vforces.at(i).rows==3){
            //aflowlib_libraries specifies precision of 8
            vvs.push_back(xvecDouble2vecString(vforces.at(i),8));
          } else {
            odd_xvec_count=true;
            break;
          }
        }
        if(!odd_xvec_count&&vvs.size()){
          for(uint i=0;i<vvs.size();i++){
            vs.push_back("["+aurostd::joinWDelimiter(vvs.at(i),",")+"]");
          }
          if(vs.size()){
            ss_helper << aurostd::joinWDelimiter(vs,",");
          }
        }
        vs.clear();
        vvs.clear();
      }
      if(!ss_helper.str().empty()){ //CO20180216 - !empty() is better for strings than !size()
        sscontent_json << "\"forces\":[" << ss_helper.str() << "]"; ss_helper.str("");
      } else {
        if(PRINT_NULL) sscontent_json << "\"forces\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vpositions_cartesian.size()) {
        ss_helper.str("");
        vs.clear();
        vvs.clear();
        odd_xvec_count=false;
        for(uint i=0;i<vpositions_cartesian.size();i++){
          if(vpositions_cartesian.at(i).rows==3){
            //aflowlib_libraries specifies precision of 8
            vvs.push_back(xvecDouble2vecString(vpositions_cartesian.at(i),8));
          } else {
            odd_xvec_count=true;
            break;
          }
        }
        if(!odd_xvec_count&&vvs.size()){
          for(uint i=0;i<vvs.size();i++){
            vs.push_back("["+aurostd::joinWDelimiter(vvs.at(i),",")+"]");
          }
          if(vs.size()){
            ss_helper << aurostd::joinWDelimiter(vs,",");
          }
        }
        vs.clear();
        vvs.clear();
      }
      if(!ss_helper.str().empty()){ //CO20180216 - !empty() is better for strings than !size()
        sscontent_json << "\"positions_cartesian\":[" << ss_helper.str() << "]"; ss_helper.str("");
      } else {
        if(PRINT_NULL) sscontent_json << "\"positions_cartesian\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vpositions_fractional.size()) {
        ss_helper.str("");
        vs.clear();
        vvs.clear();
        odd_xvec_count=false;
        for(uint i=0;i<vpositions_fractional.size();i++){
          if(vpositions_fractional.at(i).rows==3){
            //aflowlib_libraries specifies precision of 8
            vvs.push_back(xvecDouble2vecString(vpositions_fractional.at(i),8));
          } else {
            odd_xvec_count=true;
            break;
          }
        }
        if(!odd_xvec_count&&vvs.size()){
          for(uint i=0;i<vvs.size();i++){
            vs.push_back("["+aurostd::joinWDelimiter(vvs.at(i),",")+"]");
          }
          if(vs.size()){
            ss_helper << aurostd::joinWDelimiter(vs,",");
          }
        }
        vs.clear();
        vvs.clear();
      }
      if(!ss_helper.str().empty()){ //CO20180216 - !empty() is better for strings than !size()
        sscontent_json << "\"positions_fractional\":[" << ss_helper.str() << "]"; ss_helper.str("");
      } else {
        if(PRINT_NULL) sscontent_json << "\"positions_fractional\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_orig.size()) {
        sscontent_json << "\"Bravais_lattice_orig\":\"" << Bravais_lattice_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(lattice_variation_orig.size()) {
        sscontent_json << "\"lattice_variation_orig\":\"" << lattice_variation_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"lattice_variation_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(lattice_system_orig.size()) {
        sscontent_json << "\"lattice_system_orig\":\"" << lattice_system_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"lattice_system_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Pearson_symbol_orig.size()) {
        sscontent_json << "\"Pearson_symbol_orig\":\"" << Pearson_symbol_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Pearson_symbol_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_relax.size()) {
        sscontent_json << "\"Bravais_lattice_relax\":\"" << Bravais_lattice_relax << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_relax\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(lattice_variation_relax.size()) {
        sscontent_json << "\"lattice_variation_relax\":\"" << lattice_variation_relax << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"lattice_variation_relax\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(lattice_system_relax.size()) {
        sscontent_json << "\"lattice_system_relax\":\"" << lattice_system_relax << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"lattice_system_relax\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(Pearson_symbol_relax.size()) {
        sscontent_json << "\"Pearson_symbol_relax\":\"" << Pearson_symbol_relax << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Pearson_symbol_relax\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //DX20190124 - added original symmetry info - START
      // SYMMETRY
      //////////////////////////////////////////////////////////////////////////
      if(crystal_family_orig.size()){
        sscontent_json << "\"crystal_family_orig\":\"" << crystal_family_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"crystal_family_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(crystal_system_orig.size()){
        sscontent_json << "\"crystal_system_orig\":\"" << crystal_system_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"crystal_system_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(crystal_class_orig.size()){
        sscontent_json << "\"crystal_class_orig\":\"" << crystal_class_orig << "\"";
      } else {
        if(PRINT_NULL){ sscontent_json << "\"crystal_class_orig\":null";}
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(point_group_Hermann_Mauguin_orig.size()){
        sscontent_json << "\"point_group_Hermann_Mauguin_orig\":\"" << point_group_Hermann_Mauguin_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_Hermann_Mauguin_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(point_group_Schoenflies_orig.size()){
        sscontent_json << "\"point_group_Schoenflies_orig\":\"" << point_group_Schoenflies_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_Schoenflies_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(point_group_orbifold_orig.size()){
        sscontent_json << "\"point_group_orbifold_orig\":\"" << point_group_orbifold_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_orbifold_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(point_group_type_orig.size()){
        sscontent_json << "\"point_group_type_orig\":\"" << point_group_type_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_type_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(point_group_order_orig!=AUROSTD_NAN){
        sscontent_json << "\"point_group_order_orig\":" << point_group_order_orig; //DX20190124 - changed to number, not string 
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_order_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(point_group_structure_orig.size()){
        sscontent_json << "\"point_group_structure_orig\":\"" << point_group_structure_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_structure_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_lattice_type_orig.size()){
        sscontent_json << "\"Bravais_lattice_lattice_type_orig\":\"" << Bravais_lattice_lattice_type_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_lattice_type_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_lattice_variation_type_orig.size()){
        sscontent_json << "\"Bravais_lattice_lattice_variation_type_orig\":\"" << Bravais_lattice_lattice_variation_type_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_lattice_variation_type_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_lattice_system_orig.size()){
        sscontent_json << "\"Bravais_lattice_lattice_system_orig\":\"" << Bravais_lattice_lattice_system_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_lattice_system_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_superlattice_lattice_type_orig.size()){
        sscontent_json << "\"Bravais_superlattice_lattice_type_orig\":\"" << Bravais_superlattice_lattice_type_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_superlattice_lattice_type_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_superlattice_lattice_variation_type_orig.size()){
        sscontent_json << "\"Bravais_superlattice_lattice_variation_type_orig\":\"" << Bravais_superlattice_lattice_variation_type_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_superlattice_lattice_variation_type_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_superlattice_lattice_system_orig.size()){
        sscontent_json << "\"Bravais_superlattice_lattice_system_orig\":\"" << Bravais_superlattice_lattice_system_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_superlattice_lattice_system_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Pearson_symbol_superlattice_orig.size()){
        sscontent_json << "\"Pearson_symbol_superlattice_orig\":\"" << Pearson_symbol_superlattice_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Pearson_symbol_superlattice_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(vreciprocal_geometry_orig.size()) {
        //aflowlib_libraries specifies precision of 7
        sscontent_json << "\"reciprocal_geometry_orig\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vreciprocal_geometry_orig,7),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_geometry_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(reciprocal_volume_cell_orig!=AUROSTD_NAN) {
        sscontent_json << "\"reciprocal_volume_cell_orig\":" << reciprocal_volume_cell_orig;
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_volume_cell_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(reciprocal_lattice_type_orig.size()){
        sscontent_json << "\"reciprocal_lattice_type_orig\":\"" << reciprocal_lattice_type_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_lattice_type_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(reciprocal_lattice_variation_type_orig.size()){
        sscontent_json << "\"reciprocal_lattice_variation_type_orig\":\"" << reciprocal_lattice_variation_type_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_lattice_variation_type_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Wyckoff_letters_orig.size()){
        vs.clear();
        aurostd::string2tokens(Wyckoff_letters_orig,vs,";");
        vcontent_tmp.clear();
        for(uint w=0;w<vs.size();w++){
          vs2.clear();
          aurostd::string2tokens(vs[w],vs2,",");
          vcontent_tmp.push_back("["+aurostd::joinWDelimiter(aurostd::wrapVecEntries(vs2,"\""),",")+"]");
        }
        sscontent_json << "\"Wyckoff_letters_orig\":[" << aurostd::joinWDelimiter(vcontent_tmp,",") << "]";
        //DX20190129 [OBSOLETE] sscontent_json << "\"Wyckoff_letters_orig\":\"" << Wyckoff_letters_orig << "\"";
        vs.clear(); vs2.clear(); vcontent_tmp.clear();
      } else {
        if(PRINT_NULL) sscontent_json << "\"Wyckoff_letters_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Wyckoff_multiplicities_orig.size()){
        vs.clear();
        aurostd::string2tokens(Wyckoff_multiplicities_orig,vs,";");
        vcontent_tmp.clear();
        for(uint w=0;w<vs.size();w++){
          vs2.clear();
          aurostd::string2tokens(vs[w],vs2,",");
          vcontent_tmp.push_back("["+aurostd::joinWDelimiter(vs2,",")+"]");
        }
        sscontent_json << "\"Wyckoff_multiplicities_orig\":[" << aurostd::joinWDelimiter(vcontent_tmp,",") << "]";
        //DX20190129 [OBSOLETE] sscontent_json << "\"Wyckoff_multiplicities_orig\":\"" << Wyckoff_multiplicities_orig << "\"";
        vs.clear(); vs2.clear(); vcontent_tmp.clear();
      } else {
        if(PRINT_NULL) sscontent_json << "\"Wyckoff_multiplicities_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Wyckoff_site_symmetries_orig.size()){
        vs.clear();
        aurostd::string2tokens(Wyckoff_site_symmetries_orig,vs,";");
        vcontent_tmp.clear();
        for(uint w=0;w<vs.size();w++){
          vs2.clear();
          aurostd::string2tokens(vs[w],vs2,",");
          vcontent_tmp.push_back("["+aurostd::joinWDelimiter(aurostd::wrapVecEntries(vs2,"\""),",")+"]");
        }
        sscontent_json << "\"Wyckoff_site_symmetries_orig\":[" << aurostd::joinWDelimiter(vcontent_tmp,",") << "]";
        //DX20190129 [OBSOLETE] sscontent_json << "\"Wyckoff_site_symmetries_orig\":\"" << Wyckoff_site_symmetries_orig << "\"";
        vs.clear(); vs2.clear(); vcontent_tmp.clear();
      } else {
        if(PRINT_NULL) sscontent_json << "\"Wyckoff_site_symmetries_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);


      //DX20180823 - added more symmetry info - START
      // SYMMETRY
      //////////////////////////////////////////////////////////////////////////
      if(crystal_family.size()){
        sscontent_json << "\"crystal_family\":\"" << crystal_family << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"crystal_family\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(crystal_system.size()){
        sscontent_json << "\"crystal_system\":\"" << crystal_system << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"crystal_system\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(crystal_class.size()){
        sscontent_json << "\"crystal_class\":\"" << crystal_class << "\"";
      } else {
        if(PRINT_NULL){ sscontent_json << "\"crystal_class\":null";}
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(point_group_Hermann_Mauguin.size()){
        sscontent_json << "\"point_group_Hermann_Mauguin\":\"" << point_group_Hermann_Mauguin << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_Hermann_Mauguin\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(point_group_Schoenflies.size()){
        sscontent_json << "\"point_group_Schoenflies\":\"" << point_group_Schoenflies << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_Schoenflies\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(point_group_orbifold.size()){
        sscontent_json << "\"point_group_orbifold\":\"" << point_group_orbifold << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_orbifold\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(point_group_type.size()){
        sscontent_json << "\"point_group_type\":\"" << point_group_type << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_type\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(point_group_order!=AUROSTD_NAN){
        sscontent_json << "\"point_group_order\":" << point_group_order; //DX20190124 - changed to number, not string 
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_order\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(point_group_structure.size()){
        sscontent_json << "\"point_group_structure\":\"" << point_group_structure << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"point_group_structure\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_lattice_type.size()){
        sscontent_json << "\"Bravais_lattice_lattice_type\":\"" << Bravais_lattice_lattice_type << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_lattice_type\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_lattice_variation_type.size()){
        sscontent_json << "\"Bravais_lattice_lattice_variation_type\":\"" << Bravais_lattice_lattice_variation_type << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_lattice_variation_type\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_lattice_lattice_system.size()){
        sscontent_json << "\"Bravais_lattice_lattice_system\":\"" << Bravais_lattice_lattice_system << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_lattice_lattice_system\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_superlattice_lattice_type.size()){
        sscontent_json << "\"Bravais_superlattice_lattice_type\":\"" << Bravais_superlattice_lattice_type << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_superlattice_lattice_type\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_superlattice_lattice_variation_type.size()){
        sscontent_json << "\"Bravais_superlattice_lattice_variation_type\":\"" << Bravais_superlattice_lattice_variation_type << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_superlattice_lattice_variation_type\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Bravais_superlattice_lattice_system.size()){
        sscontent_json << "\"Bravais_superlattice_lattice_system\":\"" << Bravais_superlattice_lattice_system << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Bravais_superlattice_lattice_system\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Pearson_symbol_superlattice.size()){
        sscontent_json << "\"Pearson_symbol_superlattice\":\"" << Pearson_symbol_superlattice << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"Pearson_symbol_superlattice\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(vreciprocal_geometry_relax.size()) { //CO20220719 _relax
        //aflowlib_libraries specifies precision of 7
        sscontent_json << "\"reciprocal_geometry_relax\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vreciprocal_geometry_relax,7),",") << "]";  //CO20220719 _relax
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_geometry_relax\":null";  //CO20220719 _relax
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(reciprocal_volume_cell!=AUROSTD_NAN) {
        sscontent_json << "\"reciprocal_volume_cell\":" << reciprocal_volume_cell;
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_volume_cell\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(reciprocal_lattice_type.size()){
        sscontent_json << "\"reciprocal_lattice_type\":\"" << reciprocal_lattice_type << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_lattice_type\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(reciprocal_lattice_variation_type.size()){
        sscontent_json << "\"reciprocal_lattice_variation_type\":\"" << reciprocal_lattice_variation_type << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"reciprocal_lattice_variation_type\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Wyckoff_letters.size()){
        vs.clear();
        aurostd::string2tokens(Wyckoff_letters,vs,";");
        vcontent_tmp.clear();
        for(uint w=0;w<vs.size();w++){
          vs2.clear();
          aurostd::string2tokens(vs[w],vs2,",");
          vcontent_tmp.push_back("["+aurostd::joinWDelimiter(aurostd::wrapVecEntries(vs2,"\""),",")+"]");
        }
        sscontent_json << "\"Wyckoff_letters\":[" << aurostd::joinWDelimiter(vcontent_tmp,",") << "]";
        //DX20190129 [OBSOLETE] sscontent_json << "\"Wyckoff_letters\":\"" << Wyckoff_letters << "\"";
        vs.clear(); vs2.clear(); vcontent_tmp.clear();
      } else {
        if(PRINT_NULL) sscontent_json << "\"Wyckoff_letters\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Wyckoff_multiplicities.size()){
        vs.clear();
        aurostd::string2tokens(Wyckoff_multiplicities,vs,";");
        vcontent_tmp.clear();
        for(uint w=0;w<vs.size();w++){
          vs2.clear();
          aurostd::string2tokens(vs[w],vs2,",");
          vcontent_tmp.push_back("["+aurostd::joinWDelimiter(vs2,",")+"]");
        }
        sscontent_json << "\"Wyckoff_multiplicities\":[" << aurostd::joinWDelimiter(vcontent_tmp,",") << "]";
        //DX20190129 [OBSOLETE] sscontent_json << "\"Wyckoff_multiplicities\":\"" << Wyckoff_multiplicities << "\"";
        vs.clear(); vs2.clear(); vcontent_tmp.clear();
      } else {
        if(PRINT_NULL) sscontent_json << "\"Wyckoff_multiplicities\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(Wyckoff_site_symmetries.size()){
        vs.clear();
        aurostd::string2tokens(Wyckoff_site_symmetries,vs,";");
        vcontent_tmp.clear();
        for(uint w=0;w<vs.size();w++){
          vs2.clear();
          aurostd::string2tokens(vs[w],vs2,",");
          vcontent_tmp.push_back("["+aurostd::joinWDelimiter(aurostd::wrapVecEntries(vs2,"\""),",")+"]");
        }
        sscontent_json << "\"Wyckoff_site_symmetries\":[" << aurostd::joinWDelimiter(vcontent_tmp,",") << "]";
        //DX20190129 [OBSOLETE] sscontent_json << "\"Wyckoff_site_symmetries\":\"" << Wyckoff_site_symmetries << "\"";
        vs.clear(); vs2.clear(); vcontent_tmp.clear();
      } else {
        if(PRINT_NULL) sscontent_json << "\"Wyckoff_site_symmetries\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //DX20190208 - added anrl info - START
      // ANRL
      //////////////////////////////////////////////////////////////////////////
      if(aflow_prototype_label_orig.size()){
        sscontent_json << "\"aflow_prototype_label_orig\":\"" << aflow_prototype_label_orig << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"aflow_prototype_label_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(aflow_prototype_params_list_orig.size()){
        vs.clear(); aurostd::string2tokens(aflow_prototype_params_list_orig,vs,",");
        sscontent_json << "\"aflow_prototype_params_list_orig\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vs,"\""),",") << "]";
        vs.clear();
      } else {
        if(PRINT_NULL) sscontent_json << "\"aflow_prototype_params_list_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(aflow_prototype_params_values_orig.size()){
        vs.clear(); aurostd::string2tokens(aflow_prototype_params_values_orig,vs,",");
        sscontent_json << "\"aflow_prototype_params_values_orig\":[" << aurostd::joinWDelimiter(vs,",") << "]";
        vs.clear();
      } else {
        if(PRINT_NULL) sscontent_json << "\"aflow_prototype_params_values_orig\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(aflow_prototype_label_relax.size()){
        sscontent_json << "\"aflow_prototype_label_relax\":\"" << aflow_prototype_label_relax << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"aflow_prototype_label_relax\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(aflow_prototype_params_list_relax.size()){
        vs.clear(); aurostd::string2tokens(aflow_prototype_params_list_relax,vs,",");
        sscontent_json << "\"aflow_prototype_params_list_relax\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vs,"\""),",") << "]";
        vs.clear();
      } else {
        if(PRINT_NULL) sscontent_json << "\"aflow_prototype_params_list_relax\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////
      if(aflow_prototype_params_values_relax.size()){
        vs.clear(); aurostd::string2tokens(aflow_prototype_params_values_relax,vs,",");
        sscontent_json << "\"aflow_prototype_params_values_relax\":[" << aurostd::joinWDelimiter(vs,",") << "]";
        vs.clear();
      } else {
        if(PRINT_NULL) sscontent_json << "\"aflow_prototype_params_values_relax\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //DX20190208 - added anrl info - END

      //CO20200731
      //////////////////////////////////////////////////////////////////////////
      if(pocc_parameters.size()){
        sscontent_json << "\"pocc_parameters\":\"" << pocc_parameters << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"pocc_parameters\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      // AGL/AEL
      //////////////////////////////////////////////////////////////////////////
      if(agl_thermal_conductivity_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_thermal_conductivity_300K\":" << agl_thermal_conductivity_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_thermal_conductivity_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_debye!=AUROSTD_NAN) {
        sscontent_json << "\"agl_debye\":" << agl_debye;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_debye\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_acoustic_debye!=AUROSTD_NAN) {
        sscontent_json << "\"agl_acoustic_debye\":" << agl_acoustic_debye;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_acoustic_debye\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_gruneisen!=AUROSTD_NAN) {
        sscontent_json << "\"agl_gruneisen\":" << agl_gruneisen;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_gruneisen\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_heat_capacity_Cv_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_heat_capacity_Cv_300K\":" << agl_heat_capacity_Cv_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_heat_capacity_Cv_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_heat_capacity_Cp_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_heat_capacity_Cp_300K\":" << agl_heat_capacity_Cp_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_heat_capacity_Cp_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_thermal_expansion_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_thermal_expansion_300K\":" << agl_thermal_expansion_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_thermal_expansion_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_bulk_modulus_static_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_bulk_modulus_static_300K\":" << agl_bulk_modulus_static_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_bulk_modulus_static_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(agl_bulk_modulus_isothermal_300K!=AUROSTD_NAN) {
        sscontent_json << "\"agl_bulk_modulus_isothermal_300K\":" << agl_bulk_modulus_isothermal_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_bulk_modulus_isothermal_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(agl_poisson_ratio_source.size()) {
        sscontent_json << "\"agl_poisson_ratio_source\":\"" << agl_poisson_ratio_source << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_poisson_ratio_source\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(agl_vibrational_free_energy_300K_cell!=AUROSTD_NAN) {
        sscontent_json << "\"agl_vibrational_free_energy_300K_cell\":" << agl_vibrational_free_energy_300K_cell;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_vibrational_free_energy_300K_cell\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(agl_vibrational_free_energy_300K_atom!=AUROSTD_NAN) {
        sscontent_json << "\"agl_vibrational_free_energy_300K_atom\":" << agl_vibrational_free_energy_300K_atom;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_vibrational_free_energy_300K_atom\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(agl_vibrational_entropy_300K_cell!=AUROSTD_NAN) {
        sscontent_json << "\"agl_vibrational_entropy_300K_cell\":" << agl_vibrational_entropy_300K_cell;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_vibrational_entropy_300K_cell\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(agl_vibrational_entropy_300K_atom!=AUROSTD_NAN) {
        sscontent_json << "\"agl_vibrational_entropy_300K_atom\":" << agl_vibrational_entropy_300K_atom;
      } else {
        if(PRINT_NULL) sscontent_json << "\"agl_vibrational_entropy_300K_atom\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_poisson_ratio!=AUROSTD_NAN) {
        sscontent_json << "\"ael_poisson_ratio\":" << ael_poisson_ratio;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_poisson_ratio\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_bulk_modulus_voigt!=AUROSTD_NAN) {
        sscontent_json << "\"ael_bulk_modulus_voigt\":" << ael_bulk_modulus_voigt;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_bulk_modulus_voigt\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_bulk_modulus_reuss!=AUROSTD_NAN) {
        sscontent_json << "\"ael_bulk_modulus_reuss\":" << ael_bulk_modulus_reuss;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_bulk_modulus_reuss\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_shear_modulus_voigt!=AUROSTD_NAN) {
        sscontent_json << "\"ael_shear_modulus_voigt\":" << ael_shear_modulus_voigt;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_shear_modulus_voigt\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_shear_modulus_reuss!=AUROSTD_NAN) {
        sscontent_json << "\"ael_shear_modulus_reuss\":" << ael_shear_modulus_reuss;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_shear_modulus_reuss\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_bulk_modulus_vrh!=AUROSTD_NAN) {
        sscontent_json << "\"ael_bulk_modulus_vrh\":" << ael_bulk_modulus_vrh; //CT20190117
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_bulk_modulus_vrh\":null"; //CT20190117
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_shear_modulus_vrh!=AUROSTD_NAN) {
        sscontent_json << "\"ael_shear_modulus_vrh\":" << ael_shear_modulus_vrh;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_shear_modulus_vrh\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(ael_elastic_anisotropy!=AUROSTD_NAN) { //CO20181129
        sscontent_json << "\"ael_elastic_anisotropy\":" << ael_elastic_anisotropy; //CO20181129
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_elastic_anisotropy\":null"; //CO20181129
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_youngs_modulus_vrh!=AUROSTD_NAN) {
        sscontent_json << "\"ael_youngs_modulus_vrh\":" << ael_youngs_modulus_vrh;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_youngs_modulus_vrh\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_speed_sound_transverse!=AUROSTD_NAN) {
        sscontent_json << "\"ael_speed_sound_transverse\":" << ael_speed_sound_transverse;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_speed_sound_transverse\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_speed_sound_longitudinal!=AUROSTD_NAN) {
        sscontent_json << "\"ael_speed_sound_longitudinal\":" << ael_speed_sound_longitudinal;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_speed_sound_longitudinal\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      ////////////////////////////////////////////////////////////////////////// 

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_speed_sound_average!=AUROSTD_NAN) {
        sscontent_json << "\"ael_speed_sound_average\":" << ael_speed_sound_average;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_speed_sound_average\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_pughs_modulus_ratio!=AUROSTD_NAN) {
        sscontent_json << "\"ael_pughs_modulus_ratio\":" << ael_pughs_modulus_ratio;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_pughs_modulus_ratio\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_debye_temperature!=AUROSTD_NAN) {
        sscontent_json << "\"ael_debye_temperature\":" << ael_debye_temperature;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_debye_temperature\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_applied_pressure!=AUROSTD_NAN) {
        sscontent_json << "\"ael_applied_pressure\":" << ael_applied_pressure;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_applied_pressure\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////CT20181212
      if(ael_average_external_pressure!=AUROSTD_NAN) {
        sscontent_json << "\"ael_average_external_pressure\":" << ael_average_external_pressure;
      } else {
        if(PRINT_NULL) sscontent_json << "\"ael_average_external_pressure\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////ME20191105
      if ((ael_stiffness_tensor.rows == 6) && (ael_stiffness_tensor.cols == 6)) {
        sscontent_json << "\"ael_stiffness_tensor\":[";
        for (int i = 1; i <= 6; i++) {
          sscontent_json << "[";
          for (int j = 1; j <= 6; j++) sscontent_json << ael_stiffness_tensor[i][j] << ((j < 6)?",":"]");
          sscontent_json << ((i < 6)?",":"");
        }
        sscontent_json << "]";
      } else {
        if (PRINT_NULL) sscontent_json << "\"ael_stiffness_tensor\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //////////////////////////////////////////////////////////////////////////ME20191105
      if ((ael_compliance_tensor.rows == 6) && (ael_compliance_tensor.cols == 6)) {
        sscontent_json << "\"ael_compliance_tensor\":[";
        for (int i = 1; i <= 6; i++) {
          sscontent_json << "[";
          for (int j = 1; j <= 6; j++) sscontent_json << ael_compliance_tensor[i][j] << ((j < 6)?",":"]");
          sscontent_json << ((i < 6)?",":"");
        }
        sscontent_json << "]";
      } else {
        if (PRINT_NULL) sscontent_json << "\"ael_compliance_tensor\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);

      //ME20210927 BEGIN
      //APL
      //////////////////////////////////////////////////////////////////////////
      if (energy_free_vibrational_cell_apl_300K!=AUROSTD_NAN) {
        sscontent_json << "\"energy_free_vibrational_cell_apl_300K\":" << energy_free_vibrational_cell_apl_300K;
      } else {
        if (PRINT_NULL) sscontent_json << "\"energy_free_vibrational_cell_apl_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      if (energy_free_vibrational_atom_apl_300K!=AUROSTD_NAN) {
        sscontent_json << "\"energy_free_vibrational_atom_apl_300K\":" << energy_free_vibrational_atom_apl_300K;
      } else {
        if (PRINT_NULL) sscontent_json << "\"energy_free_vibrational_atom_apl_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      if (entropy_vibrational_cell_apl_300K!=AUROSTD_NAN) {
        sscontent_json << "\"entropy_vibrational_cell_apl_300K\":" << entropy_vibrational_cell_apl_300K;
      } else {
        if (PRINT_NULL) sscontent_json << "\"entropy_vibrational_cell_apl_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      if (entropy_vibrational_atom_apl_300K!=AUROSTD_NAN) {
        sscontent_json << "\"entropy_vibrational_atom_apl_300K\":" << entropy_vibrational_atom_apl_300K;
      } else {
        if (PRINT_NULL) sscontent_json << "\"entropy_vibrational_atom_apl_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      if (energy_internal_vibrational_cell_apl_300K!=AUROSTD_NAN) {
        sscontent_json << "\"energy_internal_vibrational_cell_apl_300K\":" << energy_internal_vibrational_cell_apl_300K;
      } else {
        if (PRINT_NULL) sscontent_json << "\"energy_internal_vibrational_cell_apl_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      if (energy_internal_vibrational_atom_apl_300K!=AUROSTD_NAN) {
        sscontent_json << "\"energy_internal_vibrational_atom_apl_300K\":" << energy_internal_vibrational_atom_apl_300K;
      } else {
        if (PRINT_NULL) sscontent_json << "\"energy_internal_vibrational_atom_apl_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      if (energy_zero_point_cell_apl!=AUROSTD_NAN) {
        sscontent_json << "\"energy_zero_point_cell_apl\":" << energy_zero_point_cell_apl;
      } else {
        if (PRINT_NULL) sscontent_json << "\"energy_zero_point_cell_apl\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      if (energy_zero_point_atom_apl!=AUROSTD_NAN) {
        sscontent_json << "\"energy_zero_point_atom_apl\":" << energy_zero_point_atom_apl;
      } else {
        if (PRINT_NULL) sscontent_json << "\"energy_zero_point_atom_apl\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      if (heat_capacity_Cv_cell_apl_300K!=AUROSTD_NAN) {
        sscontent_json << "\"heat_capacity_Cv_cell_apl_300K\":" << heat_capacity_Cv_cell_apl_300K;
      } else {
        if (PRINT_NULL) sscontent_json << "\"heat_capacity_Cv_cell_apl_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      if (heat_capacity_Cv_atom_apl_300K!=AUROSTD_NAN) {
        sscontent_json << "\"heat_capacity_Cv_atom_apl_300K\":" << heat_capacity_Cv_atom_apl_300K;
      } else {
        if (PRINT_NULL) sscontent_json << "\"heat_capacity_Cv_atom_apl_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////
      //ME20210927 END
      //AS20200901 BEGIN
      // QHA
      //////////////////////////////////////////////////////////////////////////
      if(gruneisen_qha!=AUROSTD_NAN) {
        sscontent_json << "\"gruneisen_qha\":" << gruneisen_qha;
      } else {
        if(PRINT_NULL) sscontent_json << "\"gruneisen_qha\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(gruneisen_qha_300K!=AUROSTD_NAN) {
        sscontent_json << "\"gruneisen_qha_300K\":" << gruneisen_qha_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"gruneisen_qha_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(thermal_expansion_qha_300K!=AUROSTD_NAN) {
        sscontent_json << "\"thermal_expansion_qha_300K\":" << thermal_expansion_qha_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"thermal_expansion_qha_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(modulus_bulk_qha_300K!=AUROSTD_NAN) {
        sscontent_json << "\"modulus_bulk_qha_300K\":" << modulus_bulk_qha_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"modulus_bulk_qha_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(modulus_bulk_derivative_pressure_qha_300K!=AUROSTD_NAN) {
        sscontent_json << "\"modulus_bulk_derivative_pressure_qha_300K\":" << modulus_bulk_derivative_pressure_qha_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"modulus_bulk_derivative_pressure_qha_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(heat_capacity_Cv_atom_qha_300K!=AUROSTD_NAN) {
        sscontent_json << "\"heat_capacity_Cv_atom_qha_300K\":" << heat_capacity_Cv_atom_qha_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"heat_capacity_Cv_atom_qha_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(heat_capacity_Cv_cell_qha_300K!=AUROSTD_NAN) {
        sscontent_json << "\"heat_capacity_Cv_cell_qha_300K\":" << heat_capacity_Cv_cell_qha_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"heat_capacity_Cv_cell_qha_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(heat_capacity_Cp_atom_qha_300K!=AUROSTD_NAN) {
        sscontent_json << "\"heat_capacity_Cp_atom_qha_300K\":" << heat_capacity_Cp_atom_qha_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"heat_capacity_Cp_atom_qha_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(heat_capacity_Cp_cell_qha_300K!=AUROSTD_NAN) {
        sscontent_json << "\"heat_capacity_Cp_cell_qha_300K\":" << heat_capacity_Cp_cell_qha_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"heat_capacity_Cp_cell_qha_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(volume_atom_qha_300K!=AUROSTD_NAN) {
        sscontent_json << "\"volume_atom_qha_300K\":" << volume_atom_qha_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"volume_atom_qha_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(energy_free_atom_qha_300K!=AUROSTD_NAN) {
        sscontent_json << "\"energy_free_atom_qha_300K\":" << energy_free_atom_qha_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"energy_free_atom_qha_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(energy_free_cell_qha_300K!=AUROSTD_NAN) {
        sscontent_json << "\"energy_free_cell_qha_300K\":" << energy_free_cell_qha_300K;
      } else {
        if(PRINT_NULL) sscontent_json << "\"energy_free_cell_qha_300K\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////
      //AS20200901 END

      // BADER
      //////////////////////////////////////////////////////////////////////////
      if(vbader_net_charges.size()) {
        sscontent_json << "\"bader_net_charges\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vbader_net_charges,6),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"bader_net_charges\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vbader_atomic_volumes.size()) {
        sscontent_json << "\"bader_atomic_volumes\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vbader_atomic_volumes,4),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"bader_atomic_volumes\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      // FILES
      //////////////////////////////////////////////////////////////////////////
      if(vfiles.size()) {
        aurostd::sort(vfiles);
        sscontent_json << "\"files\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vfiles,"\""),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"files\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      // CPUS
      //////////////////////////////////////////////////////////////////////////
      if(node_CPU_Model.size()) {
        sscontent_json << "\"node_CPU_Model\":\"" << node_CPU_Model << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"node_CPU_Model\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(node_CPU_Cores!=AUROSTD_NAN) {
        sscontent_json << "\"node_CPU_Cores\":" << node_CPU_Cores;
      } else {
        if(PRINT_NULL) sscontent_json << "\"node_CPU_Cores\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(node_CPU_MHz!=AUROSTD_NAN) {
        sscontent_json << "\"node_CPU_MHz\":" << node_CPU_MHz;
      } else {
        if(PRINT_NULL) sscontent_json << "\"node_CPU_MHz\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(node_RAM_GB!=INF) {
        sscontent_json << "\"node_RAM_GB\":" << node_RAM_GB;
      } else {
        if(PRINT_NULL) sscontent_json << "\"node_RAM_GB\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      // VERSION/DATE
      //////////////////////////////////////////////////////////////////////////
      if(aflow_version.size()) {
        sscontent_json << "\"aflow_version\":\"" << aflow_version << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"aflow_version\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(catalog.size()) {
        sscontent_json << "\"catalog\":\"" << catalog << "\"";
      } else {
        if(PRINT_NULL) sscontent_json << "\"catalog\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      sscontent_json << "\"aflowlib_version\":\"" << string(AFLOW_VERSION) << "\"";  //CO20170613
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      if(vaflowlib_date.size()) {
        sscontent_json << "\"aflowlib_date\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vaflowlib_date,"\""),",") << "]";
      } else {
        if(PRINT_NULL) sscontent_json << "\"aflowlib_date\":null";
      }
      vcontent_json.push_back(sscontent_json.str()); aurostd::StringstreamClean(sscontent_json);
      //////////////////////////////////////////////////////////////////////////

      sss << "{" << aurostd::joinWDelimiter(vcontent_json,eendl)  << "}";
      vcontent_json.clear();

      sss << endl;
    } // json

    return sss.str();
  }

  //  bool aflowlib2file(
  string _aflowlib_entry::aflowlib2file(string file,string mode) {
    string aflowlib_out=aflowlib2string(mode);
    aurostd::string2file(aflowlib_out,file);
    return aflowlib_out;
  }
}

//CO20200624 - CCE corrections
namespace aflowlib {
  double _aflowlib_entry::enthalpyFormationCell(int T) const { //CO20200624
    if(!XHOST.vflag_control.flag("NEGLECT_CCE")){ //CO20210115
      if(T==300 && enthalpy_formation_cce_300K_cell!=AUROSTD_NAN) return enthalpy_formation_cce_300K_cell;
      if(T==0 && enthalpy_formation_cce_0K_cell!=AUROSTD_NAN) return enthalpy_formation_cce_0K_cell;
    }
    return enthalpy_formation_cell;
  }
  double _aflowlib_entry::enthalpyFormationAtom(int T) const { //CO20200624
    if(!XHOST.vflag_control.flag("NEGLECT_CCE")){ //CO20210115
      if(T==300 && enthalpy_formation_cce_300K_atom!=AUROSTD_NAN) return enthalpy_formation_cce_300K_atom;
      if(T==0 && enthalpy_formation_cce_0K_atom!=AUROSTD_NAN) return enthalpy_formation_cce_0K_atom;
    }
    return enthalpy_formation_atom;
  }
}

//CO20171202 - apennsy fixes!
namespace aflowlib {
  void _aflowlib_entry::correctBadDatabase(bool verbose,ostream& oss){
    ofstream FileMESSAGE;
    return correctBadDatabase(FileMESSAGE,verbose,oss);
  }

  void _aflowlib_entry::correctBadDatabase(ofstream& FileMESSAGE,bool verbose,ostream& oss){
    //CO20180828 - LIB2 also contains unaries //so far we only know of bad binaries
    //APENNSY neglect - LIB2 only //CO20180828 - LIB2 also contains unaries  //binaries only
    stringstream message;
    if(vspecies_pp.size()==1 || vspecies_pp.size()==2) {
      string pseudoA="",pseudoB="";
      pseudoA = vspecies_pp[0];
      if(vspecies_pp.size()==2){pseudoB = vspecies_pp[1];}
      //[OBSOLETE CO20180828]string pseudoA = vspecies_pp[0];
      //[OBSOLETE CO20180828]string pseudoB = vspecies_pp[1];
      // tiny corrections
      //gamma_IrV
      if(pseudoA=="Cd" && pseudoB=="Pt" && prototype=="181") {
        enthalpy_formation_atom -= 0.0013;
        enthalpy_formation_cell = natoms * enthalpy_formation_atom;
        if(verbose){
          message << "Fixing enthalpy_formation of " << pseudoA << pseudoB << ":" << prototype;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        }
      }
      //gamma_IrV
      if(pseudoA=="Ir" && pseudoB=="V_sv" && prototype=="291") {
        enthalpy_formation_cell -= 0.001;
        enthalpy_formation_atom -= 0.0005;
        enthalpy_cell -= 0.001;
        enthalpy_atom -= 0.005;
        if(verbose){
          message << "Fixing enthalpy/enthalpy_formation of " << pseudoA << pseudoB << ":" << prototype;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        }
      }
      // HfPd
      if(pseudoA=="Hf_pv" && pseudoB=="Pd_pv" && prototype=="192") {
        enthalpy_formation_atom -= 0.003;
        enthalpy_formation_cell = natoms * enthalpy_formation_atom;
        enthalpy_atom = enthalpy_formation_atom;
        enthalpy_cell = natoms * enthalpy_atom;
        if(verbose){
          message << "Fixing enthalpy/enthalpy_formation of " << pseudoA << pseudoB << ":" << prototype;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        }
      }
      // sigma
      if(pseudoA=="Ir" && pseudoB=="Nb_sv" && prototype=="600.ABBAB") {
        enthalpy_formation_cell += 0.001;
        enthalpy_formation_atom += 0.0005;
        enthalpy_cell += 0.001;
        enthalpy_atom += 0.005;
        enthalpy_formation_cell += 0.001;
        enthalpy_formation_atom += 0.0005;
        enthalpy_cell += 0.001;
        enthalpy_atom += 0.005;
        if(verbose){
          message << "Fixing enthalpy/enthalpy_formation of " << pseudoA << pseudoB << ":" << prototype;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        }
      }
      // sigma
      if(pseudoA=="Os_pv" && pseudoB=="Re_pv" && prototype=="122") {
        enthalpy_formation_atom -= 0.001;
        enthalpy_formation_cell = natoms * enthalpy_formation_atom;
        enthalpy_atom = enthalpy_formation_atom;
        enthalpy_cell = natoms * enthalpy_atom;
        if(verbose){
          message << "Fixing enthalpy/enthalpy_formation of " << pseudoA << pseudoB << ":" << prototype;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
        }
      }
    }
  }
  bool _aflowlib_entry::ignoreBadDatabase() const{
    string reason;
    return ignoreBadDatabase(reason);
  }
  bool _aflowlib_entry::ignoreBadDatabase(string& reason) const{
    reason="";
    //grab bad database protos
    vector<string> _DEVIL_PROTOTYPES_;
    aurostd::string2tokens(_DEVIL_PROTOTYPES_STRING_,_DEVIL_PROTOTYPES_,",");

    //so far we only know of bad binaries
    //we need something more robust than just exact string match, case: 549 and 549.bis vs. 549.tetra
    bool match=false;
    //DEVIL
    for(uint di=0;di<_DEVIL_PROTOTYPES_.size() && !match;di++){if(pflow::prototypeMatch(prototype, _DEVIL_PROTOTYPES_[di])){match=true;}}
    if(match){
      reason=compound+":"+prototype+" is ill-calculated in the database";
      return true;
    }
    //find .old's
    if(1){
      uint prototype_size=prototype.size();
      string search_string=".old";uint search_string_size=search_string.size();
      if(prototype_size>search_string_size && prototype.substr(prototype_size-search_string_size,search_string_size)==search_string){  //look only at the end of the prototype
        reason=compound+":"+prototype+" is ill-calculated in the database";
        return true;
      }
    }
    //APENNSY neglect - LIB2 only //CO20180828 - LIB2 also contains unaries  //binaries only
    if(vspecies_pp.size()==1 || vspecies_pp.size()==2) {
      string pseudoA="",pseudoB="";
      pseudoA = vspecies_pp[0];
      if(vspecies_pp.size()==2){pseudoB = vspecies_pp[1];}
      //[OBSOLETE CO20180828]string pseudoA = vspecies_pp[0];
      //[OBSOLETE CO20180828]string pseudoB = vspecies_pp[1];
      // bad Ag is a wrong relaxation
      if((pseudoA=="Ag" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Ag" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ag is a wrong relaxation
      if((pseudoA=="Ag" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Ag" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Au is a wrong relaxation
      if((pseudoA=="Au" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Au" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Al_h pseudopotential !
      if((pseudoA=="Al_h" && pflow::prototypeMatch(prototype,"307")) || (pseudoB=="Al_h" && pflow::prototypeMatch(prototype,"308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Al_h pseudopotential !
      if((pseudoA=="Al_h" && pflow::prototypeMatch(prototype,"A7.A")) || (pseudoB=="Al_h" && pflow::prototypeMatch(prototype,"A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Al_h pseudopotential !
      if((pseudoA=="Al_h" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Al_h" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ca_sv is a wrong relaxation
      if((pseudoA=="Ca_sv" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Ca_sv" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ca_sv is a wrong relaxation
      if((pseudoA=="Ca_sv" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Ca_sv" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Cd is a wrong relaxation
      if((pseudoA=="Cd" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Cd" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Cu_pv is a wrong relaxation
      if((pseudoA=="Cu_pv" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Cu_pv" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Cu_pv is a wrong relaxation
      if((pseudoA=="Cu_pv" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Cu_pv" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Fe_pv is a wrong relaxation
      if((pseudoA=="Fe_pv" && pflow::prototypeMatch(prototype,"307")) || (pseudoB=="Fe_pv" && pflow::prototypeMatch(prototype,"308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Fe_pv is a wrong relaxation
      if((pseudoA=="Fe_pv" && pflow::prototypeMatch(prototype,"A7.A")) || (pseudoB=="Fe_pv" && pflow::prototypeMatch(prototype,"A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ge_h is a wrong relaxation
      if((pseudoA=="Ge_h" && pflow::prototypeMatch(prototype,"305")) || (pseudoB=="Ge_h" && pflow::prototypeMatch(prototype,"306"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad In_d is a wrong relaxation
      if((pseudoA=="In_d" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="In_d" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ir is a wrong relaxation
      if((pseudoA=="Ir" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Ir" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad K_sv is a wrong relaxation
      if((pseudoA=="K_sv" && pflow::prototypeMatch(prototype,"307")) || (pseudoB=="K_sv" && pflow::prototypeMatch(prototype,"308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad K_sv is a wrong relaxation
      if((pseudoA=="K_sv" && pflow::prototypeMatch(prototype,"A7.A")) || (pseudoB=="K_sv" && pflow::prototypeMatch(prototype,"A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad La is a wrong relaxation
      if((pseudoA=="La" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="La" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad La is a wrong relaxation
      if((pseudoA=="La" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="La" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Li_sv is a wrong relaxation
      if((pseudoA=="Li_sv" && pflow::prototypeMatch(prototype,"307")) || (pseudoB=="Li_sv" && pflow::prototypeMatch(prototype,"308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Li_sv is a wrong relaxation
      if((pseudoA=="Li_sv" && pflow::prototypeMatch(prototype,"A7.A")) || (pseudoB=="Li_sv" && pflow::prototypeMatch(prototype,"A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Na_pv is a wrong relaxation
      if((pseudoA=="Na_pv" && pflow::prototypeMatch(prototype,"307")) || (pseudoB=="Na_pv" && pflow::prototypeMatch(prototype,"308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Na_pv is a wrong relaxation
      if((pseudoA=="Na_pv" && pflow::prototypeMatch(prototype,"A7.A")) || (pseudoB=="Na_pv" && pflow::prototypeMatch(prototype,"A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ni_pv is a wrong relaxation
      if((pseudoA=="Ni_pv" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Ni_pv" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ni_pv is a wrong relaxation
      if((pseudoA=="Ni_pv" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Ni_pv" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pb_d is a wrong relaxation
      if((pseudoA=="Pb_d" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Pb_d" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pb_d is a wrong relaxation
      if((pseudoA=="Pb_d" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Pb_d" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pd_pv is a wrong relaxation
      if((pseudoA=="Pd_pv" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Pd_pv" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pd_pv is a wrong relaxation
      if((pseudoA=="Pd_pv" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Pd_pv" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pt is a wrong relaxation
      if((pseudoA=="Pt" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Pt" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Pt is a wrong relaxation
      if((pseudoA=="Pt" && pflow::prototypeMatch(prototype,"317")) || (pseudoB=="Pt" && pflow::prototypeMatch(prototype,"318"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Rh_pv is a wrong relaxation
      if((pseudoA=="Rh_pv" && pflow::prototypeMatch(prototype,"303")) || (pseudoB=="Rh_pv" && pflow::prototypeMatch(prototype,"304"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Si_h is a wrong relaxation
      if((pseudoA=="Si_h" && pflow::prototypeMatch(prototype,"305")) || (pseudoB=="Si_h" && pflow::prototypeMatch(prototype,"306"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Si_h is a wrong relaxation
      if((pseudoA=="Si_h" && pflow::prototypeMatch(prototype,"307")) || (pseudoB=="Si_h" && pflow::prototypeMatch(prototype,"308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Si_h is a wrong relaxation
      if((pseudoA=="Si_h" && pflow::prototypeMatch(prototype,"A7.A")) || (pseudoB=="Si_h" && pflow::prototypeMatch(prototype,"A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Si_h is a wrong relaxation
      if((pseudoA=="Si_h" && pflow::prototypeMatch(prototype,"323")) || (pseudoB=="Si_h" && pflow::prototypeMatch(prototype,"324"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ta_pv is a wrong relaxation
      if((pseudoA=="Ta_pv" && pflow::prototypeMatch(prototype,"307")) || (pseudoB=="Ta_pv" && pflow::prototypeMatch(prototype,"308"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad Ta_pv is a wrong relaxation
      if((pseudoA=="Ta_pv" && pflow::prototypeMatch(prototype,"A7.A")) || (pseudoB=="Ta_pv" && pflow::prototypeMatch(prototype,"A7.B"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // bad B_h is a wrong relaxation
      if((pseudoA=="B_h" && pflow::prototypeMatch(prototype,"317")) || (pseudoB=="B_h" && pflow::prototypeMatch(prototype,"318"))) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }

      // sigma
      if(pseudoA=="Os_pv" && pseudoB=="Re_pv" && pflow::prototypeMatch(prototype,"448")) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
      // wrong channel, bug
      if(pseudoA=="Rh_pv" && pseudoB=="Zr_sv" && pflow::prototypeMatch(prototype,"381")) {
        reason=pseudoA+pseudoB+":"+prototype+" is ill-calculated in the database";
        return true;
      }
    }
    return false;
  }
} // namespace aflowlib

namespace aflowlib {
  string _aflowlib_entry::getPathAURL(ostream& oss, bool load_from_common) const { //CO20200404
    ofstream FileMESSAGE;
    return getPathAURL(FileMESSAGE, oss, load_from_common);
  }
  string _aflowlib_entry::getPathAURL(ofstream& FileMESSAGE,ostream& oss, bool load_from_common) const {  //CO20200404
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    string path = "";
    if (aurl.empty()) {return path;}
    vector<string> tokens;
    aurostd::string2tokens(aurl, tokens, ":");
    //LIB1 presents problems here (3 colons): aflowlib.duke.edu:AFLOWDATA/LIB1_RAW/Pt:PAW_PBE:05Jan2001/A6
    if(0){
      if (tokens.size() != 2) {
        message << "Odd AURL format for entry " << auid << ": " << aurl;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_WARNING_);
        return path;
      }
    }

    //instead just erase first item, join others, assume we're okay...
    tokens.erase(tokens.begin());
    path=aurostd::joinWDelimiter(tokens,":");

    string server="",path_full="";
    if(1&&load_from_common){
      //attempt 1: try replacing _RAW with _LIB
      if(1){
        server="/www";
        path_full=path;
        aurostd::StringSubst(path_full,"_RAW","_LIB");
        path_full=server+"/"+path_full;
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " attempt 1 path=" << path_full << endl;}
        if(aurostd::IsDirectory(path_full)){return path_full;}
      }

      //attempt 2: try finding LIB directory
      if(1){
        server="/common";
        path_full=path;
        aurostd::StringSubst(path_full,"AFLOWDATA/","");
        aurostd::StringSubst(path_full,"ICSD_WEB","ICSD/LIB"); //CO20200223
        aurostd::StringSubst(path_full,"_RAW","/LIB");
        path_full=server+"/"+path_full;
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " attempt 2 path=" << path_full << endl;}
        if(aurostd::IsDirectory(path_full)){return path_full;}
      }

      //attempt 3: try no replacement (RAW)
      if(1){
        server="/www";
        path_full=server+"/"+path;
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " attempt 3 path=" << path_full << endl;}
        if(aurostd::IsDirectory(path_full)){return path_full;}
      }
    }
    if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")){server=XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER");}
    else{server="aflowlib.duke.edu";}

    path_full=server+"/"+path;
    return path_full;
  }
  /// @brief get the full directory path of the aflowlib entry
  /// @param filesystem_collection filesystem collection (LIB, RAW, WEB)
  /// @return directory path
  /// @authors
  /// @mod{SD,20221207,created}
  string _aflowlib_entry::getPathDirectory(const string& filesystem_collection) {
    string path_full = DEFAULT_ENTRY_LOADER_FS_PATH + (*this).aurl.substr(28) + "/";
    path_full = std::regex_replace(path_full, aurostd::regex_aurl2file, "$1/" + filesystem_collection + "/");
    return path_full;
  }
  vector<string> _aflowlib_entry::getSpecies() const{ //CO20221110
    if(!vspecies.empty()){return vspecies;}
    vector<string> vspecies_tmp;
    uint i=0;
    if(!vspecies_pp.empty()){
      vspecies_tmp.clear();
      for(i=0;i<vspecies_pp.size();i++){
        vspecies_tmp.push_back(vspecies_pp[i]);
        aurostd::VASP_PseudoPotential_CleanName_InPlace(vspecies_tmp.back());
      }
      return vspecies_tmp;
    }
    //use REST-API
    if(aurl.empty()){return vspecies_tmp;}
    string url(aurl);
    aurostd::StringSubst(url,":AFLOWDATA","/AFLOWDATA");
    string species_restapi="";
    if(species_restapi.empty()){aurostd::url2string(url+"/?species",species_restapi);}
    if(species_restapi.empty()){aurostd::url2string(url+"/?species_pp",species_restapi);}
    if(species_restapi.empty()){return vspecies_tmp;}
    aurostd::string2tokens(species_restapi,vspecies_tmp,",");
    for(i=0;i<vspecies_tmp.size();i++){aurostd::VASP_PseudoPotential_CleanName_InPlace(vspecies_tmp.back());}
    return vspecies_tmp;
  }
  //NOTE: getSpeciesAURL() does not work for LIB2 unaries
  vector<string> _aflowlib_entry::getSpeciesAURL(ostream& oss) const{ofstream FileMESSAGE;return getSpeciesAURL(FileMESSAGE, oss);} //CO20200404
  vector<string> _aflowlib_entry::getSpeciesAURL(ofstream& FileMESSAGE,ostream& oss) const {  //CO20200404
    bool LDEBUG=(false || XHOST.DEBUG);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " aurl=" << aurl << endl;}
    vector<string> vspecies_aurl;
    if(aurl.empty()){return vspecies_aurl;}
    string path=aurl.substr(aurl.find(":")+1,string::npos);
    vector<string> tokens;
    //[CO20221111 - slow]aurostd::string2tokens(aurl,tokens,":");
    //[CO20221111 - slow]//erase first item (aflowlib.duke.edu), join others, assume we're okay...
    //[CO20221111 - slow]tokens.erase(tokens.begin());
    //[CO20221111 - slow]string path=aurostd::joinWDelimiter(tokens,":");
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " path=" << path << endl;}

    bool isICSD=false;
    string::size_type nth_slash=2;  //AFLOWDATA/LIB2_RAW/TeZr_sv/10
    elements_string_type e_str_type=pp_string;
    if(path.find("_ICSD_")!=string::npos){nth_slash=3;e_str_type=composition_string;isICSD=true;} //AFLOWDATA/ICSD_WEB/HEX/Te2Zr1_ICSD_653213
    string::size_type n=0,start=0;
    while(n++<nth_slash){start=path.find("/",start+1);}
    string::size_type end=path.find("/",start+1);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " start=" << start << " end=" << end << endl;}
    
    string species_string=path.substr(start+1,end-start-1);
    if(isICSD){start=0;end=species_string.find("_ICSD_");species_string=species_string.substr(start,end);}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " species_string=" << species_string << endl;}
    if(species_string.empty()){
      stringstream message;
      message << "Odd AURL format for entry " << auid << ": " << aurl;pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_WARNING_);
      return vspecies_aurl;
    }

    //[CO20221111 - slow]//split by /
    //[CO20221111 - slow]aurostd::string2tokens(path,tokens,"/");
    //[CO20221111 - slow]if(tokens.size()<4){
    //[CO20221111 - slow]  message << "Odd AURL format for entry " << auid << ": " << aurl;
    //[CO20221111 - slow]  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_WARNING_);
    //[CO20221111 - slow]  return vspecies_aurl;
    //[CO20221111 - slow]}
    //[CO20221111 - slow]string species_string="";
    //[CO20221111 - slow]if(path.find("_ICSD_")!=string::npos){  //if ICSD: AFLOWDATA/ICSD_WEB/HEX/Te2Zr1_ICSD_653213
    //[CO20221111 - slow]  species_string=tokens[3];
    //[CO20221111 - slow]  string::size_type loc;loc=species_string.find("_ICSD_");
    //[CO20221111 - slow]  species_string=species_string.substr(0,loc);
    //[CO20221111 - slow]}
    //[CO20221111 - slow]else{ //otherwise: AFLOWDATA/LIB2_RAW/TeZr_sv/10
    //[CO20221111 - slow]  species_string=tokens[2];
    //[CO20221111 - slow]  //fix LIB1: Zr_sv:PAW_PBE:07Sep2000
    //[CO20221111 - slow]  string::size_type loc;loc=species_string.find(":");
    //[CO20221111 - slow]  species_string=species_string.substr(0,loc);
    //[CO20221111 - slow]}

    vector<double> vcomposition;  //save some overloads
    bool clean=true;          //clean LIB1
    bool sort_elements=false; //not necessary, already sorted
    bool keep_pp=false;       //remove pp
    vspecies_aurl=aurostd::getElements(species_string,vcomposition,e_str_type,FileMESSAGE,clean,sort_elements,keep_pp,oss);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vspecies_aurl=" << aurostd::joinWDelimiter(vspecies_aurl,",") << endl;}
    return vspecies_aurl;
  }
}

// **************************************************************************
// directory2auid
// auid2directory
// auid2present
// **************************************************************************
namespace aflowlib {
  string _aflowlib_entry::POCCdirectory2MetadataAUIDjsonfile(const string& directory,uint salt){  //CO20200624
    //CO20200624 - THIS IS HOW WE CREATE AUID FOR POCC STRUCTURES
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;

    if(aurl.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"AURL has not been calculated",_INPUT_MISSING_);}

    string system_name=KBIN::ExtractSystemName(directory);
    if(system_name.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No system name found",_FILE_CORRUPT_);}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " system_name=" << system_name << endl;}
    system_name=pocc::addDefaultPOCCTOL2string(system_name);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " system_name(with TOL)=" << system_name << endl;}

    _aflags aflags;aflags.Directory=directory;
    pocc::POccCalculator pcalc(aflags);
    pcalc.loadDataIntoCalculator();
    if(pcalc.m_ARUN_directories.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No ARUN.POCC_* runs found",_FILE_CORRUPT_);}
    vector<string> vauid_aruns;
    string aurl_arun="";
    for(uint i=0;i<pcalc.m_ARUN_directories.size();i++){
      aurl_arun=aurl+"/"+pcalc.m_ARUN_directories[i];
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " m_ARUN_directories[" << i << "]=" << pcalc.m_ARUN_directories[i] << endl;
        cerr << __AFLOW_FUNC__ << " m_ARUN_directories[" << i << "].aurl=" << aurl_arun << endl;
      }
      vauid_aruns.push_back(VASPdirectory2auid(aurostd::CleanFileName(directory+"/"+pcalc.m_ARUN_directories[i]),aurl_arun));
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_ARUN_directories[" << i << "].auid=" << vauid_aruns.back() << endl;}
    }

    stringstream sscontent_json;
    string eendl="";
    vector<string> vcontent_json;

    sscontent_json << "\"salt\":" << aurostd::utype2string(salt);vcontent_json.push_back(sscontent_json.str());aurostd::StringstreamClean(sscontent_json);
    sscontent_json << "\"aflow_type\":" << "\"aggregate\"";vcontent_json.push_back(sscontent_json.str());aurostd::StringstreamClean(sscontent_json);
    sscontent_json << "\"method\":" << "\"aflow_pocc\"";vcontent_json.push_back(sscontent_json.str());aurostd::StringstreamClean(sscontent_json);
    sscontent_json << "\"aggregate_parameters\":" << "\"" << system_name << "\"";vcontent_json.push_back(sscontent_json.str());aurostd::StringstreamClean(sscontent_json);
    sscontent_json << "\"aggregate_content\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vauid_aruns,"\""),",") << "]";vcontent_json.push_back(sscontent_json.str());aurostd::StringstreamClean(sscontent_json);

    string metadata_auid_json="";
    uint64_t crc=0;

    metadata_auid_json="{"+aurostd::joinWDelimiter(vcontent_json,",")+"}";
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " METADATA_AUID.JSON=" << endl << metadata_auid_json << endl;}
    //aurostd::string2file(metadata_auid_json,"metadata_auid.json");
    crc=aurostd::crc64(0,metadata_auid_json); // DONT TOUCH THIS
    auid="aflow:"+aurostd::crc2string(crc);

    //find conflicts
    string aurl_found="";
    if(aflowlib::auid2present(auid,aurl_found,1)){
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " conflict auid=" << auid << endl;	
      message << "CONFLICT POTENTIAL " << auid << " " << aurl_found << " " << aurl;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,_LOGGER_MESSAGE_);
      if(aurl_found!=aurl) { // avoid conflict with yourself
        salt++;
        metadata_auid_json=POCCdirectory2MetadataAUIDjsonfile(directory,salt);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " METADATA_AUID.JSON=" << endl << metadata_auid_json << endl;}
        //aurostd::string2file(metadata_auid_json,"metadata_auid.json");
        crc=aurostd::crc64(0,metadata_auid_json); // DONT TOUCH THIS
        auid="aflow:"+aurostd::crc2string(crc);
      } else {
        message << "CONFLICT TRIVIAL   " << auid << " " << aurl_found << " " << aurl;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,_LOGGER_MESSAGE_);
      }
    }

    return metadata_auid_json;
  }
  string VASPdirectory2auid(const string& directory,const string& aurl){
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;

    string auid="";
    bool conflict=TRUE; 
    while (conflict) {
      uint64_t crc=0;
      // DONT TOUCH THE AUID FLOW

      vector<string> vfiles2;
      for(uint iext=1;iext<XHOST.vext.size();iext++) {  // SKIP uncompressed
        vfiles2.push_back("OUTCAR.relax1"+XHOST.vext.at(iext));
        vfiles2.push_back("OUTCAR.relax2"+XHOST.vext.at(iext));
        vfiles2.push_back("OUTCAR.relax3"+XHOST.vext.at(iext));
        vfiles2.push_back("OUTCAR.relax4"+XHOST.vext.at(iext));
        vfiles2.push_back("OUTCAR.static"+XHOST.vext.at(iext));
        vfiles2.push_back("OUTCAR.bands"+XHOST.vext.at(iext));
        vfiles2.push_back("OUTCAR"+XHOST.vext.at(iext));
      }

      if(LDEBUG) cerr << __AFLOW_FUNC__ << " [0]" << endl;

      for(uint i=0;i<vfiles2.size();i++) 
        if(aurostd::FileExist(directory+"/"+vfiles2.at(i)))
          crc=aurostd::crc64(crc,aurostd::efile2string(directory+"/"+vfiles2.at(i))); // DONT TOUCH THIS
      auid="aflow:"+aurostd::crc2string(crc);

      if(LDEBUG) cerr << __AFLOW_FUNC__ << " [1]" << endl;

      if(LDEBUG) cerr << __AFLOW_FUNC__ << " auid=" << auid << endl;
      conflict=FALSE;
      string aurl_found="";
      if(aflowlib::auid2present(auid,aurl_found,1)) {
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " conflict auid=" << auid << endl;	
        message << "CONFLICT POTENTIAL " << auid << " " << aurl_found << " " << aurl;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,_LOGGER_MESSAGE_); //CO20200624
        if(aurl_found!=aurl) { // avoid conflict with yourself
          string salt="AUID_salt["+aurostd::utype2string<long double>(aurostd::get_useconds())+"]";
          message << "CONFLICT TRUE      " << auid << " " << aurl_found << " " << aurl << "  " << salt;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,_LOGGER_WARNING_); //CO20200624
          string file=vfiles2.at(0);
          //
          for(uint iext=0;iext<XHOST.vext.size();iext++) { aurostd::StringSubst(file,XHOST.vext.at(iext),""); }
          stringstream sss;aurostd::efile2stringstream(directory+"/"+file+DEFAULT_KZIP_EXT,sss); sss << endl << salt << endl;
          aurostd::execute("mv \""+directory+"/"+file+DEFAULT_KZIP_EXT+"\""+" \""+directory+"/"+file+".conflict_auid"+DEFAULT_KZIP_EXT+"\"");
          aurostd::stringstream2compressfile(DEFAULT_KZIP_BIN,sss,directory+"/"+file);
          //
          conflict=TRUE; // recheck
        } else {
          message << "CONFLICT TRIVIAL   " << auid << " " << aurl_found << " " << aurl;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,_LOGGER_MESSAGE_); //CO20200624
        }
      }
    }
    return auid;
  }
  bool _aflowlib_entry::directory2auid(const string& directory) { //CO20200624
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "_aflowlib_entry::directory2auid: BEGIN" << endl;
    auid=VASPdirectory2auid(directory,aurl); // now it has auid   
    vauid.clear();
    aflowlib::auid2vauid(auid,vauid);

    if(LDEBUG) cerr << "directory2auid: END" << endl;

    cout << XPID << "_aflowlib_entry::directory2auid: DIRECTORY=" << directory << endl; // DONT TOUCH THIS
    cout << XPID << "_aflowlib_entry::directory2auid: AURL_ID=" << aurostd::crc2string(aurostd::crc64(0,directory)) << endl; // DONT TOUCH THIS

    return TRUE;
  }
}

namespace aflowlib {
  bool json2aflowlib(const string& json,string key,string& value) { //SC20200415
    // return TRUE if something has been found
    value="";
    key="\""+key+"\":";
    string::size_type start,end;
    start=json.find(key);
    if(start!=string::npos) {
      start+=key.length();
      end=json.find("\":",start);
      if(end!=string::npos){
        value=json.substr(start,end-start);
        end=value.find_last_of(",");
        value=value.substr(0,end);
      } else {
        end=json.find("}",start);
        value=json.substr(start,end-start);
      }
      //    if((value[0]=='\"') && (value[value.size()-1]=='\"')) value=value.substr(1,value.size()-2);  // Remove quotes
    } else {
      value="";
    }
    // cleanup
    aurostd::StringSubst(value,"[","");  // Remove brakets
    aurostd::StringSubst(value,"]","");  // Remove brakets
    aurostd::StringSubst(value,"\"",""); // Remove quotes

    return !value.empty();
  }

  uint auid2present(string auid,string& aurl,int mode) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "aflowlib::auid2present: BEGIN mode=" << mode << endl;
    string loop="",json="";aurl="";
    if(auid=="" || auid.size()!=22) { cerr << XPID << "aflowlib::auid2present: auid.size() needs to be 22 characters long" << endl; return FALSE;}

    // [OBSOLETE]    if(mode==0) {
    // [OBSOLETE]    if(XHOST_vAUID.size()==0 || XHOST_vAURL.size()==0) init::InitGlobalObject("vLIBS","",FALSE);
    // [OBSOLETE]    if(LDEBUG) cerr << XPID << "aflowlib::auid2present: [4] XHOST_vAURL.size()=" << XHOST_vAURL.size() << "  XHOST_vAUID.size()=" << XHOST_vAUID.size() << endl;
    // [OBSOLETE]    bool found=FALSE;
    // [OBSOLETE]    for(uint j=0;j<XHOST_vAUID.size()&&!found;j++) {
    // [OBSOLETE]  	if(LDEBUG && XHOST_vAUID.at(j)==auid) cerr << "[" << auid << "] [" << XHOST_vAUID.at(j) << "] [" << XHOST_vAURL.at(j) << "]" << " [" << j << "]" << endl;
    // [OBSOLETE]  	if(XHOST_vAUID.at(j)==auid) {found=TRUE;aurl=XHOST_vAURL.at(j);}
    // [OBSOLETE]    }
    // [OBSOLETE]    if(LDEBUG) cerr << XPID << "aflowlib::auid2present: END  auid=" << auid << "  aurl=" << aurl << endl;
    // [OBSOLETE]    return found;
    // [OBSOLETE]   }
    if(mode==1) { // PICK THIS ONE DEFAULT
      //  bool aflowlib::auid2present(string auid,string& aurl)
      string jsonl_file=XHOST_LIBRARY_JSONL+"/"; for(uint i=0;i<8;i++) jsonl_file+=auid.at(i); jsonl_file+=".jsonl";
      bool found=FALSE;
      jsonl_file=aurostd::CleanFileName(jsonl_file);
      for(uint i=0;i<XHOST.vext.size()&&aurl.empty()&&!found;i++) {
        if(LDEBUG) cerr << XPID << "aflowlib::auid2present: TESTING=" << jsonl_file << XHOST.vext.at(i) << endl; 
        //	cout << XPID << "aflowlib::auid2present: TESTING=" << jsonl_file << XHOST.vext.at(i) << endl; 
        if(aurostd::FileExist(jsonl_file+XHOST.vext.at(i))) {
          found=TRUE;
          if(LDEBUG) cerr << XPID << "aflowlib::auid2present: FOUND=" << jsonl_file << XHOST.vext.at(i) << endl; 
          //  cout << XPID << "aflowlib::auid2present: FOUND=" << jsonl_file << XHOST.vext.at(i) << endl; 
          json=aurostd::execute2string(XHOST.vcat.at(i)+" "+jsonl_file+XHOST.vext.at(i)+" | grep "+auid);
          aflowlib::json2aflowlib(json,"aurl",aurl);
          aflowlib::json2aflowlib(json,"loop",loop);
        }
      }
      if(LDEBUG) cerr << XPID << "aflowlib::auid2present: END  auid=" << auid << "  aurl=" << aurl << "  loop=" << loop << "  json.size()=" << json.size() << endl;
      cout << XPID << "aflowlib::auid2present: auid=" << auid << "  aurl=" << aurl << "  loop=" << loop << "  json.size()=" << json.size() << endl;
      return json.size();
    }
    if(mode==2) { // not that faster and does not keep an outside vAUID table so it does not see the TRIVIAL CONFLICTS
      //  bool aflowlib::auid2present(string auid,string& aurl)
      string jsonl_file=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_AUID)+"/"+aflowlib::auid2directory(auid)+"/RAW/aflowlib.json";
      bool found=FALSE;
      jsonl_file=aurostd::CleanFileName(jsonl_file);
      for(uint i=0;i<XHOST.vext.size()&&!found;i++) {
        if(LDEBUG) cerr << XPID << "aflowlib::auid2present: TESTING=" << jsonl_file << XHOST.vext.at(i) << endl; 
        if(aurostd::FileExist(jsonl_file+XHOST.vext.at(i))) {
          found=TRUE;
          if(LDEBUG) cerr << XPID << "aflowlib::auid2present: FOUND=" << jsonl_file << XHOST.vext.at(i) << endl;
          json=aurostd::execute2string(XHOST.vcat.at(i)+" "+jsonl_file+XHOST.vext.at(i));
          aflowlib::json2aflowlib(json,"aurl",aurl);
          aflowlib::json2aflowlib(json,"loop",loop);
        }
      }
      if(LDEBUG) cerr << XPID << "aflowlib::auid2present: END  auid=" << auid << "  aurl=" << aurl << "  loop=" << loop << "  json.size()=" << json.size() << endl;
      cout << XPID << "aflowlib::auid2present: auid=" << auid << "  aurl=" << aurl << "  loop=" << loop << "  json.size()=" << json.size() << endl;
      return json.size();
    }

    return FALSE;
  }

  // [OBSOLETE]   time ./aflow --beep --force --showPID --lib2raw=/common/LIB3/LIB/Cu_pvHgSn/TFCC004.CAB
  // [OBSOLETE]   time aflow --beep --force --showPID --lib2raw=/common/LIB3/LIB/Cu_pvHgSn/TFCC004.CAB  
  // [OBSOLETE]   time ./aflow --beep --force --showPID --lib2raw=/common/LIB3/LIB/AgCdCo/TFCC001.ABC   
  // [OBSOLETE]   time aflow --beep --force --showPID --lib2raw=/common/LIB3/LIB/AgCdCo/TFCC001.ABC  
  // [OBSOLETE]   time ./aflow --beep --force --showPID --lib2raw="/common/LIB4/LIB/AgCdCoZr_sv:PAW_PBE/ABCD_cF16_216_c_d_b_a.ABCD"
  // [OBSOLETE]   time aflow --beep --force --showPID --lib2raw="/common/LIB4/LIB/AgCdCoZr_sv:PAW_PBE/ABCD_cF16_216_c_d_b_a.ABCD"


  uint auid2vauid(const string auid, deque<string>& vauid) {                // splits the auid into vauid
    vauid.clear();
    //    vauid.push_back(auid.substr(0,6)); for(uint i=6;i<=20;i+=2) vauid.push_back(auid.substr(i,2));  // splitting aflow:/ab/cd..
    vauid.push_back(auid.substr(0,8)); for(uint i=8;i<=20;i+=2) vauid.push_back(auid.substr(i,2));  // splitting aflow:ab/cd..
    return vauid.size();
  }

  string auid2directory(const string auid) {                                // gives AUID directory from existence of vauid
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "aflowlib::auid2directory: BEGIN" << endl;
    string directory="";  //CO20200624
    deque<string> vauid;
    aflowlib::auid2vauid(auid,vauid);
    if(vauid.size()>0) {
      directory=vauid.at(0);
      for(uint i=1;i<vauid.size();i++) {
        directory+="/"+vauid.at(i);
      }
    }
    if(LDEBUG) cerr << XPID << "aflowlib::auid2directory: END" << endl;
    return directory;
  }

}

// **************************************************************************
// Operate on CIFs for JMOL for one structure
// **************************************************************************
namespace aflowlib {
  bool cif2data(string file,double& a,double& b,double& c,double& alpha,double& beta,double& gamma) {
    vector<string> vline,tokens;
    aurostd::file2vectorstring(file,vline);
    for(uint i=0;i<vline.size();i++) {
      aurostd::string2tokens(vline.at(i),tokens," ");
      if(aurostd::substring2bool(vline.at(i),"_cell_length_a") && tokens.size()>1) a=aurostd::string2utype<double>(tokens.at(1));
      if(aurostd::substring2bool(vline.at(i),"_cell_length_b") && tokens.size()>1) b=aurostd::string2utype<double>(tokens.at(1));
      if(aurostd::substring2bool(vline.at(i),"_cell_length_c") && tokens.size()>1) c=aurostd::string2utype<double>(tokens.at(1));
      if(aurostd::substring2bool(vline.at(i),"_cell_angle_alpha") && tokens.size()>1) alpha=aurostd::string2utype<double>(tokens.at(1));
      if(aurostd::substring2bool(vline.at(i),"_cell_angle_beta") && tokens.size()>1) beta=aurostd::string2utype<double>(tokens.at(1));
      if(aurostd::substring2bool(vline.at(i),"_cell_angle_gamma") && tokens.size()>1) gamma=aurostd::string2utype<double>(tokens.at(1));
    }
    return TRUE;
  }

  bool cif2oss(string file,string label,string sgnumber,ostream& oss) {
    double aCIF,bCIF,cCIF,alphaCIF,betaCIF,gammaCIF;
    aflowlib::cif2data(file,aCIF,bCIF,cCIF,alphaCIF,betaCIF,gammaCIF);
    oss << " font echo 14;color echo white; ";
    oss << " set echo 3% 97%; echo \\\"AFLOW-JSmol consortium (AFLOW V" << string(AFLOW_VERSION) << ") | entry="<< label << "  |  | ";
    if(aurostd::string2utype<int>(sgnumber)>0) {
      oss <<"  Spacegroup = "<< GetSpaceGroupName(aurostd::string2utype<int>(sgnumber)) << " (#" <<  sgnumber << ")   |";
    }
    oss << " a=" << aCIF <<"\u212B, b=" << bCIF <<"\u212B, c=" << cCIF <<"\u212B    |" << " \u03B1=" << alphaCIF <<"\u00B0, \u03B2=" << betaCIF <<"\u00B0, \u03B3=" << gammaCIF <<"\u00B0   \\\"; ";

    return TRUE;
  }
}

//BEGIN JJPR
// **************************************************************************
// GET BADER jvxl file
// **************************************************************************

namespace aflowlib {
  bool iso2oss(string file,string label, string element,string cutoff,int index,ostream& oss) {
    vector<string> kk;
    kk.resize(9);
    kk[0]="red";
    kk[1]="green";
    kk[2]="yellow";
    kk[3]="blue";
    kk[4]="orange";
    kk[5]="white";
    kk[6]="purple";
    kk[7]="brown";
    kk[8]="pink";
    oss <<"  ISOSURFACE  " << element << " \\\""<< file+"/"+label+"_Bader_"+cutoff+"_"+element+".jvxl\\\"  ;isosurface mesh;  color isosurface " << kk[index] << " translucent;";
    oss << "x = load(\\\""<< file+"/"+label+"_abader.out" << "\\\"); charges = x.split(\\\"=\\\")[2].split(\\\"(\\\")[1].split(\\\",\\\"); {*}.label = charges; label \%[label];";
    // FOR LOCAL TEST: oss <<"  ISOSURFACE  " << element << " "<< "Bader_"+cutoff+"_"+element+".jvxl  ;isosurface mesh;  color isosurface " << kk[index] << " translucent;";
    return TRUE;
  }
}
//END JJPR

// **************************************************************************
// GET SPACE GROUP for one structure
// **************************************************************************
namespace aflowlib {
  uint SGtoNSG(string sgroup) {
    string::size_type idx1;
    string strsub("#");
    idx1=sgroup.find(strsub);
    if(idx1!=string::npos)  
      return (int) atoi(sgroup.substr(sgroup.find(strsub)+strsub.length()).c_str());
    else return 0;
  }

  void _aflowlib_entry::GetSGROUP(string aflowlibentry) {
    string message = "";
    vector<string> vaflowlib_entry;
    aurostd::string2tokens(aflowlibentry,vaflowlib_entry,"|");

    bool VERBOSE_LOCAL=(FALSE || XHOST.DEBUG);
    if(vaflowlib_entry.size()==0) {
      string message = DEFAULT_FILE_AFLOWLIB_ENTRY_OUT + " file not found";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    vsgroup.clear();vNsgroup.clear();

    string data_aurl="";vector<string> tokens;
    // CHECK FOR HOLES
    if(XGNDSTATE_HOLES==0)  // SAFE NO HOLES IN THE XMATRIX
      if(vaflowlib_entry.size()==0) {
        message = DEFAULT_FILE_AFLOWLIB_ENTRY_OUT + " file not found";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
      }
    if(XGNDSTATE_HOLES==1)  // ALLOW HOLES WITH FAKE VALUES
      if(vaflowlib_entry.size()==0) {
        //cerr << "FOUND SG HOLE = " << alloy_dir << "/" << params.structures[2].at(str_number).name << "   === " << str_number << " " << structures_number[str_number]<< endl;
        vsgroup.clear();vsgroup.push_back(NOSG);vsgroup.push_back(NOSG);vsgroup.push_back(NOSG);
        vNsgroup.clear();vNsgroup.push_back(0);vNsgroup.push_back(0);vNsgroup.push_back(0);
        return;
      }

    for(uint i=0;i<vaflowlib_entry.size();i++) {
      if(aurostd::substring2bool(vaflowlib_entry.at(i),"aurl=")) data_aurl=aurostd::substring2string(vaflowlib_entry.at(i),"arl=");
      if(aurostd::substring2bool(vaflowlib_entry.at(i),"sg=")) {
        aurostd::string2tokens(vaflowlib_entry.at(i),tokens,",");
        if(tokens.size()==0) {
          message = "geometry not enough tokens";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
        }
        if(tokens.size()==3) { // ok
          vsgroup.clear();
          aurostd::StringSubst(tokens.at(0),"sg=","");aurostd::StringSubst(tokens.at(0)," ","");aurostd::StringSubst(tokens.at(0),"#"," #");vsgroup.push_back(tokens.at(0));
          aurostd::StringSubst(tokens.at(1)," ","");aurostd::StringSubst(tokens.at(1),"#"," #");vsgroup.push_back(tokens.at(1));
          aurostd::StringSubst(tokens.at(2)," ","");aurostd::StringSubst(tokens.at(2),"#"," #");vsgroup.push_back(tokens.at(2));
          if(VERBOSE_LOCAL) cout << "[5] "<<"["<<tokens.at(0) << "," << tokens.at(1) << "," << tokens.at(2) <<"]" << endl;
        } else {
          vsgroup.clear();vsgroup.push_back(NOSG);vsgroup.push_back(NOSG);vsgroup.push_back(NOSG);
          vNsgroup.clear();vNsgroup.push_back(0);vNsgroup.push_back(0);vNsgroup.push_back(0);
        }
      }
    }
    // done, now makes the numbers
    for(uint ii=0;ii<vsgroup.size();ii++) vNsgroup.push_back(SGtoNSG(vsgroup.at(ii)));
    // if(NsgroupPRE==0) {sgroupPRE=sgroupMID;NsgroupPRE=SGtoNSG(sgroupPRE);}	
    // if(NsgroupPRE==0) {sgroupPRE=sgroupPOST;NsgroupPRE=SGtoNSG(sgroupPRE);}	
    for(uint i=vsgroup.size()-2;i<vsgroup.size();i++) {
      if(vNsgroup.at(i)==0) {
        vsgroup.at(i)=vsgroup.at(i-1);
        vNsgroup.at(i)=SGtoNSG(vsgroup.at(i));
      }
    }

    if(vNsgroup.size()!=3) for(uint ii=vNsgroup.size();ii<3;ii++) vNsgroup.push_back(0);
    if(vsgroup.size()!=3) for(uint ii=vsgroup.size();ii<3;ii++) vsgroup.push_back("NNN #0");

    return;
  }
}


// ***************************************************************************
// aflowlib::AflowlibLocator
// ***************************************************************************
namespace aflowlib { // move to web interface
  bool AflowlibLocator(const string& in,string& out,const string& mode) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string message = "";
    if(LDEBUG) cerr << XPID << "aflowlib::AflowlibLocator: BEGIN" << endl;

    if(mode!="AFLOWLIB_AUID2AURL" && mode!="AFLOWLIB_AURL2AUID" && mode!="AFLOWLIB_AUID2LOOP" && mode!="AFLOWLIB_AURL2LOOP") {
      message = "wrong mode=" + mode;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    if(XHOST_vAUID.size()==0) { init::InitGlobalObject("vLIBS"); }
    if(XHOST_vAURL.size()==0) { init::InitGlobalObject("vLIBS"); }
    if(XHOST_vLOOP.size()==0) { init::InitGlobalObject("vLIBS"); }
    if(XHOST_vAUID.size()!=XHOST_vAURL.size() || XHOST_vAUID.size()!=XHOST_vLOOP.size()) {
      message = "XHOST_vAUID.size()!=XHOST_vAURL.size() || XHOST_vAUID.size()!=XHOST_vLOOP.size()\n";
      message += "    XHOST_vAUID.size()=" + aurostd::utype2string<uint>(XHOST_vAUID.size()) + '\n';
      message += "    XHOST_vAURL.size()=" + aurostd::utype2string<uint>(XHOST_vAURL.size()) + '\n';
      message += "    XHOST_vLOOP.size()=" + aurostd::utype2string<uint>(XHOST_vLOOP.size());
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    out="";
    for(uint i=0;i<XHOST_vAUID.size()&&out.empty();i++) {
      if(mode=="AFLOWLIB_AUID2AURL" && XHOST_vAUID.at(i)==in) out=XHOST_vAURL.at(i);
      if(mode=="AFLOWLIB_AURL2AUID" && XHOST_vAURL.at(i)==in) out=XHOST_vAUID.at(i);
      if(mode=="AFLOWLIB_AUID2LOOP" && XHOST_vAUID.at(i)==in) out=XHOST_vLOOP.at(i);
      if(mode=="AFLOWLIB_AURL2LOOP" && XHOST_vAURL.at(i)==in) out=XHOST_vLOOP.at(i);
      //     cerr << i << endl;
    }
    if(LDEBUG) cerr << XPID << "aflowlib::AflowlibLocator: END" << endl;
    return !out.empty();
  }
} // namespace aflowlib

namespace aflowlib {
  string AflowlibLocator(string options, string mode) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "aflowlib::AflowlibLocator: BEGIN" << endl;
    if(mode!="AFLOWLIB_AUID2AURL" && mode!="AFLOWLIB_AURL2AUID" && mode!="AFLOWLIB_AUID2LOOP" && mode!="AFLOWLIB_AURL2LOOP") {
      string message = "wrong mode=" + mode;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()==0) {
      if(mode=="AFLOWLIB_AUID2AURL") init::ErrorOption(options,"aflowlib::AflowlibLocator","aflow --aflowlib_auid2aurl=auid1,auid2....");
      if(mode=="AFLOWLIB_AURL2AUID") init::ErrorOption(options,"aflowlib::AflowlibLocator","aflow --aflowlib_aurl2auid=aurl1,aurl2....");
      if(mode=="AFLOWLIB_AUID2LOOP") init::ErrorOption(options,"aflowlib::AflowlibLocator","aflow --aflowlib_auid2loop=auid1,auid2....");
      if(mode=="AFLOWLIB_AURL2LOOP") init::ErrorOption(options,"aflowlib::AflowlibLocator","aflow --aflowlib_aurl2loop=aurl1,aurl2....");
    } 
    // move on
    stringstream output;
    string locator;
    for(uint i=0;i<tokens.size();i++) {
      if(aflowlib::AflowlibLocator(tokens.at(i),locator,mode)) {
        output << locator << endl;
      } else {
        output << tokens.at(i) << " not found" << endl;
      }
    }

    if(LDEBUG) cerr << XPID << "aflowlib::AflowlibLocator: END" << endl;
    return output.str();
  }
} // namespace aflowlib

//AFLUX integration
//FR+CO20180329
//Obsolete with aurostd::xhttp //HE20220407
//namespace aflowlib {
//  bool APIget::establish(){
//    struct hostent * host = gethostbyname( Domain.c_str() );
//
//    //[CO20181226 - OBSOLETE]PORT=80;  //CO20180401
//
//    if ( (host == NULL) || (host->h_addr == NULL) ) {
//      cerr << "Error retrieving DNS information." << endl;
//      return false;
//    }
//
//    bzero(&client, sizeof(client));
//    client.sin_family = AF_INET;
//    client.sin_port = htons( PORT );
//    memcpy(&client.sin_addr, host->h_addr, host->h_length);
//
//    sock = socket(AF_INET, SOCK_STREAM, 0);
//
//    if (sock < 0) {
//      cerr << "Error creating socket." << endl;
//      return false;
//    }
//
//    if ( connect(sock, (struct sockaddr *)&client, sizeof(client)) < 0 ) {
//      close(sock);
//      cerr << "Could not connect" << endl;
//      return false;
//    }
//
//    stringstream ss;
//    ss << "GET " << API_Path << Summons << " HTTP/1.0\r\n" ;
//    //    cerr << "GET " << API_Path << Summons << " HTTP/1.0\r\n" ;
//    ss << "HOST: " << Domain << "\r\n";
//    ss << "Connection: close\r\n";
//    ss << "\r\n";
//    string request = ss.str();
//
//    if (send(sock, request.c_str(), request.length(), 0) != (int)request.length()) {
//      cerr << "Error sending request." << endl;
//      return false;
//    }
//    return true;
//  }
//  void APIget::reset( string a_Summons, string a_API_Path, string a_Domain ) {
//    if( a_Summons == "#" ) {
//      Summons = "";
//      //DX20210615 [OBSOLETE - old path] API_Path = "/search/API/?";
//      //DX20210615 [OBSOLETE - old domain] Domain = "aflowlib.duke.edu";
//      API_Path = "/API/aflux/?"; //DX20210615 - new path
//      Domain = "aflow.org"; //DX20210615 - new domain
//    } else {
//      Summons = a_Summons;
//      if( ! a_API_Path.empty() ) API_Path = a_API_Path;
//      if( ! a_Domain.empty() ) Domain = a_Domain;
//    }
//  }
//  ostream& operator<<( ostream& output, APIget& a ) {
//    char cur;
//    bool responsedata = false;
//    bool waslinefeed = false;
//    if( a.establish() ) {
//      while ( ! responsedata ) { //discard headers
//        read(a.sock, &cur, 1);
//        //cerr << cur << ":" << (int)cur << endl;
//        if( waslinefeed  && cur == '\r') responsedata = true;
//        if( cur == '\n' ) waslinefeed = true;
//        else waslinefeed = false;
//      };
//      read(a.sock, &cur, 1); //discard final \n in header \r\n\r\n
//      while ( read(a.sock, &cur, 1) > 0 ) output << cur; //cout << cur;
//      close(a.sock);
//    }
//    return output;
//  }
//}

//DX+FR20190206 - AFLUX functionality via command line - START
// ***************************************************************************
namespace aflowlib {
  string AFLUXCall(const aurostd::xoption& vpflow){

    // Performs AFLUX call based on summons input from command line
    if(vpflow.flag("AFLUX::USAGE")) {
      string usage="aflow --aflux=<summons>";
      string options="";
      //[CO20200624 - OBSOLETE]stringstream ss_usage;
      //[CO20200624 - OBSOLETE]init::ErrorOption(ss_usage,vpflow.getattachedscheme("AFLUX"),"aflowlib::AFLUXCall()",aurostd::liststring2string(usage,options));
      //[CO20200624 - OBSOLETE]return ss_usage.str();
      init::MessageOption(vpflow.getattachedscheme("AFLUX"),"aflowlib::AFLUXCall()",aurostd::liststring2string(usage,options));
      return "";
    }

    string summons = "";
    if(vpflow.flag("AFLUX::SUMMONS")) { //CO20200520 - AFLUX::SUMMONS
      summons=vpflow.getattachedscheme("AFLUX::SUMMONS"); //CO20200520 - AFLUX::SUMMONS
      // check if string is enclosed in double or single quotes
      // (since bash throws error for unprotected parentheses)
      if(!summons.empty() && ((summons[0] == '\"' && summons[summons.size()-1] == '\"') || (summons[0] == '\'' && summons[summons.size()-1] == '\''))){
        summons.erase(summons.begin()); summons.erase(summons.begin()+summons.size()-1);
      }
    }
    return AFLUXCall(summons);
  }
  string AFLUXCall(const vector<string>& matchbook){

    // Performs AFLUX call based on vector of matchbook entries
    string summons = aurostd::joinWDelimiter(matchbook,",");
    return AFLUXCall(summons);
  }
  string AFLUXCall(const string& summons){
    // Performs AFLUX call based on summons input
    // switched to aurostd::xhttp //HE20220407
    bool LDEBUG=(false || XHOST.DEBUG);
    string url = "http://aflow.org/API/aflux/?" + aurostd::httpPercentEncodingFull(summons);
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << ": Summons = " << summons << endl;
      cerr << __AFLOW_FUNC__ << ": URL = " << url << endl;
      cerr << __AFLOW_FUNC__ << ": Performing call ... please be patient ..." << endl;
    }
    //ME20220426 - add error handling
    int status_code = 0;
    string response = aurostd::httpGet(url, status_code);
    if (status_code < 200 || status_code >= 400) {
      string message = "Bad status code for AFLUX request (" + aurostd::utype2string<int>(status_code) + ").";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_HTTP_);
    }
    if ((response.find("Lux Fail") != string::npos)
      || (response.find("DB Fail") != string::npos)
      || (response.find("Count Fail") != string::npos)) {
      string message = "Bad response: " + response + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    return response;
  }
}

namespace aflowlib {
  vector<vector<std::pair<string,string> > > getPropertiesFromAFLUXResponse(const string& response){

    // Puts list of keyword-value pairs into a vector corresponding to each entry
    // Assumes the response format to be "format(aflow)", i.e., "|" delimiter
    // Here, pair.first=<keyword> and pair.second=<value>
    // In order to be general, all keywords and values are stored as a string 


    vector<vector<std::pair<string,string> > > properties_response;

    vector<string> entries,fields,key_value;
    aurostd::string2vectorstring(response,entries); //CO20200520
    //[CO20200520 - OBSOLETE]aurostd::string2tokens(response,entries,"\n");

    // for each entry in response
    for(uint e=0;e<entries.size();e++){
      // split into key-value pairs
      aurostd::string2tokens(entries[e],fields,"|");

      // properties for a particular entry
      vector<std::pair<string,string> > property_pairs;
      for(uint i=0;i<fields.size();i++){
        aurostd::string2tokens(fields[i],key_value,"=");
        if(key_value.size()<2 && !aurostd::substring2bool(fields[i],"example") && !aurostd::substring2bool(fields[i],"description")){
          string message = "Cannot find key-value pair splitting on \"=\" for the following field: \"" + fields[i] + "\".";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        }
        std::pair<string,string> property; 
        property.first = aurostd::RemoveWhiteSpaces(key_value[0]);  // key 
        // ME20220419 - some entries have = inside a value,
        // so remove first element and join the rest
        key_value.erase(key_value.begin());
        property.second = aurostd::RemoveWhiteSpaces(aurostd::joinWDelimiter(key_value, "=")); // value
        property_pairs.push_back(property);
      }
      properties_response.push_back(property_pairs);
    }
    return properties_response;
  }
}

//DX+FR20190206 - AFLUX functionality via command line - END

namespace aflowlib {  //CO20201220
  // ***************************************************************************
  // mergeEntries()
  // ***************************************************************************
  //simple helper function for loading multiple libraries together, will take
  //combinations of nested vectors and convert them to other nested vectors.
  //naries is output, entries_new is input
  //these assume vector<aflowlib::_aflowlib_entry> entries_new is all of the same
  //type, e.g., binaries of MnPd (e.g., MnPd2, Mn2Pd, etc., similar structure to LIBS)
  //also assumes ordered vector<vector<vector<aflowlib::_aflowlib_entry> > >& naries
  //outer most - unary, binary, etc.
  //next outer - species match, Mn, Pd, MnPd, etc.
  //inner most - entries
  bool mergeEntries(vector<vector<vector<_aflowlib_entry> > >& naries,const vector<vector<vector<_aflowlib_entry> > >& entries_new) {
    for (uint i = 0; i < entries_new.size(); i++) {
      if(!mergeEntries(naries, entries_new[i], true)) { return false; }  //structured data
    }
    return true;
  }
  bool mergeEntries(vector<vector<vector<_aflowlib_entry> > >& naries,const vector<vector<_aflowlib_entry> >& entries_new,bool entries_new_same_type) {
    for (uint i = 0; i < entries_new.size(); i++) {
      if(naries.size() <= i + 1) {
        naries.push_back(vector<vector<_aflowlib_entry> >(0));
      }
      if(!mergeEntries(naries[i], entries_new[i], entries_new_same_type, true)) {  //triple vector<> naries implies this structure
        return false;
      }
    }
    return true;
  }
  bool mergeEntries(vector<vector<vector<_aflowlib_entry> > >& naries,const vector<_aflowlib_entry>& entries_new,bool entries_new_same_type) {
    if(!entries_new.size()) { return true; }
    int index_layer1=-1, index_layer2=-1;
    if(entries_new_same_type) {
      if(!mergeEntries(naries, entries_new[0], index_layer1, index_layer2)) { return false; }
      for (uint i = 1; i < entries_new.size(); i++) { naries[index_layer1][index_layer2].push_back(entries_new[i]); }
    } else {
      for (uint i = 0; i < entries_new.size(); i++) {
        if(!entries_new[i].vspecies.size()) { return false; }  //what the heck is this?
        if(!mergeEntries(naries, entries_new[i], index_layer1, index_layer2)) { return false; }
      }
    }
    return true;
  }
  bool mergeEntries(vector<vector<vector<_aflowlib_entry> > >& naries,const _aflowlib_entry& entry_new) {
    int index_layer1 = -1;
    int index_layer2 = -1;
    return mergeEntries(naries, entry_new, index_layer1, index_layer2);
  }
  bool mergeEntries(vector<vector<vector<_aflowlib_entry> > >& naries,const _aflowlib_entry& entry_new,int& index_layer1,int& index_layer2) {
    if(!entry_new.vspecies.size()) { return false; }    //what the heck is this?
    while (naries.size() < entry_new.vspecies.size()) {  //assumes entries_new all have the same vspecies.size()
      naries.push_back(vector<vector<_aflowlib_entry> >(0));
    }
    index_layer1 = entry_new.vspecies.size() - 1;
    return mergeEntries(naries[index_layer1], entry_new, index_layer2, true);  //triple vector<> naries implies this structure
  }
  //naries takes on two forms depending on sort_by_species
  //if sort_by_species==true, then naries is truly a single nary (unary, binary, etc.) with
  //the second layer consisting of different species
  //otherwise, naries is the total entries (similar to naries above), where second layer
  //is unaries, binary, etc. (no layer with different species)
  //if sort_by_species and we are coming from LIBs (entries_new_same_type==true), then we don't need to check every entry, we already
  //know they have the same type (binary of same species)
  bool mergeEntries(vector<vector<_aflowlib_entry> >& naries,const vector<vector<vector<_aflowlib_entry> > >& entries_new,bool sort_by_species) {
    for (uint i = 0; i < entries_new.size(); i++) {
      if(!mergeEntries(naries, entries_new[i], true, sort_by_species)) { return false; }
    }
    return true;
  }
  bool mergeEntries(vector<vector<_aflowlib_entry> >& naries,const vector<vector<_aflowlib_entry> >& entries_new,bool entries_new_same_type,bool sort_by_species) {
    for (uint i = 0; i < entries_new.size(); i++) {
      if(!mergeEntries(naries, entries_new[i], entries_new_same_type, sort_by_species)) { return false; }
    }
    return true;
  }
  bool mergeEntries(vector<vector<_aflowlib_entry> >& naries,const vector<_aflowlib_entry>& entries_new,bool entries_new_same_type,bool sort_by_species) {
    if(!entries_new.size()) { return true; }
    int index;
    if(entries_new_same_type) {
      if(!mergeEntries(naries, entries_new[0], index, sort_by_species)) { return false; }
      for (uint i = 1; i < entries_new.size(); i++) { naries[index].push_back(entries_new[i]); }
    } else {
      for (uint i = 0; i < entries_new.size(); i++) {
        if(!mergeEntries(naries, entries_new[i], index, sort_by_species)) { return false; }
      }
    }
    return true;
  }
  //naries takes on two forms depending on sort_by_species
  //if sort_by_species==true, then naries is truly a single nary (unary, binary, etc.) with
  //the second layer consisting of different species
  //otherwise, naries is the total entries (similar to naries above), where second layer
  //is unaries, binary, etc. (no layer with different species)
  bool mergeEntries(vector<vector<_aflowlib_entry> >& naries,const _aflowlib_entry& entry_new,bool sort_by_species) {
    int index = -1;
    return mergeEntries(naries, entry_new, index, sort_by_species);
  }
  bool mergeEntries(vector<vector<_aflowlib_entry> >& naries,const _aflowlib_entry& entry_new,int& index,bool sort_by_species) {
    if(!entry_new.vspecies.size()) { return false; }
    index = -1;
    if(sort_by_species) {
      //all of naries is unary, binary, etc., now just need to index species
      for (uint i = 0; i < naries.size() && index == -1; i++) {
        //test of stupidity
        if(naries[i][0].vspecies.size() != entry_new.vspecies.size()) { return false; }
        if(naries[i][0].vspecies == entry_new.vspecies) { index = i; }
      }
      if(index == -1) {
        naries.push_back(vector<_aflowlib_entry>(0));
        index = naries.size() - 1;
      }
    } else {
      //just need to create space for unary, binary, etc.
      while (naries.size() < entry_new.vspecies.size()) {
        naries.push_back(vector<_aflowlib_entry>(0));
      }
      index = entry_new.vspecies.size() - 1;
    }
    naries[index].push_back(entry_new);
    return true;
  }
  bool mergeEntries(vector<_aflowlib_entry>& naries,const vector<vector<vector<_aflowlib_entry> > >& entries_new) {
    for (uint i = 0; i < entries_new.size(); i++) {
      for (uint j = 0; j < entries_new[i].size(); j++) {
        if(!mergeEntries(naries, entries_new[i][j])) { return false; }
      }
    }
    return true;
  }
  bool mergeEntries(vector<_aflowlib_entry>& naries,const vector<vector<_aflowlib_entry> >& entries_new) {
    for (uint i = 0; i < entries_new.size(); i++) {
      if(!mergeEntries(naries, entries_new[i])) { return false; }
    }
    return true;
  }
  bool mergeEntries(vector<_aflowlib_entry>& naries,const vector<_aflowlib_entry>& entries_new) {
    for (uint i = 0; i < entries_new.size(); i++) {
      if(!mergeEntries(naries, entries_new[i])) { return false; }
    }
    return true;
  }
  bool mergeEntries(vector<_aflowlib_entry>& naries,const _aflowlib_entry& entry_new) {
    naries.push_back(entry_new);
    return true;
  }
}  // namespace pflow

//DX20200929 - START
namespace aflowlib {
  string getSpaceGroupMatchbook(const vector<uint>& space_groups, uint relaxation_step){

    vector<string> vsummons(space_groups.size());

    for(uint i=0;i<space_groups.size();i++){
      vsummons[i] = getSpaceGroupMatchbook(space_groups[i], relaxation_step, false); //false - signals more than one space group
    }
    if(relaxation_step==_COMPARE_DATABASE_GEOMETRY_ORIGINAL_){ return "spacegroup_orig(" + aurostd::joinWDelimiter(vsummons,":") + ")"; } //DX20210615 - relaxation-step specific keyword
    else if(relaxation_step==_COMPARE_DATABASE_GEOMETRY_RELAX1_){ return "sg2(" + aurostd::joinWDelimiter(vsummons,":") + ")"; } //DX20210615 - relaxation-step specific keyword
    else if(relaxation_step==_COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_){ return "spacegroup_relax(" + aurostd::joinWDelimiter(vsummons,":") + ")"; } //DX20210615 - relaxation-step specific keyword
    return ""; //DX20210615 - no other relaxations supported, return empty so the compiler will not complain
  }
}

namespace aflowlib {
  string getSpaceGroupMatchbook(uint space_group_number, uint relaxation_step, bool only_one_sg){

    // Formats the space group summons for an AFLUX matchbook
    // This also grabs the relative enantiomorphs, since they
    // are simply mirror images of one another (structurally the same)
    // The summons are formatted differently depending on the relaxation type
    // (i.e. placement of comma(s) or lack thereof)
    // May need to be reformatted later with AFLOW+AFLUX integration

    // percent encode special characters first
    string squote_encoded = aurostd::PercentEncodeASCII('\''); // single quote
    string hashtag_encoded = aurostd::PercentEncodeASCII('#'); // hashtag/octothorpe
    string space_encoded = aurostd::PercentEncodeASCII(' '); // single space

    string space_group_summons = "";
    // check if enantiomorphic space group
    uint enantiomorph_space_group_number = SYM::getEnantiomorphSpaceGroupNumber(space_group_number);
    if(space_group_number == enantiomorph_space_group_number){
      // relaxed: need to match last in string, i.e., "*,<sg_symbol> <sg_number>" (comma necessary or we may grab the orig symmetry)
      if(relaxation_step==_COMPARE_DATABASE_GEOMETRY_ORIGINAL_){
        //DX20210615 [OBOSLETE] space_group_summons = "%27" + GetSpaceGroupName(space_group_number) + "%20%23" + aurostd::utype2string<int>(space_group_number) + ",%27*";
        space_group_summons = aurostd::utype2string<int>(space_group_number); //DX20210615
      }
      else if(relaxation_step==_COMPARE_DATABASE_GEOMETRY_RELAX1_){
        //DX20210615 [OBSOLETE] space_group_summons = "*%27," + GetSpaceGroupName(space_group_number) + "%20%23" + aurostd::utype2string<int>(space_group_number) + ",%27*";
        space_group_summons = squote_encoded + ",\"" + GetSpaceGroupName(space_group_number) + space_encoded + hashtag_encoded + aurostd::utype2string<int>(space_group_number) + "\"," + squote_encoded; //DX20210615
      }
      else if(relaxation_step==_COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_){
        //DX20210615 [OBSOLETE] space_group_summons = "*%27," + GetSpaceGroupName(space_group_number) + "%20%23" + aurostd::utype2string<int>(space_group_number) + "%27";
        space_group_summons = aurostd::utype2string<int>(space_group_number); //DX20210615
      }
      else{
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Unexpected relaxation step input: " + aurostd::utype2string(_COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_), _FILE_NOT_FOUND_);
      }
    }
    else { // need to get enantiomorph too
      if(relaxation_step==_COMPARE_DATABASE_GEOMETRY_ORIGINAL_){
        //DX20210615 [OBSOLETE] space_group_summons = "%27" + GetSpaceGroupName(space_group_number) + "%20%23" + aurostd::utype2string<int>(space_group_number) + ",%27*";
        //DX20210615 [OBSOLETE] space_group_summons += ":%27" + GetSpaceGroupName(enantiomorph_space_group_number) + "%20%23" + aurostd::utype2string<int>(enantiomorph_space_group_number) + ",%27*";
        space_group_summons = aurostd::utype2string<int>(space_group_number); //DX20210615
        space_group_summons += ":" + aurostd::utype2string<int>(enantiomorph_space_group_number); //DX20210615
      }
      else if(relaxation_step==_COMPARE_DATABASE_GEOMETRY_RELAX1_){
        //DX20210615 [OBSOLETE] space_group_summons = "*%27," + GetSpaceGroupName(space_group_number) + "%20%23" + aurostd::utype2string<int>(space_group_number) + ",%27*";
        //DX20210615 [OBSOLETE] space_group_summons += ":*%27," + GetSpaceGroupName(enantiomorph_space_group_number) + "%20%23" + aurostd::utype2string<int>(enantiomorph_space_group_number) + ",%27*";
        space_group_summons = squote_encoded + ",\"" + GetSpaceGroupName(space_group_number) + space_encoded + hashtag_encoded + aurostd::utype2string<int>(space_group_number) + "\"," + squote_encoded; //DX20210615
        space_group_summons += ":" + squote_encoded + ",\"" + GetSpaceGroupName(enantiomorph_space_group_number) + space_encoded + hashtag_encoded + aurostd::utype2string<int>(enantiomorph_space_group_number) + "\"," + squote_encoded; //DX20210615
      }
      else if(relaxation_step==_COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_){
        //DX20210615 [OBSOLETE] space_group_summons = "*%27," + GetSpaceGroupName(space_group_number) + "%20%23" + aurostd::utype2string<int>(space_group_number) + "%27";
        //DX20210615 [OBSOLETE] space_group_summons += ":*%27," + GetSpaceGroupName(enantiomorph_space_group_number) + "%20%23" + aurostd::utype2string<int>(enantiomorph_space_group_number) + "%27";
        space_group_summons = aurostd::utype2string<int>(space_group_number); //DX20210615
        space_group_summons += ":" + aurostd::utype2string<int>(enantiomorph_space_group_number); //DX20210615
      }
      else{
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Unexpected relaxation step input: " + aurostd::utype2string(_COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_), _FILE_NOT_FOUND_);
      }
    }

    // ---------------------------------------------------------------------------
    // if there is only one space group in the query, put summons in sg2();
    // otherwise this is done outside this function
    if(only_one_sg) { 
      //DX20210615 [OBSOLETE] space_group_summons = "sg2(" + space_group_summons + ")";
      if(relaxation_step==_COMPARE_DATABASE_GEOMETRY_ORIGINAL_){ space_group_summons = "spacegroup_orig(" + space_group_summons + ")"; } //DX20210615 - relaxation-step specific keyword
      else if(relaxation_step==_COMPARE_DATABASE_GEOMETRY_RELAX1_){ space_group_summons = "sg2(" + space_group_summons + ")"; } //DX20210615 - relaxation-step specific keyword
      else if(relaxation_step==_COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_){ space_group_summons = "spacegroup_relax(" + space_group_summons + ")"; } //DX20210615 - relaxation-step specific keyword
    }

    return space_group_summons;
  }
}
//DX20200929 - END

namespace aflowlib {
  uint WEB_Aflowlib_Entry(string options,ostream& oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cout << XPID << "aflowlib::WEB_Aflowlib_Entry: begin<br>" << endl;

    stringstream num_prec;
    vector<string> voptions;
    aurostd::string2tokens(options,voptions,",");
    string errormsg = "";
    if(voptions.size()==0) {
      errormsg = "--aflowlib= has no arguments.";
      //init::ErrorOption(options,"aflowlib::WEB_Aflowlib_Entry","aflow --aflowlib=entry");  //CO20200624 - soft patch for FR+web // OBSOLETE - use errormsg instead
      //[CO+ME20200731 - ErrorOption throws an error, so no need to return]return 0; //CO20200624 - 0 is error here
    } 

    // move on
    for(uint ioption=0;ioption<voptions.size();ioption++) {
      //  oss << voptions.at(ioption) << endl;
      string option=voptions.at(ioption); //aurostd::args2attachedstring(argv,"--aflowlib=",(string) "nan");
      if(option.at(option.size()-1)=='/'|| option.at(option.size()-1)=='.') option.erase(option.end()-1,option.end()-0); //  some demoronization
      if(option.at(0)=='/'|| option.at(0)=='.') option.erase(option.begin(),option.begin()+1); //  some demoronization
      string directory="";
      string directory_RAW="",directory_LIB="",directory_WEB="";
      string directory_AUID_LIB="",directory_AUID_RAW="",directory_AUID_WEB="";
      string url_WEB;
      string label="";
      // string line_gif="<br><img border=0 width=60% height=2 src=http://materials.duke.edu/auro/images/line.gif><br><br>";
      // string art058_link=" https://doi.org/10.1016/j.commatsci.2010.05.010";
      // string art064_link=" https://doi.org/10.1021/co200012w";
      // string icsd_link=" https://www.fiz-karlsruhe.com/icsd.html";
      // string aflow_ael_readme=" http://materials.duke.edu/AFLOW/README_AFLOW_AEL.TXT"; //CO20180817
      // string art096_link=" https://doi.org/10.1103/PhysRevB.90.174107";
      // string art100_link=" https://www.nature.com/articles/sdata20159";
      // string aflow_agl_readme=" http://materials.duke.edu/AFLOW/README_AFLOW_AGL.TXT"; //CO20180817
      // string art115_link=" https://doi.org/10.1103/PhysRevMaterials.1.015401"; //CO20180817
      // string aflow_sym_readme=" http://materials.duke.edu/AFLOW/README_AFLOW_SYM.TXT"; //CO20180817
      // string art135_link=" https://doi.org/10.1107/S2053273318003066"; //CO20180817

      //DX20180817
      // string bravais_lattice_orig_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#bravais_lattice_orig";
      // string bravais_lattice_relax_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#bravais_lattice_relax";
      // string lattice_system_orig_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_system_orig";
      // string lattice_variation_orig_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_variation_orig";
      // string lattice_system_relax_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_system_relax";
      // string lattice_variation_relax_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#lattice_variation_relax";
      // string Pearson_symbol_orig_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#pearson_symbol_orig";
      // string Pearson_symbol_relax_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#pearson_symbol_relax";
      // string sg_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#sg";
      // string sg2_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#sg2";
      // string spacegroup_orig_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#spacegroup_orig";
      // string spacegroup_relax_wiki_link=" http://aflowlib.duke.edu/aflowwiki/doku.php?id=documentation:all_keywords&#spacegroup_relax";

      aflowlib::_aflowlib_entry aentry;

      xoption vflags;
      vflags.flag("FLAG::PREAMBLE",TRUE);
      vflags.flag("FLAG::CALCULATION",TRUE);
      vflags.flag("FLAG::JMOL",TRUE);
      vflags.flag("FLAG::EDATA_ORIG",FALSE);
      vflags.flag("FLAG::EDATA_RELAX",TRUE);
      vflags.flag("FLAG::THERMODYNAMICS",TRUE);
      vflags.flag("FLAG::MAGNETIC",TRUE);
      vflags.flag("FLAG::ELECTRONIC",FALSE);     // will setup later
      vflags.flag("FLAG::SCINTILLATION",TRUE);   // will setup later
      vflags.flag("FLAG::AGL",FALSE);            // will setup later
      vflags.flag("FLAG::AEL",FALSE);            // will setup later
      vflags.flag("FLAG::BADER",FALSE);          // will setup later
      vflags.flag("FLAG::POCC",FALSE);           // will setup later  //CO20201220

      // check if ICSD inside (anyway)
      vector<string> vline,tokens;

      string html_TAB=" target=\"_blank\"";

      if(LDEBUG) cout << "WEB_Aflowlib_Entry: [1]<br>" << endl;
      if(LDEBUG) cout << "WEB_Aflowlib_Entry: [4]<br>" << endl;
      if(LDEBUG) cout << "WEB_Aflowlib_Entry: option=" << option << endl;

      // START SEARCH

      vflags.flag("FLAG::FOUND",FALSE);
      string catalog="",auid="",strtmp="";

      // **********************************************************************************************************
      // TRY AUID
      // **********************************************************************************************************
      // ME20200707 - Also check for AUIDs without the aflow: prefix
      string option_orig = option;
      if ((option.size() == 16) && aurostd::_ishex(option)) option = "aflow:" + option;
      if(!vflags.flag("FLAG::FOUND") && aurostd::substring2bool(aurostd::tolower(option),"aflow:")) { // CHECK AUID
        if(LDEBUG) cout << "WEB_Aflowlib_Entry: option=" << option << endl;
        string auid=aurostd::tolower(option);
        if(auid.size()!=22) {
          errormsg="aflowlib::WEB_Aflowlib_Entry:_error_on_size_of_auid="+auid;
        } else {
          // NEW
          directory_AUID_LIB=init::AFLOW_Projects_Directories("AUID")+"/"+auid.substr(0,8); for(uint i=8;i<=20;i+=2) directory_AUID_LIB+="/"+auid.substr(i,2);  // splitting aflow:ab/cd..
          directory_AUID_WEB=directory_AUID_LIB+"/WEB";
          directory_AUID_RAW=directory_AUID_LIB+"/RAW";
          directory_AUID_LIB=directory_AUID_LIB+"/LIB";
          if(LDEBUG) cout << "WEB_Aflowlib_Entry: directory_AUID_LIB=" << directory_AUID_LIB << endl;
          if(LDEBUG) cout << "WEB_Aflowlib_Entry: directory_AUID_RAW=" << directory_AUID_RAW << endl;
          if(LDEBUG) cout << "WEB_Aflowlib_Entry: directory_AUID_WEB=" << directory_AUID_WEB << endl;
          directory="";
          if(!aurostd::FileExist(directory_AUID_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            errormsg="aflowlib::WEB_Aflowlib_Entry:_entry_does_not_exist="+directory_AUID_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT;
          } else {
            _aflowlib_entry entry_tmp(string(directory_AUID_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT));
            directory=entry_tmp.aurl;
            auid=entry_tmp.aurl;
            aurostd::StringSubst(directory,"aflowlib.duke.edu:","");
            aurostd::StringSubst(directory,"materials.duke.edu:","");
            aurostd::StringSubst(directory,"AFLOWDATA/ICSD_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/ICSD_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB0_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB0_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB1_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB1_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB2_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB2_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB3_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB3_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB4_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB4_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB5_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB5_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB6_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB6_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB7_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB7_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB8_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB8_WEB/","");
            aurostd::StringSubst(directory,"AFLOWDATA/LIB9_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/LIB9_WEB/","");
            vflags.flag("FLAG::AUID",TRUE);
            vflags.flag("FLAG::FOUND",TRUE);
            catalog=entry_tmp.catalog;
            label=directory;
            for(uint ilat=0;ilat<BRAVAIS_LATTICES.size();ilat++) aurostd::StringSubst(label,BRAVAIS_LATTICES[ilat]+"/",""); //HE20220420 switch to global bravais lattices list
            aurostd::StringSubst(label,"/",".");
          }
        }
      }
      //ME20200707 - Restore
      option = option_orig;

      // **********************************************************************************************************
      // TRY DIRECTORY
      // **********************************************************************************************************
      if(!vflags.flag("FLAG::FOUND")) { // tests with proto name
        string dir2test;
        dir2test=option;

        vector<string> vdir2test;
        for(uint i0=0;i0<dir2test.size();i0++) { // allow all identical indices...
          for(uint i1=i0;i1<dir2test.size();i1++) { // allow all identical indices
            for(uint i2=i1;i2<dir2test.size();i2++) { // allow all identical indices
              //	    for(uint i3=i2;i3<dir2test.size();i3++) { // allow all identical indices
              dir2test=option;
              if(dir2test.at(i0)=='.') dir2test.at(i0)='/';
              if(dir2test.at(i1)=='.') dir2test.at(i1)='/';
              if(dir2test.at(i2)=='.') dir2test.at(i2)='/';
              //      if(dir2test.at(i3)=='.') dir2test.at(i3)='/';
              bool found=false;
              for(uint i=0;i<vdir2test.size()&&!found;i++) if(dir2test==vdir2test.at(i)) found=true;
              if(!found) {
                vdir2test.push_back(dir2test);
              }
            }
            //	  }
          }
        }   
        if(LDEBUG) cout << "WEB_Aflowlib_Entry: testing dir2test=" << dir2test << endl;
        for(uint i=0;i<vdir2test.size()&&!vflags.flag("FLAG::FOUND");i++) { // allow i=j=k so that 1 OR 2 OR 3 dots are also tested
          dir2test=vdir2test.at(i);
          //		cout << "testing(" << i << ") = " << dir2test << endl;
          for(uint ilat=0;ilat<BRAVAIS_LATTICES.size()&&!vflags.flag("FLAG::FOUND");ilat++) { //HE20220420 switch to global bravais lattices list
            if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("ICSD")+"/RAW/"+BRAVAIS_LATTICES[ilat]+"/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
              catalog="ICSD";directory=BRAVAIS_LATTICES[ilat]+"/"+dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);
            }
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB0")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB0";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB1")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB1";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB2")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB2";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB3")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB3";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB4")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB4";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB5")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB5";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB6")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB6";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB7")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB7";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB8")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB8";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
          if(!vflags.flag("FLAG::FOUND") && aurostd::FileExist(init::AFLOW_Projects_Directories("LIB9")+"/RAW/"+dir2test+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
            catalog="LIB9";directory=dir2test;label=option;vflags.flag("FLAG::FOUND",TRUE);		
          }
        }
      }

      // **********************************************************************************************************
      // TRY ICSD LINK
      // **********************************************************************************************************
      if(!vflags.flag("FLAG::FOUND")) { // icsd link
        //ME20200707 - Remove possible icsd: prefix
        option_orig = option;
        if (option.find("icsd") == 0) option = option.substr(5);
        string directory_ICSD2LINK=init::AFLOW_Projects_Directories("AUID")+"/icsd:/"+option;
        aurostd::StringSubst(directory_ICSD2LINK,"ICSD:","icsd:");
        aurostd::StringSubst(directory_ICSD2LINK,"icsd:icsd:","icsd:");    
        //	cerr << directory_ICSD2LINK << endl;
        if(aurostd::FileExist(directory_ICSD2LINK+"/RAW/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
          _aflowlib_entry entry_tmp(string(directory_ICSD2LINK+"/RAW/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT));
          directory=entry_tmp.aurl;
          auid=entry_tmp.aurl;
          aurostd::StringSubst(directory,"aflowlib.duke.edu:","");
          aurostd::StringSubst(directory,"materials.duke.edu:","");
          aurostd::StringSubst(directory,"AFLOWDATA/ICSD_RAW/","");aurostd::StringSubst(directory,"AFLOWDATA/ICSD_WEB/","");
          vflags.flag("FLAG::ICSD",TRUE);
          vflags.flag("FLAG::FOUND",TRUE);
          catalog=entry_tmp.catalog;
          //ME20200923 - Use the last subdirectory or XHOST.label will be inconsistent.
          //For example, using --aflowlib=7bed936e9d5a44ca results in XHOST.label = Ni1Sn1Ti1_ICSD_174568
          //whereas using --aflowlib=174568 results in XHOST.label=FCC.Ni1Sn1Ti1_ICSD_174568
          //[OBSOLETE] label=directory;
          vector<string> tokens;
          aurostd::string2tokens(directory, tokens, "/");
          label = tokens.back();
          //	  cerr << directory_ICSD2LINK+"/RAW/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << endl;
        }
        //ME20200707 - Restore
        option = option_orig;
      }

      // **********************************************************************************************************
      // CONSIDERED FOUND
      // **********************************************************************************************************    
      if(vflags.flag("FLAG::FOUND")) {
        if(catalog=="ICSD") {
          vflags.flag("FLAG::ICSD",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("ICSD")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("ICSD")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("ICSD")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/ICSD_WEB/"+directory;
        }

        if(catalog=="LIB0") {
          vflags.flag("FLAG::LIB0",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB0")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB0")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB0")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB0_RAW/"+directory;
        }
        if(vflags.flag("FLAG::FOUND") && catalog=="LIB1") {
          vflags.flag("FLAG::LIB1",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB1")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB1")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB1")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB1_RAW/"+directory;
        }
        if(catalog=="LIB2") {
          vflags.flag("FLAG::LIB2",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB2")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB2")+"/RAW/"+directory; // June 2016
          directory_RAW=init::AFLOW_Projects_Directories("LIB2")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB2_RAW/"+directory; // May 2014
        }
        if(catalog=="LIB3") {
          vflags.flag("FLAG::LIB3",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB3")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB3")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB3")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB3_WEB/"+directory;
        }
        if(catalog=="LIB4") {
          vflags.flag("FLAG::LIB4",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB4")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB4")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB4")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB4_WEB/"+directory;
        }
        if(catalog=="LIB5") {
          vflags.flag("FLAG::LIB5",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB5")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB5")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB5")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB5_WEB/"+directory;
        }
        if(catalog=="LIB6") {
          vflags.flag("FLAG::LIB6",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB6")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB6")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB6")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB6_WEB/"+directory;
        }
        if(catalog=="LIB7") {
          vflags.flag("FLAG::LIB7",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB7")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB7")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB7")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB7_WEB/"+directory;
        }
        if(catalog=="LIB8") {
          vflags.flag("FLAG::LIB8",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB8")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB8")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB8")+"/RAW/"+directory;
          url_WEB="/AFLOWDATA/LIB8_WEB/"+directory;
        }
        if(catalog=="LIB9") {
          vflags.flag("FLAG::LIB9",TRUE);
          directory_LIB=init::AFLOW_Projects_Directories("LIB9")+"/LIB/"+directory;
          directory_WEB=init::AFLOW_Projects_Directories("LIB9")+"/WEB/"+directory;
          directory_RAW=init::AFLOW_Projects_Directories("LIB9")+"/RAW/"+directory;	
          url_WEB="/AFLOWDATA/LIB9_WEB/"+directory;
        }
      } else {
        errormsg="aflowlib::WEB_Aflowlib_Entry:_entry_does_not_exist="+option;
      }

      // got it  ?
      if(!aurostd::FileExist(directory_RAW+"/"+_AFLOWIN_)) directory_RAW="";

      if(!directory.empty()) { // play with aentry.entry
        aurostd::StringSubst(label,"/",".");
        aentry.file2aflowlib(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,oss);  //   oss << aentry.entry << endl;
        directory_AUID_LIB=init::AFLOW_Projects_Directories("AUID")+"/"+aflowlib::auid2directory(aentry.auid);
        directory_AUID_WEB=directory_AUID_LIB+"/WEB";
        directory_AUID_RAW=directory_AUID_LIB+"/RAW";
        directory_AUID_LIB=directory_AUID_LIB+"/LIB";    
        aurostd::string2tokens(aentry.sg2,tokens,"#");if(tokens.size()>0) aentry.sg2=tokens.at(tokens.size()-1);
        if(aentry.vfiles_WEB.size()==0) aentry.vfiles_WEB=aentry.vfiles;
      }

      // ME20210209 - Show original structure when symmetry changed
      if (vflags.flag("FLAG::FOUND")) {
        bool same_symmetry = ((aentry.spacegroup_orig == aentry.spacegroup_relax)
            && (aentry.Wyckoff_multiplicities_orig == aentry.Wyckoff_multiplicities)
            && (aentry.Wyckoff_site_symmetries_orig == aentry.Wyckoff_site_symmetries_orig));
        vflags.flag("FLAG::EDATA_ORIG", !same_symmetry);
      }

      if(vflags.flag("FLAG::FOUND")) {
        // check AGL/AEL
        vflags.flag("FLAG::ELECTRONIC",aurostd::substring2bool(aentry.vloop,"bands"));
        vflags.flag("FLAG::SCINTILLATION",aurostd::substring2bool(aentry.vloop,"bands"));
        vflags.flag("FLAG::AGL",aurostd::substring2bool(aentry.vloop,"agl"));
        vflags.flag("FLAG::AEL",aurostd::substring2bool(aentry.vloop,"ael"));
        vflags.flag("FLAG::BADER",aurostd::substring2bool(aentry.vloop,"bader"));
        vflags.flag("FLAG::POCC",aurostd::substring2bool(aentry.vloop,"pocc")); //CO20201220
        if(vflags.flag("FLAG::POCC")){vflags.flag("FLAG::JMOL",FALSE);} //CO20201220 - no POSCAR/CIF to plot for PARTCAR (yet)
      }
      stringstream aflowlib_json;
      if(vflags.flag("FLAG::FOUND")) {
        strtmp=aurostd::efile2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON);
        aurostd::StringSubst(strtmp,"}\n",""); // remove trailing bracket add it at the end
        aflowlib_json << strtmp;
      } else {
        aflowlib_json << "{DUMMY"; // will remove at the end
      }
      stringstream aflowlib_out;
      if(vflags.flag("FLAG::FOUND")) {
        strtmp=aurostd::efile2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT);
        aurostd::StringSubst(strtmp,"\n",""); // remove trailing bracket add it at the end
        aflowlib_out << strtmp;
      } else {
        aflowlib_out << "DUMMY"; // will remove at the end
      }

      // adding pieces 
      //   aurostd::StringSubst(aflowlib_json,"}\n",""); // remove trailing bracket add it at the end

      // XHOST.hostname
      aflowlib_json << "," << "\"XHOST.hostname\":" << "\"" << XHOST.hostname << "\"";
      aflowlib_out << " | " << "XHOST.hostname=" << XHOST.hostname;
      // option
      aflowlib_json << "," << "\"XHOST.option\":" << "\"" << option << "\"";
      aflowlib_out << " | " << "XHOST.option=" << option;
      // label
      aflowlib_json << "," << "\"XHOST.label\":" << "\"" << label << "\"";
      aflowlib_out << " | " << "XHOST.label=" << label;
      // directory
      aflowlib_json << "," << "\"XHOST.directory\":" << "\"" << directory << "\"";
      aflowlib_out << " | " << "XHOST.directory=" << directory;

      if(vflags.flag("FLAG::FOUND")) {
        // directory_LIB
        aflowlib_json << "," << "\"XHOST.directory_LIB\":" << "\"" << directory_LIB << "\"";
        aflowlib_out << " | " << "XHOST.directory_LIB=" << directory_LIB;
        // directory_RAW
        aflowlib_json << "," << "\"XHOST.directory_RAW\":" << "\"" << directory_RAW << "\"";
        aflowlib_out << " | " << "XHOST.directory_RAW=" << directory_RAW;
        // directory_WEB
        aflowlib_json << "," << "\"XHOST.directory_WEB\":" << "\"" << directory_WEB << "\"";
        aflowlib_out << " | " << "XHOST.directory_WEB=" << directory_WEB;
        // directory_AUID_LIB
        aflowlib_json << "," << "\"XHOST.directory_AUID_LIB\":" << "\"" << directory_AUID_LIB << "\"";
        aflowlib_out << " | " << "XHOST.directory_AUID_LIB=" << directory_AUID_LIB;
        // directory_AUID_RAW
        aflowlib_json << "," << "\"XHOST.directory_AUID_RAW\":" << "\"" << directory_AUID_RAW << "\"";
        aflowlib_out << " | " << "XHOST.directory_AUID_RAW=" << directory_AUID_RAW;
        // directory_AUID_WEB
        aflowlib_json << "," << "\"XHOST.directory_AUID_WEB\":" << "\"" << directory_AUID_WEB << "\"";
        aflowlib_out << " | " << "XHOST.directory_AUID_WEB=" << directory_AUID_WEB;
      }

      if(vflags.flag("FLAG::FOUND")) {
        // XHOST.FLAG::ICSD
        aflowlib_json << "," << "\"XHOST.FLAG::ICSD\":" << (vflags.flag("FLAG::ICSD")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::ICSD=" << (vflags.flag("FLAG::ICSD")?"1":"0");
        // XHOST.FLAG::LIB0
        aflowlib_json << "," << "\"XHOST.FLAG::LIB0\":" << (vflags.flag("FLAG::LIB0")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB0=" << (vflags.flag("FLAG::LIB0")?"1":"0");
        // XHOST.FLAG::LIB1
        aflowlib_json << "," << "\"XHOST.FLAG::LIB1\":" << (vflags.flag("FLAG::LIB1")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB1=" << (vflags.flag("FLAG::LIB1")?"1":"0");
        // XHOST.FLAG::LIB2
        aflowlib_json << "," << "\"XHOST.FLAG::LIB2\":" << (vflags.flag("FLAG::LIB2")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB2=" << (vflags.flag("FLAG::LIB2")?"1":"0");
        // XHOST.FLAG::LIB3
        aflowlib_json << "," << "\"XHOST.FLAG::LIB3\":" << (vflags.flag("FLAG::LIB3")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB3=" << (vflags.flag("FLAG::LIB3")?"1":"0");
        // XHOST.FLAG::LIB4
        aflowlib_json << "," << "\"XHOST.FLAG::LIB4\":" << (vflags.flag("FLAG::LIB4")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB4=" << (vflags.flag("FLAG::LIB4")?"1":"0");
        // XHOST.FLAG::LIB5
        aflowlib_json << "," << "\"XHOST.FLAG::LIB5\":" << (vflags.flag("FLAG::LIB5")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB5=" << (vflags.flag("FLAG::LIB5")?"1":"0");
        // XHOST.FLAG::LIB6
        aflowlib_json << "," << "\"XHOST.FLAG::LIB6\":" << (vflags.flag("FLAG::LIB6")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB6=" << (vflags.flag("FLAG::LIB6")?"1":"0");
        // XHOST.FLAG::LIB7
        aflowlib_json << "," << "\"XHOST.FLAG::LIB7\":" << (vflags.flag("FLAG::LIB7")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB7=" << (vflags.flag("FLAG::LIB7")?"1":"0");
        // XHOST.FLAG::LIB8
        aflowlib_json << "," << "\"XHOST.FLAG::LIB8\":" << (vflags.flag("FLAG::LIB8")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB8=" << (vflags.flag("FLAG::LIB8")?"1":"0");
        // XHOST.FLAG::LIB9
        aflowlib_json << "," << "\"XHOST.FLAG::LIB9\":" << (vflags.flag("FLAG::LIB9")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::LIB9=" << (vflags.flag("FLAG::LIB9")?"1":"0");
        // XHOST.FLAG::AUID
        aflowlib_json << "," << "\"XHOST.FLAG::AUID\":" << (vflags.flag("FLAG::AUID")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::AUID=" << (vflags.flag("FLAG::AUID")?"1":"0");
      }
      // XHOST.FLAG::FOUND
      aflowlib_json << "," << "\"XHOST.FLAG::FOUND\":" << (vflags.flag("FLAG::FOUND")?"true":"false");
      aflowlib_out << " | " << "XHOST.FLAG::FOUND=" << (vflags.flag("FLAG::FOUND")?"1":"0");
      //   if(!vflags.flag("FLAG::FOUND"))
      {
        // errormsg
        aflowlib_json << "," << "\"XHOST.errormsg\":" << "\"" << errormsg << "\"";
        aflowlib_out << " | " << "XHOST.errormsg=" << errormsg;
      }

      if(vflags.flag("FLAG::FOUND")) {
        // XHOST.FLAG::PREAMBLE
        aflowlib_json << "," << "\"XHOST.FLAG::PREAMBLE\":" << (vflags.flag("FLAG::PREAMBLE")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::PREAMBLE=" << (vflags.flag("FLAG::PREAMBLE")?"1":"0");
        // XHOST.FLAG::CALCULATION
        aflowlib_json << "," << "\"XHOST.FLAG::CALCULATION\":" << (vflags.flag("FLAG::CALCULATION")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::CALCULATION=" << (vflags.flag("FLAG::CALCULATION")?"1":"0");
        // XHOST.FLAG::JMOL
        aflowlib_json << "," << "\"XHOST.FLAG::JMOL\":" << (vflags.flag("FLAG::JMOL")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::JMOL=" << (vflags.flag("FLAG::JMOL")?"1":"0");
        // XHOST.FLAG::EDATA_ORIG
        aflowlib_json << "," << "\"XHOST.FLAG::EDATA_ORIG\":" << (vflags.flag("FLAG::EDATA_ORIG")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::EDATA_ORIG=" << (vflags.flag("FLAG::EDATA_ORIG")?"1":"0");
        // XHOST.FLAG::EDATA_RELAX
        aflowlib_json << "," << "\"XHOST.FLAG::EDATA_RELAX\":" << (vflags.flag("FLAG::EDATA_RELAX")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::EDATA_RELAX=" << (vflags.flag("FLAG::EDATA_RELAX")?"1":"0");
        // XHOST.FLAG::THERMODYNAMICS
        aflowlib_json << "," << "\"XHOST.FLAG::THERMODYNAMICS\":" << (vflags.flag("FLAG::THERMODYNAMICS")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::THERMODYNAMICS=" << (vflags.flag("FLAG::THERMODYNAMICS")?"1":"0");
        // XHOST.FLAG::MAGNETIC
        aflowlib_json << "," << "\"XHOST.FLAG::MAGNETIC\":" << (vflags.flag("FLAG::MAGNETIC")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::MAGNETIC=" << (vflags.flag("FLAG::MAGNETIC")?"1":"0");
        // XHOST.FLAG::ELECTRONIC
        aflowlib_json << "," << "\"XHOST.FLAG::ELECTRONIC\":" << (vflags.flag("FLAG::ELECTRONIC")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::ELECTRONIC=" << (vflags.flag("FLAG::ELECTRONIC")?"1":"0");
        // XHOST.FLAG::SCINTILLATION
        aflowlib_json << "," << "\"XHOST.FLAG::SCINTILLATION\":" << (vflags.flag("FLAG::SCINTILLATION")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::SCINTILLATION=" << (vflags.flag("FLAG::SCINTILLATION")?"1":"0");
        // XHOST.FLAG::AGL
        aflowlib_json << "," << "\"XHOST.FLAG::AGL\":" << (vflags.flag("FLAG::AGL")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::AGL=" << (vflags.flag("FLAG::AGL")?"1":"0");
        // XHOST.FLAG::AEL
        aflowlib_json << "," << "\"XHOST.FLAG::AEL\":" << (vflags.flag("FLAG::AEL")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::AEL=" << (vflags.flag("FLAG::AEL")?"1":"0");
        // XHOST.FLAG::BADER
        aflowlib_json << "," << "\"XHOST.FLAG::BADER\":" << (vflags.flag("FLAG::BADER")?"true":"false");
        aflowlib_out << " | " << "XHOST.FLAG::BADER=" << (vflags.flag("FLAG::BADER")?"1":"0");

        //ME20191004 START
        // Grab compressed files
        if(XHOST.vflag_control.flag("PRINT_MODE::JSON") || !XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
          string content;
          // fgroup for JMOL applet
          if (aurostd::EFileExist(directory_RAW + "/aflow.fgroup.bands.json")) {
            content = aurostd::efile2string(directory_RAW + "/aflow.fgroup.bands.json");
          } else if (aurostd::EFileExist(directory_RAW + "/aflow.fgroup.relax.json")) {
            content = aurostd::efile2string(directory_RAW + "/aflow.fgroup.relax.json");
          }
          aflowlib_json << ", \"fgroup\":" << (content.empty()?"null":content);

          content = "";
          if (vflags.flag("FLAG::ELECTRONIC")) {
            //ME20200616 - Made less dependent on file name conventions
            //string system_name = KBIN::ExtractSystemName(directory_LIB);
            vector<string> vfiles;
            aurostd::DirectoryLS(directory_RAW, vfiles);
            uint nfiles = vfiles.size();
            for (uint f = 0; f < nfiles; f++) {
              if (vfiles[f].find("_bandsdata.json") != string::npos) content = aurostd::efile2string(directory_RAW + "/" + vfiles[f]);
            }
          }
          aflowlib_json << ", \"bandsdata\":" << (content.empty()?"null":content);
        }
        //ME20191004 STOP
      }

      //ME20191217 START
      // additional web output
      aflowlib_json << "," << "\"aflow_version\":\"" << AFLOW_VERSION << "\"";
      aflowlib_out << "|" << "aflow_version=" << AFLOW_VERSION;
      aflowlib_json << "," << "\"aflow_build_date\":\"" << TODAY << "\""; //CO20200624 - more semantic
      aflowlib_out << "|" << "aflow_build_date=" << TODAY;  //CO20200624 - more semantic
      //ME20191217 STOP

      // XHOST.machine_type
      aflowlib_json << "," << "\"XHOST.machine_type\":" << "\"" << XHOST.machine_type << "\"";
      aflowlib_out << " | " << "XHOST.machine_type=" << XHOST.machine_type;
      // XHOST.user
      aflowlib_json << "," << "\"XHOST.user\":" << "\"" << XHOST.user << "\"";
      aflowlib_out << " | " << "XHOST.user=" << XHOST.user;
      // XHOST.group
      aflowlib_json << "," << "\"XHOST.group\":" << "\"" << XHOST.group << "\"";
      aflowlib_out << " | " << "XHOST.group=" << XHOST.group;
      // XHOST.shell
      aflowlib_json << "," << "\"XHOST.shell\":" << "\"" << XHOST.shell << "\"";
      aflowlib_out << " | " << "XHOST.shell=" << XHOST.shell;
      // XHOST.progname
      aflowlib_json << "," << "\"XHOST.progname\":" << "\"" << XHOST.progname << "\"";
      aflowlib_out << " | " << "XHOST.progname=" << XHOST.progname;
      // XHOST.generator
      aflowlib_json << "," << "\"XHOST.generator\":" << "\"" << "aflowlib::WEB_Aflowlib_Entry" << "\"";
      aflowlib_out << " | " << "XHOST.generator=" << "aflowlib::WEB_Aflowlib_Entry";

      // wrap up
      // TXT
      if(XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
        oss << aflowlib_out.str() << endl;
      }
      // JSON
      if(XHOST.vflag_control.flag("PRINT_MODE::JSON") || !XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
        aflowlib_json << "}";
        strtmp=aflowlib_json.str();
        aurostd::StringSubst(strtmp,"{DUMMY,","{");
        //    oss << "[" << option << "]" << endl;
        oss << strtmp << endl;
      }
    }
    return voptions.size();
  }
}


#endif  // _AURO_IMPLEMENTATIONS_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

