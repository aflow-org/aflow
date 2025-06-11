
#ifndef AFLOW_DEFS_H
#define AFLOW_DEFS_H

#include <cmath>
#include <string>

// #define  _AFLOW_TEMP_PRESERVE_  // to preseve /tmp files for debug

//[CO20200502 - moved to aurostd.h]#define NNN   -123456
//[CO20200502 - moved to aurostd.h]#define GCC_VERSION (__GNUC__ * 10000  + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#define _ANRL_NOWEB_ //DX
// hard-coded prototype generator (ANRL/ subdirectory required) //DX20200623
// to revert to the hard-coded prototypes, do the following sequence:
// 1) set USE_HARDCODED_PROTOTYPES to true in aflow_makefile.cpp
// 2) compile
// 3) run aflow --makefile
// 4) set USE_HARDCODED_PROTOTYPES (below) to true
// 5) recompile
#define USE_HARDCODED_PROTOTYPES false

// toggle symbolic math
// (for now it is coupled with USE_HARDCODED_PROTOTYPES, although it does not have to be)
#define USE_SYMBOLIC_SOURCE !(USE_HARDCODED_PROTOTYPES) // true

//ZERO PRECISION DEFINITIONS - TIGHT (DEFAULT) AND LOOSE
#define _ZERO_PRECISION_ 10
#define _ZERO_TOL_ std::pow(10,-_ZERO_PRECISION_) //DX
#define _ZERO_PRECISION_LOOSE_ 3
#define _ZERO_TOL_LOOSE_ std::pow(10,-_ZERO_PRECISION_LOOSE_) //DX
#define _DOUBLE_PRECISION_ 8
#define _DOUBLE_TOL_ std::pow(10,-_DOUBLE_PRECISION_)
#define _FLOAT_PRECISION_ 6
#define _FLOAT_TOL_ std::pow(10,-_FLOAT_PRECISION_)  //ME20200519 - tolerance for float precision
//PRECISION and TOLERANCE definitions
#define _DOUBLE_WRITE_PRECISION_MAX_ 14  //CO20180509 - used for xstrctures
#define _DOUBLE_WRITE_PRECISION_ 12  //CO20180509 - used in writing doubles in qmvasp
#define _AFLOWLIB_STOICH_PRECISION_ _DOUBLE_PRECISION_ //[CO20200731 - too many different precisions... just use default]9  //CO20200731
#define _AFLOWLIB_DATA_DOUBLE_PREC_ _DOUBLE_PRECISION_ //[CO20200731 - too many different precisions... just use default]6 //CO20200731
#define _AFLOWLIB_DATA_GEOMETRY_PREC_ _DOUBLE_PRECISION_ //[CO20200731 - too many different precisions... just use default]7 //CO20200731
#define _AFLOW_POCC_PRECISION_ _DOUBLE_PRECISION_ //8 //must be less than _DOUBLE_WRITE_PRECISION_MAX_, which is currently set to 14
#define _AFLOW_POCC_ZERO_TOL_ std::pow(10,-_AFLOW_POCC_PRECISION_)
#define _XPROTO_TOO_CLOSE_ERROR_ 0.60 // was 0.75
#define _XPROTO_ZERO_VOL_ _FLOAT_TOL_  //CO20190218

#define _AFLOW_MAX_ARGV_ 1024 //CO20211104 - moved from aflowlib_libraries.cpp

//MESSAGE defaults - CO20200502
#define _AFLOW_MESSAGE_DEFAULTS_ "user,host,pid,time" //tid //CO20200624 - only depends on XHOST (not aflags)

//XRD
#define XRAY_RADIATION_COPPER_Kalpha 1.5418   //Angstroms     //CO20190622


//moved from avasp.cpp for broader access (chull.cpp)
#define SPECIE_TRANSITION_METALS std::string("Ag,Au,Cd,Co,Cr_pv,Cu_pv,Fe_pv,Hf_pv,Hg,Ir,La,Mn_pv,Mo_pv,Nb_sv,Ni_pv,Os_pv,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sc_sv,Ta_pv,Tc_pv,Ti_sv,V_sv,W_pv,Y_sv,Zn,Zr_sv")
#define SPECIE_RAW_LIB2U SPECIE_TRANSITION_METALS
#define SPECIE_RAW_LIB2 std::string("Ag,Al,As,Au,B_h,Ba_sv,Be_sv,Bi_d,Br,Ca_sv,Cd,Cl,Co,Cr_pv,Cu_pv,Fe_pv,Ga_h,Ge_h,Hf_pv,Hg,In_d,Ir,K_sv,La,Li_sv,Mg_pv,Mn_pv,Mo_pv,Na_pv,Nb_sv,Ni_pv,Os_pv,P,Pb_d,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sb,Sc_sv,Se,Si,Sn,Sr_sv,Ta_pv,Tc_pv,Te,Ti_sv,Tl_d,V_sv,W_pv,Y_sv,Zn,Zr_sv")

#define SPECIE_RAW_LIB3 std::string("Ag,Al,As,Au,B_h,Ba_sv,Be_sv,Bi_d,Br,Ca_sv,Cd,Cl,Co,Cr_pv,Cu_pv,Fe_pv,Ga_h,Ge_h,Hf_pv,Hg,In_d,Ir,K_sv,La,Li_sv,Mg_pv,Mn_pv,Mo_pv,Na_sv,Nb_sv,Ni_pv,Os_pv,P,Pb_d,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sb,Sc_sv,Se,Si,Sn,Sr_sv,Ta_pv,Tc_pv,Te,Ti_sv,Tl_d,V_sv,W_pv,Y_sv,Zn,Zr_sv")
//#define SPECIE_RAW_LIB3 std::string("Ag,Au,Cd,Co,Cr_pv,Cu_pv,Fe_pv,Hf_pv,Hg,Ir,La,Mn_pv,Mo_pv,Nb_sv,Ni_pv,Os_pv,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sc_sv,Ta_pv,Tc_pv,Ti_sv,V_sv,W_pv,Y_sv,Zn,Zr_sv")
//#define SPECIE_RAW_LIB3 std::string("Ag,Al,As,Au,B_h,Bi_d,Cd,Co,Cr_pv,Cu_pv,Fe_pv,Ga_h,Ge_h,Hf_pv,Hg,In_d,Ir,La,Mg_pv,Mn_pv,Mo_pv,Nb_sv,Ni_pv,Os_pv,P,Pb_d,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sb,Sc_sv,Se,Si,Sn,Ta_pv,Te,Tc_pv,Ti_sv,V_sv,W_pv,Y_sv,Zn,Zr_sv")

//Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,La,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg - not in order!
#define SPECIE_RAW_LIB4 SPECIE_TRANSITION_METALS


//XELEMENT_PROPERTIES_ALL (define early)
#define _AFLOW_XELEMENT_PROPERTIES_ALL_ "name,symbol,Z,period,group,series,block,mass,volume_molar,volume,area_molar_Miedema,valence_std,valence_iupac,valence_PT,valence_s,valence_p,valence_d,valence_f,density_PT,crystal,crystal_structure_PT,spacegroup,spacegroup_number,variance_parameter_mass,lattice_constants,lattice_angles,phase,radius_Saxena,radius_PT,radius_covalent_PT,radius_covalent,radius_VanDerWaals_PT,radii_Ghosh08,radii_Slatter,radii_Pyykko,conductivity_electrical,electronegativity_Pauling,hardness_chemical_Ghosh,electronegativity_Pearson,electronegativity_Ghosh,electronegativity_Allen,oxidation_states,oxidation_states_preferred,electron_affinity_PT,energies_ionization,work_function_Miedema,density_line_electron_WS_Miedema,energy_surface_0K_Miedema,chemical_scale_Pettifor,Mendeleev_number,temperature_boiling,temperature_melting,enthalpy_fusion,enthalpy_vaporization,enthalpy_atomization_WE,energy_cohesive,specific_heat_PT,critical_pressure,critical_temperature_PT,thermal_expansion,conductivity_thermal,hardness_mechanical_Brinell,hardness_mechanical_Mohs,hardness_mechanical_Vickers,hardness_chemical_Pearson,hardness_chemical_Putz,hardness_chemical_RB,modulus_shear,modulus_Young,modulus_bulk,Poisson_ratio_PT,modulus_bulk_x_volume_molar_Miedema,magnetic_type_PT,susceptibility_magnetic_mass,susceptibility_magnetic_volume,susceptibility_magnetic_molar,temperature_Curie,refractive_index,color_PT,HHIP,HHIR,xray_scatt" //CO20201111
#define _ENERGIES_IONIZATION_MAX_AFLOWMACHL_ 5

//MONITOR_VASP
#define VERBOSE_MONITOR_VASP false
#define AFLOW_MEMORY_TAG "AFLOW ERROR: AFLOW_MEMORY"

// definitions for MULTHREADS
//#define MAX_ALLOCATABLE_PTHREADS     1024
#define MAX_ALLOCATABLE_PTHREADS     256
#define PTHREADS_DEFAULT 8


#define _AFLOWIN_DEFAULT_     std::string("aflow.in")  //CO20210302
#define _AFLOWIN_AEL_DEFAULT_ std::string("ael_aflow.in")  //CO20210302
#define _AFLOWIN_AGL_DEFAULT_ std::string("agl_aflow.in")  //CO20210302
#define _AFLOWIN_QHA_DEFAULT_ std::string("aflow_qha.in")  //CO20210302 - moved from APL/aflow_qha.cpp

#define _AFLOWIN_AEL_VARIANTS_ std::string("ael_aflow.in,aflow_ael.in")  //CO20210302
#define _AFLOWIN_AGL_VARIANTS_ std::string("agl_aflow.in,aflow_agl.in")  //CO20210302

#define _AFLOWLOCK_DEFAULT_     std::string("LOCK")  //CO20210302
#define _AFLOWLOCK_AEL_DEFAULT_ std::string("ael.LOCK")  //CO20210302
#define _AFLOWLOCK_AGL_DEFAULT_ std::string("agl.LOCK")  //CO20210302
#define _AFLOWLOCK_QHA_DEFAULT_ std::string("LOCK.qha")  //CO20210302 - moved from APL/aflow_qha.cpp

#define _AFLOWLOCK_AEL_VARIANTS_ std::string("ael.LOCK,LOCK.ael")  //CO20210302
#define _AFLOWLOCK_AGL_VARIANTS_ std::string("agl.LOCK,LOCK.agl")  //CO20210302


// --------------------------------------------------------------------------
// definitions for aflow
// aflow2 default definitions
#define AFLOW_MATERIALS_SERVER_DEFAULT        std::string("materials.duke.edu")
#define AFLOW_WEB_SERVER_DEFAULT              std::string("nietzsche.mems.duke.edu")
#define AFLOWLIB_SERVER_DEFAULT               std::string("aflowlib.duke.edu")
#define AFLOWLIB_MATERIALS_SERVER             std::string("aflow.org")
#define AFLOWLIB_CONSORTIUM_STRING            std::string("aflow.org")
#define _XENTRY_ std::string("index.php")

#define DEFAULT_KBIN_ALIEN_BIN        std::string("ls -las")
#define DEFAULT_KBIN_MATLAB_BIN       std::string("/usr/local/bin/matlab -nodesktop -nosplash -nodisplay ")

#define QSUB_COMMAND_DEFAULT          "qsub"
#define QSUB_PARAMS_DEFAULT           " "

#define KBIN_SYMMETRY_SGROUP_RADIUS_DEFAULT 3.0
#define KBIN_SYMMETRY_SGROUP_MAX_NUMBER 1000000

#define KBIN_SUBDIRECTORIES           std::string("ARUN.")


#define ALIEN_INPUT_FILE_NAME_DEFAULT  "./input"
#define ALIEN_EXTERNAL_INPUT_DEFAULT   "../input_external"
#define ALIEN_OUTPUT_FILE_NAME_DEFAULT  "./output"

// aflow1 definitions (soon to be obsolete)
#define _MPI_NP_STRINGS_ "MPI_NP","mpi_np","-MPI_NP","-mpi_np"
#define _MPI_NCPUS_DEF_ 4
#define VASP_OPTIONS_MPI_DEFAULT         ""
#define VASPLS_BIN_POSTFIX_DEFAULT       "LS"
#define GRND_BIN_DEFAULT                 "./grnd_intel"

#define _VASP_POSCAR_MODE_EXPLICIT_START_ "[VASP_POSCAR_MODE_EXPLICIT]START"  //SD20220501
#define _VASP_POSCAR_MODE_EXPLICIT_STOP_ "[VASP_POSCAR_MODE_EXPLICIT]STOP"  //SD20220501
#define _VASP_POSCAR_MODE_EXPLICIT_START_P_ _VASP_POSCAR_MODE_EXPLICIT_START_ "."  //CO20200624
#define _VASP_POSCAR_MODE_EXPLICIT_STOP_P_ _VASP_POSCAR_MODE_EXPLICIT_STOP_ "."  //CO20200624

#define _AFLOWLIB_ENTRY_SEPARATOR_       std::string(" | ")

// --------------------------------------------------------------------------
// definition for frozsl files
#define _FROZSL_VASPSETUP_FILE_ "./aflow.frozsl_vaspsetup_file"
// --------------------------------------------------------------------------
// definitions for WEB PHP
#define AFLOW_PHP_APOOL_REFERENCES       std::string("19,20,49,50,51,53,54,55,56,57,59,61,62,63,65,66,67,70,71,74,75,76,81,87,99")

// --------------------------------------------------------------------------
// Definitions
//#define DEFAULT_AFLOW_FIND_PARAMETERS "-follow"
#define DEFAULT_AFLOW_FIND_PARAMETERS_NORMAL     std::string("-follow")
#define DEFAULT_AFLOW_FIND_PARAMETERS_NOLEAF     std::string("-noleaf -follow")
#define BUFFER_MAXLEN 1024

// --------------------------------------------------------------------------
// include all prototypes for aflow
#ifndef SWAP
#define SWAP(a,b)      {temp=(a);(a)=(b);(b)=temp;}
#endif
#define RCYCLIC(a,b,c) {temp=(c);(b)=(a);(c)=(b);a=temp;}
#define LCYCLIC(a,b,c) {temp=(a);(a)=(b);(b)=(c);c=temp;}

#define NANOPARTICLE_RADIUS_DEFAULT   10.0
#define NANOPARTICLE_DISTANCE_DEFAULT 10.0

//BANDGAP  //CO20191110
#define _METALGAP_ -1.0*AUROSTD_NAN
#define _METALEDGE_ -1.0

//moved from avasp
#define _AFLOWINPAD_ 60

//DX20180131 - add symmetry definitions - START
// symmetry
#define SG_SETTING_1    1
#define SG_SETTING_2    2
#define SG_SETTING_ANRL 3
//DX20180131 - add symmetry definitions - END

//DX20191122 START
// atom environment modes
#define ATOM_ENVIRONMENT_MODE_1    1 // minimum coordination shell - element split
#define ATOM_ENVIRONMENT_MODE_2    2 // [FUTURE] out to a given radius
#define ATOM_ENVIRONMENT_MODE_3    3 // [FUTURE] largest gap in radial distribution function (GFA)
//DX20191122 END

// --------------------------------------------------------------------------

// Structures for flags and properties to share FAST !
// STRUCTURES
#define AFLOWIN_SEPARATION_LINE  std::string("[AFLOW] ************************************************************************************************************************** ")
#define SEPARATION_LINE_DASH std::string("------------------------------------------------------------------------------------------------") //DX+CO20210429 - generic dash-line separator (used between symmetry operators)
#define SEPARATION_LINE_DASH_SHORT std::string("---------------------------------------------------------------------------") //DX+CO20210429 - generic dash-line separator, short (used in symmetry log output)

#define PRINT_NULL_JSON false //DX20210430 - add global flag to print "null" for empty JSON values


#define _COORDS_FRACTIONAL_ 0
#define _COORDS_CARTESIAN_  1
#define _UPDATE_LATTICE_VECTORS_TO_ABCANGLES_   2
#define _UPDATE_LATTICE_ABCANGLES_TO_VECTORS_   3

#define _PGROUP_ 1             // for point group lattice
#define _PGROUPK_ 6            // for point group klattice
#define _PGROUP_XTAL_ 7        // for point group crystal
#define _PGROUPK_XTAL_ 8       // for point group kcrystal
#define _PGROUPK_PATTERSON_ 9  // for point group Patterson //DX20200129
#define _FGROUP_ 2             // for factor group
#define _SGROUP_ 3             // for space group
#define _AGROUP_ 4             // for site positions point group
#define _IATOMS_ 5             // for equivalent atoms

#define NUM_ELEMENTS (103+1)  // up to Uranium


#define MAX_TITLE_SIZE 512

#define IOAFLOW_AUTO   0
#define IOVASP_AUTO    1
#define IOVASP_POSCAR  2
#define IOVASP_ABCCAR  3
#define IOVASP_WYCKCAR 4
#define IOQE_AUTO      5
#define IOQE_GEOM      6
#define IOABINIT_AUTO  7
#define IOABINIT_GEOM  8
#define IOAIMS_AUTO    9
#define IOAIMS_GEOM   10
#define IOCIF         11 //DX20180723
#define IOELK_AUTO    12 //DX20200310
#define IOELK_GEOM    13 //DX20200310
#define IOATAT_STR    14 //SD20220114
#define IOAFLUX_QRY   15 //HE20220210
#define IOLMP_DATA    16 //SD20240111

#define NOSG std::string("NNN #0")

#define _EQUIV_FPOS_EPS_    2.0e-5    // NOV 2009 Israel  used to be 1.0e-6 too small for ICSD
#define _pocc_no_sublattice_ -1

#define RADIANS 0
#define DEGREES  1
#define _calculate_symmetry_default_sgroup_radius_   2.0

// aflow_xproto.cpp
#define _HTQC_PROJECT_STRING_ "HTQC Project"
#define _TERNARY_PROJECT_STRING_ "HTQC^3 Project"
#define _ICSD_STRING_ "(icsd library)"
#define _ICSD_PROJECT_STRING_ "ICSD Project"
#define _ICSD_AFLOWLIB_STRING_ "(icsd_aflowlib library)"

// for HTQC
#define STRUCTURE_MODE_NONE             0
#define STRUCTURE_MODE_RAW              1
#define STRUCTURE_MODE_ABC              2
#define STRUCTURE_MODE_WYC              3
#define STRUCTURE_MODE_ICSD             4
#define STRUCTURE_MODE_HTQC_ICSD        5
#define STRUCTURE_MODE_USE              6
#define STRUCTURE_MODE_REMOVE           7
#define STRUCTURE_MODE_SPECIES          8
#define STRUCTURE_MODE_SWAP_AB          9
#define STRUCTURE_MODE_SWAP_BC         10
#define STRUCTURE_MODE_SWAP_AC         11
#define STRUCTURE_MODE_SWAP_XY         12
#define STRUCTURE_MODE_PRIM            13
#define STRUCTURE_MODE_CONVENTIONAL    14
#define STRUCTURE_MODE_VOLUME          15
#define LIBRARY_MODE_ICSD               0
#define LIBRARY_MODE_ICSD_AFLOWLIB      1
#define LIBRARY_MODE_HTQC               2
#define LIBRARY_MODE_HTQC_ICSD          3
#define LIBRARY_MODE_HTQC_ICSD_AFLOWLIB 4
#define LIBRARY_MODE_LIB0               5
#define LIBRARY_MODE_LIB3               6
#define LIBRARY_MODE_LIB4               7
#define LIBRARY_MODE_LIB5               8
#define LIBRARY_MODE_LIB6               9
#define LIBRARY_MODE_LIB7               10
#define LIBRARY_MODE_LIB8               11
#define LIBRARY_MODE_LIB9               12
#define LIBRARY_MODE_PROTOTYPE          13
#define LIBRARY_MODE_XSTRUCTURE         14
#define LIBRARY_MODE_AUID               15
#define LIBRARY_MODE_ARUN               16  //ME20181226

//DX20180710 - updated - #define DOI_ANRL " [ANRL doi: arXiv:1607.02532]"
#define DOI_ANRL " [ANRL doi: 10.1016/j.commatsci.2017.01.017 (part 1), doi: 10.1016/j.commatsci.2018.10.043 (part 2)]" //DX20180710 - updated //DX20190214 updated part 2 doi
#define DOI_POCC " [POCC doi: 10.1021/acs.chemmater.6b01449]"

#define _AVASP_PSEUDOPOTENTIAL_AUTO_ std::string("AUTO")
#define _AVASP_PSEUDOPOTENTIAL_DELIMITER_ std::string(":")
#define _AVASP_PSEUDOPOTENTIAL_POTENTIAL_TYPE_ std::string("TYPE") //CO20191020
#define _AVASP_PSEUDOPOTENTIAL_POTENTIAL_COMPLETE_ std::string("COMPLETE")

// aflow_mix.cpp  aflow_nomix.cpp   aflow_mix_pauling.cpp
#define MISCIBILITY_SYSTEM_NOT_STUDIED  3
#define MISCIBILITY_SYSTEM_SOLUTION     2
#define MISCIBILITY_SYSTEM_MISCIBLE     1
#define MISCIBILITY_SYSTEM_NOMIX        0
#define MISCIBILITY_SYSTEM_UNKNOWN     -1
#define MISCIBILITY_SYSTEM_CUTOFF     200
#define MIEDEMA_MIX_SLOPE 3.069              // Miedema Rule Table 1a Physica 100B (1980) 1-28

#define DEFAULT_TOTAL_LAYERS 10 //CO20190601 - move to .aflow.rc eventually
#define DEFAULT_V3_ANGLE_DEVIATION 5.0  //CO20190803 - move to .aflow.rc eventually

#define _VAR_THRESHOLD_STD_ 0.001
#define _Y_CORR_THRESHOLD_STD_ 0.0
#define _SELF_CORR_THRESHOLD_STD_ 0.95

// aflow_xprototype.h stuff by DAVID
#define _AFLOW_PROTOTYPE_ENCYCLOPEDIA_ std::string("http://aflow.org/CrystalDatabase/")

#endif //AFLOW_DEFS_H
