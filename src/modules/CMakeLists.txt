
set(ael_source
        AEL/aflow_ael_elastic_fit.cpp
        AEL/aflow_ael_elasticity.cpp
        AEL/aflow_ael_get_stress.cpp
)

set(agl_source
        AGL/aflow_agl_debye.cpp
        AGL/aflow_agl_electronic.cpp
        AGL/aflow_agl_eqn_state.cpp
        AGL/aflow_agl_get_ev.cpp
        AGL/aflow_agl_hugoniot.cpp
        AGL/aflow_agl_plot_write.cpp
        AGL/aflow_agl_polynomial.cpp
        AGL/aflow_agl_rungibbs.cpp
        AGL/aflow_agl_thermal.cpp
)

set(apl_source
        APL/aflow_aapl_cluster.cpp
        APL/aflow_aapl_ifcs.cpp
        APL/aflow_aapl_tcond.cpp
        APL/aflow_apl_atomic_disp.cpp
        APL/aflow_apl_doscalc.cpp
        APL/aflow_apl_fccalc.cpp
        APL/aflow_apl_kphonons.cpp
        APL/aflow_apl_pdisc.cpp
        APL/aflow_apl_phoncalc.cpp
        APL/aflow_apl_qmesh.cpp
        APL/aflow_apl_supercell.cpp
        APL/aflow_apl_thermalpc.cpp
)

set(cce_source
        CCE/aflow_cce.cpp
)

set(compare_source
        COMPARE/aflow_compare_structure.cpp
        COMPARE/aflow_compare_structure_function.cpp
)

set(cumulants_source
        CUMULANTS/aflow_cumulants.cpp
)

set(gfa_source
        GFA/aflow_gfa.cpp
)

set(hull_source
        HULL/aflow_nhull.cpp
        HULL/aflow_nhull_entry.cpp
        HULL/aflow_nhull_latex.cpp
        HULL/aflow_nhull_stats.cpp
        HULL/aflow_nhull_util.cpp
        HULL/aflow_nhull_facet.cpp
        HULL/aflow_nhull_facet.h
)

set(pocc_source
        POCC/aflow_pocc.cpp
        POCC/aflow_pocc_ael_agl.cpp
        POCC/aflow_pocc_apl.cpp
        POCC/aflow_pocc_old.cpp
        POCC/aflow_pocc_qha.cpp
        POCC/aflow_pocc_thermo.cpp
)

set(prototypes_source
        PROTOTYPES/aflow_anrl.cpp
        PROTOTYPES/aflow_xproto.cpp
        PROTOTYPES/aflow_xprototype.cpp
)

set(qca_source
        QCA/aflow_qca.cpp
)

set(qha_source
        QHA/aflow_qha.cpp
)

set(sym_source
        SYM/aflow_spacegroup.cpp
        SYM/aflow_symmetry.cpp
        SYM/aflow_symmetry_spacegroup.cpp
        SYM/aflow_symmetry_spacegroup_functions.cpp
        SYM/aflow_symmetry_spacegroup_ITC_library.cpp
        SYM/aflow_wyckoff.cpp
        HULL/aflow_nhull_entry.cpp
        HULL/aflow_nhull_entry.h
        HULL/aflow_nhull_facet.cpp
        HULL/aflow_nhull_facet.h
)

add_library(modules OBJECT
        ${ael_source}
        ${agl_source}
        ${apl_source}
        ${cce_source}
        ${compare_source}
        ${cumulants_source}
        ${gfa_source}
        ${hull_source}
        ${pocc_source}
        ${prototypes_source}
        ${qca_source}
        ${qha_source}
        ${sym_source}
)

target_include_directories(modules
        PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_BINARY_DIR}
)

target_link_libraries(modules PRIVATE config)

# reduce optimize for certain data files to avoid extreme long build times in RELEASE builds
SET_SOURCE_FILES_PROPERTIES(SYM/aflow_wyckoff.cpp PROPERTIES COMPILE_FLAGS -O0)
SET_SOURCE_FILES_PROPERTIES(SYM/aflow_symmetry_spacegroup_ITC_library.cpp PROPERTIES COMPILE_FLAGS -O0)
