cmake_minimum_required(VERSION 3.20)

# set the project name and version
project(aflow VERSION 4.0.4)
# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
include(FetchContent)

# set definitions
if( CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    list( APPEND CMAKE_CXX_FLAGS "-fp-model=precise")
endif()
# find Git and if available set GIT_HASH variable
find_package(Git QUIET)
set(AFLOW_VERSION "${CMAKE_PROJECT_VERSION}")
set(AFLOW_PACK_VERSION "${CMAKE_PROJECT_VERSION}")
if(GIT_FOUND)
    execute_process(
            COMMAND git log -1 --pretty=format:%h
            OUTPUT_VARIABLE GIT_HASH
            RESULT_VARIABLE GIT_CHECK
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
    )
    if(${GIT_CHECK} EQUAL "0")
        set(AFLOW_VERSION "${CMAKE_PROJECT_VERSION}-${GIT_HASH}")
        set(AFLOW_PACK_VERSION "${AFLOW_VERSION}")
    endif ()
endif()
add_compile_definitions(AFLOW_VERSION="${AFLOW_VERSION}")
add_compile_definitions(AFLOW_VERSION_SHORT="${CMAKE_PROJECT_VERSION}")
message("Preparing AFLOW ${AFLOW_VERSION}")

string(TIMESTAMP TODAY_TMP "%Y-%m-%d")
add_compile_definitions(TODAY="${TODAY_TMP}")
add_compile_definitions(AFLOW_MULTITHREADS_ENABLE)
if(APPLE)
    add_compile_definitions(_MACOSX_)
elseif(UNIX)
    add_compile_definitions(_UNIX_)
endif()

# ensure other external libraries are available
find_package(Python REQUIRED)
# add thread support
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
# add libarchive
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # Homebrew ships libarchive keg only, include dirs have to be set manually
    execute_process(
            COMMAND brew --prefix libarchive
            OUTPUT_VARIABLE LIBARCHIVE_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE STATUS
    )
    if(STATUS AND NOT STATUS EQUAL 0)
        message("libarchive not installed through brew")
    else()
        set(LibArchive_ROOT "${LIBARCHIVE_PREFIX}")
    endif()

    execute_process(
            COMMAND brew --prefix openssl
            OUTPUT_VARIABLE OPENSSL_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE STATUS
    )
    if(STATUS AND NOT STATUS EQUAL 0)
        message("openSSL not installed through brew")
    else()
        set(openSSL_ROOT "${OPENSSL_PREFIX}")
    endif()

    execute_process(
            COMMAND brew --prefix curl
            OUTPUT_VARIABLE CURL_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE STATUS
    )
    if(STATUS AND NOT STATUS EQUAL 0)
        message("curl not installed through brew")
    else()
        set(CURL_ROOT "${CURL_PREFIX}")
    endif()

endif()

find_package(LibArchive REQUIRED)
include_directories(${LibArchive_INCLUDE_DIR})
find_package(OpenSSL REQUIRED)
include_directories(${OPENSSL_INCLUDE_DIR})
find_package(CURL REQUIRED)
include_directories(${CURL_INCLUDE_DIR})


set(aflowlib_sources
        aflowlib/aflowlib_database.cpp
        aflowlib/aflowlib_entry_loader.cpp
        aflowlib/aflowlib_libraries.cpp
        aflowlib/aflowlib_libraries_scrubber.cpp
        aflowlib/aflowlib_web_interface.cpp
        aflowlib/aflowlib_web_outreach.cpp
)

set(contrib_sources
        contrib/aflow_contrib_cormac.cpp
        contrib/aflow_contrib_shidong_cluster_expansion.cpp
        contrib/aflow_contrib_wahyu.cpp
)

set(flow_sources
        flow/aflow_avasp.cpp
        flow/aflow_iaims.cpp
        flow/aflow_ialien.cpp
        flow/aflow_ifrozsl.cpp
        flow/aflow_ivasp.cpp
        flow/aflow_kaims.cpp
        flow/aflow_kalien.cpp
        flow/aflow_kbin.cpp
        flow/aflow_kvasp.cpp
        flow/aflow_modules.cpp
        flow/aflow_oaims.cpp
        flow/aflow_ovasp.cpp
        flow/aflow_pflow_funcs.cpp
        flow/aflow_pflow_print.cpp
        flow/aflow_xclasses.cpp
        flow/aflow_bader.cpp
        flow/aflow_estructure.cpp
        flow/aflow_machl.cpp
        flow/aflow_mix.cpp
        flow/aflow_mix_experiments.cpp
        flow/aflow_nomix.cpp
        flow/aflow_plotter.cpp
        flow/aflow_xpseudopotentials.cpp
)

set(interface_sources
        interfaces/aflow_matlab.cpp
        interfaces/aflow_qsub.cpp
        interfaces/aflow_sqlite.cpp
        interfaces/aflow_symbolic.cpp
        interfaces/aflow_voro.cpp
        interfaces/aflow_pthreads.cpp
)

set(structure_sources
        structure/aflow_defects.cpp
        structure/aflow_lattice.cpp
        structure/aflow_surface.cpp
        structure/aflow_xatom.cpp
        structure/aflow_xelement.cpp
        structure/aflow_xstructure.cpp
        structure/aflow_xstructure_parser.cpp
)

set(main_sources
        aflow.cpp
        aflow_aconvasp.cpp
        aflow_aconvasp_main.cpp
        aflow_aflowrc.cpp
        aflow_init.cpp
        aflow_xhost.cpp
        aflow_xthread.cpp
)

# main aflow executable target
add_executable(aflow
        ${main_sources}
        ${aflowlib_sources}
        ${contrib_sources}
        ${flow_sources}
        ${interface_sources}
        ${structure_sources}
)
target_link_options(aflow PRIVATE "LINKER:-v")
target_include_directories(aflow PUBLIC ${CMAKE_SOURCE_DIR})

# set up and link internal libraries
add_subdirectory(AUROSTD)
target_link_libraries(aflow PUBLIC aurostd)
add_subdirectory(modules)
target_link_libraries(aflow PUBLIC modules)

# set up and link packaged external libraries
add_subdirectory(extern/SQLITE)
target_link_libraries(aflow PUBLIC sqlite)
add_subdirectory(extern/SYMBOLICCPLUSPLUS)
target_link_libraries(aflow PUBLIC symbolic)

# set up and link internal tests
add_subdirectory(test)
target_link_libraries(aflow PUBLIC aflow_tests)

# link external libraries
target_link_libraries(aflow PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(aflow PRIVATE ${LibArchive_LIBRARIES})
target_link_libraries(aflow PRIVATE ${OPENSSL_CRYPTO_LIBRARIES})
target_link_libraries(aflow PRIVATE ${CURL_LIBRARIES})
if(CMAKE_DL_LIBS)
    target_link_libraries(aflow PUBLIC ${CMAKE_DL_LIBS})
endif()

# reduce optimize for certain data files to avoid extreme long build times in RELEASE builds
SET_SOURCE_FILES_PROPERTIES(flow/aflow_xpseudopotentials.cpp PROPERTIES COMPILE_FLAGS -O0)

# installation targets
install(TARGETS aflow RUNTIME)

# cpack config
set(CPACK_SOURCE_IGNORE_FILES "/CVS/;/\\.svn/;\\.swp$;\\.#;/#;.*~;cscope.*;cmake-;CMakeUser;\\.xz$;")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Automatic-FLOW for materials discovery")
set(CPACK_PACKAGE_CONTACT "Stefano Curtarolo <stefano@duke.edu>")
set(CPACK_PACKAGE_VERSION ${AFLOW_PACK_VERSION})
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/../LICENSE")

if(APPLE)
    set(CPACK_GENERATOR "TGZ;STGZ")
else ()
    set(CPACK_GENERATOR "TGZ;STGZ;DEB")
    find_program(RPM_EXECUTABLE rpm)
    if (RPM_EXECUTABLE)
        set(CPACK_GENERATOR "${CPACK_GENERATOR};RPM")
    endif ()
endif ()

# comes after all cpack config/settings
include(CPack)

# set up tests
enable_testing()

### Test set internal unittest
    # aurpstd
    add_test(ut.aurostd.xscalar aflow --unit_test xscalar)
    add_test(ut.aurostd.xvector aflow --unit_test xvector)
    add_test(ut.aurostd.xmatrix aflow --unit_test xmatrix)
    add_test(ut.aurostd.xfit aflow --unit_test xfit)
    add_test(ut.aurostd.xfile aflow --unit_test xfile)
    add_test(ut.aurostd.aurostd_main aflow --unit_test aurostd_main)
    add_test(ut.aurostd.xparser aflow --unit_test xparser)
    # json serialization
    add_test(ut.serialization.json aflow --unit_test json_serialization)
    # database
    add_test(ut.database.schema aflow --unit_test schema)
    add_test(ut.database.entry_loader aflow --unit_test entry_loader)
    # aflowlib internal - requires internal database access and debian packages: gnuplot texlive-binaries imagemagick texlive-full
    add_test(ut.aflowlib.lib2raw aflow --unit_test lib2raw)
    # structure
    add_test(ut.structure.atomic_environment aflow --unit_test atomic_environment)
    add_test(ut.structure.xstructure aflow --unit_test xstructure)
    add_test(ut.structure.xstructure_parser aflow --unit_test xstructure_parser)
    add_test(ut.structure.xstructure_kpath aflow --unit_test xstructure_kpath)
        # structure_gen
    add_test(ut.structure_gen.ceramgen aflow --unit_test ceramgen)
    add_test(ut.structure_gen.proto aflow --unit_test proto)
    # surface
    add_test(ut.surface.triangle_area aflow --unit_test surface_TriangleArea)
    add_test(ut.surface.planeget_abcd aflow --unit_test surface_PlaneGetABCD)
    add_test(ut.surface.planedistance aflow --unit_test surface_PlaneDistance)
    add_test(ut.surface.planeget_proj aflow --unit_test surface_PlaneGetProjection)
    add_test(ut.surface.planeget_hkl aflow --unit_test surface_PlaneGetHKL)
    add_test(ut.surface.planeget_vertices aflow --unit_test surface_PlaneGetVVV)
    add_test(ut.surface.planeget_vertices_v2 aflow --unit_test surface_PlaneGetVVV_V2)
    add_test(ut.surface.planeget_atomicdensity aflow --unit_test surface_GetPlaneDensityAtoms)

### Test set chull
    # hull distance
    add_test(chull.distance.AlCo aflow --chull --alloy=AlCo --d2h=aflow:ef70a469857d3c47 --screen_only)
    set_property(TEST chull.distance.AlCo PROPERTY PASS_REGULAR_EXPRESSION "{\"aflow:ef70a469857d3c47\":242")
    add_test(chull.distance.MnPd aflow --chull --alloy=MnPd --d2h=aflow:9001c322296fc162 --screen_only)
    set_property(TEST chull.distance.MnPd PROPERTY PASS_REGULAR_EXPRESSION "{\"aflow:9001c322296fc162\":62")
    # stability criterion
    add_test(chull.stability.TeZr aflow --chull --alloy=TeZr --sc=aflow:570e80d4daa0a9d4 --screen_only)
    set_property(TEST chull.stability.TeZr PROPERTY PASS_REGULAR_EXPRESSION "{\"aflow:570e80d4daa0a9d4\":48")
    add_test(chull.stability.MnPd aflow --chull --alloy=MnPd --sc=aflow:4fc97e549ca17d4e --screen_only)
    set_property(TEST chull.stability.MnPd PROPERTY PASS_REGULAR_EXPRESSION "{\"aflow:4fc97e549ca17d4e\":13")
    add_test(chull.stability.MnPdPt aflow --chull --alloy=MnPdPt --sc=aflow:fb9eaa58604ce774 --screen_only)
    set_property(TEST chull.stability.MnPdPt PROPERTY PASS_REGULAR_EXPRESSION "{\"aflow:fb9eaa58604ce774\":32")
    # N+1 enthalpy gain
    add_test(chull.n+1_enthalpy.MnPdPt aflow --chull --alloy=MnPdPt --n+1_enthalpy_gain=aflow:fb9eaa58604ce774 --screen_only)
    set_property(TEST chull.n+1_enthalpy.MnPdPt PROPERTY PASS_REGULAR_EXPRESSION "{\"aflow:fb9eaa58604ce774\":32")
    # hull energy at stoichiometry
    add_test(chull.stoichiometry.MnPd.0.25 aflow --chull --alloy=MnPd --hull_energy=0.25 --screen_only)
    set_property(TEST chull.stoichiometry.MnPd.0.25 PROPERTY PASS_REGULAR_EXPRESSION "{\"hull_energy\\[0\\.2500,0\\.0000\\]\": -235")
    add_test(chull.stoichiometry.CuZr.0.5 aflow --chull --alloy=CuZr --hull_energy=0.5 --screen_only)
    set_property(TEST chull.stoichiometry.CuZr.0.5 PROPERTY PASS_REGULAR_EXPRESSION "{\"hull_energy\\[0\\.5000,0\\.0000\\]\": -158")

### Test set pocc
    # Mg_(x)Zn_(1-x)O and ZnS_(1-x)Se_(x)
    add_test(pocc.count.201:Mg:O:Zn  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=201:Mg:O:Zn --pocc_params=P0-1xB_P1-0.75xA-0.25xC | ${CMAKE_CURRENT_BINARY_DIR}/aflow --pocc_count_unique")
    set_property(TEST pocc.count.201:Mg:O:Zn PROPERTY PASS_REGULAR_EXPRESSION "COMPLETE.*Calculated 7 unique supercells")
    add_test(pocc.count.218:S:Se:Zn  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=218:S:Se:Zn --pocc_params=S0-1xC_S1-0.75xB-0.25xA | ${CMAKE_CURRENT_BINARY_DIR}/aflow --pocc_count_unique")
    set_property(TEST pocc.count.218:S:Se:Zn PROPERTY PASS_REGULAR_EXPRESSION "COMPLETE.*Calculated 7 unique supercells")
    add_test(pocc.file_count.AB_cF8_225_a_b.AB ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 2 --arguments="--aflow_proto=AB_cF8_225_a_b.AB:C:Hf:N:Zr --pocc_params=S0-1xA_S1-0.5xB-0.5xD,S0-0.5xA-0.5xC_S1-0.5xB-0.5xD --generate_aflowin_only")
    add_test(pocc.file_count.bcc,fcc:Cu:Fe ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 22 --arguments="--aflow_proto=bcc,fcc:Cu:Fe --pocc_params=P0-0.1xA-0.9xB,P0-0.2xA-0.8xB,P0-0.25xA-0.75xB,P0-0.333xA-0.667xB,P0-0.4xA-0.6xB,P0-0.5xA-0.5xB,P0-0.6xA-0.4xB,P0-0.667xA-0.333xB,P0-0.75xA-0.25xB,P0-0.8xA-0.2xB,P0-0.9xA-0.1xB --generate_aflowin_only")
    # bcc Fe_(x)Cu_(1-x)
    add_test(pocc.file_count.bcc.Fe0.1Cu0.9 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 19 --arguments="--aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.1xA-0.9xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.bcc.Fe0.2Cu0.8 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 6 --arguments="--aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.2xA-0.8xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.bcc.Fe0.25Cu0.75 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 8 --arguments="--aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.25xA-0.75xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.bcc.Fe0.333Cu0.667 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 4 --arguments="--aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.333xA-0.667xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.bcc.Fe0.4Cu0.6 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 10 --arguments="--aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.4xA-0.6xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.bcc.Fe0.5Cu0.5 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 3 --arguments="--aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.5xA-0.5xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.bcc.Fe0.6Cu0.4 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 10 --arguments="--aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.6xA-0.4xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.bcc.Fe0.667Cu0.333 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 4 --arguments="--aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.667xA-0.333xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.bcc.Fe0.75Cu0.25 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 8 --arguments="--aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.75xA-0.25xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.bcc.Fe0.8Cu0.2 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 6 --arguments="--aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.8xA-0.2xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.bcc.Fe0.9Cu0.1 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 19 --arguments="--aflow_proto=bcc:Cu:Fe --pocc_params=P0-0.9xA-0.1xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    # fcc Fe_(x)_Cu_(1-x)
    add_test(pocc.file_count.fcc.Fe0.1Cu0.9 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 19 --arguments="--aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.1xA-0.9xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.fcc.Fe0.2Cu0.8 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 6 --arguments="--aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.2xA-0.8xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.fcc.Fe0.25Cu0.75 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 8 --arguments="--aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.25xA-0.75xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.fcc.Fe0.333Cu0.667 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 4 --arguments="--aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.333xA-0.667xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.fcc.Fe0.4Cu0.6 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 10 --arguments="--aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.4xA-0.6xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.fcc.Fe0.5Cu0.5 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 3 --arguments="--aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.5xA-0.5xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.fcc.Fe0.6Cu0.4 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 10 --arguments="--aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.6xA-0.4xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.fcc.Fe0.667Cu0.333 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 4 --arguments="--aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.667xA-0.333xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.fcc.Fe0.75Cu0.25 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 8 --arguments="--aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.75xA-0.25xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.fcc.Fe0.8Cu0.2 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 6 --arguments="--aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.8xA-0.2xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    add_test(pocc.file_count.fcc.Fe0.9Cu0.1 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 19 --arguments="--aflow_proto=fcc:Cu:Fe --pocc_params=P0-0.9xA-0.1xB --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    # HNF optimization
    add_test(pocc.HNF.cementite  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB3_oP16_62_c_cd-001.AB:C:Mo:Nb:Ta:V:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ${CMAKE_CURRENT_BINARY_DIR}/aflow --hnf")
    set_property(TEST pocc.HNF.cementite PROPERTY PASS_REGULAR_EXPRESSION "Optimized HNF value = 5")
    add_test(pocc.HNF.C2HfNb  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb --pocc_params=S0-1xA_S1-0.5xB-0.5xC | ${CMAKE_CURRENT_BINARY_DIR}/aflow --hnf")
    set_property(TEST pocc.HNF.C2HfNb PROPERTY PASS_REGULAR_EXPRESSION "Optimized HNF value = 2")
    # 5-metal carbides
    add_test(pocc.file_count.5metal_carbide ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 56 --arguments="--aflow_proto=AB_cF8_225_a_b.AB:C:Hf,Mo,Nb,Ta,Ti,V,W,Zr:Hf,Mo,Nb,Ta,Ti,V,W,Zr:Hf,Mo,Nb,Ta,Ti,V,W,Zr:Hf,Mo,Nb,Ta,Ti,V,W,Zr:Hf,Mo,Nb,Ta,Ti,V,W,Zr --pocc_params=P0-1xA_P1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF --generate_aflowin_only")
    add_test(pocc.count.carbide.C5HfMoNbTaTi  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Mo:Nb:Ta:Ti --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ${CMAKE_CURRENT_BINARY_DIR}/aflow --pocc_count_unique")
    set_property(TEST pocc.count.carbide.C5HfMoNbTaTi PROPERTY PASS_REGULAR_EXPRESSION "COMPLETE.*Calculated 49 unique supercells")
    add_test(pocc.count.carbide.C5HfNbTiWZr  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Nb:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ${CMAKE_CURRENT_BINARY_DIR}/aflow --pocc_count_unique")
    set_property(TEST pocc.count.carbide.C5HfNbTiWZr PROPERTY PASS_REGULAR_EXPRESSION "COMPLETE.*Calculated 49 unique supercells")
    add_test(pocc.count.carbide.C5HfTaTiWZr  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB_cF8_225_a_b.AB:C:Hf:Ta:Ti:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ${CMAKE_CURRENT_BINARY_DIR}/aflow --pocc_count_unique")
    set_property(TEST pocc.count.carbide.C5HfTaTiWZr PROPERTY PASS_REGULAR_EXPRESSION "COMPLETE.*Calculated 49 unique supercells")
    add_test(pocc.count.carbide.C5MoNbTaTiW  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:Ti:W --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ${CMAKE_CURRENT_BINARY_DIR}/aflow --pocc_count_unique")
    set_property(TEST pocc.count.carbide.C5MoNbTaTiW PROPERTY PASS_REGULAR_EXPRESSION "COMPLETE.*Calculated 49 unique supercells")
    add_test(pocc.count.carbide.C5MoNbTaWZr  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:W:Zr --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ${CMAKE_CURRENT_BINARY_DIR}/aflow --pocc_count_unique")
    set_property(TEST pocc.count.carbide.C5MoNbTaWZr PROPERTY PASS_REGULAR_EXPRESSION "COMPLETE.*Calculated 49 unique supercells")
    add_test(pocc.file_count.carbide.C5MoNbTaVW ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_files --executable ${CMAKE_BINARY_DIR}/aflow
            --expected 50 --arguments="--aflow_proto=AB_cF8_225_a_b.AB:C:Mo:Nb:Ta:V:W --pocc_params=P0-1xA_P1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF --generate_aflowin_only" --arguments2="--run --generate_aflowin_only")
    # 5-metal borides
    add_test(pocc.count.boride.B10CrHfMoNbTa  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB2_hP3_191_a_d-001.BA:B:Cr:Hf:Mo:Nb:Ta --pocc_params=S0-1xA_S1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF | ${CMAKE_CURRENT_BINARY_DIR}/aflow --pocc_count_unique")
    set_property(TEST pocc.count.boride.B10CrHfMoNbTa PROPERTY PASS_REGULAR_EXPRESSION "COMPLETE.*Calculated 84 unique supercells")
    # 3-metal carbonitrides
    add_test(pocc.count.carbonitride.C3N3Hf2Ta2Ti2  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB_cF8_225_a_b.AB:C:Hf:N:Ta:Ti --pocc_params=S0-0.5xA-0.5xC_S1-0.333xB-0.333xD-0.333xE | ${CMAKE_CURRENT_BINARY_DIR}/aflow --pocc_count_unique")
    set_property(TEST pocc.count.carbonitride.C3N3Hf2Ta2Ti2 PROPERTY PASS_REGULAR_EXPRESSION "COMPLETE.*Calculated 1355 unique supercells")
    add_test(pocc.count.carbonitride.C3N3Ta2Ti2Zr2  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB_cF8_225_a_b.AB:C:N:Ta:Ti:Zr --pocc_params=S0-0.5xA-0.5xB_S1-0.333xC-0.333xD-0.333xE | ${CMAKE_CURRENT_BINARY_DIR}/aflow --pocc_count_unique")
    set_property(TEST pocc.count.carbonitride.C3N3Ta2Ti2Zr2 PROPERTY PASS_REGULAR_EXPRESSION "COMPLETE.*Calculated 1355 unique supercells")
    add_test(pocc.count.carbonitride.C3N3Nb2Ta2Ti2  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB_cF8_225_a_b.AB:C:N:Nb:Ta:Ti --pocc_params=S0-0.5xA-0.5xB_S1-0.333xC-0.333xD-0.333xE | ${CMAKE_CURRENT_BINARY_DIR}/aflow --pocc_count_unique")
    set_property(TEST pocc.count.carbonitride.C3N3Nb2Ta2Ti2 PROPERTY PASS_REGULAR_EXPRESSION "COMPLETE.*Calculated 1355 unique supercells")
    add_test(pocc.count.carbonitride.C3N3Hf2Nb2Zr2  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB_cF8_225_a_b.AB:C:Hf:N:Nb:Zr --pocc_params=S0-0.5xA-0.5xC_S1-0.333xB-0.333xD-0.333xE | ${CMAKE_CURRENT_BINARY_DIR}/aflow --pocc_count_unique")
    set_property(TEST pocc.count.carbonitride.C3N3Hf2Nb2Zr2 PROPERTY PASS_REGULAR_EXPRESSION "COMPLETE.*Calculated 1355 unique supercells")
    # 5-metal perovskite oxides
    add_test(pocc.count.perovskite.Ba5HfNbSnTiZrO15  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB3C_cP5_221_a_c_b.BAC:Ba:Hf:Nb:O:Sn:Ti:Zr --pocc_params=S0-1xD_S1-1xA_S2-0.2xB-0.2xC-0.2xE-0.2xF-0.2xG | ${CMAKE_CURRENT_BINARY_DIR}/aflow --pocc_count_unique")
    set_property(TEST pocc.count.perovskite.Ba5HfNbSnTiZrO15 PROPERTY PASS_REGULAR_EXPRESSION "COMPLETE.*Calculated 49 unique supercells")
    add_test(pocc.count.perovskite.Ba5HfSnTiYZrO15  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB3C_cP5_221_a_c_b.BAC:Ba:Hf:O:Sn:Ti:Y:Zr --pocc_params=S0-1xC_S1-1xA_S2-0.2xB-0.2xD-0.2xE-0.2xF-0.2xG | ${CMAKE_CURRENT_BINARY_DIR}/aflow --pocc_count_unique")
    set_property(TEST pocc.count.perovskite.Ba5HfSnTiYZrO15 PROPERTY PASS_REGULAR_EXPRESSION "COMPLETE.*Calculated 49 unique supercells")

### Test set SG
    add_test(sg.AB_cP2_221_b_a  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB_cP2_221_b_a:Cl:Cs  | ${CMAKE_CURRENT_BINARY_DIR}/aflow --aflowSG")
    set_property(TEST sg.AB_cP2_221_b_a PROPERTY PASS_REGULAR_EXPRESSION "Pm-3m #221")

### Test set AFLOWSYM
    add_test(sym.ABC2_mP8_10_ac_eh_mn-001 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py count_words --executable ${CMAKE_BINARY_DIR}/aflow
            --expected "pgroup:8,pgroupk:8,fgroup:4,pgroup_xtal:4,pgroupk_xtal:4,agroup:24,sgroup:108" --shell_string="${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=ABC2_mP8_10_ac_eh_mn-001:Au:Ag:Te | ${CMAKE_CURRENT_BINARY_DIR}/aflow --aflowSYM --screen_only")

### Test set edata
    add_test(edata.AB2_tI6_139_a_e  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB2_tI6_139_a_e:Os:Rh --params=2.7095,4.27466,0.662829 | ${CMAKE_CURRENT_BINARY_DIR}/aflow --edata")
    set_property(TEST edata.AB2_tI6_139_a_e PROPERTY PASS_REGULAR_EXPRESSION "Real space Bravais Lattice Primitive.*=.*BCT")

### Test set xtal finder
    # ideal case
    add_test(xtal_finder.ideal.A2B3_hR10_167_c_e  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=A2B3_hR10_167_c_e:Al:O --params=5.22337,2.72958,0.64784,0.9439 | ${CMAKE_CURRENT_BINARY_DIR}/aflow --prototype")
    set_property(TEST xtal_finder.ideal.A2B3_hR10_167_c_e PROPERTY PASS_REGULAR_EXPRESSION "AFLOW label.*A2B3_hR10_167_c_e")
    add_test(xtal_finder.ideal.AB18C8_cF108_225_a_eh_f  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB18C8_cF108_225_a_eh_f --params=10.3439,0.325,0.65833,0.16 | ${CMAKE_CURRENT_BINARY_DIR}/aflow --prototype")
    set_property(TEST xtal_finder.ideal.AB18C8_cF108_225_a_eh_f PROPERTY PASS_REGULAR_EXPRESSION "AFLOW label.*AB18C8_cF108_225_a_eh_f")
    add_test(xtal_finder.ideal.ABC4_mP12_13_e_a_2g  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=ABC4_mP12_13_e_a_2g --params=8.6244,0.500336,1.63361,145.35,0.5182,0.2986,0.0278,0.0003,0.2821,0.4045,0.2366 | ${CMAKE_CURRENT_BINARY_DIR}/aflow --prototype")
    set_property(TEST xtal_finder.ideal.ABC4_mP12_13_e_a_2g PROPERTY PASS_REGULAR_EXPRESSION "AFLOW label.*ABC4_mP12_13_e_a_2g")
    # compare_materials
    add_test(xtal_finder.compare.AB_cF8_225_a_b ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py temp_folder_run
            --shell_string="${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB_cF8_225_a_b:Cl:Na > AB_cF8_225_a_b.poscar\; cat AB_cF8_225_a_b.poscar | ${CMAKE_CURRENT_BINARY_DIR}/aflow --supercell=4,1,1 > AB_cF8_225_a_b_411_supercell.poscar\; ${CMAKE_CURRENT_BINARY_DIR}/aflow --compare_materials=AB_cF8_225_a_b.poscar,AB_cF8_225_a_b_411_supercell.poscar --quiet")
    set_property(TEST xtal_finder.compare.AB_cF8_225_a_b PROPERTY PASS_REGULAR_EXPRESSION ": MATCH")
    add_test(xtal_finder.compare.AB_cF8_225_a_b.reverse ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py temp_folder_run
            --shell_string="${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB_cF8_225_a_b:Cl:Na > AB_cF8_225_a_b.poscar\; cat AB_cF8_225_a_b.poscar | ${CMAKE_CURRENT_BINARY_DIR}/aflow --supercell=4,1,1 > AB_cF8_225_a_b_411_supercell.poscar\; ${CMAKE_CURRENT_BINARY_DIR}/aflow --compare_materials=AB_cF8_225_a_b_411_supercell.poscar,AB_cF8_225_a_b.poscar --quiet")
    set_property(TEST xtal_finder.compare.AB_cF8_225_a_b.reverse PROPERTY PASS_REGULAR_EXPRESSION ": MATCH")
    # atom decorations
    add_test(xtal_finder.decoration.ABCD_cF16_216_c_d_b_a  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=ABCD_cF16_216_c_d_b_a | ${CMAKE_CURRENT_BINARY_DIR}/aflow --unique_atom_decorations --print_misfit --quiet")
    set_property(TEST xtal_finder.decoration.ABCD_cF16_216_c_d_b_a PROPERTY PASS_REGULAR_EXPRESSION "Unique atom decorations.*3")
    add_test(xtal_finder.decoration.ABC_hP3_156_a_c_b  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=ABC_hP3_156_a_c_b --params=4.346,1.57271,4.44089e-16,0.4227,0.7214  | ${CMAKE_CURRENT_BINARY_DIR}/aflow --unique_atom_decorations --print_misfit --quiet")
    set_property(TEST xtal_finder.decoration.ABC_hP3_156_a_c_b PROPERTY PASS_REGULAR_EXPRESSION "Unique atom decorations.*3")
    # compare2proto
    add_test(xtal_finder.compare2proto.A2B3_hR10_167_c_e  sh -c "${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=A2B3_hR10_167_c_e:Al:O --params=-1,2.72957758313,0.36,0.5561 | ${CMAKE_CURRENT_BINARY_DIR}/aflow --compare2protos --print=txt --catalog=anrl --quiet")
    set_property(TEST xtal_finder.compare2proto.A2B3_hR10_167_c_e PROPERTY PASS_REGULAR_EXPRESSION "duplicate_compounds=4")
    # magnetic configurations
    add_test(xtal_finder.compare.magnetic.antiferro.antiferro ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py temp_folder_run
            --shell_string="${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=A_cI2_229_a:Cr | ${CMAKE_CURRENT_BINARY_DIR}/aflow --sconv >A_cI2_229_a_sconv.poscar\; ${CMAKE_CURRENT_BINARY_DIR}/aflow --compare_materials=A_cI2_229_a_sconv.poscar,A_cI2_229_a_sconv.poscar --magmoms=-1,1:1,-1 --quiet")
    set_property(TEST xtal_finder.compare.magnetic.antiferro.antiferro PROPERTY PASS_REGULAR_EXPRESSION ": MATCH")
    add_test(xtal_finder.compare.magnetic.ferro.antiferro ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py temp_folder_run
            --shell_string="${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=A_cI2_229_a:Cr | ${CMAKE_CURRENT_BINARY_DIR}/aflow --sconv >A_cI2_229_a_sconv.poscar\; ${CMAKE_CURRENT_BINARY_DIR}/aflow --compare_materials=A_cI2_229_a_sconv.poscar,A_cI2_229_a_sconv.poscar --magmoms=1,1:1,-1 --quiet")
    set_property(TEST xtal_finder.compare.magnetic.ferro.antiferro PROPERTY PASS_REGULAR_EXPRESSION "UNMATCHABLE")
    add_test(xtal_finder.compare.magnetic.ferro_up.ferro_down ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py temp_folder_run
            --shell_string="${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=A_cI2_229_a:Cr | ${CMAKE_CURRENT_BINARY_DIR}/aflow --sconv >A_cI2_229_a_sconv.poscar\; ${CMAKE_CURRENT_BINARY_DIR}/aflow --compare_materials=A_cI2_229_a_sconv.poscar,A_cI2_229_a_sconv.poscar --magmoms=1,1:-1,-1 --quiet")
    set_property(TEST xtal_finder.compare.magnetic.ferro_up.ferro_down PROPERTY PASS_REGULAR_EXPRESSION "UNMATCHABLE")

### Test set CCE
    # prototype based test
    add_test(cce.proto.MgO ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py temp_folder_run
            --shell_string="${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB_cF8_225_a_b:Mg:O --params=4.253 > POSCAR_proto\; ${CMAKE_CURRENT_BINARY_DIR}/aflow --cce=POSCAR_proto --cce_test --dft_formation_energies=-5.434,-6.220,-6.249 --functionals=PBE,LDA,SCAN --oxidation_numbers=2,-2")
    set_property(TEST cce.proto.MgO PROPERTY PASS_REGULAR_EXPRESSION "0.801,0.763,0.015,-0.025,-0.014,-0.054,-6.235,-6.197,-6.235,-6.195,-6.235,-6.195")
    add_test(cce.proto.Mg2SiO4 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py temp_folder_run
            --shell_string="${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=A2B4C_oP28_62_ac_2cd_c:Mg:O:Si --params=10.193,0.586382811734,0.466202295693,0.2774,-0.0085,0.0913,0.7657,0.4474,0.2215,0.094,0.4262,0.1628,0.0331,0.2777 > POSCAR_proto\; ${CMAKE_CURRENT_BINARY_DIR}/aflow --cce=POSCAR_proto --cce_test --dft_formation_energies=-79.46,-90.112,-90.444 --functionals=PBE,LDA,SCAN --oxidation_numbers=2,2,2,2,2,2,2,2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,4,4,4,4")
    set_property(TEST cce.proto.Mg2SiO4 PROPERTY PASS_REGULAR_EXPRESSION "10.456,9.912,-0.252,-0.824,-0.444,-1.020,-89.916,-89.372,-89.860,-89.288,-90.000,-89.424")
    add_test(cce.proto.BaTiO3 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py temp_folder_run
            --shell_string="${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB3C_cP5_221_a_c_b:Ba:O:Ti --params=4.033 > POSCAR_proto\; ${CMAKE_CURRENT_BINARY_DIR}/aflow --cce=POSCAR_proto --cce_test --dft_formation_energies=-15.450,-17.641,-17.702 --functionals=PBE,LDA,SCAN")
    set_property(TEST cce.proto.BaTiO3 PROPERTY PASS_REGULAR_EXPRESSION "2.063,1.981,-0.461,-0.553,-0.692,-0.783,-17.513,-17.431,-17.180,-17.088,-17.010,-16.919")
    add_test(cce.proto.Na2WO4 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py temp_folder_run
            --shell_string="${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=A2BC4_cF56_227_d_a_e:Na:W:O --params=9.0433,0.260612 > POSCAR_proto\; ${CMAKE_CURRENT_BINARY_DIR}/aflow --cce=POSCAR_proto --cce_test --dft_formation_energies=-30.456,-33.592,-33.872 --functionals=PBE,LDA,SCAN --oxidation_numbers=1,1,1,1,-2,-2,-2,-2,-2,-2,-2,-2,6,6")
    set_property(TEST cce.proto.Na2WO4 PROPERTY PASS_REGULAR_EXPRESSION "2.014,1.843,-1.765,-1.963,-1.359,-1.561,-32.470,-32.299,-31.827,-31.629,-32.513,-32.311")
    add_test(cce.proto.Fe3O4 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py temp_folder_run
            --shell_string="${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=A3B4_cF56_227_ad_e:Fe:O --params=8.3886,0.745448 > POSCAR_proto\; ${CMAKE_CURRENT_BINARY_DIR}/aflow --cce=POSCAR_proto --cce_test --dft_formation_energies=-18.194,-21.236,-25.404 --functionals=PBE,LDA,SCAN")
    set_property(TEST cce.proto.Fe3O4 PROPERTY PASS_REGULAR_EXPRESSION "5.381,5.246,1.832,1.654,-1.084,-1.263,-23.575,-23.440,-23.068,-22.890,-24.320,-24.141")
    add_test(cce.proto.Li2O ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py temp_folder_run
            --shell_string="${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB2_cF12_225_a_c:O:Li --params=4.5709 > POSCAR_proto\; ${CMAKE_CURRENT_BINARY_DIR}/aflow --cce_test --cce=POSCAR_proto --dft_formation_energies=-5.584,-6.320,-6.291 --functionals=PBE,LDA,SCAN")
    set_property(TEST cce.proto.Li2O PROPERTY PASS_REGULAR_EXPRESSION "0.613,0.563,-0.123,-0.178,-0.094,-0.149,-6.197,-6.147,-6.197,-6.142,-6.197,-6.142")
    add_test(cce.proto.BaTiO3 ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py temp_folder_run
            --shell_string="${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB3C_hR10_148_c_f_c:Fe:O:Ti --params=5.08264,2.77119,0.648852,0.852345,0.436666,0.768665,0.0556432 > POSCAR_proto_Ba1O3Ti1_ICSD_187292 > POSCAR_proto\; ${CMAKE_CURRENT_BINARY_DIR}/aflow --cce=POSCAR_proto --cce_test --dft_formation_energies=-22.384,-25.198,-27.066 --functionals=PBE,LDA,SCAN")
    set_property(TEST cce.proto.BaTiO3 PROPERTY PASS_REGULAR_EXPRESSION "3.402,3.240,0.414,0.258,-1.258,-1.442,-25.786,-25.624,-25.612,-25.456,-25.808,-25.624")
    add_test(cce.proto.HfN ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/external_unit_tests.py temp_folder_run
            --shell_string="${CMAKE_CURRENT_BINARY_DIR}/aflow --proto=AB_cF8_225_a_b:Hf:N --params=4.4938 > POSCAR_proto_Ba1O3Ti1_ICSD_187292 > POSCAR_proto\; ${CMAKE_CURRENT_BINARY_DIR}/aflow --cce=POSCAR_proto --cce_test --dft_formation_energies=-3.513,-4.429,-3.958,-3.516,0 --functionals=PBE,LDA,SCAN,PBE+U:ICSD,exp")
    set_property(TEST cce.proto.HfN PROPERTY PASS_REGULAR_EXPRESSION "0.359,0.338,-0.557,-0.581,-0.086,-0.109,0.356,0.335,-3.872,-3.851,-3.872,-3.848,-3.872,-3.849,-3.872,-3.851,-3.872")

# Label test sets
set_property(TEST ut.aurostd.xscalar ut.aurostd.xvector ut.aurostd.xmatrix ut.aurostd.xfit ut.aurostd.aurostd_main
        ut.aurostd.xparser ut.database.schema ut.structure_gen.ceramgen pocc.file_count.AB_cF8_225_a_b.AB
        sym.ABC2_mP8_10_ac_eh_mn-001 ut.structure.atomic_environment
        PROPERTY LABELS quick)
