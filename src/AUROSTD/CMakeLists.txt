# set definitions

add_compile_definitions(AFLOW_VERSION="${AFLOW_VERSION}")
add_compile_definitions(AFLOW_VERSION_SHORT="${CMAKE_PROJECT_VERSION}")
add_compile_definitions(TODAY="${TODAY_TMP}")
add_compile_definitions(CMAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}/")

# building apack to avoid dependency on xz being installed
add_executable(apack apack.c)
target_link_libraries(apack PRIVATE ${LibArchive_LIBRARIES})

# building aurostd
add_library(aurostd OBJECT aurostd.h
                           aurostd_argv.cpp
                           aurostd_hash.cpp
                           aurostd_data.cpp
                           aurostd_time.cpp
                           aurostd_xfile.cpp
                           aurostd_main.cpp
                           aurostd_xcombos.cpp
                           aurostd_xcomplex.cpp
                           aurostd_xerror.cpp
                           aurostd_xfit.cpp
                           aurostd_xhttp.cpp
                           aurostd_xmatrix.cpp
                           aurostd_xoption.cpp
                           aurostd_xparser.cpp
                           aurostd_xparser_json.cpp
                           aurostd_xserialization.cpp
                           aurostd_xrandom.cpp
                           aurostd_xscalar.cpp
                           aurostd_xtensor.cpp
                           aurostd_xvector.cpp)

target_include_directories(aurostd PUBLIC ${CMAKE_SOURCE_DIR})

# set DATA folder
set(DATA_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/../DATA/)
set(DATA_GROUPS README FINDSYM FROZSL PHVASP TEST SCRIPTS IMAGES)

# commands to generate compressed/archived variants
foreach (DATA_GROUP ${DATA_GROUPS})
    # build a list of files to pack
    file(GLOB_RECURSE FILE_LIST_${DATA_GROUP} LIST_DIRECTORIES false RELATIVE ${DATA_FOLDER}/${DATA_GROUP}
            ${DATA_FOLDER}/${DATA_GROUP}/*)

    # remove hidden files from the list
    foreach( FILE ${FILE_LIST_${DATA_GROUP}})
        if (${FILE} MATCHES "^\\..*" OR ${FILE} MATCHES ".*/\\..*")
            list(REMOVE_ITEM FILE_LIST_TEST ${FILE})
        endif()
    endforeach()

    add_custom_command(
            OUTPUT ${DATA_FOLDER}/${DATA_GROUP}.tar.xz
            COMMAND apack ${DATA_FOLDER}/${DATA_GROUP}.tar.xz ${FILE_LIST_${DATA_GROUP}}
            WORKING_DIRECTORY ${DATA_FOLDER}/${DATA_GROUP}/
            DEPENDS apack
            USES_TERMINAL)
endforeach ()

add_custom_command(
        OUTPUT ${DATA_FOLDER}/PROTO.tar.xz
        COMMAND ${Python_EXECUTABLE} ${DATA_FOLDER}/PROTO/prepare_proto.py && ${CMAKE_BINARY_DIR}/AUROSTD/apack ${DATA_FOLDER}/PROTO.tar.xz lookup.json data.json && rm ${DATA_FOLDER}/PROTO/data.json ${DATA_FOLDER}/PROTO/lookup.json
        DEPENDS apack
        WORKING_DIRECTORY ${DATA_FOLDER}/PROTO
        USES_TERMINAL)

# set dependency
foreach (DATA_GROUP ${DATA_GROUPS})
    set_property(SOURCE aurostd_data.cpp APPEND PROPERTY OBJECT_DEPENDS ${DATA_FOLDER}/${DATA_GROUP}.tar.xz)
endforeach ()
set_property(SOURCE aurostd_data.cpp APPEND PROPERTY OBJECT_DEPENDS ${DATA_FOLDER}/PROTO.tar.xz)
