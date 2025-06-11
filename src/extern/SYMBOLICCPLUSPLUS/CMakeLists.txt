# set definitions

add_compile_definitions(AFLOW_VERSION="${AFLOW_VERSION}")
add_compile_definitions(AFLOW_VERSION_SHORT="${CMAKE_PROJECT_VERSION}")
add_compile_definitions(TODAY="${TODAY_TMP}")

add_library(symbolic OBJECT symbolic_main.cpp)

