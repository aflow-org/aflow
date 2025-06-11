
add_library(aflow_tests OBJECT aflow_unit_test.cpp)

target_include_directories(aflow_tests PUBLIC ${CMAKE_SOURCE_DIR})
