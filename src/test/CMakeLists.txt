
add_library(aflow_tests OBJECT
        aflow_test.cpp
        aflow_unit_test.cpp
)

target_include_directories(aflow_tests
        PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_BINARY_DIR}
)

target_link_libraries(aflow_tests PRIVATE config)