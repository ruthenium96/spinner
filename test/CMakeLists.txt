add_executable(julytest
        ${CMAKE_SOURCE_DIR}/test/main.cpp 
        ${CMAKE_SOURCE_DIR}/test/PerformanceTest.cpp 
        ${CMAKE_SOURCE_DIR}/common/LexicographicIndexConverter.cpp
        ${CMAKE_SOURCE_DIR}/common/LexicographicIndexConverter.h
        ${CMAKE_SOURCE_DIR}/entities/space/Space.cpp
        ${CMAKE_SOURCE_DIR}/test/group_test.cpp
        ${CMAKE_SOURCE_DIR}/test/symmetrizer_tests.cpp
        ${CMAKE_SOURCE_DIR}/test/non_abelian_simplifier_tests.cpp
        ${CMAKE_SOURCE_DIR}/test/matrix_builder.cpp)

target_link_libraries(julytest gtest gtest_main groups components space common operator matrix)