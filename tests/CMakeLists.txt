mark_as_advanced(
        BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
        gmock_build_tests gtest_build_samples gtest_build_tests
        gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

set_target_properties(gtest PROPERTIES FOLDER libs)
set_target_properties(gtest_main PROPERTIES FOLDER libs)
set_target_properties(gmock PROPERTIES FOLDER libs)
set_target_properties(gmock_main PROPERTIES FOLDER libs)

add_executable(july_unittests
        main_unittest.cpp
        blm_unittest.cpp
        matrices_unittest.cpp
        )
target_link_libraries(july_unittests gtest gmock gtest_main blm wignerSymbols)

set_target_properties(july_unittests PROPERTIES
        FOLDER tests
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO)

target_include_directories(july_unittests PRIVATE "${PROJECT_SOURCE_DIR}")