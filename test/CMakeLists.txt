add_executable(julytest
        ${CMAKE_SOURCE_DIR}/test/main.cpp 
        ${CMAKE_SOURCE_DIR}/test/PerformanceTest.cpp 
        ${CMAKE_SOURCE_DIR}/components/c2_symmetrizer.cpp 
        ${CMAKE_SOURCE_DIR}/common/LexicographicIndexWorker.cpp
        ${CMAKE_SOURCE_DIR}/common/LexicographicIndexWorker.h
        ${CMAKE_SOURCE_DIR}/components/tz_sorter.cpp 
        ${CMAKE_SOURCE_DIR}/entities/Space.cpp
        )

target_link_libraries(julytest gtest gtest_main)