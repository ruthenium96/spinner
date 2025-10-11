add_library(common
        runner/Runner.cpp
        runner/Executer.cpp runner/Executer.h
        PrintingFunctions.cpp PrintingFunctions.h
        UncertainValue.cpp UncertainValue.h
        physical_optimization/OptimizationList.cpp physical_optimization/OptimizationList.h
        runner/ConsistentModelOptimizationList.h runner/ConsistentModelOptimizationList.cpp
        MumxHash.h)

add_library(logger Logger.h)
set_target_properties(logger PROPERTIES LINKER_LANGUAGE CXX)

add_library(index_converter
        index_converter/lexicographic/IndexConverter.cpp index_converter/lexicographic/IndexConverter.h
        index_converter/lexicographic/IndexPermutator.cpp index_converter/lexicographic/IndexPermutator.h
        index_converter/AbstractIndexConverter.cpp index_converter/AbstractIndexConverter.h
        index_converter/s_squared/IndexConverter.cpp index_converter/s_squared/IndexConverter.h
        index_converter/s_squared/IndexPermutator.cpp index_converter/s_squared/IndexPermutator.h
        index_converter/s_squared/Level.cpp index_converter/s_squared/Level.h
        index_converter/s_squared/OrderOfSummation.cpp index_converter/s_squared/OrderOfSummation.h)

find_package(spdlog)
Include(FetchContent)
FetchContent_Declare(
   spdlog
   GIT_REPOSITORY https://github.com/gabime/spdlog.git
   GIT_TAG v1.9.2
)
FetchContent_MakeAvailable(spdlog)

target_link_libraries(logger spdlog::spdlog magic_enum)

target_link_libraries(index_converter logger)

target_link_libraries(common index_converter model space magnetic_susceptibility nonlinear_solver eigendecompositor
        logger magic_enum)
