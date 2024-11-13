add_library(common
        lexicographic/IndexConverter.cpp lexicographic/IndexConverter.h
        Logger.h
        runner/Runner.cpp
        runner/Executer.cpp runner/Executer.h
        PrintingFunctions.cpp PrintingFunctions.h
        physical_optimization/OptimizationList.cpp physical_optimization/OptimizationList.h
        runner/ConsistentModelOptimizationList.h runner/ConsistentModelOptimizationList.cpp
        MumxHash.h)

find_package(spdlog)
#todo: download if necessary
# Include(FetchContent)
#FetchContent_Declare(
#    spdlog
#    GIT_REPOSITORY https://github.com/gabime/spdlog.git
#    GIT_TAG v1.9.2
#)
#FetchContent_MakeAvailable(spdlog)


target_link_libraries(common model space magnetic_susceptibility nonlinear_solver eigendecompositor
        spdlog::spdlog magic_enum)
