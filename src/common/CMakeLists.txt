add_library(common
        lexicographic/IndexConverter.cpp lexicographic/IndexConverter.h
        Logger.h
        runner/Runner.cpp
        runner/Executer.cpp runner/Executer.h
        Quantity.cpp
        PrintingFunctions.cpp PrintingFunctions.h
        physical_optimization/OptimizationList.cpp physical_optimization/OptimizationList.h
        runner/ConsistentModelOptimizationList.h runner/ConsistentModelOptimizationList.cpp
        MumxHash.h)

target_link_libraries(common model space magnetic_susceptibility nonlinear_solver eigendecompositor)
