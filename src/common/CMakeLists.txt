add_library(common
        lexicographic/IndexConverter.cpp lexicographic/IndexConverter.h
        Logger.h
        runner/Runner.cpp
        Quantity.cpp
        physical_optimization/OptimizationList.cpp physical_optimization/OptimizationList.h
        runner/ConsistentModelOptimizationList.h runner/ConsistentModelOptimizationList.cpp)

target_link_libraries(common model matrix spectrum space magnetic_susceptibility nonlinear_solver)