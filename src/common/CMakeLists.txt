add_library(common
        lexicographic/IndexConverter.cpp lexicographic/IndexConverter.h
        lexicographic/SparseMatrix.h
        lexicographic/SparseMatrix-arma-sp-mat.cpp
        Logger.h
        runner/Runner.cpp
        Quantity.cpp
        physical_optimization/OptimizationList.cpp physical_optimization/OptimizationList.h
        runner/ConsistentModelOptimizationList.h runner/ConsistentModelOptimizationList.cpp
        runner/ModelOptimizationListConsistence.h runner/ModelOptimizationListConsistence.cpp)

target_link_libraries(common model space magnetic_susceptibility stlbfgs)