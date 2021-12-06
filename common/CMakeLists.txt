add_library(common
        lexicographic/IndexConverter.cpp lexicographic/IndexConverter.h
        lexicographic/SparseMatrix.h
        lexicographic/SparseMatrix-arma-sp-mat.cpp
        symbols/Symbols.h
        symbols/Symbols.cpp
        Logger.h
        runner/Runner.cpp
        Quantity.cpp)

target_link_libraries(common components magnetic_susceptibility stlbfgs)