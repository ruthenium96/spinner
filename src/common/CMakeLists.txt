add_library(common
        lexicographic/IndexConverter.cpp lexicographic/IndexConverter.h
        lexicographic/SparseMatrix.h
        lexicographic/SparseMatrix-arma-sp-mat.cpp
        symbols/Symbols.h
        symbols/Symbols.cpp
        symbols/SymbolName.cpp
        symbols/SymbolName.h
        Logger.h
        runner/Runner.cpp
        model/Model.cpp
        model/Model.h
        Quantity.cpp)

target_link_libraries(common components magnetic_susceptibility stlbfgs)