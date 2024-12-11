add_library(model
        ModelInput.h ModelInput.cpp
        operators/Operator.h operators/Operator.cpp
        operators/terms/ConstantTerm.h operators/terms/ConstantTerm.cpp
        operators/terms/lexicographic/ScalarProductTerm.h operators/terms/lexicographic/ScalarProductTerm.cpp
        operators/terms/lexicographic/SzSzTwoCenterTerm.cpp operators/terms/lexicographic/SzSzTwoCenterTerm.h
        operators/terms/lexicographic/SzSzOneCenterTerm.cpp operators/terms/lexicographic/SzSzOneCenterTerm.h
        operators/terms/LocalSSquaredOneCenterTerm.cpp operators/terms/LocalSSquaredOneCenterTerm.h
        operators/terms/ito/ScalarProductTerm.cpp operators/terms/ito/ScalarProductTerm.h
        operators/terms/ito/TotalSSquaredTerm.cpp operators/terms/ito/TotalSSquaredTerm.h
        operators/terms/Term.h
        symbols/SymbolicWorker.h symbols/SymbolicWorker.cpp
        symbols/NumericalWorker.cpp symbols/NumericalWorker.h
        symbols/SymbolName.h symbols/SymbolName.cpp
        Model.cpp Model.h)

target_link_libraries(model group index_converter)