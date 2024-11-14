add_library(model
        ModelInput.h ModelInput.cpp
        operators/Operator.h operators/Operator.cpp
        operators/terms/ConstantTerm.h operators/terms/ConstantTerm.cpp
        operators/terms/ScalarProductTerm.h operators/terms/ScalarProductTerm.cpp
        operators/terms/SzSzTwoCenterTerm.cpp operators/terms/SzSzTwoCenterTerm.h
        operators/terms/SzSzOneCenterTerm.cpp operators/terms/SzSzOneCenterTerm.h
        operators/terms/LocalSSquaredOneCenterTerm.cpp operators/terms/LocalSSquaredOneCenterTerm.h
        operators/terms/ScalarProductITOTerm.cpp operators/terms/ScalarProductITOTerm.h
        operators/terms/Term.h
        symbols/SymbolicWorker.h symbols/SymbolicWorker.cpp
        symbols/NumericalWorker.cpp symbols/NumericalWorker.h
        symbols/SymbolName.h symbols/SymbolName.cpp
        Model.cpp Model.h)

target_link_libraries(model group)