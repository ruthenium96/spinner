add_library(model
        Model.h Model.cpp
        operators/Operator.h operators/Operator.cpp
        operators/terms/ConstantTerm.h operators/terms/ConstantTerm.cpp
        operators/terms/ScalarProductTerm.h operators/terms/ScalarProductTerm.cpp
        operators/terms/SzSzTwoCenterTerm.cpp operators/terms/SzSzTwoCenterTerm.h
        operators/terms/SzSzOneCenterTerm.cpp operators/terms/SzSzOneCenterTerm.h
        operators/terms/Term.h
        symbols/Symbols.h symbols/Symbols.cpp
        symbols/SymbolName.h symbols/SymbolName.cpp
        )

target_link_libraries(model group)