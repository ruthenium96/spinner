add_library(model
        Model.h Model.cpp
        operators/Operator.h operators/Operator.cpp
        operators/ConstantTerm.h operators/ConstantTerm.cpp
        operators/ScalarProductTerm.h operators/ScalarProductTerm.cpp
        operators/Term.h
        symbols/Symbols.h symbols/Symbols.cpp
        symbols/SymbolName.h symbols/SymbolName.cpp
        )

target_link_libraries(model group)