add_library(model
        Model.h Model.cpp
        operators/Operator.h operators/Operator.cpp
        operators/terms/ConstantTerm.h operators/terms/ConstantTerm.cpp
        operators/terms/ScalarProductTerm.h operators/terms/ScalarProductTerm.cpp
        operators/terms/ZProjectionsProductTerm.cpp operators/terms/ZProjectionsProductTerm.h
        operators/terms/Term.h
        symbols/Symbols.h symbols/Symbols.cpp
        symbols/SymbolName.h symbols/SymbolName.cpp
        )

target_link_libraries(model group)