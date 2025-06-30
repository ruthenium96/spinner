add_library(eigendecompositor
        ExactEigendecompositor.cpp
        FTLMEigendecompositor.cpp FTLMEigendecompositor.h
        ImplicitQuantityEigendecompositor.cpp ImplicitQuantityEigendecompositor.h
        OneSymbolInHamiltonianEigendecompositor.cpp OneSymbolInHamiltonianEigendecompositor.h
        AbstractEigendecompositor.cpp
        FlattenedSpectra.cpp FlattenedSpectra.h
        ExplicitQuantitiesEigendecompositor.cpp ExplicitQuantitiesEigendecompositor.h
        EigendecompositorConstructor.cpp EigendecompositorConstructor.h)

target_link_libraries(eigendecompositor model space matrix spectrum)
