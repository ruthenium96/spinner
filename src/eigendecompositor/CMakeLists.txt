add_library(eigendecompositor
        ExactEigendecompositor.cpp
        ImplicitSSquareEigendecompositor.cpp ImplicitSSquareEigendecompositor.h)

target_link_libraries(eigendecompositor model space matrix spectrum)
