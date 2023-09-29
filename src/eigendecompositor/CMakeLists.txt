add_library(eigendecompositor
        ExactEigendecompositor.cpp
        )

target_link_libraries(eigendecompositor model space matrix spectrum)
