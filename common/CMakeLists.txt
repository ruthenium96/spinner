add_library(common
        LexicographicIndexConverter.cpp LexicographicIndexConverter.h
        LexicographicSparseMatrix.h
        LexicographicSparseMatrix-arma-sp-mat.cpp
        Logger.h
        Runner.cpp
        Quantity.cpp)

target_link_libraries(common components)