add_library(matrix
        matrix/Matrix.cpp matrix/Matrix.h
        matrix/Submatrix.cpp matrix/Submatrix.h)

add_library(spectrum
        spectrum/Spectrum.cpp spectrum/Spectrum.h spectrum/Subspectrum.cpp spectrum/Subspectrum.h)

add_subdirectory(data_structures)

add_library(magnetic_susceptibility
        magnetic_susceptibility/AbstractMuSquaredWorker.h
        magnetic_susceptibility/assistant/EnsembleAverager.cpp
        magnetic_susceptibility/assistant/EnsembleAverager.h
        magnetic_susceptibility/assistant/ExperimentalValuesWorker.cpp
        magnetic_susceptibility/assistant/ExperimentalValuesWorker.h
        magnetic_susceptibility/UniqueGOnlySSquaredMuSquaredWorker.cpp
        magnetic_susceptibility/UniqueGOnlySSquaredMuSquaredWorker.h
        magnetic_susceptibility/CurieWeissMuSquaredWorker.cpp
        magnetic_susceptibility/CurieWeissMuSquaredWorker.h magnetic_susceptibility/BasicMuSquaredWorker.cpp magnetic_susceptibility/BasicMuSquaredWorker.h magnetic_susceptibility/AbstractMuSquaredWorker.cpp)