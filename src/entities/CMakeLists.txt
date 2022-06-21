add_library(matrix
        matrix/Matrix.cpp matrix/Matrix.h
        matrix/Submatrix.cpp matrix/Submatrix.h)

add_library(spectrum
        spectrum/Spectrum.cpp spectrum/Spectrum.h spectrum/Subspectrum.cpp spectrum/Subspectrum.h)

add_subdirectory(data_structures)

add_library(magnetic_susceptibility
        magnetic_susceptibility/MuSquaredWorker.cpp
        magnetic_susceptibility/MuSquaredWorker.h
        magnetic_susceptibility/EnsembleAverager.cpp
        magnetic_susceptibility/EnsembleAverager.h
        magnetic_susceptibility/ExperimentalValuesWorker.cpp
        magnetic_susceptibility/ExperimentalValuesWorker.h
        magnetic_susceptibility/UniqueGOnlySSquaredMuSquaredWorker.cpp
        magnetic_susceptibility/UniqueGOnlySSquaredMuSquaredWorker.h)