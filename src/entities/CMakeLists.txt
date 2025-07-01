add_library(matrix
        matrix/Matrix.cpp matrix/Matrix.h
        matrix/Submatrix.cpp matrix/Submatrix.h)

target_link_libraries(matrix model)

add_library(spectrum
        spectrum/Spectrum.cpp spectrum/Spectrum.h spectrum/Subspectrum.cpp spectrum/Subspectrum.h)

add_subdirectory(data_structures)

add_library(magnetic_susceptibility
        magnetic_susceptibility/MagneticSusceptibilityController.h
        magnetic_susceptibility/assistant/CommonEnsembleAverager.cpp
        magnetic_susceptibility/assistant/CommonEnsembleAverager.h
        magnetic_susceptibility/assistant/ExperimentalValuesWorker.cpp
        magnetic_susceptibility/assistant/ExperimentalValuesWorker.h
        magnetic_susceptibility/assistant/FTLMEnsembleAverager.cpp
        magnetic_susceptibility/assistant/FTLMEnsembleAverager.h
        magnetic_susceptibility/assistant/AbstractEnsembleAverager.cpp
        magnetic_susceptibility/assistant/AbstractEnsembleAverager.h
        magnetic_susceptibility/worker/UniqueGWorker.cpp
        magnetic_susceptibility/worker/UniqueGWorker.h
        magnetic_susceptibility/worker/CurieWeissWorker.cpp
        magnetic_susceptibility/worker/CurieWeissWorker.h
        magnetic_susceptibility/worker/BasicWorker.cpp
        magnetic_susceptibility/worker/BasicWorker.h
        magnetic_susceptibility/MagneticSusceptibilityController.cpp
        magnetic_susceptibility/worker/AbstractWorker.h
        magnetic_susceptibility/worker/DifferentGWorker.cpp
        magnetic_susceptibility/worker/DifferentGWorker.h
        magnetic_susceptibility/assistant/WeightingScheme.cpp
        magnetic_susceptibility/assistant/WeightingScheme.h
        magnetic_susceptibility/worker/WorkerConstructor.cpp 
        magnetic_susceptibility/worker/WorkerConstructor.h)

target_link_libraries(magnetic_susceptibility common)