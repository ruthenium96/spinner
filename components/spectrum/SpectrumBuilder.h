#ifndef JULY_SPECTRUMBUILDER_H
#define JULY_SPECTRUMBUILDER_H

#include "entities/matrix/Matrix.h"
#include "entities/spectrum/Spectrum.h"

class SpectrumBuilder {
public:
    Spectrum apply(const Matrix& hamiltonian_matrix, const std::vector<Matrix>& non_hamiltonian_matrices);
};


#endif //JULY_SPECTRUMBUILDER_H
