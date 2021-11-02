#ifndef JULY_SPECTRUMBUILDER_H
#define JULY_SPECTRUMBUILDER_H

#include "entities/matrix/Matrix.h"
#include "entities/spectrum/Spectrum.h"

class SpectrumBuilder {
public:
    Spectrum apply_to_energy(const Matrix& hamiltonian_matrix);
    Spectrum apply_to_non_energy(const Matrix& non_hamiltonian_matrix);
private:
    std::vector<DenseMatrix> unitary_transformation_matrices;
};


#endif //JULY_SPECTRUMBUILDER_H
