#ifndef JULY_SPECTRUMBUILDER_H
#define JULY_SPECTRUMBUILDER_H

#include "entities/matrix/Matrix.h"
#include "entities/spectrum/Spectrum.h"

class SpectrumBuilder {
public:
    Spectrum apply_to_energy(const Matrix& hamiltonian_matrix);
    Spectrum apply_to_non_energy(const Matrix& non_hamiltonian_matrix);
    Subspectrum apply_to_subentity_energy(const Submatrix& hamiltonian_submatrix, DenseMatrix& unitary_transformation_matrix);
    Subspectrum apply_to_subentity_non_energy(const Submatrix& non_hamiltonian_submatrix, const DenseMatrix& unitary_transformation_matrix);
private:
    std::vector<DenseMatrix> unitary_transformation_matrices;
};


#endif //JULY_SPECTRUMBUILDER_H
