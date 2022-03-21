#ifndef JULY_SPECTRUM_H
#define JULY_SPECTRUM_H

#include <vector>

#include "Subspectrum.h"
#include "src/entities/matrix/Matrix.h"

struct Spectrum {
    Spectrum() = default;

    std::vector<Subspectrum> blocks;
    explicit Spectrum(std::vector<Subspectrum>&& m);

    static Spectrum energy(
        const Matrix& hamiltonian_matrix,
        std::vector<DenseMatrix>& unitary_transformation_matrices);
    static Spectrum non_energy(
        const Matrix& non_hamiltonian_matrix,
        const std::vector<DenseMatrix>& unitary_transformation_matrices);
};

std::ostream& operator<<(std::ostream& os, const Spectrum& space);

#endif  //JULY_SPECTRUM_H
