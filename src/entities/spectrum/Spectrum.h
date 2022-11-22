#ifndef SPINNER_SPECTRUM_H
#define SPINNER_SPECTRUM_H

#include <vector>

#include "Subspectrum.h"
#include "src/entities/matrix/Matrix.h"

struct Spectrum {
    Spectrum() = default;

    std::vector<Subspectrum> blocks;
    explicit Spectrum(std::vector<Subspectrum>&& m);

    static std::
        pair<Spectrum, std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseMatrix>>>
        energy(const Matrix& hamiltonian_matrix);
    static Spectrum non_energy(
        const Matrix& non_hamiltonian_matrix,
        const std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseMatrix>>&
            unitary_transformation_matrices);
};

std::ostream& operator<<(std::ostream& os, const Spectrum& space);

#endif  //SPINNER_SPECTRUM_H
