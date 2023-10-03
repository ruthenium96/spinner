#ifndef SPINNER_SUBSPECTRUM_H
#define SPINNER_SUBSPECTRUM_H

#include "src/entities/BlockProperties.h"
#include "src/entities/data_structures/AbstractDenseVector.h"
#include "src/entities/matrix/Submatrix.h"

struct Subspectrum {
    BlockProperties properties;
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector> raw_data;

    Subspectrum() = delete;
    Subspectrum(
        std::unique_ptr<quantum::linear_algebra::AbstractDenseVector> raw_data_,
        BlockProperties properties_);

    static Subspectrum energy_without_eigenvectors(const Submatrix& hamiltonian_submatrix);
    static std::
        pair<Subspectrum, std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
        energy(const Submatrix& hamiltonian_submatrix);
    static Subspectrum non_energy(
        const Submatrix& non_hamiltonian_submatrix,
        const std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>&
            unitary_transformation_matrix);

    friend std::ostream& operator<<(std::ostream& os, const Subspectrum& subspectrum);
};

#endif  //SPINNER_SUBSPECTRUM_H
