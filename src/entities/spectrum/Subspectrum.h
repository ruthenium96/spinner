#ifndef SPINNER_SUBSPECTRUM_H
#define SPINNER_SUBSPECTRUM_H

#include "src/entities/BlockProperties.h"
#include "src/entities/data_structures/AbstractMatrix.h"
#include "src/entities/data_structures/AbstractVector.h"
#include "src/entities/matrix/Submatrix.h"

struct Subspectrum {
    BlockProperties properties;
    std::unique_ptr<quantum::linear_algebra::AbstractVector> raw_data;

    Subspectrum() = delete;
    Subspectrum(
        std::unique_ptr<quantum::linear_algebra::AbstractVector> raw_data_,
        BlockProperties properties_);

    static std::pair<Subspectrum, std::unique_ptr<quantum::linear_algebra::AbstractMatrix>>
    energy(const Submatrix& hamiltonian_submatrix);
    static Subspectrum non_energy(
        const Submatrix& non_hamiltonian_submatrix,
        const std::unique_ptr<quantum::linear_algebra::AbstractMatrix>&
            unitary_transformation_matrix);

    friend std::ostream& operator<<(std::ostream& os, const Subspectrum& subspectrum);
};

#endif  //SPINNER_SUBSPECTRUM_H
