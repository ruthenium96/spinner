#ifndef SPINNER_SUBSPACE_H
#define SPINNER_SUBSPACE_H

#include <ostream>

#include "src/entities/BlockProperties.h"
#include "src/entities/data_structures/AbstractSparseMatrix.h"

namespace space {
struct Subspace {
    BlockProperties properties;
    std::unique_ptr<quantum::linear_algebra::AbstractSparseMatrix> decomposition;

    explicit Subspace(std::unique_ptr<quantum::linear_algebra::AbstractSparseMatrix>&&);
    Subspace();

    friend std::ostream& operator<<(std::ostream& os, const Subspace& subspace);
};
}  // namespace space

#endif  // SPINNER_SUBSPACE_H
