#ifndef SPINNER_SUBSPACE_H
#define SPINNER_SUBSPACE_H

#include <ostream>

#include "src/entities/BlockProperties.h"
#include "src/entities/data_structures/UnitarySparseMatrix.h"

namespace space {
struct Subspace {
    BlockProperties properties;
    UnitarySparseMatrix decomposition;

    explicit Subspace(UnitarySparseMatrix&&);
    Subspace() = default;

    friend std::ostream& operator<<(std::ostream& os, const Subspace& subspace);
};
}  // namespace space

#endif  // SPINNER_SUBSPACE_H
