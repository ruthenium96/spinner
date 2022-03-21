#ifndef JULY_SUBSPACE_H
#define JULY_SUBSPACE_H

#include <ostream>

#include "src/entities/BlockProperties.h"
#include "src/entities/data_structures/UnitarySparseMatrix.h"

struct Subspace {
    BlockProperties properties;
    UnitarySparseMatrix decomposition;

    explicit Subspace(UnitarySparseMatrix&&);
    Subspace() = default;

    friend std::ostream& operator<<(std::ostream& os, const Subspace& subspace);
};

#endif  // JULY_SUBSPACE_H
