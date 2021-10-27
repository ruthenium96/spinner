#ifndef JULY_SUBSPACE_H
#define JULY_SUBSPACE_H

#include <ostream>
#include "entities/BlockProperties.h"
#include "entities/data_structures/UnitarySpaseMatrix.h"

struct Subspace {

    BlockProperties properties;
    UnitarySpaseMatrix decomposition;

    explicit Subspace(UnitarySpaseMatrix&&);
    Subspace() = default;

    friend std::ostream &operator<<(std::ostream &os, const Subspace &subspace);
};

#endif // JULY_SUBSPACE_H
