#ifndef JULY_SUBSPACE_H
#define JULY_SUBSPACE_H

#include <ostream>
#include "entities/BlockProperties.h"
#include "entities/space/SubspaceData.h"

struct Subspace {

    BlockProperties properties;
    SubspaceData decomposition;

    explicit Subspace(SubspaceData&&);
    Subspace() = default;

    friend std::ostream &operator<<(std::ostream &os, const Subspace &subspace);
};

#endif // JULY_SUBSPACE_H
