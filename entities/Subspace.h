#ifndef JULY_SUBSPACE_H
#define JULY_SUBSPACE_H

#include <cstdint>
#include <map>
#include <vector>
#include "BlockProperties.h"

using DecompositionMap = std::map<uint32_t, double>;

struct Subspace {
    BlockProperties properties;
    std::vector<DecompositionMap> basis;
};

#endif // JULY_SUBSPACE_H
