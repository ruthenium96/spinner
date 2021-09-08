#ifndef JULY_SUBSPACE_H
#define JULY_SUBSPACE_H

#include <cstdint>
#include <map>
#include <vector>

using DecompositionMap = std::map<uint32_t, double>;

struct Subspace {
    // TODO: I guess, we can pack all subspace properties to Properties class.
    int n_proj = -1;
    int degeneracy = 1;
    std::vector<int> representation;

    std::vector<DecompositionMap> basis;
};

#endif // JULY_SUBSPACE_H
