#ifndef JULY_SUBSPACE_H
#define JULY_SUBSPACE_H

#include <map>
#include <vector>

using Index = unsigned long;
using Coefficient = double;
using Decomposition = std::map<Index, Coefficient>;

struct Subspace {
    int n_proj = -1;
    int representation = -1;

    std::vector<Decomposition> basis;
};

#endif // JULY_SUBSPACE_H
