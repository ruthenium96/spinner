#ifndef JULY_SUBSPACE_H
#define JULY_SUBSPACE_H

#include <cstdint>
#include <map>
#include <vector>

using Lex_Index = uint32_t;
using Coefficient = double;
using Decomposition = std::map<Lex_Index, Coefficient>;

struct Subspace {
    // TODO: I guess, we can pack all subspace properties to Properties class.
    int n_proj = -1;
    int representation = -1;

    std::vector<Decomposition> basis;
};

#endif // JULY_SUBSPACE_H
