#ifndef JULY_SUBSPACE_H
#define JULY_SUBSPACE_H

#include <map>
#include <vector>

typedef unsigned long Index;
typedef double Coefficient;

struct Subspace {
    int n_proj = -1;
    int representation = -1;

    std::vector<std::map<Index, Coefficient>> basis;
};

#endif // JULY_SUBSPACE_H
