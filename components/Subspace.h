#ifndef JULY_SUBSPACE_H
#define JULY_SUBSPACE_H

#include <vector>

typedef unsigned long Index;
typedef double Coefficient;

struct Decomposition_ {
    unsigned long index;
    double coeff;
};

struct Subspace {
    int n_proj = -1;
    int representation = -1;

    std::vector<std::map<Index, Coefficient>> basis;

//    std::vector<std::vector<Decomposition_>> basis;
};

#endif //JULY_SUBSPACE_H
