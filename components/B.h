#ifndef JULY_B_H
#define JULY_B_H

#include <vector>

struct Decomposition_ {
    unsigned long index;
    double coeff;
};

struct Subspace {
    int n_proj = -1;
    int representation = -1;

    std::vector<std::vector<Decomposition_>> basis;
};

#endif //JULY_B_H
