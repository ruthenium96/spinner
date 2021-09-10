#ifndef JULY_BLOCKPROPERTIES_H
#define JULY_BLOCKPROPERTIES_H

struct BlockProperties {
    int n_proj = -1;
    int dimensionality = 1;
    int degeneracy = 1;
    std::vector<int> representation;
};

#endif //JULY_BLOCKPROPERTIES_H
