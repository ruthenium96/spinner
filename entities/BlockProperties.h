#ifndef JULY_BLOCKPROPERTIES_H
#define JULY_BLOCKPROPERTIES_H

#include <ostream>
#include <string>
#include <vector>

struct BlockProperties {
    friend std::ostream &operator<<(std::ostream &os, const BlockProperties &properties);
    std::string get_representation_name() const;

    int n_proj = -1;
    int dimensionality = 1;
    int degeneracy = 1;
    std::vector<int> representation;
};

#endif //JULY_BLOCKPROPERTIES_H
