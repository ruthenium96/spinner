#ifndef JULY_GROUP_H
#define JULY_GROUP_H

#include <vector>
#include <cstdint>
#include "Group_Info.h"

using Permutation = std::vector<uint8_t>;

class Group {
public:
    explicit Group(group::GroupNames group_name, std::vector<Permutation> generators);

    std::vector<std::vector<uint8_t>> permutate(const std::vector<uint8_t>& initial) const;

    bool operator==(const Group &rhs) const;

    bool operator!=(const Group &rhs) const;

    const group::GroupInfo& info;

    std::vector<Permutation> elements_;

private:

    // There are generators and elements for _specific_ input.
    // Length of Permutation vectors -- number of spins.
    std::vector<Permutation> generators_;
};


#endif //JULY_GROUP_H
