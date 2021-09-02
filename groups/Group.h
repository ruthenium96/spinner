#ifndef JULY_GROUP_H
#define JULY_GROUP_H

#include <vector>
#include <cstdint>
#include "groups/Group_Info.h"

using Permutation = std::vector<uint8_t>;

// To create specific group, one has to define
// 1) group_size,
// 2) number_of_representations,
// 3) group_in_form_of_generators,
// 4) coefficients_of_projectors
class Group {
public:
    explicit Group(GroupNames group_name, std::vector<Permutation> generators);

    std::vector<std::vector<uint8_t>> permutate(const std::vector<uint8_t>& initial) const;

    const GroupInfo& groupInfo;

    std::vector<Permutation> elements_;

private:

    // There are generators and elements for _specific_ input.
    // Length of Permutation vectors -- number of spins.
    std::vector<Permutation> generators_;
};


#endif //JULY_GROUP_H
