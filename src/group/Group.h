#ifndef SPINNER_GROUP_H
#define SPINNER_GROUP_H

#include <set>
#include <vector>

#include "AlgebraicProperties.h"

namespace spinner::group {

using PermutationArrayType = uint8_t;
using Permutation = std::vector<PermutationArrayType>;

class Group {
  public:
    explicit Group(GroupType group_type, std::vector<Permutation> generators);

    size_t size_of_permutations() const;

    bool do_groups_commute(const Group& rhs) const;

    std::vector<std::set<size_t>> construct_orbits_of_mults() const;

    bool operator==(const Group& rhs) const;

    bool operator!=(const Group& rhs) const;

    const AlgebraicProperties properties;

    const std::vector<Permutation>& getElements() const;
    const std::vector<Permutation>& getGenerators() const;

  private:
    std::vector<Permutation> elements_;
    // There are generators and elements for _specific_ input.
    // Length of Permutation vectors -- number of spins.
    std::vector<Permutation> generators_;
};

} // namespace spinner::group

#endif  //SPINNER_GROUP_H
