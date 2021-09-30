#pragma once

#include <utility>

#include "entities/Space.h"
#include "groups/Group.h"

namespace runner {
class Runner {
  public:
    explicit Runner(std::vector<int> mults);

    void NonAbelianSimplify();

    void Symmetrize(Group new_group);
    void Symmetrize(group::GroupNames group_name, std::vector<Permutation> generators);

    void TzSort();

    const Space& getSpace() const;

    uint32_t getTotalSpaceSize() const;

  private:
    struct SpaceHistory {
        bool isTzSorted = false;
        std::vector<Group> applied_groups;
        uint32_t number_of_non_simplified_abelian_groups = 0;
    };

    // TODO: should we implement it as method of Group_Info?
    bool IsItAbelianGroup(const Group& group) const;

    const spaces::LexicographicIndexConverter converter_;
    Space space_;
    SpaceHistory history_;
};
} // namespace runner