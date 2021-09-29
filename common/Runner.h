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
    };

    const spaces::LexicographicIndexConverter converter_;
    Space space_;
    SpaceHistory history_;
};
} // namespace runner