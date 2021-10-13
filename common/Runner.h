#ifndef JULY_RUNNER_H
#define JULY_RUNNER_H

#include <utility>

#include "entities/space/Space.h"
#include "groups/Group.h"

namespace runner {
class Runner {
  public:
    explicit Runner(std::vector<int> mults);

    void NonAbelianSimplify();

    void Symmetrize(Group new_group);
    void Symmetrize(group::GroupNames group_name, std::vector<Permutation> generators);

    void TzSort();

    [[nodiscard]] const Space& getSpace() const;

    [[nodiscard]] uint32_t getTotalSpaceSize() const;

  private:
    struct SpaceHistory {
        std::vector<Group> applied_groups;
        uint32_t number_of_non_simplified_abelian_groups = 0;
        bool isTzSorted = false;
    };

    // TODO: should we implement it as method of Group_Info?
    [[nodiscard]] bool IsAbelianGroup(const Group& group) const noexcept;

    const spaces::LexicographicIndexConverter converter_;
    Space space_;
    SpaceHistory history_;
};
} // namespace runner

#endif //JULY_RUNNER_H
