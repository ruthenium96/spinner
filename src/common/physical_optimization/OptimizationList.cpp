#include "OptimizationList.h"

#include <algorithm>

namespace common::physical_optimization {

OptimizationList& OptimizationList::TzSort() {
    isTzSorted_ = true;
    return *this;
}

OptimizationList& OptimizationList::EliminatePositiveProjections() {
    if (!isTzSorted_) {
        throw std::invalid_argument("Cannot eliminate positive projections without tz-sort");
    }
    isPositiveProjectionsEliminated_ = true;
    return *this;
}
OptimizationList& OptimizationList::Symmetrize(group::Group new_group) {
    // check if user trying to use the same Group for a second time:
    if (std::count(groupsToApply_.begin(), groupsToApply_.end(), new_group)) {
        return *this;
    }
    //    // TODO: symmetrizer does not work correct after non-Abelian simplifier. Fix it.
    //    if (space_history_.isNonAbelianSimplified && !new_group.properties.is_abelian) {
    //        throw std::invalid_argument(
    //            "Symmetrization after using of non-Abelian simplifier causes bugs.");
    //    }
    //    if (!new_group.properties.is_abelian) {
    //        ++space_history_.number_of_non_simplified_abelian_groups;
    //    }
    groupsToApply_.emplace_back(std::move(new_group));

    return *this;
}

OptimizationList& OptimizationList::Symmetrize(
    group::Group::GroupTypeEnum group_name,
    std::vector<group::Permutation> generators) {
    group::Group new_group(group_name, std::move(generators));
    return Symmetrize(new_group);
}

bool OptimizationList::isTzSorted() const {
    return isTzSorted_;
}

bool OptimizationList::isPositiveProjectionsEliminated() const {
    return isPositiveProjectionsEliminated_;
}

const std::vector<group::Group>& OptimizationList::getGroupsToApply() const {
    return groupsToApply_;
}

}  // namespace common::physical_optimization