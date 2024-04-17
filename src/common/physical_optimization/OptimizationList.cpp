#include "OptimizationList.h"

#include <algorithm>

namespace common::physical_optimization {

OptimizationList& OptimizationList::TzSort() {
    if (isSSquaredTransformed_) {
        throw std::invalid_argument("Cannot TzSort *AFTER* S2-transformation");
    }
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

OptimizationList& OptimizationList::SSquaredTransform() {
    if (!isTzSorted_) {
        throw std::invalid_argument("Cannot perform S2-transformation without tz-sort");
        // actually, can, but it is inefficient
    }
    for (const auto& group : groupsToApply_) {
        if (!group.properties.is_abelian) {
            throw std::invalid_argument(
                "Currently cannot perform S2-transformation with non-Abelian symmetries");
        }
    }
    isSSquaredTransformed_ = true;
    return *this;
}

OptimizationList& OptimizationList::ITOCalculate() {
    if (!isSSquaredTransformed_) {
        throw std::invalid_argument("Cannot perform ITO-calculation without S2-transformation");
    }
    isITOCalculated_ = true;
    return *this;
}

OptimizationList& OptimizationList::Symmetrize(group::Group new_group) {
    // check if user trying to use the same Group for a second time:
    if (std::count(groupsToApply_.begin(), groupsToApply_.end(), new_group)) {
        return *this;
    }
    // check if groups can be applied at the same time
    for (const auto& old_group : groupsToApply_) {
        if (old_group.size_of_permutations() != new_group.size_of_permutations()) {
            throw std::invalid_argument(
                "Trying to apply groups with different sizes of permutations");
        }
        // we are trying to construct so-called "internal direct product" of groups.
        // it is necessary and sufficient if groups elements commute.
        if (!old_group.do_groups_commute(new_group)) {
            throw std::invalid_argument("Groups do not commute!");
        }
    }
    if (isSSquaredTransformed_) {
        throw std::invalid_argument("Cannot symmetrize *AFTER* S2-transformation");
    }
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

bool OptimizationList::isSSquaredTransformed() const {
    return isSSquaredTransformed_;
}

bool OptimizationList::isITOCalculated() const {
    return isITOCalculated_;
}

const std::vector<group::Group>& OptimizationList::getGroupsToApply() const {
    return groupsToApply_;
}

}  // namespace common::physical_optimization