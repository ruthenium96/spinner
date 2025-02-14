#include "OptimizationList.h"

#include <algorithm>

namespace common::physical_optimization {

OptimizationList::OptimizationList(BasisType basis_type) {
    basis_type_ = basis_type;
}

OptimizationList& OptimizationList::TzSort() {
    if (isSSquaredTransformed_) {
        throw std::invalid_argument("Cannot TzSort *AFTER* S2-transformation");
    }
    isTzSorted_ = true;
    return *this;
}

OptimizationList& OptimizationList::TSquaredSort() {
    isTSquaredSorted_ = true;
    return *this;
}

OptimizationList& OptimizationList::EliminatePositiveProjections() {
    if (isLexBasis() && !isTzSorted_) {
        throw std::invalid_argument("Cannot eliminate positive projections without tz-sort");
    }
    if (isITOBasis() && (!isTzSorted_ || !isTSquaredSorted_)) {
        throw std::invalid_argument("Cannot eliminate non-minimal projections without tz-sort or t2-sort");
    }
    isPositiveProjectionsEliminated_ = true;
    return *this;
}

OptimizationList& OptimizationList::SSquaredTransform() {
    if (basis_type_ == ITO) {
        throw std::invalid_argument("Cannot use S2-transformation with ITO-basis");
    }
    if (!isTzSorted_) {
        throw std::invalid_argument("Cannot perform S2-transformation without tz-sort");
        // actually, can, but it is inefficient
    }
    if (!isTSquaredSorted_) {
        throw std::invalid_argument("Cannot perform S2-transformation without t2-sort");
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
    if (basis_type_ == ITO && !new_group.properties.is_abelian) {
        throw std::invalid_argument("Currently cannot perform use ITO-basis with non-Abelian symmetries");
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

bool OptimizationList::isLexBasis() const {
    return basis_type_ == LEX;
}

bool OptimizationList::isITOBasis() const {
    return basis_type_ == ITO;
}

bool OptimizationList::isTzSorted() const {
    return isTzSorted_;
}

bool OptimizationList::isTSquaredSorted() const {
    return isTSquaredSorted_;
}

bool OptimizationList::isPositiveProjectionsEliminated() const {
    return isPositiveProjectionsEliminated_;
}

bool OptimizationList::isSSquaredTransformed() const {
    return isSSquaredTransformed_;
}

const std::vector<group::Group>& OptimizationList::getGroupsToApply() const {
    return groupsToApply_;
}

}  // namespace common::physical_optimization