#include "OptimizationList.h"

#include <algorithm>
#include <stdexcept>

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
    if (isNonMinimalProjectionsEliminated_) {
        throw std::invalid_argument("Cannot simultaneously eliminate both non-minimal and positive projections");
    }
    if (!isTzSorted_) {
        throw std::invalid_argument("Cannot eliminate positive projections without tz-sort");
    }
    isPositiveProjectionsEliminated_ = true;
    return *this;
}

OptimizationList& OptimizationList::EliminateNonMininalProjections() {
    if (isPositiveProjectionsEliminated_) {
        throw std::invalid_argument("Cannot simultaneously eliminate both positive and non-minimal projections");
    }
    if (!isTzSorted_ || !isTSquaredSorted_) {
        throw std::invalid_argument("Cannot eliminate non-minimal projections without tz-sort or t2-sort");
    }
    isNonMinimalProjectionsEliminated_ = true;
    return *this;
}

OptimizationList& OptimizationList::SSquaredTransform() {
    if (basis_type_ == ITO) {
        throw std::invalid_argument("Cannot use S2-transformation with ITO-basis");
    }
    if (!isTzSorted_) {
        throw std::invalid_argument("Cannot perform S2-transformation without tz-sort");
    }
    if (!isTSquaredSorted_) {
        throw std::invalid_argument("Cannot perform S2-transformation without t2-sort");
    }
    if (!isNonMinimalProjectionsEliminated_) {
        throw std::invalid_argument("Cannot perform S2-transformation without elimination of non-minimal projections");
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

OptimizationList& OptimizationList::FTLMApproximate(FTLMSettings ftlmSettings) {
    if (isSSquaredTransformed()) {
        throw std::invalid_argument("Cannot use FTLM with SSquaredTransformation yet");
    }
    if (!isPositiveProjectionsEliminated()) {
        throw std::invalid_argument("Positive projections elimination required for FTLM due accuracy issues");
    }
    if (ftlmSettings.krylov_subspace_size > ftlmSettings.exact_decomposition_threshold) {
        throw std::invalid_argument("Size of Krylov subspace bigger than threshold of exact decomposition");
    }
    if (ftlmSettings.number_of_seeds == 0) {
        throw std::invalid_argument("Number of seeds in FTLM equals to zero");
    }
    isFTLMApproximated_ = true;
    ftlmSettings_ = ftlmSettings;
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
    if (isNonAbelianSimplified()) {
        throw std::invalid_argument("Cannot symmetrize after non-Abelian simplification");
    }
    if (basis_type_ == ITO && !new_group.properties.is_abelian) {
        throw std::invalid_argument("Currently cannot perform use ITO-basis with non-Abelian symmetries");
    }
    if (isSSquaredTransformed_ && !new_group.properties.is_abelian) {
        throw std::invalid_argument("Currently cannot perform S2-transformation with non-Abelian symmetries");
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

OptimizationList& OptimizationList::NonAbelianSimplify() {
    if (getGroupsToApply().empty()) {
        throw std::invalid_argument("Cannot non-Abelian simplify without symmetrizing");
    }
    size_t number_of_nonabelian_groups = 0;
    for (const auto& group : getGroupsToApply()) {
        if (!group.properties.is_abelian) {
            ++number_of_nonabelian_groups;
        }
    }
    if (number_of_nonabelian_groups == 0) {
        throw std::invalid_argument("Only Abelian groups. Cannot non-Abelian simplify them");
    }
    if (number_of_nonabelian_groups > 1) {
        throw std::invalid_argument("More than one non-Abelian group. Cannot non-Abelian simplify them yet");
    }
    if (getGroupsToApply()[0].properties.is_abelian) {
        throw std::invalid_argument("Non-Abelian group need to be the first group to apply for non-Abelian simplification");
    }
    isNonAbelianSimplified_ = true;
    return *this;
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

bool OptimizationList::isNonMinimalProjectionsEliminated() const {
    return isNonMinimalProjectionsEliminated_;
}

bool OptimizationList::isSSquaredTransformed() const {
    return isSSquaredTransformed_;
}

bool OptimizationList::isFTLMApproximated() const {
    return isFTLMApproximated_;
}

const std::vector<group::Group>& OptimizationList::getGroupsToApply() const {
    return groupsToApply_;
}

bool OptimizationList::isNonAbelianSimplified() const {
    return isNonAbelianSimplified_;
}

const OptimizationList::FTLMSettings& OptimizationList::getFTLMSettings() const {
    return ftlmSettings_.value();
}

}  // namespace common::physical_optimization