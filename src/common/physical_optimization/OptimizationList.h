#ifndef SPINNER_OPTIMIZATIONLIST_H
#define SPINNER_OPTIMIZATIONLIST_H

#include <optional>
#include "src/group/Group.h"

namespace common::physical_optimization {
// This class memorizes (physical) optimizations to apply and checks their _internal_ consistence.
class OptimizationList {
  public:

    enum BasisType {LEX, ITO};
    struct FTLMSettings {
      size_t krylov_subspace_size;
      size_t exact_decomposition_threshold;
      size_t number_of_seeds;
    };

    explicit OptimizationList(BasisType basis_type = BasisType::LEX);

    OptimizationList& TzSort();
    OptimizationList& TSquaredSort();
    OptimizationList& EliminateNonMininalProjections();
    OptimizationList& EliminatePositiveProjections();
    OptimizationList& SSquaredTransform();
    OptimizationList& Symmetrize(group::Group new_group);
    OptimizationList&
    Symmetrize(group::Group::GroupTypeEnum group_name, std::vector<group::Permutation> generators);
    // TODO: OptimizationList& NonAbelianSimplify();
    OptimizationList& FTLMApproximate(FTLMSettings ftlmSettings);

    bool isLexBasis() const;
    bool isITOBasis() const;
    bool isTzSorted() const;
    bool isTSquaredSorted() const;
    bool isPositiveProjectionsEliminated() const;
    bool isNonMinimalProjectionsEliminated() const;
    bool isSSquaredTransformed() const;
    bool isFTLMApproximated() const;
    const std::vector<group::Group>& getGroupsToApply() const;
    const FTLMSettings& getFTLMSettings() const;

  private:
    bool isTzSorted_ = false;
    bool isTSquaredSorted_ = false;
    bool isPositiveProjectionsEliminated_ = false;
    bool isNonMinimalProjectionsEliminated_ = false;
    bool isSSquaredTransformed_ = false;
    std::vector<group::Group> groupsToApply_;
    BasisType basis_type_;
    // TODO: something about NonAbelianSimplifier

    bool isFTLMApproximated_ = false;
    std::optional<FTLMSettings> ftlmSettings_;
};
}  // namespace common::physical_optimization
#endif  //SPINNER_OPTIMIZATIONLIST_H
