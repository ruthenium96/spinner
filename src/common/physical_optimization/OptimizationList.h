#ifndef SPINNER_OPTIMIZATIONLIST_H
#define SPINNER_OPTIMIZATIONLIST_H

#include <src/group/Group.h>
namespace common::physical_optimization {
// This class memorizes (physical) optimizations to apply and checks their _internal_ consistence.
class OptimizationList {
  public:
    OptimizationList& TzSort();
    OptimizationList& EliminatePositiveProjections();
    OptimizationList& SSquaredTransform();
    OptimizationList& Symmetrize(group::Group new_group);
    OptimizationList&
    Symmetrize(group::Group::GroupTypeEnum group_name, std::vector<group::Permutation> generators);
    // TODO: OptimizationList& NonAbelianSimplify();

    bool isTzSorted() const;
    bool isPositiveProjectionsEliminated() const;
    bool isSSquaredTransformed() const;
    const std::vector<group::Group>& getGroupsToApply() const;

  private:
    bool isTzSorted_ = false;
    bool isPositiveProjectionsEliminated_ = false;
    bool isSSquaredTransformed_ = false;
    std::vector<group::Group> groupsToApply_;
    // TODO: something about NonAbelianSimplifier
};
}  // namespace common::physical_optimization
#endif  //SPINNER_OPTIMIZATIONLIST_H
