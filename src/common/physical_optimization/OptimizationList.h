#ifndef JULY_OPTIMIZATIONLIST_H
#define JULY_OPTIMIZATIONLIST_H

#include <src/group/Group.h>
namespace common::physical_optimization {
// This class memorizes (physical) optimizations to apply and checks their _internal_ consistence.
class OptimizationList {
  public:
    OptimizationList& TzSort();
    OptimizationList& EliminatePositiveProjections();
    OptimizationList& Symmetrize(group::Group new_group);
    OptimizationList&
    Symmetrize(group::Group::GroupTypeEnum group_name, std::vector<group::Permutation> generators);
    // TODO: OptimizationList& NonAbelianSimplify();
    // TODO: OptimizationList& SSquaredTransform();

    bool isTzSorted() const;
    bool isPositiveProjectionsEliminated() const;
    const std::vector<group::Group>& getGroupsToApply() const;

  private:
    bool isTzSorted_ = false;
    bool isPositiveProjectionsEliminated_ = false;
    std::vector<group::Group> groupsToApply_;
    // TODO: something about NonAbelianSimplifier
    // TODO: something about SSquaredTransform
};
}  // namespace common::physical_optimization
#endif  //JULY_OPTIMIZATIONLIST_H
