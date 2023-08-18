#ifndef SPINNER_GROUPADAPTER_H
#define SPINNER_GROUPADAPTER_H

#include "src/group/Group.h"
#include "src/spin_algebra/OrderOfSummation.h"
#include "src/spin_algebra/RepresentationsMultiplier.h"

namespace spin_algebra {

class GroupAdapter {
  public:
    GroupAdapter(const std::vector<group::Group>& groups, size_t number_of_mults);
    std::shared_ptr<const OrderOfSummation> getOrderOfSummations() const;
    const RepresentationsMultiplier& getRepresentationMultiplier() const;

  private:
    std::shared_ptr<const OrderOfSummation> order_of_summations_;
    RepresentationsMultiplier representationsMultiplier_;
};

}  // namespace spin_algebra

#endif  //SPINNER_GROUPADAPTER_H
