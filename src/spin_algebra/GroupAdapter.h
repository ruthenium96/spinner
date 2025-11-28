#ifndef SPINNER_GROUPADAPTER_H
#define SPINNER_GROUPADAPTER_H

#include "src/group/Group.h"
#include "src/common/index_converter/s_squared/OrderOfSummation.h"

namespace spin_algebra {

class GroupAdapter {
  public:
    GroupAdapter(const std::vector<group::Group>& groups, size_t number_of_mults);
    std::shared_ptr<const index_converter::s_squared::OrderOfSummation> getOrderOfSummations() const;

  private:
    std::shared_ptr<const index_converter::s_squared::OrderOfSummation> order_of_summations_;
};

}  // namespace spin_algebra

#endif  //SPINNER_GROUPADAPTER_H
