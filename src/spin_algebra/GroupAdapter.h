#ifndef SPINNER_GROUPADAPTER_H
#define SPINNER_GROUPADAPTER_H

#include "src/group/Group.h"
#include "src/spin_algebra/OrderOfSummation.h"

namespace spin_algebra {

class GroupAdapter {
  public:
    GroupAdapter(const std::vector<group::Group>& groups, size_t number_of_mults);
    std::shared_ptr<const OrderOfSummation> getOrderOfSummations() const;
    const std::vector<group::CayleyTable>& getAllGroupsCayleyTables() const;

  private:
    std::shared_ptr<const OrderOfSummation> order_of_summations_;
    std::vector<group::CayleyTable> all_groups_cayley_tables_;
};

}  // namespace spin_algebra

#endif  //SPINNER_GROUPADAPTER_H
