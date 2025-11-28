#include "GroupAdapter.h"

namespace spin_algebra {

GroupAdapter::GroupAdapter(const std::vector<group::Group>& groups, size_t number_of_mults) {
    // TODO: fix it somehow based on the information from groups:
    const size_t number_of_summation = number_of_mults - 1;

    std::vector<std::vector<std::set<size_t>>> all_groups_orbits_of_mults;
    for (const auto& group : groups) {
        all_groups_orbits_of_mults.emplace_back(group.construct_orbits_of_mults());
    }

    order_of_summations_ = index_converter::s_squared::OrderOfSummation::constructFromOrbits(
        all_groups_orbits_of_mults,
        number_of_mults,
        number_of_summation);
}

std::shared_ptr<const index_converter::s_squared::OrderOfSummation> GroupAdapter::getOrderOfSummations() const {
    return order_of_summations_;
}

}  // namespace spin_algebra