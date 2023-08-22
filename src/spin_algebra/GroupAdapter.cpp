#include "GroupAdapter.h"

namespace spin_algebra {

GroupAdapter::GroupAdapter(const std::vector<group::Group>& groups, size_t number_of_mults) {
    // TODO: fix it somehow based on the information from groups:
    const size_t number_of_summation = number_of_mults - 1;

    std::vector<std::vector<std::set<size_t>>> all_groups_orbits_of_mults;
    std::vector<group::CayleyTable> all_groups_cayley_tables;
    for (const auto& group : groups) {
        all_groups_orbits_of_mults.emplace_back(group.construct_orbits_of_mults());
        all_groups_cayley_tables.emplace_back(group.properties.cayley_table);
    }

    representationsMultiplier_ = RepresentationsMultiplier(all_groups_cayley_tables);

    order_of_summations_ = spin_algebra::OrderOfSummation::constructFromOrbits(
        all_groups_orbits_of_mults,
        number_of_mults,
        number_of_summation);
}

std::shared_ptr<const OrderOfSummation> GroupAdapter::getOrderOfSummations() const {
    return order_of_summations_;
}

const RepresentationsMultiplier& GroupAdapter::getRepresentationMultiplier() const {
    return representationsMultiplier_;
}
}  // namespace spin_algebra