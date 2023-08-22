#include "OrderOfSummation.h"

#include <algorithm>
#include <optional>

namespace {
template<typename T>
bool has_duplicates(std::vector<T> vector) {
    std::sort(vector.begin(), vector.end());
    return std::adjacent_find(vector.begin(), vector.end()) != vector.end();
}
}  // namespace

namespace spin_algebra {

std::shared_ptr<const OrderOfSummation> spin_algebra::OrderOfSummation::constructFromOrbits(
    const std::vector<std::vector<std::set<size_t>>>& all_groups_orbits_of_mults,
    size_t number_of_mults,
    size_t number_of_summation) {
    size_t performed_summations = 0;

    std::vector<std::optional<size_t>> pos_of_results_of_sum(
        number_of_mults + number_of_summation - 1,
        std::nullopt);

    auto order_of_summations = std::make_shared<spin_algebra::OrderOfSummation>();

    for (size_t number_of_group = 0; number_of_group < all_groups_orbits_of_mults.size();
         ++number_of_group) {
        const auto& current_orbits = all_groups_orbits_of_mults.at(number_of_group);
        for (const auto& orbit : current_orbits) {
            if (orbit.size() == 1) {
                continue;
            }
            OrderOfSummation::AdditionInstruction instruction;
            instruction.number_of_group = number_of_group;
            instruction.position_of_sum = number_of_mults + performed_summations;
            for (auto pos : orbit) {
                while (pos_of_results_of_sum[pos].has_value()) {
                    pos = pos_of_results_of_sum[pos].value();
                }
                instruction.positions_of_summands.push_back(pos);
            }

            if (has_duplicates(instruction.positions_of_summands)) {
                // analogue of orbit.size() == 1:
                continue;
            }

            for (const auto pos : instruction.positions_of_summands) {
                pos_of_results_of_sum[pos] = instruction.position_of_sum;
            }
            order_of_summations->instructions_.push_back(instruction);

            performed_summations++;
            if (performed_summations == number_of_summation) {
                break;
            }
        }
    }

    while (performed_summations < number_of_summation) {
        OrderOfSummation::AdditionInstruction instruction;
        instruction.position_of_sum = number_of_mults + performed_summations;

        for (size_t i = 0; i < number_of_mults + number_of_summation - 1; ++i) {
            if (pos_of_results_of_sum[i] == std::nullopt) {
                instruction.positions_of_summands.push_back(i);
                if (instruction.positions_of_summands.size() == 2) {
                    break;
                }
            }
        }

        for (const auto pos : instruction.positions_of_summands) {
            pos_of_results_of_sum[pos] = instruction.position_of_sum;
        }
        order_of_summations->instructions_.push_back(instruction);
        performed_summations++;
    }

    return order_of_summations;
}

std::vector<OrderOfSummation::AdditionInstruction>::const_iterator OrderOfSummation::begin() const {
    return instructions_.begin();
}

std::vector<OrderOfSummation::AdditionInstruction>::const_iterator OrderOfSummation::end() const {
    return instructions_.end();
}

const OrderOfSummation::AdditionInstruction& OrderOfSummation::at(size_t i) const {
    return instructions_.at(i);
}

size_t OrderOfSummation::size() const {
    return instructions_.size();
}

}  // namespace spin_algebra