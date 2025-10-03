#include "IndexPermutator.h"
#include <cassert>
#include <stdexcept>
#include "src/group/Group.h"
#include "src/spin_algebra/Multiplicity.h"

namespace index_converter::s_squared {

IndexPermutator::IndexPermutator(std::shared_ptr<const IndexConverter> converter, const group::Group& group) :
	converter_(converter) {
    if (group.getElements().size() > 2) {
        throw std::invalid_argument("We cannot use non-abelian groups in s_squared::IndexConverter");
    }

	initial_size_of_permutation_ = group.size_of_permutations();

	extended_permutations_.reserve(group.getElements().size());
	for (int g = 0; g < group.getElements().size(); ++g) {
		extended_permutations_.push_back(extendPermutation(group.getElements()[g]));
	}
}

uint32_t IndexPermutator::get_total_space_size() const {
	return converter_->get_total_space_size();
}

inline void IndexPermutator::construct_level_and_sign(
    const group::Permutation& extended_g,
    const Level& level,
    Level& permutated_level, 
    int8_t& sign) const {
    for (int position = initial_size_of_permutation_; position < extended_g.size(); ++position) {
        spin_algebra::Multiplicity multiplicity_to_push = level.getMultiplicity(extended_g[position]);
        permutated_level.setMultiplicity(position, multiplicity_to_push);
    }

    std::map<unsigned long, spin_algebra::Multiplicity> visited_multiplicities;

    for (const auto& instruction : converter_->getOrderOfSummation()->getInstructions()) {
        const auto pos_one = instruction.positions_of_summands[0];
        const auto pos_two = instruction.positions_of_summands[1];
        const auto pos_res = instruction.position_of_sum;

        if (pos_one == extended_g[pos_one] || pos_two == extended_g[pos_two]) {
            continue;
        }

        auto mult_one = permutated_level.getMultiplicity(pos_one);
        auto mult_two = permutated_level.getMultiplicity(pos_two);
        auto mult_res = permutated_level.getMultiplicity(pos_res);
        visited_multiplicities[pos_one] = mult_one;
        visited_multiplicities[pos_two] = mult_two;
        visited_multiplicities[pos_res] = mult_res;
    }
    size_t acc = 0;
    for (const auto& [key, value] : visited_multiplicities)
    {
        acc += value;
    }
    if (acc % 4 != visited_multiplicities.size() % 4) {
        sign *= -1;
    }
}

std::vector<IndexWithSign> IndexPermutator::convert_index_to_permutated_indexes(uint32_t index) const {
	auto order_of_summation = converter_->getOrderOfSummation();

    const auto state = converter_->convert_index_to_state(index);
    const auto& level = state.first;
    auto unshifted_projection = state.second;

    auto empty_level = Level(level.getInitialMultiplicities(), 
                                     order_of_summation->size());

    std::vector<Level> permutated_levels(extended_permutations_.size(), empty_level);
    std::vector<int8_t> signs(extended_permutations_.size(), +1);

    for (int g_position = 0; g_position < extended_permutations_.size(); ++g_position) {
        construct_level_and_sign(extended_permutations_[g_position], level, permutated_levels[g_position], signs[g_position]);
    }

    std::vector<IndexWithSign> answer;
    answer.resize(permutated_levels.size());

    for (uint8_t g = 0; g < permutated_levels.size(); ++g) {
        answer[g].index = converter_->convert_state_to_index(permutated_levels[g], unshifted_projection).value();
        answer[g].sign = signs[g];
    }
    return answer;
}


group::Permutation IndexPermutator::extendPermutation(const group::Permutation& permutation) const {
    std::vector<std::optional<group::PermutationArrayType>> mb_extended_permutation;
	auto order_of_summation = converter_->getOrderOfSummation();
    mb_extended_permutation.resize(permutation.size() + order_of_summation->size(), std::nullopt);
    for (int pos = 0; pos < permutation.size(); ++pos) {
        mb_extended_permutation[pos] = permutation[pos];
    }
    int initialized_values = permutation.size();

    while(initialized_values < mb_extended_permutation.size()) {
        std::vector<OrderOfSummation::AdditionInstruction> permuted_order;
        permuted_order.reserve(order_of_summation->size());
        for (const auto& instruction : order_of_summation->getInstructions()) {
            auto pos_one = instruction.positions_of_summands[0];
            auto pos_two = instruction.positions_of_summands[1];
            if (mb_extended_permutation[pos_one].has_value() && 
                mb_extended_permutation[pos_two].has_value()) {
                auto perm_pos_one = mb_extended_permutation[pos_one].value();
                auto perm_pos_two = mb_extended_permutation[pos_two].value();
                OrderOfSummation::AdditionInstruction perm_instruction;
                perm_instruction.positions_of_summands = {perm_pos_one, perm_pos_two};
                perm_instruction.position_of_sum = instruction.position_of_sum;
                permuted_order.push_back(perm_instruction);
            }
        }

        for (const auto& perm_instruction : permuted_order) {
            auto perm_pos_one = perm_instruction.positions_of_summands[0];
            auto perm_pos_two = perm_instruction.positions_of_summands[1];
            auto perm_pos_sum = perm_instruction.position_of_sum;
            for (const auto& instruction : order_of_summation->getInstructions()) {
                auto pos_one = instruction.positions_of_summands[0];
                auto pos_two = instruction.positions_of_summands[1];
                auto pos_sum = instruction.position_of_sum;

                if ((perm_pos_one == pos_one && perm_pos_two == pos_two) ||
                (perm_pos_two == pos_one && perm_pos_one == pos_two)) {
                    if (!mb_extended_permutation[pos_sum].has_value()) {
                        mb_extended_permutation[pos_sum] = perm_pos_sum;
                        ++initialized_values;
                    }
                }
            }
        }
    }
    group::Permutation extended_permutation;
    extended_permutation.reserve(mb_extended_permutation.size());
    for (const auto& optional : mb_extended_permutation) {
        extended_permutation.push_back(optional.value());
    }

    return extended_permutation;
}

} // namespace index_converter::s_squared