#include "IndexPermutator.h"
#include <cassert>

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
        const group::Permutation& extended_g = extended_permutations_[g_position];

        for (int position = initial_size_of_permutation_; position < extended_g.size(); ++position) {
            spin_algebra::Multiplicity multiplicity_to_push = level.getMultiplicity(extended_g[position]);
            permutated_levels[g_position].setMultiplicity(position, multiplicity_to_push);
        }

        for (size_t position_from = 0; position_from < extended_g.size(); ++position_from) {
            auto position_to = extended_g[position_from];
            if (position_to == position_from) {
                continue;
            }
            // this code works only for P2, so we need to explicitly check:
            assert(extended_g[position_to] == position_from);
            for (const auto& instruction : order_of_summation->getInstructions()) {
                auto pos_one = instruction.positions_of_summands[0];
                auto pos_two = instruction.positions_of_summands[1];
                auto pos_res = instruction.position_of_sum;
                auto mult_one = permutated_levels[g_position].getMultiplicity(pos_one);
                auto mult_two = permutated_levels[g_position].getMultiplicity(pos_two);
                auto mult_res = permutated_levels[g_position].getMultiplicity(pos_res);
                if (mult_one == mult_two) {
                    auto maxMultiplicity = mult_one + mult_two - 1;
                    if ((maxMultiplicity - mult_res) % 4 == 0) {
                        signs[g_position] *= 1;
                    } else {
                        signs[g_position] *= -1;
                    }
                } else {
					signs[g_position] *= 1;
				}
            }
        }
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