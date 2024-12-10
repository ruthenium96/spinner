#include "IndexConverter.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <numeric>
#include <optional>

#include "src/common/index_converter/AbstractIndexConverter.h"
#include "src/spin_algebra/Multiplicity.h"
#include "src/spin_algebra/MultiplicityDirectSum.h"
#include "OrderOfSummation.h"
#include "src/common/PrintingFunctions.h"

namespace index_converter::s_squared {

IndexConverter::IndexConverter(
	const std::vector<spin_algebra::Multiplicity>& mults,
    const std::shared_ptr<const OrderOfSummation>& order_of_summation) :
	AbstractIndexConverter(mults), order_of_summation_(order_of_summation) {
	auto number_of_initial_mults_ = get_mults().size();
    auto number_of_all_mults_ = number_of_initial_mults_ + order_of_summation_->size();

    auto empty_level = Level(
        std::make_shared<std::vector<spin_algebra::Multiplicity>>(get_mults()),
        order_of_summation_->size());

    std::vector<Level> result_of_summation = {empty_level};

    for (const auto& instruction : *order_of_summation_) {
        // we currently do not support cases with three or more summands:
        assert(instruction.positions_of_summands.size() == 2);
        size_t pos_one = instruction.positions_of_summands[0];
        size_t pos_two = instruction.positions_of_summands[1];
        size_t pos_sum = instruction.position_of_sum;
        std::vector<Level> temp_result;
        for (const auto& level : result_of_summation) {
            spin_algebra::MultiplicityDirectSum mult_one(level.getMultiplicity(pos_one));
            spin_algebra::MultiplicityDirectSum mult_two(level.getMultiplicity(pos_two));
            auto mult_direct_sum = mult_one * mult_two;
            for (auto mult_sum : mult_direct_sum.getMultiplicities()) {
                temp_result.push_back(level);
                temp_result.back().setMultiplicity(pos_sum, mult_sum);
            }
        }
        std::swap(result_of_summation, temp_result);
        temp_result.clear();
    }

    auto max_mult = std::accumulate(get_mults().begin(), get_mults().end(), 1 - get_mults().size());

    s_squared_levels_.resize((max_mult - 1) / 2 + 1);

    for (auto& level : result_of_summation) {
        auto index_of_total_spin_block = (level.total() - 1) / 2;

        s_squared_levels_[index_of_total_spin_block].emplace_back(std::move(level));
    }

    for (auto& block : s_squared_levels_) {
        std::sort(block.begin(), block.end());
    }

    cumulative_sum_.resize(s_squared_levels_.size() + 1);
    cumulative_sum_[0] = 0;

    size_t i = 0;
    for (const auto& block : s_squared_levels_) {
        auto multiplicity_of_block = block[0].total();
        cumulative_sum_[i + 1] = cumulative_sum_[i] + block.size() * multiplicity_of_block;
        ++i;
    }

    common::sSquaredIndexConverterPrint(*this);
}

spin_algebra::Multiplicity IndexConverter::convert_index_to_total_multiplicity(uint32_t index) const {
    auto block_iterator = std::upper_bound(cumulative_sum_.begin(),
                                        cumulative_sum_.end(),
                                        index) - 1;
    auto index_in_block = index - *block_iterator;
    auto index_of_block = std::distance(cumulative_sum_.begin(), block_iterator);
    return s_squared_levels_[index_of_block][0].total();
}

uint8_t IndexConverter::convert_index_to_tz_projection(uint32_t index) const {
    auto block_iterator = std::upper_bound(cumulative_sum_.begin(),
                                        cumulative_sum_.end(),
                                        index) - 1;
    auto index_in_block = index - *block_iterator;
    auto index_of_block = std::distance(cumulative_sum_.begin(), block_iterator);
    auto multiplicity_of_block = s_squared_levels_[index_of_block][0].total();

    auto shift = (get_max_ntz_proj() - multiplicity_of_block) / 2;

    return index_in_block % multiplicity_of_block + shift;
}

std::pair<const Level&, uint8_t> IndexConverter::convert_index_to_state(uint32_t index) const {
    auto block_iterator = std::upper_bound(cumulative_sum_.begin(),
                                    cumulative_sum_.end(),
                                    index) - 1;
    auto index_in_block = index - *block_iterator;
    auto index_of_block = std::distance(cumulative_sum_.begin(), block_iterator);
    auto multiplicity_of_block = s_squared_levels_[index_of_block][0].total();

    auto number_of_state = index_in_block / multiplicity_of_block;
    auto unshifted_projection = index_in_block % multiplicity_of_block;

    return {s_squared_levels_[index_of_block][number_of_state], unshifted_projection};
}

std::optional<uint32_t> IndexConverter::convert_state_to_index(const Level& state, uint8_t projection) const {
    auto total_multiplicity = state.total();

    auto index_of_block = (total_multiplicity - 1) / 2;

    const auto& block = s_squared_levels_[index_of_block];

    auto iterator_in_block = std::lower_bound(block.begin(), block.end(), state);
    if (iterator_in_block == block.end() || std::is_neq(*iterator_in_block <=> state)) {
        return std::nullopt;
    }
    
    auto number_of_state = std::distance(block.begin(), iterator_in_block);

    auto index_in_block = number_of_state * total_multiplicity + projection;

    return index_in_block + cumulative_sum_[index_of_block];
}

std::shared_ptr<const OrderOfSummation> IndexConverter::getOrderOfSummation() const {
    return order_of_summation_;
}

}