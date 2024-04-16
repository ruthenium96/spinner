#include "SSquaredConverter.h"

#include <cassert>
#include <numeric>

namespace spin_algebra {
SSquaredConverter::SSquaredConverter(
    const std::vector<Multiplicity>& multiplicities_to_sum,
    const std::shared_ptr<const OrderOfSummation>& order_of_summation,
    const RepresentationsMultiplier& representationsMultiplier):
    order_of_summation_(order_of_summation) {

    number_of_initial_mults_ = multiplicities_to_sum.size();
    number_of_all_mults_ = number_of_initial_mults_ + order_of_summation_->size();

    auto empty_history = SSquaredState(
        std::make_shared<std::vector<Multiplicity>>(multiplicities_to_sum),
        order_of_summation_->size(),
        representationsMultiplier.getNumberOfRepresentations());

    std::vector<SSquaredState> result_of_summation = {empty_history};

    for (const auto& instruction : *order_of_summation_) {
        // we currently do not support cases with three or more summands:
        assert(instruction.positions_of_summands.size() == 2);
        size_t pos_one = instruction.positions_of_summands[0];
        size_t pos_two = instruction.positions_of_summands[1];
        size_t pos_sum = instruction.position_of_sum;
        std::vector<SSquaredState> temp_result;
        for (const auto& history : result_of_summation) {
            MultiplicityDirectSum mult_one(history.getMultiplicity(pos_one));
            MultiplicityDirectSum mult_two(history.getMultiplicity(pos_two));
            auto mult_direct_sum = mult_one * mult_two;
            for (auto mult_sum : mult_direct_sum.getMultiplicities()) {
                const auto& representations_one = history.getRepresentations(pos_one);
                const auto& representations_two = history.getRepresentations(pos_two);
                auto representations_sum = representationsMultiplier.multiplyRepresentations(
                    representations_one,
                    representations_two,
                    instruction.number_of_group,
                    mult_one.getMultiplicities().at(0),
                    mult_two.getMultiplicities().at(0),
                    mult_sum);
                temp_result.push_back(history);
                temp_result.back().setMultiplicity(pos_sum, mult_sum);
                temp_result.back().setRepresentations(pos_sum, std::move(representations_sum));
            }
        }
        std::swap(result_of_summation, temp_result);
        temp_result.clear();
    }

    for (auto& history : result_of_summation) {
        auto final_properties = history.back();

        auto property_iterator = properties_indexes_.find(final_properties);
        if (property_iterator == properties_indexes_.end()) {
            states_.emplace_back();
            properties_indexes_[final_properties] = states_.size() - 1;
        }
        size_t index = properties_indexes_[final_properties];

        states_[index].emplace_back(std::move(history));
    }

    cumulative_sum_.resize(states_.size() + 1);
    cumulative_sum_.back() = 1;

    size_t i = 0;
    for (const auto& block : states_) {
        cumulative_sum_[i + 1] = cumulative_sum_[i] + block.size();
        ++i;
    }
}

std::optional<std::reference_wrapper<const std::vector<SSquaredState>>>
SSquaredConverter::block_with_property(SSquaredState::Properties properties) const {
    const auto property_iterator = properties_indexes_.find(properties);
    if (property_iterator == properties_indexes_.end()) {
        return std::nullopt;
    } else {
        return states_.at(property_iterator->second);
    }
}

std::vector<uint8_t>
SSquaredConverter::constructRanksOfTZero(uint32_t center_a, uint32_t center_b) const {
// todo: cache results, there is just O(N^3) memory required
    std::vector<uint8_t> answer(number_of_all_mults_, 255);

    for (size_t i = 0; i < number_of_initial_mults_; ++i) {
        if (i == center_a || i == center_b) {
            answer[i] = 1;
        } else {
            answer[i] = 0;
        }
    }

    for (const auto& instruction : *order_of_summation_) {
        size_t first_spin = instruction.positions_of_summands.at(0);
        size_t second_spin = instruction.positions_of_summands.at(1);
        size_t result_spin = instruction.position_of_sum;

        uint8_t result_rank = (answer[first_spin] + answer[second_spin]) % 2;
        answer[result_spin] = result_rank;
    }

    return answer;
}

const SSquaredState& SSquaredConverter::at(size_t number) const {
    auto it =
        std::lower_bound(cumulative_sum_.begin(), cumulative_sum_.end(), number);
    auto start_of_block = *it;
    auto number_of_block = std::distance(cumulative_sum_.begin(), it);
    auto number_in_block = number - start_of_block;
    assert(number_in_block >= 0);

    return states_.at(number_of_block).at(number_in_block);
}

size_t SSquaredConverter::number_in_block(size_t number) const {
    auto it =
        std::lower_bound(cumulative_sum_.begin(), cumulative_sum_.end(), number);
    auto start_of_block = *it;
    auto number_in_block = number - start_of_block;

    return number_in_block;
}

std::optional<std::vector<size_t>>
SSquaredConverter::indexes_with_property(SSquaredState::Properties properties) const {
    const auto property_iterator = properties_indexes_.find(properties);
    if (property_iterator == properties_indexes_.end()) {
        return std::nullopt;
    } else {
        auto number_of_block = property_iterator->second;
        auto start_of_block = cumulative_sum_.at(number_of_block);
        auto end_of_block = cumulative_sum_.at(number_of_block + 1);

        std::vector<size_t> answer(end_of_block - start_of_block);
        std::iota(answer.begin(), answer.end(), start_of_block);

        return answer;
    }
}

const std::vector<SSquaredState>& SSquaredConverter::block_with_number(size_t number) const {
    auto it =
        std::lower_bound(cumulative_sum_.begin(), cumulative_sum_.end(), number);
    auto number_of_block = std::distance(cumulative_sum_.begin(), it);

    return states_.at(number_of_block);
}

std::shared_ptr<const OrderOfSummation> SSquaredConverter::getOrderOfSummation() const {
    return order_of_summation_;
}

double SSquaredConverter::total_CG_coefficient(
    const SSquaredState& s_squared_state,
    const std::vector<double>& projections) const {
    double c = 1;

    for (const auto& instruction : *order_of_summation_) {
        size_t pos_one = instruction.positions_of_summands[0];
        size_t pos_two = instruction.positions_of_summands[1];
        size_t pos_sum = instruction.position_of_sum;

        double spin_one = s_squared_state.getSpin(pos_one);
        double spin_two = s_squared_state.getSpin(pos_two);
        double spin_sum = s_squared_state.getSpin(pos_sum);

        double proj_one = projections[pos_one];
        double proj_two = projections[pos_two];
        c *= clebshGordanCalculator_
                 .clebsh_gordan_coefficient(spin_one, spin_two, spin_sum, proj_one, proj_two);
        if (c == 0.0) {
            return c;
        }
    }

    return c;
}
}  // namespace spin_algebra