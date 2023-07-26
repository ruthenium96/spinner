#include "SSquaredState.h"

#include <cassert>
#include <iostream>
#include <numeric>
#include <utility>

namespace spin_algebra {

SSquaredState::SSquaredState(
    std::shared_ptr<std::vector<Multiplicity>> initialMultiplicities,
    size_t number_of_summations) {
    initialMultiplicities_ = std::move(initialMultiplicities);
    intermediateMultiplicities_.resize(number_of_summations, 0);
}

Multiplicity SSquaredState::getMultiplicity(size_t number) const {
    size_t initial_multiplicities_size = initialMultiplicities_->size();
    if (number < initial_multiplicities_size) {
        return initialMultiplicities_->at(number);
    } else {
        return intermediateMultiplicities_.at(number - initial_multiplicities_size);
    }
}

void SSquaredState::setMultiplicity(size_t number, Multiplicity multiplicity) {
    size_t initial_multiplicities_size = initialMultiplicities_->size();
    if (number < initial_multiplicities_size) {
        throw std::invalid_argument("Cannot set initial multiplicity");
    } else {
        size_t position_to_change = number - initial_multiplicities_size;
        if (intermediateMultiplicities_.at(position_to_change) == 0) {
            intermediateMultiplicities_[position_to_change] = multiplicity;
        } else {
            throw std::invalid_argument("Position was already set");
        }
    }
}

size_t SSquaredState::getSize() const {
    return initialMultiplicities_->size() + intermediateMultiplicities_.size();
}

Multiplicity SSquaredState::back() const {
    if (!intermediateMultiplicities_.empty()) {
        return intermediateMultiplicities_.back();
    } else {
        return initialMultiplicities_->back();
    }
}

std::map<Multiplicity, std::vector<SSquaredState>> SSquaredState::addAllMultiplicitiesAndSort(
    const std::vector<Multiplicity>& multiplicities_to_sum,
    const std::shared_ptr<const std::vector<AdditionInstruction>>& order_of_summation) {
    auto empty_history = SSquaredState(
        std::make_shared<std::vector<Multiplicity>>(multiplicities_to_sum),
        order_of_summation->size());

    std::vector<SSquaredState> result_of_summation = {empty_history};

    for (const auto& instruction : *order_of_summation) {
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
                temp_result.push_back(history);
                temp_result.back().setMultiplicity(pos_sum, mult_sum);
            }
        }
        std::swap(result_of_summation, temp_result);
        temp_result.clear();
    }

    auto result_of_sum_and_sort = std::map<Multiplicity, std::vector<SSquaredState>>();

    for (auto history : result_of_summation) {
        // TODO: implement .final() instead of .back()
        //  currently, position_of_sum of the last addition should point to
        //  the last element of intermediateMultiplicities_
        auto final_multiplicity = history.back();

        if (result_of_sum_and_sort.count(final_multiplicity) == 0) {
            result_of_sum_and_sort[final_multiplicity] = std::vector<SSquaredState>();
        }
        result_of_sum_and_sort[final_multiplicity].emplace_back(std::move(history));
    }

    return result_of_sum_and_sort;
}

}  // namespace spin_algebra