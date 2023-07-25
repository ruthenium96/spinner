#include "SSquaredState.h"

#include <cassert>
#include <iostream>
#include <numeric>
#include <utility>

namespace spin_algebra {

SSquaredState::SSquaredState(std::shared_ptr<std::vector<Multiplicity>> initialMultiplicities) {
    initialMultiplicities_ = std::move(initialMultiplicities);
}

void SSquaredState::pushNewIntermediateMultiplicity(Multiplicity multiplicity) {
    intermediateMultiplicities_.push_back(multiplicity);
}

Multiplicity SSquaredState::getMultiplicity(size_t number) const {
    size_t initial_multiplicities_size = initialMultiplicities_->size();
    if (number < initial_multiplicities_size) {
        return initialMultiplicities_->at(number);
    } else {
        return intermediateMultiplicities_.at(number - initial_multiplicities_size);
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
    auto empty_history =
        SSquaredState(std::make_shared<std::vector<Multiplicity>>(multiplicities_to_sum));

    std::vector<SSquaredState> result_of_summation = {empty_history};

    for (const auto& instruction : *order_of_summation) {
        // we currently do not support cases with three or more summands:
        assert(instruction.size() == 2);
        size_t pos_one = instruction[0];
        size_t pos_two = instruction[1];
        std::vector<SSquaredState> temp_result;
        for (const auto& history : result_of_summation) {
            MultiplicityDirectSum mult_one(history.getMultiplicity(pos_one));
            MultiplicityDirectSum mult_two(history.getMultiplicity(pos_two));
            auto mult_direct_sum = mult_one * mult_two;
            for (auto mult_sum : mult_direct_sum.getMultiplicities()) {
                temp_result.push_back(history);
                temp_result.back().pushNewIntermediateMultiplicity(mult_sum);
            }
        }
        std::swap(result_of_summation, temp_result);
        temp_result.clear();
    }

    auto result_of_sum_and_sort = std::map<Multiplicity, std::vector<SSquaredState>>();

    for (auto history : result_of_summation) {
        auto final_multiplicity = history.back();

        if (result_of_sum_and_sort.count(final_multiplicity) == 0) {
            result_of_sum_and_sort[final_multiplicity] = std::vector<SSquaredState>();
        }
        result_of_sum_and_sort[final_multiplicity].emplace_back(std::move(history));
    }

    return result_of_sum_and_sort;
}

}  // namespace spin_algebra