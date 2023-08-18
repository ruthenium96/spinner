#include "SSquaredState.h"

#include <cassert>
#include <iostream>
#include <numeric>
#include <utility>

namespace spin_algebra {

SSquaredState::SSquaredState(
    std::shared_ptr<std::vector<Multiplicity>> initialMultiplicities,
    size_t number_of_summations,
    size_t number_of_representations) {
    initialMultiplicities_ = std::move(initialMultiplicities);
    intermediateMultiplicities_.resize(number_of_summations, 0);
    intermediateRepresentations_ = std::vector<std::vector<std::optional<uint8_t>>>(
        getSize(),
        std::vector<std::optional<uint8_t>>(number_of_representations, std::nullopt));
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

SSquaredState::Properties SSquaredState::back() const {
    Properties properties;
    properties.multiplicity = intermediateMultiplicities_.back();
    for (const auto& mb_representaton : intermediateRepresentations_.back()) {
        properties.representations.push_back(mb_representaton.value());
    }
    return properties;
}

std::map<SSquaredState::Properties, std::vector<SSquaredState>>
SSquaredState::addAllMultiplicitiesAndSort(
    const std::vector<Multiplicity>& multiplicities_to_sum,
    const std::shared_ptr<const OrderOfSummation>& order_of_summation,
    const spin_algebra::RepresentationsMultiplier& representationsMultiplier) {
    auto empty_history = SSquaredState(
        std::make_shared<std::vector<Multiplicity>>(multiplicities_to_sum),
        order_of_summation->size(),
        representationsMultiplier.getNumberOfRepresentations());

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
                auto representations_one = history.getRepresentations(pos_one);
                auto representations_two = history.getRepresentations(pos_two);
                auto representations_sum = representationsMultiplier.multiplyRepresentations(
                    representations_one,
                    representations_two,
                    instruction.number_of_group,
                    mult_one.getMultiplicities().at(0),
                    mult_two.getMultiplicities().at(0),
                    mult_sum);
                temp_result.push_back(history);
                temp_result.back().setMultiplicity(pos_sum, mult_sum);
                temp_result.back().setRepresentations(pos_sum, representations_sum);
            }
        }
        std::swap(result_of_summation, temp_result);
        temp_result.clear();
    }

    auto result_of_sum_and_sort = std::map<Properties, std::vector<SSquaredState>>();

    for (auto history : result_of_summation) {
        auto final_properties = history.back();

        if (result_of_sum_and_sort.count(final_properties) == 0) {
            result_of_sum_and_sort[final_properties] = std::vector<SSquaredState>();
        }
        result_of_sum_and_sort[final_properties].emplace_back(std::move(history));
    }

    return result_of_sum_and_sort;
}

const std::vector<std::optional<uint8_t>>& SSquaredState::getRepresentations(size_t number) const {
    return intermediateRepresentations_.at(number);
}

void SSquaredState::setRepresentations(
    size_t number,
    std::vector<std::optional<uint8_t>> representation) {
    intermediateRepresentations_[number] = std::move(representation);
}

}  // namespace spin_algebra