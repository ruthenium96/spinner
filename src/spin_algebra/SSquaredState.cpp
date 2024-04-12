#include "SSquaredState.h"

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

const std::vector<std::optional<uint8_t>>& SSquaredState::getRepresentations(size_t number) const {
    return intermediateRepresentations_.at(number);
}

void SSquaredState::setRepresentations(
    size_t number,
    std::vector<std::optional<uint8_t>> representation) {
    intermediateRepresentations_[number] = std::move(representation);
}

}  // namespace spin_algebra