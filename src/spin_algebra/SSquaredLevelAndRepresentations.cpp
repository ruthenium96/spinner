#include "SSquaredLevelAndRepresentations.h"

#include <iostream>
#include <numeric>
#include <utility>
#include "src/common/index_converter/s_squared/Level.h"

namespace spin_algebra {

SSquaredLevelAndRepresentations::SSquaredLevelAndRepresentations(
    std::shared_ptr<std::vector<Multiplicity>> initialMultiplicities,
    size_t number_of_summations,
    size_t number_of_representations) : 
    index_converter::s_squared::Level(initialMultiplicities, number_of_summations) {
    intermediateRepresentations_ = std::vector<std::vector<std::optional<uint8_t>>>(
        getSize(),
        std::vector<std::optional<uint8_t>>(number_of_representations, std::nullopt));
}

SSquaredLevelAndRepresentations::Properties SSquaredLevelAndRepresentations::back() const {
    Properties properties;
    properties.multiplicity = total();
    for (const auto& mb_representaton : intermediateRepresentations_.back()) {
        properties.representations.push_back(mb_representaton.value());
    }
    return properties;
}

const std::vector<std::optional<uint8_t>>& SSquaredLevelAndRepresentations::getRepresentations(size_t number) const {
    return intermediateRepresentations_.at(number);
}

void SSquaredLevelAndRepresentations::setRepresentations(
    size_t number,
    std::vector<std::optional<uint8_t>> representation) {
    intermediateRepresentations_[number] = std::move(representation);
}

double SSquaredLevelAndRepresentations::getSpin(size_t number) const {
    return ((double)getMultiplicity(number) - 1.0) / 2.0;
}

}  // namespace spin_algebra