#include "RepresentationsMultiplier.h"

#include <utility>

namespace spin_algebra {
RepresentationsMultiplier::RepresentationsMultiplier(
    std::vector<std::map<std::pair<uint8_t, uint8_t>, std::set<uint8_t>>> all_cayley_tables) :
    all_cayley_tables_(std::move(all_cayley_tables)) {}

std::vector<std::optional<uint8_t>> RepresentationsMultiplier::multiplyRepresentations(
    const std::vector<std::optional<uint8_t>>& representations_one,
    const std::vector<std::optional<uint8_t>>& representations_two,
    const std::optional<size_t>& mb_number_of_group,
    Multiplicity mult_one,
    Multiplicity mult_two,
    Multiplicity mult_sum) const {
    auto representations_sum = std::vector<std::optional<uint8_t>>(representations_one.size());
    for (size_t number_of_group = 0; number_of_group < representations_one.size();
         ++number_of_group) {
        if (mb_number_of_group != number_of_group) {
            const auto& cayley_table = all_cayley_tables_[number_of_group];
            representations_sum[number_of_group] = outsideOrbitSummation(
                cayley_table,
                representations_one[number_of_group],
                representations_two[number_of_group]);
        } else {
            representations_sum[number_of_group] =
                insideOrbitSummation(mult_one, mult_two, mult_sum);
        }
    }
    return representations_sum;
}

uint32_t RepresentationsMultiplier::outsideOrbitSummation(
    const std::map<std::pair<uint8_t, uint8_t>, std::set<uint8_t>>& cayley_table,
    std::optional<uint32_t> single_representation_one,
    std::optional<uint32_t> single_representation_two) {
    // outside orbit summation
    if (single_representation_one.has_value() && single_representation_two.has_value()) {
        auto result_representation =
            cayley_table.at({single_representation_one.value(), single_representation_two.value()});
        // TODO: the first (and currently the only) representation
        return *result_representation.begin();
    } else if (single_representation_one.has_value()) {
        return single_representation_one.value();
    } else if (single_representation_two.has_value()) {
        return single_representation_two.value();
    } else {
        return 0;  // all-symmetric representation
    }
}

uint32_t RepresentationsMultiplier::insideOrbitSummation(
    Multiplicity mult_one,
    Multiplicity mult_two,
    Multiplicity mult_sum) {
    // inside orbit summation
    // TODO: It is correct only for P2 group:
    if (mult_one == mult_two) {
        // M1 = 2S1 + 1
        // M2 = 2S2 + 1
        // Mm = 2(S1+S2) + 1 = (2S1+1) + (2S2+1) - 1 = M1 + M2 - 1
        Multiplicity maxMultiplicity = mult_one + mult_two - 1;
        if ((maxMultiplicity - mult_sum) % 4 == 0) {
            return 0;
        } else {
            return 1;
        }
    } else if (mult_one > mult_two) {
        return 0;
    } else if (mult_one < mult_two) {
        return 1;
    }
}

size_t RepresentationsMultiplier::getNumberOfRepresentations() const {
    return all_cayley_tables_.size();
}
}  // namespace spin_algebra