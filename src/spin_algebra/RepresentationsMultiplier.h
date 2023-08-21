#ifndef SPINNER_REPRESENTATIONSMULTIPLIER_H
#define SPINNER_REPRESENTATIONSMULTIPLIER_H

#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <vector>

#include "src/spin_algebra/Multiplicity.h"

namespace spin_algebra {

class RepresentationsMultiplier {
  public:
    RepresentationsMultiplier() = default;
    explicit RepresentationsMultiplier(
        std::vector<std::map<std::pair<uint8_t, uint8_t>, std::set<uint8_t>>> all_cayley_tables);
    std::vector<std::optional<uint8_t>> multiplyRepresentations(
        const std::vector<std::optional<uint8_t>>& representation_one,
        const std::vector<std::optional<uint8_t>>& representation_two,
        const std::optional<size_t>& mb_number_of_group,
        Multiplicity mult_one,
        Multiplicity mult_two,
        Multiplicity mult_sum) const;
    size_t getNumberOfRepresentations() const;

  private:
    std::vector<std::map<std::pair<uint8_t, uint8_t>, std::set<uint8_t>>> all_cayley_tables_;
    // TODO: move it to group module?
    static uint32_t outsideOrbitSummation(
        const std::map<std::pair<uint8_t, uint8_t>, std::set<uint8_t>>& cayley_table,
        std::optional<uint32_t> single_representation_one,
        std::optional<uint32_t> single_representation_two);
    static uint32_t
    insideOrbitSummation(Multiplicity mult_one, Multiplicity mult_two, Multiplicity mult_sum);
};

}  // namespace spin_algebra

#endif  //SPINNER_REPRESENTATIONSMULTIPLIER_H
