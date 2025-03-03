#ifndef SPINNER_SSQUAREDLEVELANDREPRESENTATIONS_H
#define SPINNER_SSQUAREDLEVELANDREPRESENTATIONS_H

#include <memory>
#include <optional>
#include <vector>

#include "src/common/index_converter/s_squared/Level.h"

namespace spin_algebra {

class SSquaredLevelAndRepresentations : public index_converter::s_squared::Level {
  private:
    // first vector over different centers, second -- over different groups
    std::vector<std::vector<std::optional<uint8_t>>> intermediateRepresentations_;
  public:
    SSquaredLevelAndRepresentations(
        std::shared_ptr<std::vector<Multiplicity>> initialMultiplicities,
        size_t number_of_summations,
        size_t number_of_representations = 0);
    void setRepresentations(size_t number, std::vector<std::optional<uint8_t>> representation);

    struct Properties {
        Multiplicity multiplicity;
        std::vector<uint8_t> representations;
        std::strong_ordering operator<=>(const Properties& rhs) const = default;
    };

    double getSpin(size_t number) const;
    const std::vector<std::optional<uint8_t>>& getRepresentations(size_t number) const;
    Properties back() const;
};

}  // namespace spin_algebra

#endif  //SPINNER_SSQUAREDLEVELANDREPRESENTATIONS_H
