#ifndef SPINNER_SSQUAREDSTATE_H
#define SPINNER_SSQUAREDSTATE_H

#include <memory>
#include <vector>

#include "MultiplicityDirectSum.h"
#include "src/common/index_converter/s_squared/OrderOfSummation.h"
#include "RepresentationsMultiplier.h"

namespace spin_algebra {

class SSquaredState {
  private:
    std::shared_ptr<std::vector<Multiplicity>> initialMultiplicities_;
    std::vector<Multiplicity> intermediateMultiplicities_;
    // first vector over different centers, second -- over different groups
    std::vector<std::vector<std::optional<uint8_t>>> intermediateRepresentations_;
  public:
    SSquaredState(
        std::shared_ptr<std::vector<Multiplicity>> initialMultiplicities,
        size_t number_of_summations,
        size_t number_of_representations = 0);
    void setMultiplicity(size_t number, Multiplicity multiplicity);
    void setRepresentations(size_t number, std::vector<std::optional<uint8_t>> representation);

    struct Properties {
        Multiplicity multiplicity;
        std::vector<uint8_t> representations;
        std::strong_ordering operator<=>(const Properties& rhs) const = default;
    };

    Multiplicity getMultiplicity(size_t number) const;
    double getSpin(size_t number) const;
    const std::vector<std::optional<uint8_t>>& getRepresentations(size_t number) const;
    size_t getSize() const;
    Properties back() const;
};

}  // namespace spin_algebra

#endif  //SPINNER_SSQUAREDSTATE_H
