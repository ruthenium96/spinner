#ifndef SPINNER_SSQUAREDSTATE_H
#define SPINNER_SSQUAREDSTATE_H

#include <map>
#include <memory>
#include <vector>

#include "MultiplicityDirectSum.h"
#include "OrderOfSummation.h"
#include "RepresentationsMultiplier.h"

namespace spin_algebra {

class SSquaredState {
  private:
    std::shared_ptr<std::vector<Multiplicity>> initialMultiplicities_;
    std::vector<Multiplicity> intermediateMultiplicities_;
    // first vector over different centers, second -- over different groups
    std::vector<std::vector<std::optional<uint8_t>>> intermediateRepresentations_;
    SSquaredState(
        std::shared_ptr<std::vector<Multiplicity>> initialMultiplicities,
        size_t number_of_summations,
        size_t number_of_representations = 0);
    void setMultiplicity(size_t number, Multiplicity multiplicity);
    void setRepresentations(size_t number, std::vector<std::optional<uint8_t>> representation);

  public:
    struct Properties {
        Multiplicity multiplicity;
        std::vector<uint8_t> representations;
        std::strong_ordering operator<=>(const Properties& rhs) const = default;
    };

    Multiplicity getMultiplicity(size_t number) const;
    const std::vector<std::optional<uint8_t>>& getRepresentations(size_t number) const;
    size_t getSize() const;
    Properties back() const;
    static std::map<Properties, std::vector<SSquaredState>> addAllMultiplicitiesAndSort(
        const std::vector<Multiplicity>& multiplicities_to_sum,
        const std::shared_ptr<const OrderOfSummation>& order_of_summation,
        const spin_algebra::RepresentationsMultiplier& representationsMultiplier);
};

}  // namespace spin_algebra

#endif  //SPINNER_SSQUAREDSTATE_H
