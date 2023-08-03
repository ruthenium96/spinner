#ifndef SPINNER_SSQUAREDSTATE_H
#define SPINNER_SSQUAREDSTATE_H

#include <map>
#include <memory>
#include <vector>

#include "MultiplicityDirectSum.h"
#include "OrderOfSummation.h"

namespace spin_algebra {

class SSquaredState {
  private:
    std::shared_ptr<std::vector<Multiplicity>> initialMultiplicities_;
    std::vector<Multiplicity> intermediateMultiplicities_;
    // TODO: std::vector of representations?
    SSquaredState(
        std::shared_ptr<std::vector<Multiplicity>> initialMultiplicities,
        size_t number_of_summations);
    void setMultiplicity(size_t number, Multiplicity multiplicity);

  public:
    Multiplicity getMultiplicity(size_t number) const;
    size_t getSize() const;
    Multiplicity back() const;
    static std::map<Multiplicity, std::vector<SSquaredState>> addAllMultiplicitiesAndSort(
        const std::vector<Multiplicity>& multiplicities_to_sum,
        const std::shared_ptr<const OrderOfSummation>& order_of_summation);
};

}  // namespace spin_algebra

#endif  //SPINNER_SSQUAREDSTATE_H
