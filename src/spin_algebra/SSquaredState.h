#ifndef SPINNER_SSQUAREDSTATE_H
#define SPINNER_SSQUAREDSTATE_H

#include <map>
#include <memory>
#include <vector>

#include "MultiplicityDirectSum.h"

namespace spin_algebra {

struct AdditionInstruction {
    std::vector<uint64_t> positions_of_summands;
    uint64_t position_of_sum;
};

// TODO: refactor this class.
class SSquaredState {
  private:
    std::shared_ptr<std::vector<Multiplicity>> initialMultiplicities_;
    std::vector<Multiplicity> intermediateMultiplicities_;
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
        const std::shared_ptr<const std::vector<AdditionInstruction>>& order_of_summation);
};

}  // namespace spin_algebra

#endif  //SPINNER_SSQUAREDSTATE_H
