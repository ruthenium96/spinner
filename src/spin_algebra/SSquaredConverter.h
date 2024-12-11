#ifndef SPINNER_SSQUAREDCONVERTER_H
#define SPINNER_SSQUAREDCONVERTER_H

#include <map>

#include "ClebshGordanCalculator.h"
#include "SSquaredState.h"
#include "src/common/index_converter/s_squared/State.h"

namespace spin_algebra {

class SSquaredConverter {
  public:
    SSquaredConverter(const std::vector<Multiplicity>& multiplicities_to_sum,
                      const std::shared_ptr<const index_converter::s_squared::OrderOfSummation>& order_of_summation,
                      const spin_algebra::RepresentationsMultiplier& representationsMultiplier);

    const SSquaredState& at(size_t number) const;
    size_t number_in_block(size_t number) const;
    const std::vector<SSquaredState>& block_with_number(size_t number) const;

    std::shared_ptr<const index_converter::s_squared::OrderOfSummation> getOrderOfSummation() const;

    double total_CG_coefficient(
        const spin_algebra::SSquaredState& s_squared_state,
        const std::vector<double>& projections) const;

    std::optional<std::vector<size_t>> indexes_with_property(SSquaredState::Properties) const;
    std::optional<std::reference_wrapper<const std::vector<SSquaredState>>>
        block_with_property(SSquaredState::Properties) const;
  private:
    std::shared_ptr<const index_converter::s_squared::OrderOfSummation> order_of_summation_;
    std::vector<std::vector<SSquaredState>> states_;
    std::vector<size_t> cumulative_sum_;
    std::map<SSquaredState::Properties, size_t> properties_indexes_;
    size_t number_of_initial_mults_;
    size_t number_of_all_mults_;

    spin_algebra::ClebshGordanCalculator clebshGordanCalculator_;
};

}  // namespace spin_algebra

#endif  //SPINNER_SSQUAREDCONVERTER_H
