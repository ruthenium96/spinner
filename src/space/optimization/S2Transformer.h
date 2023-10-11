#ifndef SPINNER_S2TRANSFORMER_H
#define SPINNER_S2TRANSFORMER_H

#include <numeric>
#include <utility>

#include "src/common/lexicographic/IndexConverter.h"
#include "src/entities/matrix/Matrix.h"
#include "src/space/Space.h"
#include "src/spin_algebra/ClebshGordanCalculator.h"
#include "src/spin_algebra/RepresentationsMultiplier.h"
#include "src/spin_algebra/SSquaredState.h"

namespace space::optimization {
class S2Transformer {
  public:
    S2Transformer(
        lexicographic::IndexConverter converter,
        quantum::linear_algebra::FactoriesList factories,
        std::shared_ptr<const spin_algebra::OrderOfSummation> order_of_summation,
        const spin_algebra::RepresentationsMultiplier& representationsMultiplier);

    space::Space apply(space::Space&& space) const;

  private:
    // TODO: move these functions to another class?
    std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>
    constructTransformationMatrix(
        const std::vector<spin_algebra::SSquaredState>& history,
        const space::Subspace& subspace) const;
    double total_CG_coefficient(
        const spin_algebra::SSquaredState& s_squared_state,
        const std::vector<double>& projections) const;
    std::vector<double> construct_projections(uint32_t lex_index) const;

    const lexicographic::IndexConverter converter_;
    const quantum::linear_algebra::FactoriesList factories_;

    std::shared_ptr<const spin_algebra::OrderOfSummation> order_of_summation_;
    std::map<spin_algebra::SSquaredState::Properties, std::vector<spin_algebra::SSquaredState>>
        sorted_s_squared_states_;

    spin_algebra::ClebshGordanCalculator clebshGordanCalculator_;
};
}  // namespace space::optimization
#endif  //SPINNER_S2TRANSFORMER_H
