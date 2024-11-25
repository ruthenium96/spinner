#ifndef SPINNER_S2TRANSFORMER_H
#define SPINNER_S2TRANSFORMER_H

#include <numeric>
#include <utility>

#include "src/common/index_converter/lexicographic/IndexConverter.h"
#include "src/entities/matrix/Matrix.h"
#include "src/space/Space.h"
#include "src/spin_algebra/ClebshGordanCalculator.h"
#include "src/spin_algebra/RepresentationsMultiplier.h"
#include "src/spin_algebra/SSquaredConverter.h"
#include "src/spin_algebra/SSquaredState.h"

namespace space::optimization {
class S2Transformer {
  public:
    S2Transformer(
        std::shared_ptr<const index_converter::lexicographic::IndexConverter> converter,
        quantum::linear_algebra::FactoriesList factories,
        std::shared_ptr<const spin_algebra::SSquaredConverter> ssquared_converter,
        bool ito = false);

    space::Space apply(space::Space&& space) const;

  private:
    // TODO: move these functions to another class?
    std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>
    constructTransformationMatrix(
        const std::vector<spin_algebra::SSquaredState>& history,
        const space::Subspace& subspace) const;
    std::vector<double> construct_projections(uint32_t lex_index) const;

    std::shared_ptr<const index_converter::lexicographic::IndexConverter> converter_;
    const quantum::linear_algebra::FactoriesList factories_;
    bool ito_;

    std::shared_ptr<const spin_algebra::SSquaredConverter> ssquared_converter_;
};
}  // namespace space::optimization
#endif  //SPINNER_S2TRANSFORMER_H
