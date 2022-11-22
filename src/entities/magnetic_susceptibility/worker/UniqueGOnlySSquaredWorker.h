#ifndef SPINNER_UNIQUEGONLYSSQUAREDWORKER_H
#define SPINNER_UNIQUEGONLYSSQUAREDWORKER_H

#include "BasicWorker.h"

namespace magnetic_susceptibility::worker {

// mu^2 = mu_B^2 * g_{iso}^2 * <S^2_{total}>
class UniqueGOnlySSquaredWorker: public BasicWorker {
  public:
    UniqueGOnlySSquaredWorker(
        std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& energy,
        std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& degeneracy,
        std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& s_squared,
        double g_unique);

    double calculateTheoreticalMuSquared(double temperature) const override;
    std::vector<ValueAtTemperature> calculateDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        std::
            map<common::QuantityEnum, std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>
                values_derivatives_map) const override;

  private:
    double g_unique_;
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector> s_squared_;
};

}  // namespace magnetic_susceptibility

#endif  //SPINNER_UNIQUEGONLYSSQUAREDWORKER_H
