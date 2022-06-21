#ifndef SPINNER_UNIQUEGONLYSSQUAREDMUSQUAREDWORKER_H
#define SPINNER_UNIQUEGONLYSSQUAREDMUSQUAREDWORKER_H

#include "MuSquaredWorker.h"

namespace magnetic_susceptibility {

// mu^2 = mu_B^2 * g_{iso}^2 * <S^2_{total}>
class UniqueGOnlySSquaredMuSquaredWorker: public MuSquaredWorker {
  public:
    UniqueGOnlySSquaredMuSquaredWorker(
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& energy,
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& degeneracy,
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& s_squared,
        double g_unique);

    double calculateTheoreticalMuSquared(double temperature) const override;

  private:
    double g_unique_;
    std::unique_ptr<quantum::linear_algebra::AbstractVector> s_squared_;
    std::vector<ValueAtTemperature> calculateDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& derivative_value) const override;
    std::vector<ValueAtTemperature>
    calculateDerivative(model::symbols::SymbolTypeEnum symbol_type) const override;
};

}  // namespace magnetic_susceptibility

#endif  //SPINNER_UNIQUEGONLYSSQUAREDMUSQUAREDWORKER_H
