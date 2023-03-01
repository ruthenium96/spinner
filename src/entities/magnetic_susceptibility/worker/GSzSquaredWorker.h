#ifndef SPINNER_GSZSQUAREDWORKER_H
#define SPINNER_GSZSQUAREDWORKER_H

#include "BasicWorker.h"

namespace magnetic_susceptibility::worker {

// mu^2 = mu_B^2 * <(\sum_i g_i S_{iz})^2>
class GSzSquaredWorker: public BasicWorker {
  public:
    GSzSquaredWorker(
        std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& energy,
        std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& degeneracy,
        std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& g_sz_squared);

    double calculateTheoreticalMuSquared(double temperature) const override;
    std::vector<ValueAtTemperature> calculateDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        std::
            map<common::QuantityEnum, std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>
                values_derivatives_map) const override;

  private:
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector> g_sz_squared_;
};
}  // namespace magnetic_susceptibility::worker
#endif  //SPINNER_GSZSQUAREDWORKER_H
