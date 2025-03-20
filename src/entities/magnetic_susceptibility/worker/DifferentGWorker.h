#ifndef SPINNER_GSZSQUAREDWORKER_H
#define SPINNER_GSZSQUAREDWORKER_H

#include "BasicWorker.h"

namespace magnetic_susceptibility::worker {

// mu^2 = 3 * mu_B^2 * <(\sum_i g_i S_{iz})^2> or 
// mu^2 = 3 * mu_B^2 * <(\sum_i g_i^2 T^{(0)}_0 + \sum_i \sum_j g_i g_j T^{(0)}_0)>
class DifferentGWorker: public BasicWorker {
  public:
    DifferentGWorker(
      std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra,
      common::QuantityEnum quantity_enum_for_averaging, 
      double quantity_factor);

    double calculateTheoreticalMuSquared(double temperature) const override;
    std::vector<ValueAtTemperature> calculateDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        std::
            map<common::QuantityEnum, std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>
                values_derivatives_map) const override;

  private:
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra_;
    common::QuantityEnum quantity_enum_for_averaging_;
    double quantity_factor_;
};
}  // namespace magnetic_susceptibility::worker
#endif  //SPINNER_GSZSQUAREDWORKER_H
