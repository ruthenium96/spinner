#ifndef SPINNER_UNIQUEGWORKER_H
#define SPINNER_UNIQUEGWORKER_H

#include "BasicWorker.h"
#include "src/common/Quantity.h"
#include "src/model/symbols/SymbolName.h"

namespace magnetic_susceptibility::worker {

// mu^2 = mu_B^2 * g_{iso}^2 * <S^2_{total}> or
// mu^2 = 3 * mu_B^2 * g_{iso}^2 * <M^2_{total}>
class UniqueGWorker: public BasicWorker {
  public:
    UniqueGWorker(
        std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra,
        double g_unique, 
        common::QuantityEnum quantity_enum_for_averaging, 
        double quantity_factor);

    double calculateTheoreticalMuSquared(double temperature) const override;
    std::vector<ValueAtTemperature> calculateDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        model::symbols::SymbolName symbol_name) const override;

  private:
    double g_unique_;
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra_;
    double quantity_factor_;
    common::QuantityEnum quantity_enum_for_averaging_;
};

}  // namespace magnetic_susceptibility

#endif  //SPINNER_UNIQUEGWORKER_H
