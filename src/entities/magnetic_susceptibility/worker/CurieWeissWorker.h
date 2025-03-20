#ifndef SPINNER_CURIEWEISSWORKER_H
#define SPINNER_CURIEWEISSWORKER_H

#include <optional>

#include "src/entities/magnetic_susceptibility/worker/AbstractWorker.h"
#include "src/model/symbols/SymbolName.h"

namespace magnetic_susceptibility::worker {

// It is a decorator. Takes concrete worker with mu^2(T) and returns T*mu^2(T)/(T-\Theta)
class CurieWeissWorker: public AbstractWorker {
  public:
    CurieWeissWorker(
        std::unique_ptr<AbstractWorker>&& muSquaredWorker,
        std::shared_ptr<const double> Theta);
    double calculateTheoreticalMuSquared(double temperature) const override;

    std::shared_ptr<ExperimentalValuesWorker> getExperimentalValuesWorker() override;
    std::shared_ptr<const ExperimentalValuesWorker> getExperimentalValuesWorker() const override;
    void setExperimentalValuesWorker(
        const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker) override;
    std::vector<ValueAtTemperature> calculateDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        model::symbols::SymbolName symbol_name) const override;         

  protected:
    std::unique_ptr<AbstractWorker> worker_;
    std::shared_ptr<const double> Theta_;
};

}  // namespace magnetic_susceptibility

#endif  //SPINNER_CURIEWEISSWORKER_H
