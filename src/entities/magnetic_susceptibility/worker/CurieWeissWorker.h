#ifndef SPINNER_CURIEWEISSWORKER_H
#define SPINNER_CURIEWEISSWORKER_H

#include <optional>

#include "src/entities/magnetic_susceptibility/worker/AbstractWorker.h"

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

  protected:
    std::unique_ptr<AbstractWorker> worker_;
    std::shared_ptr<const double> Theta_;

    std::vector<ValueAtTemperature> calculateDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& derivative_value) const override;
    std::vector<ValueAtTemperature>
    calculateDerivative(model::symbols::SymbolTypeEnum symbol_type) const override;
};

}  // namespace magnetic_susceptibility

#endif  //SPINNER_CURIEWEISSWORKER_H
