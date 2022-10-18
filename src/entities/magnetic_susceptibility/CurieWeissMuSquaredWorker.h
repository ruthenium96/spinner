#ifndef SPINNER_CURIEWEISSMUSQUAREDWORKER_H
#define SPINNER_CURIEWEISSMUSQUAREDWORKER_H

#include <optional>

#include "AbstractMuSquaredWorker.h"

namespace magnetic_susceptibility {

// It is a decorator. Takes concrete worker with mu^2(T) and returns T*mu^2(T)/(T-\Theta)
class CurieWeissMuSquaredWorker: public AbstractMuSquaredWorker {
  public:
    CurieWeissMuSquaredWorker(
        std::unique_ptr<AbstractMuSquaredWorker>&& muSquaredWorker,
        std::shared_ptr<const double> Theta);
    double calculateTheoreticalMuSquared(double temperature) const override;

    std::shared_ptr<ExperimentalValuesWorker> getExperimentalValuesWorker() override;
    std::shared_ptr<const ExperimentalValuesWorker> getExperimentalValuesWorker() const override;
    void setExperimentalValuesWorker(
        const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker) override;

  protected:
    std::unique_ptr<AbstractMuSquaredWorker> muSquaredWorker_;
    std::shared_ptr<const double> Theta_;

    std::vector<ValueAtTemperature> calculateDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& derivative_value) const override;
    std::vector<ValueAtTemperature>
    calculateDerivative(model::symbols::SymbolTypeEnum symbol_type) const override;
};

}  // namespace magnetic_susceptibility

#endif  //SPINNER_CURIEWEISSMUSQUAREDWORKER_H
