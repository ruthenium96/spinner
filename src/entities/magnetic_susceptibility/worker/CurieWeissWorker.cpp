#include "CurieWeissWorker.h"

#include <utility>

namespace magnetic_susceptibility::worker {

CurieWeissWorker::CurieWeissWorker(
    std::unique_ptr<AbstractWorker>&& worker,
    std::shared_ptr<const double> Theta) :
    worker_(std::move(worker)),
    Theta_(std::move(Theta)) {}

double CurieWeissWorker::calculateTheoreticalMuSquared(double temperature) const {
    double theoreticalMuSquared = worker_->calculateTheoreticalMuSquared(temperature);
    theoreticalMuSquared = temperature * theoreticalMuSquared / (temperature - *Theta_);
    return theoreticalMuSquared;
}

std::shared_ptr<ExperimentalValuesWorker> CurieWeissWorker::getExperimentalValuesWorker() {
    return worker_->getExperimentalValuesWorker();
}

std::shared_ptr<const ExperimentalValuesWorker>
CurieWeissWorker::getExperimentalValuesWorker() const {
    return worker_->getExperimentalValuesWorker();
}

void CurieWeissWorker::setExperimentalValuesWorker(
    const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker) {
    worker_->setExperimentalValuesWorker(experimental_values_worker);
}

std::vector<ValueAtTemperature> CurieWeissWorker::calculateDerivative(
    model::symbols::SymbolTypeEnum symbol_type,
    std::map<common::QuantityEnum, std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>
        values_derivatives_map) const {
    if (symbol_type == model::symbols::Theta) {
        std::vector<double> temperatures = getExperimentalValuesWorker()->getTemperatures();
        std::vector<ValueAtTemperature> theoretical_mu_squared_values(temperatures.size());
        for (size_t i = 0; i < temperatures.size(); ++i) {
            double value = temperatures[i] * worker_->calculateTheoreticalMuSquared(temperatures[i])
                / ((temperatures[i] - *Theta_) * (temperatures[i] - *Theta_));
            theoretical_mu_squared_values[i] = {temperatures[i], value};
        }
        return theoretical_mu_squared_values;
    } else {
        return worker_->calculateDerivative(symbol_type, std::move(values_derivatives_map));
    }
}

}  // namespace magnetic_susceptibility::worker