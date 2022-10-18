#include "CurieWeissMuSquaredWorker.h"

#include <utility>

namespace magnetic_susceptibility {

CurieWeissMuSquaredWorker::CurieWeissMuSquaredWorker(
    std::unique_ptr<AbstractMuSquaredWorker>&& muSquaredWorker,
    std::shared_ptr<const double> Theta) :
    muSquaredWorker_(std::move(muSquaredWorker)),
    Theta_(std::move(Theta)) {}

double CurieWeissMuSquaredWorker::calculateTheoreticalMuSquared(double temperature) const {
    double theoreticalMuSquared = muSquaredWorker_->calculateTheoreticalMuSquared(temperature);
    theoreticalMuSquared = temperature * theoreticalMuSquared / (temperature - *Theta_);
    return theoreticalMuSquared;
}

std::shared_ptr<ExperimentalValuesWorker> CurieWeissMuSquaredWorker::getExperimentalValuesWorker() {
    return muSquaredWorker_->getExperimentalValuesWorker();
}

std::shared_ptr<const ExperimentalValuesWorker>
CurieWeissMuSquaredWorker::getExperimentalValuesWorker() const {
    return muSquaredWorker_->getExperimentalValuesWorker();
}

void CurieWeissMuSquaredWorker::setExperimentalValuesWorker(
    const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker) {
    muSquaredWorker_->setExperimentalValuesWorker(experimental_values_worker);
}

std::vector<ValueAtTemperature> CurieWeissMuSquaredWorker::calculateDerivative(
    model::symbols::SymbolTypeEnum symbol_type,
    std::unique_ptr<quantum::linear_algebra::AbstractVector>&& derivative_value) const {
    return muSquaredWorker_->calculateDerivative(symbol_type, std::move(derivative_value));
}

std::vector<ValueAtTemperature>
CurieWeissMuSquaredWorker::calculateDerivative(model::symbols::SymbolTypeEnum symbol_type) const {
    if (symbol_type == model::symbols::Theta) {
        std::vector<double> temperatures = getExperimentalValuesWorker()->getTemperatures();
        std::vector<ValueAtTemperature> theoretical_mu_squared_values(temperatures.size());
        for (size_t i = 0; i < temperatures.size(); ++i) {
            double value = temperatures[i]
                * muSquaredWorker_->calculateTheoreticalMuSquared(temperatures[i])
                / ((temperatures[i] - *Theta_) * (temperatures[i] - *Theta_));
            theoretical_mu_squared_values[i] = {temperatures[i], value};
        }
        return theoretical_mu_squared_values;
    } else {
        return muSquaredWorker_->calculateDerivative(symbol_type);
    }
}

}  // namespace magnetic_susceptibility