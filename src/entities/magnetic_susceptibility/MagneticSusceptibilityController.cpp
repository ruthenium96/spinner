#include "MagneticSusceptibilityController.h"

namespace magnetic_susceptibility {

MagneticSusceptibilityController::MagneticSusceptibilityController(
    std::unique_ptr<worker::AbstractWorker>&& worker) :
    worker_(std::move(worker)) {}

double MagneticSusceptibilityController::calculateResidualError() const {
    return worker_->getExperimentalValuesWorker()->calculateResidualError();
}

std::vector<ValueAtTemperature> MagneticSusceptibilityController::getTheoreticalValues() const {
    return worker_->getExperimentalValuesWorker()->getTheoreticalValues();
}

double MagneticSusceptibilityController::multiplyExperimentalAndTheoreticalDerivatives(
    std::vector<ValueAtTemperature> theoretical_derivative) const {
    std::vector<ValueAtTemperature> experimental_derivative =
        worker_->getExperimentalValuesWorker()->calculateDerivative();
    if (experimental_derivative.size() != theoretical_derivative.size()) {
        throw std::invalid_argument(
            "experimental_derivative.size() != theoretical_derivative.size()");
    }
    // just dot product:
    double derivative = 0;
    for (size_t i = 0; i < experimental_derivative.size(); ++i) {
        if (experimental_derivative[i].temperature != theoretical_derivative[i].temperature) {
            throw std::invalid_argument(
                "experimental_derivative[i].temperature != theoretical_derivative[i].temperature");
        }
        derivative += experimental_derivative[i].value * theoretical_derivative[i].value;
    }
    return derivative;
}

void MagneticSusceptibilityController::initializeExperimentalValues(
    const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker) {
    worker_->setExperimentalValuesWorker(experimental_values_worker);
    std::vector<double> temperatures = worker_->getExperimentalValuesWorker()->getTemperatures();
    std::vector<ValueAtTemperature> theoretical_mu_squared_values(temperatures.size());
    for (size_t i = 0; i < temperatures.size(); ++i) {
        theoretical_mu_squared_values[i] = {
            temperatures[i],
            worker_->calculateTheoreticalMuSquared(temperatures[i])};
    }
    worker_->getExperimentalValuesWorker()->setTheoreticalMuSquared(theoretical_mu_squared_values);
}

double MagneticSusceptibilityController::calculateTheoreticalMuSquared(double temperature) const {
    return worker_->calculateTheoreticalMuSquared(temperature);
}

double MagneticSusceptibilityController::calculateTotalDerivative(
    model::symbols::SymbolTypeEnum symbol_type,
    std::map<common::QuantityEnum, std::unique_ptr<quantum::linear_algebra::AbstractVector>>
        values_derivatives_map) const {
    std::vector<ValueAtTemperature> theoretical_derivative =
        worker_->calculateDerivative(symbol_type, std::move(values_derivatives_map));
    return multiplyExperimentalAndTheoreticalDerivatives(std::move(theoretical_derivative));
}

}  // namespace magnetic_susceptibility