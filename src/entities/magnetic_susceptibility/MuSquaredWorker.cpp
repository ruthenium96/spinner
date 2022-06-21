#include "MuSquaredWorker.h"

#include <algorithm>

namespace magnetic_susceptibility {

MuSquaredWorker::MuSquaredWorker(
    std::unique_ptr<quantum::linear_algebra::AbstractVector>&& energy,
    std::unique_ptr<quantum::linear_algebra::AbstractVector>&& degeneracy) :
    ensemble_averager_(std::move(energy), std::move(degeneracy)) {}

double MuSquaredWorker::calculateResidualError() const {
    return experimental_values_worker_.value()->calculateResidualError();
}

double MuSquaredWorker::calculateTotalDerivative(
    model::symbols::SymbolTypeEnum symbol_type,
    std::unique_ptr<quantum::linear_algebra::AbstractVector>&& derivative_value) const {
    std::vector<ValueAtTemperature> theoretical_derivative =
        calculateDerivative(symbol_type, std::move(derivative_value));
    return multiplyExperimentalAndTheoreticalDerivatives(std::move(theoretical_derivative));
}

double MuSquaredWorker::calculateTotalDerivative(model::symbols::SymbolTypeEnum symbol_type) const {
    std::vector<ValueAtTemperature> theoretical_derivative = calculateDerivative(symbol_type);
    return multiplyExperimentalAndTheoreticalDerivatives(std::move(theoretical_derivative));
}

double MuSquaredWorker::multiplyExperimentalAndTheoreticalDerivatives(
    std::vector<ValueAtTemperature> theoretical_derivative) const {
    std::vector<ValueAtTemperature> experimental_derivative =
        experimental_values_worker_.value()->calculateDerivative();
    if (experimental_derivative.size() != theoretical_derivative.size()) {
        throw std::invalid_argument(
            "experimantal_derivative.size() != theoretical_derivative.size()");
    }
    double derivative = 0;
    for (size_t i = 0; i < experimental_derivative.size(); ++i) {
        if (experimental_derivative[i].temperature != theoretical_derivative[i].temperature) {
            throw std::invalid_argument(
                "experimantal_derivative[i].temperature != theoretical_derivative[i].temperature");
        }
        derivative += experimental_derivative[i].value * theoretical_derivative[i].value;
    }
    return derivative;
}

std::vector<ValueAtTemperature> MuSquaredWorker::getTheoreticalValues() const {
    return experimental_values_worker_.value()->getTheoreticalValues();
}

void MuSquaredWorker::initializeExperimentalValues(
    const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker) {
    experimental_values_worker_ = experimental_values_worker;
    std::vector<double> temperatures = experimental_values_worker_.value()->getTemperatures();
    std::vector<ValueAtTemperature> theoretical_mu_squared_values(temperatures.size());
    for (size_t i = 0; i < temperatures.size(); ++i) {
        theoretical_mu_squared_values[i] = {
            temperatures[i],
            calculateTheoreticalMuSquared(temperatures[i])};
    }
    experimental_values_worker_.value()->setTheoreticalMuSquared(theoretical_mu_squared_values);
}

}  // namespace magnetic_susceptibility