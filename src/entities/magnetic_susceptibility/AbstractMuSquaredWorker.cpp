#include "AbstractMuSquaredWorker.h"

namespace magnetic_susceptibility {

double AbstractMuSquaredWorker::calculateResidualError() const {
    return getExperimentalValuesWorker()->calculateResidualError();
}

std::vector<ValueAtTemperature> AbstractMuSquaredWorker::getTheoreticalValues() const {
    return getExperimentalValuesWorker()->getTheoreticalValues();
}

double AbstractMuSquaredWorker::multiplyExperimentalAndTheoreticalDerivatives(
    std::vector<ValueAtTemperature> theoretical_derivative) const {
    std::vector<ValueAtTemperature> experimental_derivative =
        getExperimentalValuesWorker()->calculateDerivative();
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

double AbstractMuSquaredWorker::calculateTotalDerivative(
    model::symbols::SymbolTypeEnum symbol_type,
    std::unique_ptr<quantum::linear_algebra::AbstractVector>&& derivative_value) const {
    std::vector<ValueAtTemperature> theoretical_derivative =
        calculateDerivative(symbol_type, std::move(derivative_value));
    return multiplyExperimentalAndTheoreticalDerivatives(std::move(theoretical_derivative));
}

double AbstractMuSquaredWorker::calculateTotalDerivative(
    model::symbols::SymbolTypeEnum symbol_type) const {
    std::vector<ValueAtTemperature> theoretical_derivative = calculateDerivative(symbol_type);
    return multiplyExperimentalAndTheoreticalDerivatives(std::move(theoretical_derivative));
}

void AbstractMuSquaredWorker::initializeExperimentalValues(
    const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker) {
    setExperimentalValuesWorker(experimental_values_worker);
    std::vector<double> temperatures = getExperimentalValuesWorker()->getTemperatures();
    std::vector<ValueAtTemperature> theoretical_mu_squared_values(temperatures.size());
    for (size_t i = 0; i < temperatures.size(); ++i) {
        theoretical_mu_squared_values[i] = {
            temperatures[i],
            calculateTheoreticalMuSquared(temperatures[i])};
    }
    getExperimentalValuesWorker()->setTheoreticalMuSquared(theoretical_mu_squared_values);
}

}  // namespace magnetic_susceptibility