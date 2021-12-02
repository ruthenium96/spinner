#include "MuSquaredWorker.h"

#include <algorithm>

namespace magnetic_susceptibility {

MuSquaredWorker::MuSquaredWorker(DenseVector&& energy, DenseVector&& degeneracy) :
    ensemble_averager_(std::move(energy), std::move(degeneracy)) {}

void MuSquaredWorker::initializeExperimentalValues(
    std::vector<ValueAtTemperature> experimental_values,
    ExperimentalValuesEnum experimental_quantity_type,
    double number_of_centers_ratio) {
    if (experimental_values_worker_.has_value()) {
        throw std::invalid_argument("Experimantal values have been already initialized");
    }

    experimental_values_worker_ = ExperimentalValuesWorker(
        std::move(experimental_values),
        experimental_quantity_type,
        number_of_centers_ratio);

    std::vector<double> temperatures = experimental_values_worker_->getTemperatures();
    std::vector<ValueAtTemperature> theoretical_mu_squared_values(temperatures.size());
    for (size_t i = 0; i < temperatures.size(); ++i) {
        theoretical_mu_squared_values[i] = {
            temperatures[i],
            theory_at_temperature(temperatures[i])};
    }
    experimental_values_worker_->setTheoreticalMuSquared(theoretical_mu_squared_values);
}

double MuSquaredWorker::calculateResidualError() const {
    return experimental_values_worker_->calculateResidualError();
}

double MuSquaredWorker::calculateTotalDerivative(
    symbols::SymbolTypeEnum symbol_type,
    DenseVector&& derivative_value) const {
    std::vector<ValueAtTemperature> experimantal_derivative =
        experimental_values_worker_->calculateDerivative();
    std::vector<ValueAtTemperature> theoretical_derivative =
        calculateDerivative(symbol_type, std::move(derivative_value));
    if (experimantal_derivative.size() != theoretical_derivative.size()) {
        throw std::invalid_argument(
            "experimantal_derivative.size() != theoretical_derivative.size()");
    }
    double derivative = 0;
    for (size_t i = 0; i < experimantal_derivative.size(); ++i) {
        if (experimantal_derivative[i].temperature != theoretical_derivative[i].temperature) {
            throw std::invalid_argument(
                "experimantal_derivative[i].temperature != theoretical_derivative[i].temperature");
        }
        derivative += experimantal_derivative[i].value * theoretical_derivative[i].value;
    }
    return derivative;
}

}  // namespace magnetic_susceptibility