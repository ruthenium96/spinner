#include "ExperimentalValuesWorker.h"

#include <stdexcept>

namespace magnetic_susceptibility {

ExperimentalValuesWorker::ExperimentalValuesWorker(
    std::vector<ValueAtTemperature> experimental_values,
    ExperimentalValuesEnum experimental_values_type,
    double number_of_centers_ratio) {
    experimental_values_type_ = experimental_values_type;

    if (experimental_values_type == mu_squared_in_bohr_magnetons_squared) {
        // DO NOTHING
    } else if (experimental_values_type == mu_in_bohr_magnetons) {
        for (ValueAtTemperature& v : experimental_values) {
            v.value = v.value * v.value;
        }
    } else if (experimental_values_type == chiT_in_cm_cubed_kelvin_per_mol) {
        for (ValueAtTemperature& v : experimental_values) {
            v.value /= mu_squared_in_bohr_magnetons_squared_to_chiT_in_cm_cubed_kelvin_per_mol;
        }
    }

    for (ValueAtTemperature& v : experimental_values) {
        v.value *= number_of_centers_ratio;
    }

    // TODO: sort by temperature
    experimental_mu_squared_ = std::move(experimental_values);
}

double ExperimentalValuesWorker::calculateResidualError() const {
    double residual_error = 0;
    for (size_t i = 0; i < theoretical_mu_squared_.size(); ++i) {
        double diff = theoretical_mu_squared_[i].value - experimental_mu_squared_[i].value;
        residual_error += diff * diff;
    }
    return residual_error;
}

std::vector<ValueAtTemperature> ExperimentalValuesWorker::calculateDerivative() const {
    std::vector<ValueAtTemperature> derivative(experimental_mu_squared_.size());
    for (size_t i = 0; i < theoretical_mu_squared_.size(); ++i) {
        double diff = theoretical_mu_squared_[i].value - experimental_mu_squared_[i].value;
        derivative[i] = {theoretical_mu_squared_[i].temperature, 2 * diff};
    }
    return derivative;
}

std::vector<double> ExperimentalValuesWorker::getTemperatures() const {
    std::vector<double> temperatures(experimental_mu_squared_.size());
    for (size_t i = 0; i < temperatures.size(); ++i) {
        temperatures[i] = experimental_mu_squared_[i].temperature;
    }
    return temperatures;
}

void ExperimentalValuesWorker::setTheoreticalMuSquared(
    std::vector<ValueAtTemperature> theoretical_mu_squared) {
    if (!theoretical_mu_squared_.empty()) {
        throw std::invalid_argument("Theoretical values have been already set.");
    }
    theoretical_mu_squared_ = std::move(theoretical_mu_squared);
}
}  // namespace magnetic_susceptibility