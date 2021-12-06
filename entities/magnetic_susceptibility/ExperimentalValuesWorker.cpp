#include "ExperimentalValuesWorker.h"

#include <cmath>
#include <stdexcept>

namespace magnetic_susceptibility {

ExperimentalValuesWorker::ExperimentalValuesWorker(
    std::vector<ValueAtTemperature> experimental_values,
    ExperimentalValuesEnum experimental_values_type,
    double number_of_centers_ratio) {
    experimental_values_type_ = experimental_values_type;
    number_of_centers_ratio_ = number_of_centers_ratio;

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
        v.value *= number_of_centers_ratio_;
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
    theoretical_mu_squared_ = std::move(theoretical_mu_squared);
}
std::vector<ValueAtTemperature> ExperimentalValuesWorker::getTheoreticalValues() const {
    std::vector<ValueAtTemperature> theoretical_values = theoretical_mu_squared_;

    for (ValueAtTemperature& v : theoretical_values) {
        v.value /= number_of_centers_ratio_;
    }

    if (experimental_values_type_ == mu_squared_in_bohr_magnetons_squared) {
        // DO NOTHING
    } else if (experimental_values_type_ == mu_in_bohr_magnetons) {
        for (ValueAtTemperature& v : theoretical_values) {
            v.value = sqrt(v.value);
        }
    } else if (experimental_values_type_ == chiT_in_cm_cubed_kelvin_per_mol) {
        for (ValueAtTemperature& v : theoretical_values) {
            v.value *= mu_squared_in_bohr_magnetons_squared_to_chiT_in_cm_cubed_kelvin_per_mol;
        }
    }
    return theoretical_values;
}
}  // namespace magnetic_susceptibility