#include "ChiT.h"

namespace magnetic_susceptibility {
ChiT::ChiT(double g_factor, DenseVector energy, DenseVector degeneracy, DenseVector s_squared) :
    g_factor_(g_factor),
    energy_(std::move(energy)),
    s_squared_(std::move(s_squared)),
    degeneracy_(std::move(degeneracy)) {}

double ChiT::at_temperature(double temperature) {
    double s_squared_ensemble_average =
        ensemble_average(temperature, energy_, degeneracy_, s_squared_);
    return s_squared_ensemble_average * bohr_magneton * bohr_magneton * g_factor_ * g_factor_ / 3;
}

double ChiT::ensemble_average(
    double temperature,
    const DenseVector& energy,
    const DenseVector& degeneracy,
    const DenseVector& value) {
    DenseVector divided_and_wise_exped_energy = energy.divide_and_wise_exp(-1 * temperature);
    double partition_function = divided_and_wise_exped_energy.dot(degeneracy);
    double value_numerator =
        value.element_wise_multiplication(degeneracy).dot(divided_and_wise_exped_energy);
    return value_numerator / partition_function;
}
}  // namespace magnetic_susceptibility