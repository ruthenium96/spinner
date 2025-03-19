#include "EnsembleAverager.h"

namespace magnetic_susceptibility {

EnsembleAverager::EnsembleAverager(
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& energy,
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& degeneracy) :
    energy_(std::move(energy)),
    degeneracy_(std::move(degeneracy)) {}

double EnsembleAverager::ensemble_average(
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& value,
    double temperature) const {
    auto divided_and_wise_exped_energy = energy_->multiply_by(-1.0 / temperature);
    divided_and_wise_exped_energy->wise_exp();
    double partition_function = divided_and_wise_exped_energy->dot(degeneracy_);
    double value_numerator =
        value->element_wise_multiplication(degeneracy_)->dot(divided_and_wise_exped_energy);
    return value_numerator / partition_function;
}
}  // namespace magnetic_susceptibility