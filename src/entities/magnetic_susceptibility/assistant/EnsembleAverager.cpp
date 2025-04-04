#include "EnsembleAverager.h"
#include "src/common/Quantity.h"

namespace magnetic_susceptibility {

EnsembleAverager::EnsembleAverager(
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra) :
    flattenedSpectra_(flattenedSpectra) {}

double EnsembleAverager::ensemble_average(
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& value,
    double temperature) const {
    if (flattenedSpectra_->getFlattenSpectrum(common::squared_back_projection).has_value()) {
        return ensemble_average_ftlm(value, temperature);
    } else {
        return ensemble_average_common(value, temperature);
    }
}

double EnsembleAverager::ensemble_average_ftlm(
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& value,
    double temperature) const {
    const auto& energy_vector = flattenedSpectra_->getFlattenSpectrum(common::Energy).value().get();
    const auto& degeneracy_vector = flattenedSpectra_->getDegeneracyValues();
    const auto& squared_back_projection_vector = flattenedSpectra_->getFlattenSpectrum(common::squared_back_projection).value().get();
    auto divided_and_wise_exped_energy = energy_vector->multiply_by(-1.0 / temperature);
    divided_and_wise_exped_energy->wise_exp();
    double partition_function = divided_and_wise_exped_energy->triple_dot(degeneracy_vector, squared_back_projection_vector);
    double value_numerator = value->triple_dot(degeneracy_vector, divided_and_wise_exped_energy);
    return value_numerator / partition_function;
}

double EnsembleAverager::ensemble_average_common(
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& value,
    double temperature) const {
    const auto& energy_vector = flattenedSpectra_->getFlattenSpectrum(common::Energy).value().get();
    const auto& degeneracy_vector = flattenedSpectra_->getDegeneracyValues();
    auto divided_and_wise_exped_energy = energy_vector->multiply_by(-1.0 / temperature);
    divided_and_wise_exped_energy->wise_exp();
    double partition_function = divided_and_wise_exped_energy->dot(degeneracy_vector);
    double value_numerator = value->triple_dot(degeneracy_vector, divided_and_wise_exped_energy);
    return value_numerator / partition_function;
}
}  // namespace magnetic_susceptibility