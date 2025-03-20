#include "EnsembleAverager.h"
#include "src/common/Quantity.h"

namespace magnetic_susceptibility {

EnsembleAverager::EnsembleAverager(
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra) :
    flattenedSpectra_(flattenedSpectra) {}

double EnsembleAverager::ensemble_average(
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& value,
    double temperature) const {
    const auto& energy_vector = flattenedSpectra_->getFlattenSpectrum(common::Energy).value().get();
    const auto& degeneracy_vector = flattenedSpectra_->getDegeneracyValues();
    auto divided_and_wise_exped_energy = energy_vector->multiply_by(-1.0 / temperature);
    divided_and_wise_exped_energy->wise_exp();
    double partition_function = divided_and_wise_exped_energy->dot(degeneracy_vector);
    double value_numerator =
        value->element_wise_multiplication(degeneracy_vector)->dot(divided_and_wise_exped_energy);
    return value_numerator / partition_function;
}
}  // namespace magnetic_susceptibility