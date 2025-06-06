#include "CommonEnsembleAverager.h"

#include "src/common/OneOrMany.h"
#include "src/common/Quantity.h"

namespace magnetic_susceptibility {

CommonEnsembleAverager::CommonEnsembleAverager(
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra) :
    flattenedSpectra_(flattenedSpectra) {}

double CommonEnsembleAverager::ensemble_average(
    OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>> value,
    double temperature) const {
    const auto& energy_vector = getOneRef(flattenedSpectra_->getFlattenSpectrum(common::Energy).value()).get();
    const auto& degeneracy_vector = flattenedSpectra_->getDegeneracyValues();
    const auto& value_vector = getOneRef(value).get();
    auto divided_and_wise_exped_energy = energy_vector->multiply_by(-1.0 / temperature);
    divided_and_wise_exped_energy->wise_exp();
    double partition_function = divided_and_wise_exped_energy->dot(degeneracy_vector);
    double value_numerator = value_vector->triple_dot(degeneracy_vector, divided_and_wise_exped_energy);
    return value_numerator / partition_function;
}
}  // namespace magnetic_susceptibility