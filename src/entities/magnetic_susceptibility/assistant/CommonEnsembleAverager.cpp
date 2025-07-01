#include "CommonEnsembleAverager.h"

#include "src/common/OneOrMany.h"
#include "src/common/Quantity.h"
#include "src/common/UncertainValue.h"

namespace magnetic_susceptibility {

CommonEnsembleAverager::CommonEnsembleAverager(
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra) :
    flattenedSpectra_(flattenedSpectra) {}

common::UncertainValue CommonEnsembleAverager::ensemble_average(
    OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>> value,
    double temperature) const {
    auto [value_numerator, partition_function] = 
        ensemble_average_numerator_denominator(getOneRef(value).get(), temperature);
    return common::UncertainValue(value_numerator / partition_function);
}

std::pair<double, double> CommonEnsembleAverager::ensemble_average_numerator_denominator(
    std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>> value,
    double temperature) const {
    const auto& energy_vector = getOneRef(flattenedSpectra_->getFlattenSpectrum(common::Energy).value()).get();
    const auto& weights_vector = getOneRef(flattenedSpectra_->getWeights());
    const auto& value_vector = value.get();
    return calculate_averaged_value_and_partition_function(
        temperature,
        energy_vector,
        weights_vector,
        value_vector);
}

}  // namespace magnetic_susceptibility