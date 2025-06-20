#include "FTLMEnsembleAverager.h"

#include <cmath>
#include <numeric>

#include "src/common/OneOrMany.h"
#include "src/common/Quantity.h"
#include "src/common/UncertainValue.h"

namespace {

double average_over_seeds(const std::vector<double>& values) {
    return std::accumulate(values.begin(), values.end(), 0.0) / (double)values.size();
}

} // namespace

namespace magnetic_susceptibility {

FTLMEnsembleAverager::FTLMEnsembleAverager(
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra) :
    flattenedSpectra_(flattenedSpectra) {}

common::UncertainValue FTLMEnsembleAverager::ensemble_average(
    OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>> values,
    double temperature) const {
    auto fraction = ensemble_average_numerator_denominator(getManyRef(values), temperature);
    size_t number_of_seeds = fraction.second.size();
    double averaged_partition_function = average_over_seeds(fraction.second);
    double value = average_over_seeds(fraction.first) / averaged_partition_function;
    double uncertainty = value / std::sqrt(averaged_partition_function * number_of_seeds);
    return common::UncertainValue(value, uncertainty, common::UncertaintySources::FTLM);
}

std::pair<std::vector<double>, std::vector<double>> FTLMEnsembleAverager::ensemble_average_numerator_denominator(
    const std::vector<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>& values,
    double temperature) const {
    auto energy_vectors = getManyRef(flattenedSpectra_->getFlattenSpectrum(common::Energy).value());
    const auto& degeneracy_vector = flattenedSpectra_->getDegeneracyValues();
    auto squared_back_projection_vectors = getManyRef(flattenedSpectra_->getFlattenSpectrum(common::squared_back_projection).value());
    const auto& value_vectors = values;

    std::vector<double> partition_functions(energy_vectors.size());
    std::vector<double> value_numerators(energy_vectors.size());

    for (int i = 0; i < energy_vectors.size(); ++i) {
        const auto& energy_vector = energy_vectors[i].get();
        const auto& squared_back_projection_vector = squared_back_projection_vectors[i].get();
        const auto& value_vector = value_vectors[i].get();
        auto divided_and_wise_exped_energy = energy_vector->multiply_by(-1.0 / temperature);
        divided_and_wise_exped_energy->wise_exp();
        double partition_function = divided_and_wise_exped_energy->triple_dot(degeneracy_vector, squared_back_projection_vector);
        double value_numerator = value_vector->triple_dot(degeneracy_vector, divided_and_wise_exped_energy);
        partition_functions.push_back(partition_function);
        value_numerators.push_back(value_numerator);
    }

    return {std::move(value_numerators), std::move(partition_functions)};
}

}  // namespace magnetic_susceptibility