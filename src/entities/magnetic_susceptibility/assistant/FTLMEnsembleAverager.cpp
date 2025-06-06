#include "FTLMEnsembleAverager.h"

#include <cmath>
#include <numeric>

#include "src/common/OneOrMany.h"
#include "src/common/Quantity.h"

namespace {

double average_over_seeds(const std::vector<double>& values) {
    return std::accumulate(values.begin(), values.end(), 0.0) / (double)values.size();
}

} // namespace

namespace magnetic_susceptibility {

FTLMEnsembleAverager::FTLMEnsembleAverager(
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra) :
    flattenedSpectra_(flattenedSpectra) {}

double FTLMEnsembleAverager::ensemble_average(
    OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>> values,
    double temperature) const {
        auto energy_vectors = getManyRef(flattenedSpectra_->getFlattenSpectrum(common::Energy).value());
    const auto& degeneracy_vector = flattenedSpectra_->getDegeneracyValues();
    auto squared_back_projection_vectors = getManyRef(flattenedSpectra_->getFlattenSpectrum(common::squared_back_projection).value());
    auto value_vectors = getManyRef(values);

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

    double averaged_partition_function = average_over_seeds(partition_functions);
    double averaged_value_numenator = average_over_seeds(value_numerators);

    return averaged_value_numenator / averaged_partition_function;
}
}  // namespace magnetic_susceptibility