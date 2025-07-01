#include "FTLMEnsembleAverager.h"

#include <cmath>
#include <functional>
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
    auto fraction = ensemble_average_numerator_denominator(values, temperature);
    size_t number_of_seeds = fraction.second.size();
    double averaged_partition_function = average_over_seeds(fraction.second);
    double value = average_over_seeds(fraction.first) / averaged_partition_function;
    double uncertainty = std::abs(value) / std::sqrt(averaged_partition_function * number_of_seeds); 
    return common::UncertainValue(value, uncertainty, common::UncertaintySources::FTLM);
}

std::pair<std::vector<double>, std::vector<double>> FTLMEnsembleAverager::ensemble_average_numerator_denominator(
    const OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>& values,
    double temperature) const {
    auto energy_vectors = flattenedSpectra_->getFlattenSpectrum(common::Energy).value();
    auto weights_vectors = flattenedSpectra_->getWeights();
    const auto& value_vectors = values;

    OneOrMany<std::pair<double, double>> results = transform_one_or_many(
        std::function([temperature](
            std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>> energy_vector, 
            std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>> weights_vector, 
            std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>> value_vector){
            auto divided_and_wise_exped_energy = energy_vector.get()->multiply_by(-1.0 / temperature);
            divided_and_wise_exped_energy->wise_exp();
            double partition_function = divided_and_wise_exped_energy->dot(weights_vector.get());
            double value_numerator = value_vector.get()->triple_dot(weights_vector.get(), divided_and_wise_exped_energy);
            return std::pair{partition_function, value_numerator};   
        }),
        energy_vectors,
        weights_vectors,
        value_vectors
    );

   std::vector<double> partition_functions;
   std::vector<double> value_numerators;

    for (const auto& pair : getManyRef(results)) {
        partition_functions.push_back(pair.first);
        value_numerators.push_back(pair.second);
    }

    return {std::move(value_numerators), std::move(partition_functions)};
}

}  // namespace magnetic_susceptibility