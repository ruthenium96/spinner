#include "AbstractEnsembleAverager.h"

namespace magnetic_susceptibility {

std::pair<double, double> AbstractEnsembleAverager::calculate_averaged_value_and_partition_function(
    double temperature,
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& energy_vector, 
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& weights_vector, 
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& value_vector) const {
    auto divided_and_wise_exped_energy = energy_vector->multiply_by(-1.0 / temperature);
    divided_and_wise_exped_energy->wise_exp();
    double partition_function = divided_and_wise_exped_energy->dot(weights_vector);
    double value_numerator = value_vector.get()->triple_dot(weights_vector, divided_and_wise_exped_energy);
    return std::pair{value_numerator, partition_function};   
}
}