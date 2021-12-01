#include "ChiT.h"

#include <algorithm>

namespace magnetic_susceptibility {

MuSquaredWorker::MuSquaredWorker(DenseVector&& energy, DenseVector&& degeneracy) :
    ensemble_averager_(std::move(energy), std::move(degeneracy)) {}

void MuSquaredWorker::initializeExperimentalValues(
    std::vector<ValueAtTemperature> experimental_values,
    ExperimentalValuesEnum experimental_quantity_type,
    double number_of_centers_ratio) {
    if (experimental_values_worker_.has_value()) {
        throw std::invalid_argument("Experimantal values have been already initialized");
    }

    experimental_values_worker_ = ExperimentalValuesWorker(
        std::move(experimental_values),
        experimental_quantity_type,
        number_of_centers_ratio);
}

}  // namespace magnetic_susceptibility