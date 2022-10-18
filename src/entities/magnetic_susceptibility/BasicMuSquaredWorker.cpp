#include "BasicMuSquaredWorker.h"

namespace magnetic_susceptibility {
BasicMuSquaredWorker::BasicMuSquaredWorker(
    std::unique_ptr<quantum::linear_algebra::AbstractVector>&& energy,
    std::unique_ptr<quantum::linear_algebra::AbstractVector>&& degeneracy) :
    ensemble_averager_(std::move(energy), std::move(degeneracy)) {}

std::shared_ptr<ExperimentalValuesWorker> BasicMuSquaredWorker::getExperimentalValuesWorker() {
    return experimental_values_worker_.value();
}

std::shared_ptr<const ExperimentalValuesWorker>
BasicMuSquaredWorker::getExperimentalValuesWorker() const {
    return experimental_values_worker_.value();
}

void BasicMuSquaredWorker::setExperimentalValuesWorker(
    const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker) {
    experimental_values_worker_ = experimental_values_worker;
}

}  // namespace magnetic_susceptibility