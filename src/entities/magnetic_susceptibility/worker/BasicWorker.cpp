#include "BasicWorker.h"

namespace magnetic_susceptibility::worker {
BasicWorker::BasicWorker(
    std::unique_ptr<quantum::linear_algebra::AbstractVector>&& energy,
    std::unique_ptr<quantum::linear_algebra::AbstractVector>&& degeneracy) :
    ensemble_averager_(std::move(energy), std::move(degeneracy)) {}

std::shared_ptr<ExperimentalValuesWorker> BasicWorker::getExperimentalValuesWorker() {
    return experimental_values_worker_.value();
}

std::shared_ptr<const ExperimentalValuesWorker> BasicWorker::getExperimentalValuesWorker() const {
    return experimental_values_worker_.value();
}

void BasicWorker::setExperimentalValuesWorker(
    const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker) {
    experimental_values_worker_ = experimental_values_worker;
}

}  // namespace magnetic_susceptibility