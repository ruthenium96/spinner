#include "BasicWorker.h"

namespace magnetic_susceptibility::worker {
BasicWorker::BasicWorker(
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra) :
    ensemble_averager_(flattenedSpectra) {}

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