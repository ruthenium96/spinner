#include "BasicWorker.h"
#include <memory>
#include "src/entities/magnetic_susceptibility/assistant/CommonEnsembleAverager.h"
#include "src/entities/magnetic_susceptibility/assistant/FTLMEnsembleAverager.h"

namespace magnetic_susceptibility::worker {
BasicWorker::BasicWorker(
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra) {
    if (holdsMany(flattenedSpectra->getFlattenSpectrum(common::Energy).value())) {
        ensemble_averager_ = std::make_shared<FTLMEnsembleAverager>(flattenedSpectra);
    } else {
        ensemble_averager_ = std::make_shared<CommonEnsembleAverager>(flattenedSpectra);
    }
}

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