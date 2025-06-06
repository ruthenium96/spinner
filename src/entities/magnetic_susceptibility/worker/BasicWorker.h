#ifndef SPINNER_BASICWORKER_H
#define SPINNER_BASICWORKER_H

#include <memory>
#include <optional>

#include "src/eigendecompositor/FlattenedSpectra.h"
#include "src/entities/magnetic_susceptibility/assistant/AbstractEnsembleAverager.h"
#include "src/entities/magnetic_susceptibility/assistant/ExperimentalValuesWorker.h"
#include "src/entities/magnetic_susceptibility/worker/AbstractWorker.h"

namespace magnetic_susceptibility::worker {
class BasicWorker: public AbstractWorker {
  public:
    BasicWorker(std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra);
    std::shared_ptr<ExperimentalValuesWorker> getExperimentalValuesWorker() override;
    std::shared_ptr<const ExperimentalValuesWorker> getExperimentalValuesWorker() const override;
    void setExperimentalValuesWorker(
        const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker) override;

  protected:
    std::shared_ptr<AbstractEnsembleAverager> ensemble_averager_;
    std::optional<std::shared_ptr<ExperimentalValuesWorker>> experimental_values_worker_;
};
}  // namespace magnetic_susceptibility
#endif  //SPINNER_BASICWORKER_H
