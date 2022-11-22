#ifndef SPINNER_BASICWORKER_H
#define SPINNER_BASICWORKER_H

#include <optional>

#include "src/entities/magnetic_susceptibility/assistant/EnsembleAverager.h"
#include "src/entities/magnetic_susceptibility/assistant/ExperimentalValuesWorker.h"
#include "src/entities/magnetic_susceptibility/worker/AbstractWorker.h"

namespace magnetic_susceptibility::worker {
class BasicWorker: public AbstractWorker {
  public:
    BasicWorker(
        std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& energy,
        std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& degeneracy);
    std::shared_ptr<ExperimentalValuesWorker> getExperimentalValuesWorker() override;
    std::shared_ptr<const ExperimentalValuesWorker> getExperimentalValuesWorker() const override;
    void setExperimentalValuesWorker(
        const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker) override;

  protected:
    const EnsembleAverager ensemble_averager_;
    std::optional<std::shared_ptr<ExperimentalValuesWorker>> experimental_values_worker_;
};
}  // namespace magnetic_susceptibility
#endif  //SPINNER_BASICWORKER_H
