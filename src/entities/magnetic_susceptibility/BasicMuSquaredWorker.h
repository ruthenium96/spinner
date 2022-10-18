#ifndef SPINNER_BASICMUSQUAREDWORKER_H
#define SPINNER_BASICMUSQUAREDWORKER_H

#include <optional>

#include "AbstractMuSquaredWorker.h"
#include "src/entities/magnetic_susceptibility/assistant/EnsembleAverager.h"
#include "src/entities/magnetic_susceptibility/assistant/ExperimentalValuesWorker.h"

namespace magnetic_susceptibility {
class BasicMuSquaredWorker: public AbstractMuSquaredWorker {
  public:
    BasicMuSquaredWorker(
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& energy,
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& degeneracy);
    std::shared_ptr<ExperimentalValuesWorker> getExperimentalValuesWorker() override;
    std::shared_ptr<const ExperimentalValuesWorker> getExperimentalValuesWorker() const override;
    void setExperimentalValuesWorker(
        const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker) override;

  protected:
    const EnsembleAverager ensemble_averager_;
    std::optional<std::shared_ptr<ExperimentalValuesWorker>> experimental_values_worker_;
};
}  // namespace magnetic_susceptibility
#endif  //SPINNER_BASICMUSQUAREDWORKER_H
