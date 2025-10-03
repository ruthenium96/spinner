#ifndef SPINNER_FTLMENSEMBLEAVERAGER_H
#define SPINNER_FTLMENSEMBLEAVERAGER_H

#include <memory>

#include "AbstractEnsembleAverager.h"
#include "src/common/UncertainValue.h"
#include "src/eigendecompositor/FlattenedSpectra.h"
#include "src/entities/data_structures/AbstractDenseVector.h"

// Calculates ensemble-averaged values using Boltzmann distribution.
namespace magnetic_susceptibility {
class FTLMEnsembleAverager : public AbstractEnsembleAverager {
  public:
    FTLMEnsembleAverager(std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra);
    common::UncertainValue ensemble_average(
        OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>> values,
        double temperature) const override;
    std::pair<std::vector<double>, std::vector<double>> ensemble_average_numerator_denominator(
        const OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>& values,
        double temperature) const;  
  private:
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra_;
};
}  // namespace magnetic_susceptibility

#endif  //SPINNER_FTLMENSEMBLEAVERAGER_H
