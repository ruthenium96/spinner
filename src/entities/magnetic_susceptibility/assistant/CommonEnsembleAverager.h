#ifndef SPINNER_COMMONENSEMBLEAVERAGER_H
#define SPINNER_COMMONENSEMBLEAVERAGER_H

#include <memory>

#include "AbstractEnsembleAverager.h"
#include "src/eigendecompositor/FlattenedSpectra.h"
#include "src/entities/data_structures/AbstractDenseVector.h"

// Calculates ensemble-averaged values using Boltzmann distribution.
namespace magnetic_susceptibility {
class CommonEnsembleAverager : public AbstractEnsembleAverager {
  public:
    CommonEnsembleAverager(std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra);
    double ensemble_average(
        OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>> value,
        double temperature) const override;
  private:
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra_;
  
};
}  // namespace magnetic_susceptibility

#endif  //SPINNER_COMMONENSEMBLEAVERAGER_H
