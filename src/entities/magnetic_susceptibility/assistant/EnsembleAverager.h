#ifndef SPINNER_ENSEMBLEAVERAGER_H
#define SPINNER_ENSEMBLEAVERAGER_H

#include <cmath>
#include <memory>

#include "src/eigendecompositor/FlattenedSpectra.h"
#include "src/entities/data_structures/AbstractDenseVector.h"

// Calculates ensemble-averaged values using Boltzmann distribution.
namespace magnetic_susceptibility {
class EnsembleAverager {
  public:
    EnsembleAverager(std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra);
    double ensemble_average(
        const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& value,
        double temperature) const;

  private:
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra_;
};
}  // namespace magnetic_susceptibility

#endif  //SPINNER_ENSEMBLEAVERAGER_H
