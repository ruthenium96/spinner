#ifndef SPINNER_ENSEMBLEAVERAGER_H
#define SPINNER_ENSEMBLEAVERAGER_H

#include <cmath>

#include "src/entities/data_structures/AbstractDenseVector.h"

// Calculates ensemble-averaged values using Boltzmann distribution.
namespace magnetic_susceptibility {
class EnsembleAverager {
  public:
    EnsembleAverager(
        std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& energy,
        std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& degeneracy);
    double ensemble_average(
        const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& value,
        double temperature) const;

  private:
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector> energy_;
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector> degeneracy_;
};
}  // namespace magnetic_susceptibility

#endif  //SPINNER_ENSEMBLEAVERAGER_H
