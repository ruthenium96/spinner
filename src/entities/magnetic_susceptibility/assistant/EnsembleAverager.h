#ifndef SPINNER_ENSEMBLEAVERAGER_H
#define SPINNER_ENSEMBLEAVERAGER_H

#include <cmath>

#include "src/entities/data_structures/AbstractVector.h"

// Calculates ensemble-averaged values using Boltzmann distribution.
namespace magnetic_susceptibility {
class EnsembleAverager {
  public:
    EnsembleAverager(
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& energy,
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& degeneracy);
    double ensemble_average(
        const std::unique_ptr<quantum::linear_algebra::AbstractVector>& value,
        double temperature) const;

  private:
    std::unique_ptr<quantum::linear_algebra::AbstractVector> energy_;
    const std::unique_ptr<quantum::linear_algebra::AbstractVector> degeneracy_;
    mutable std::unique_ptr<quantum::linear_algebra::AbstractVector> divided_and_wise_exped_energy;
    mutable double partition_function = NAN;
    mutable double last_temperature = NAN;
};
}  // namespace magnetic_susceptibility

#endif  //SPINNER_ENSEMBLEAVERAGER_H
