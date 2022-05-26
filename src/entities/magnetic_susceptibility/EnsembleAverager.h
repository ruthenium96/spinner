#ifndef SPINNER_ENSEMBLEAVERAGER_H
#define SPINNER_ENSEMBLEAVERAGER_H

#include <cmath>

#include "src/entities/data_structures/DenseMatrix.h"
namespace magnetic_susceptibility {
class EnsembleAverager {
  public:
    EnsembleAverager(DenseVector&& energy, DenseVector&& degeneracy);
    double ensemble_average(const DenseVector& value, double temperature) const;

  private:
    const DenseVector energy_;
    const DenseVector degeneracy_;
    mutable DenseVector divided_and_wise_exped_energy;
    mutable double partition_function = NAN;
    mutable double last_temperature = NAN;
};
}  // namespace magnetic_susceptibility

#endif  //SPINNER_ENSEMBLEAVERAGER_H
