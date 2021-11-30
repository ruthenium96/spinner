#ifndef JULY_ENSEMBLEAVERAGER_H
#define JULY_ENSEMBLEAVERAGER_H

#include <cmath>

#include "entities/data_structures/DenseMatrix.h"
namespace magnetic_susceptibility {
class EnsembleAverager {
  public:
    EnsembleAverager(DenseVector&& energy, DenseVector&& degeneracy);
    double ensemble_average(const DenseVector& value, double temperature);

  private:
    const DenseVector energy_;
    const DenseVector degeneracy_;
    DenseVector divided_and_wise_exped_energy;
    double partition_function = NAN;
    double last_temperature = NAN;
};
}  // namespace magnetic_susceptibility

#endif  //JULY_ENSEMBLEAVERAGER_H
