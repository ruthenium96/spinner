#ifndef SPINNER_ABSTRACTENSEMBLEAVERAGER_H
#define SPINNER_ABSTRACTENSEMBLEAVERAGER_H

#include <memory>

#include "src/common/OneOrMany.h"
#include "src/common/UncertainValue.h"
#include "src/linalg_structures/AbstractDenseVector.h"

namespace spinner::magnetic_susceptibility {
class AbstractEnsembleAverager {
  public:
    virtual common::UncertainValue ensemble_average(
        OneOrMany<std::reference_wrapper<const std::unique_ptr<linalg_structures::AbstractDenseVector>>> value,
        double temperature) const = 0;
  protected:
      std::pair<double, double> calculate_averaged_value_and_partition_function(
        double temperature,
        const std::unique_ptr<linalg_structures::AbstractDenseVector>& energy_vector, 
        const std::unique_ptr<linalg_structures::AbstractDenseVector>& weights_vector, 
        const std::unique_ptr<linalg_structures::AbstractDenseVector>& value_vector) const;
};
} // namespace spinner::magnetic_susceptibility

#endif //SPINNER_ABSTRACTENSEMBLEAVERAGER_H