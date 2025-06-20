#ifndef SPINNER_ABSTRACTENSEMBLEAVERAGER_H
#define SPINNER_ABSTRACTENSEMBLEAVERAGER_H

#include <memory>

#include "src/common/OneOrMany.h"
#include "src/common/UncertainValue.h"
#include "src/entities/data_structures/AbstractDenseVector.h"

namespace magnetic_susceptibility {
class AbstractEnsembleAverager {
  public:
    virtual common::UncertainValue ensemble_average(
        OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>> value,
        double temperature) const = 0;
};
}

#endif //SPINNER_ABSTRACTENSEMBLEAVERAGER_H