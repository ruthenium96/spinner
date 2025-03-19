#ifndef SPINNER_ALLSPARSESEMIUNITARYMATRIXFACTORIES_H
#define SPINNER_ALLSPARSESEMIUNITARYMATRIXFACTORIES_H

#include <memory>
#include <vector>

#include "src/entities/data_structures/AbstractFactories.h"

std::vector<std::shared_ptr<quantum::linear_algebra::AbstractSparseTransformFactory>>
constructAllSparseSemiunitaryMatrixFactories();

#endif  //SPINNER_ALLSPARSESEMIUNITARYMATRIXFACTORIES_H
