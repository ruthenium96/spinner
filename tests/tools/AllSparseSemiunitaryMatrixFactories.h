#ifndef SPINNER_ALLSPARSESEMIUNITARYMATRIXFACTORIES_H
#define SPINNER_ALLSPARSESEMIUNITARYMATRIXFACTORIES_H

#include "src/entities/data_structures/arma/ArmaFactories.h"
#include "src/entities/data_structures/eigen/EigenFactories.h"
#include "src/entities/data_structures/hashmap/StdFactories.h"

std::vector<std::shared_ptr<quantum::linear_algebra::AbstractSparseTransformFactory>>
constructAllSparseSemiunitaryMatrixFactories();

#endif  //SPINNER_ALLSPARSESEMIUNITARYMATRIXFACTORIES_H
