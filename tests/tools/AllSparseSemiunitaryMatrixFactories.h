#ifndef SPINNER_ALLSPARSESEMIUNITARYMATRIXFACTORIES_H
#define SPINNER_ALLSPARSESEMIUNITARYMATRIXFACTORIES_H

#include <memory>
#include <vector>

#include "src/linalg_structures/AbstractFactories.h"

std::vector<std::shared_ptr<spinner::linalg_structures::AbstractSparseTransformFactory>>
constructAllSparseSemiunitaryMatrixFactories();

#endif  //SPINNER_ALLSPARSESEMIUNITARYMATRIXFACTORIES_H
