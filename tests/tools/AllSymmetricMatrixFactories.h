#ifndef SPINNER_ALLSYMMETRICMATRIXFACTORIES_H
#define SPINNER_ALLSYMMETRICMATRIXFACTORIES_H

#include "src/entities/data_structures/arma/ArmaFactories.h"
#include "src/entities/data_structures/eigen/EigenFactories.h"

std::vector<std::shared_ptr<quantum::linear_algebra::AbstractSymmetricMatrixFactory>>
constructAllSymmetricMatrixFactories();

#endif  //SPINNER_ALLSYMMETRICMATRIXFACTORIES_H
