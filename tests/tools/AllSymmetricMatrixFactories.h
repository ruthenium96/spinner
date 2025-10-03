#ifndef SPINNER_ALLSYMMETRICMATRIXFACTORIES_H
#define SPINNER_ALLSYMMETRICMATRIXFACTORIES_H

#include <vector>

#include "src/entities/data_structures/AbstractFactories.h"

std::vector<std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>>
constructAllDenseTransformAndDiagonalizeFactories();

#endif  //SPINNER_ALLSYMMETRICMATRIXFACTORIES_H
