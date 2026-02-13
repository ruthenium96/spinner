#ifndef SPINNER_ALLSYMMETRICMATRIXFACTORIES_H
#define SPINNER_ALLSYMMETRICMATRIXFACTORIES_H

#include <vector>

#include "src/linalg_structures/AbstractFactories.h"

std::vector<std::shared_ptr<spinner::linalg_structures::AbstractDenseTransformAndDiagonalizeFactory>>
constructAllDenseTransformAndDiagonalizeFactories();

#endif  //SPINNER_ALLSYMMETRICMATRIXFACTORIES_H
