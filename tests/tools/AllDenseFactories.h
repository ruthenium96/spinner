#ifndef SPINNER_ALLDENSEFACTORIES_H
#define SPINNER_ALLDENSEFACTORIES_H

#include "src/entities/data_structures/arma/ArmaDenseFactory.h"
#include "src/entities/data_structures/eigen/EigenDenseFactory.h"

std::vector<std::shared_ptr<quantum::linear_algebra::AbstractDenseFactory>>
constructAllDenseFactories();

#endif  //SPINNER_ALLDENSEFACTORIES_H
