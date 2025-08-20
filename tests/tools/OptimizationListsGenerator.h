#ifndef SPINNER_OPTIMIZATIONLISTSGENERATOR_H
#define SPINNER_OPTIMIZATIONLISTSGENERATOR_H

#include <vector>

#include "src/common/physical_optimization/OptimizationList.h"
#include "src/group/Group.h"

std::vector<common::physical_optimization::OptimizationList> generate_all_optimization_lists(const std::vector<group::Group> groups);

#endif  //SPINNER_OPTIMIZATIONLISTSGENERATOR_H