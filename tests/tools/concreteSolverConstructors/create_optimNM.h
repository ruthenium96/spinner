#ifndef SPINNER_CREATE_OPTIMNM_H
#define SPINNER_CREATE_OPTIMNM_H

#include "src/nonlinear_solver/optimNM/optimNMAdapter.h"
#include "tests/integration_tests/fitting_theoretical_curves/abstractSolver_tests.h"
#include "tests/unit_tests/nonlinear_solver/abstractSolver_tests.h"

template<>
std::shared_ptr<spinner::nonlinear_solver::AbstractNonlinearSolver>
createConcreteSolver<spinner::nonlinear_solver::optimNMAdapter>();

typedef testing::Types<spinner::nonlinear_solver::optimNMAdapter> optimNM;

#endif  //SPINNER_CREATE_OPTIMNM_H