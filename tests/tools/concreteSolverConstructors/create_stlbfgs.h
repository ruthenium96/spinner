#ifndef SPINNER_CREATE_STLBFGS_H
#define SPINNER_CREATE_STLBFGS_H

#include "src/nonlinear_solver/stlbfgs/stlbfgsAdapter.h"
#include "tests/integration_tests/fitting_theoretical_curves/abstractSolver_tests.h"
#include "tests/unit_tests/nonlinear_solver/abstractSolver_tests.h"

template<>
std::shared_ptr<spinner::nonlinear_solver::AbstractNonlinearSolver>
createConcreteSolver<spinner::nonlinear_solver::stlbfgsAdapter>();

typedef testing::Types<spinner::nonlinear_solver::stlbfgsAdapter> stlbfgs;

#endif  //SPINNER_CREATE_STLBFGS_H