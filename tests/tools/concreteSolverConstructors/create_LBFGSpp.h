#ifndef SPINNER_CREATE_LBFGSPP_H
#define SPINNER_CREATE_LBFGSPP_H

#include "src/nonlinear_solver/LBFGSpp/LBFGSppAdapter.h"
#include "tests/integration_tests/fitting_theoretical_curves/abstractSolver_tests.h"
#include "tests/unit_tests/nonlinear_solver/abstractSolver_tests.h"

template<>
std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>
createConcreteSolver<nonlinear_solver::LBFGSppAdapter>();

typedef testing::Types<nonlinear_solver::LBFGSppAdapter> LBFGSpp;

#endif  //SPINNER_CREATE_LBFGSPP_H