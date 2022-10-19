#include "abstractSolver_tests.h"
#include "src/nonlinear_solver/stlbfgs/stlbfgsAdapter.h"

template<>
std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>
createConcreteSolver<nonlinear_solver::stlbfgsAdapter>() {
    return std::make_shared<nonlinear_solver::stlbfgsAdapter>();
};

typedef testing::Types<nonlinear_solver::stlbfgsAdapter> stlbfgs;
INSTANTIATE_TYPED_TEST_SUITE_P(stlbfgsSolverTests, fitting_magnetic_susceptibility, stlbfgs);
