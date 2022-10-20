#include "abstractSolver_tests.h"
#include "src/nonlinear_solver/optimNM/optimNMAdapter.h"

template<>
std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>
createConcreteSolver<nonlinear_solver::optimNMAdapter>() {
    return std::make_shared<nonlinear_solver::optimNMAdapter>();
};

typedef testing::Types<nonlinear_solver::optimNMAdapter> optimNM;
INSTANTIATE_TYPED_TEST_SUITE_P(optimNMSolverTests, fitting_magnetic_susceptibility, optimNM);
