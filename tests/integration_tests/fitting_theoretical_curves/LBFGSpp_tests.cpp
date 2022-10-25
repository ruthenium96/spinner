#include "abstractSolver_tests.h"
#include "src/nonlinear_solver/LBFGSpp/LBFGSppAdapter.h"

template<>
std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>
createConcreteSolver<nonlinear_solver::LBFGSppAdapter>() {
    return std::make_shared<nonlinear_solver::LBFGSppAdapter>();
};

typedef testing::Types<nonlinear_solver::LBFGSppAdapter> LBFGSpp;
INSTANTIATE_TYPED_TEST_SUITE_P(LBFGSppSolverTests, fitting_magnetic_susceptibility_simple, LBFGSpp);
