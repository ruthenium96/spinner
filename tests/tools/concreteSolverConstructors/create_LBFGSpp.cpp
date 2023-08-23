#include "create_LBFGSpp.h"

template<>
std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>
createConcreteSolver<nonlinear_solver::LBFGSppAdapter>() {
    return std::make_shared<nonlinear_solver::LBFGSppAdapter>();
};