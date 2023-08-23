#include "create_stlbfgs.h"

template<>
std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>
createConcreteSolver<nonlinear_solver::stlbfgsAdapter>() {
    return std::make_shared<nonlinear_solver::stlbfgsAdapter>();
};
