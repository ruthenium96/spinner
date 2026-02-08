#include "create_stlbfgs.h"

template<>
std::shared_ptr<spinner::nonlinear_solver::AbstractNonlinearSolver>
createConcreteSolver<spinner::nonlinear_solver::stlbfgsAdapter>() {
    return std::make_shared<spinner::nonlinear_solver::stlbfgsAdapter>();
};
