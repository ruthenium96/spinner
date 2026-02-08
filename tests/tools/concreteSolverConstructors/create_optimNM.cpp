#include "create_optimNM.h"

template<>
std::shared_ptr<spinner::nonlinear_solver::AbstractNonlinearSolver>
createConcreteSolver<spinner::nonlinear_solver::optimNMAdapter>() {
    return std::make_shared<spinner::nonlinear_solver::optimNMAdapter>();
};
