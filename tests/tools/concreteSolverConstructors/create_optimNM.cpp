#include "create_optimNM.h"

template<>
std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>
createConcreteSolver<nonlinear_solver::optimNMAdapter>() {
    return std::make_shared<nonlinear_solver::optimNMAdapter>();
};
