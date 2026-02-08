#ifndef SPINNER_CREATECONCRETESOLVER_H
#define SPINNER_CREATECONCRETESOLVER_H

#include "src/nonlinear_solver/AbstractNonlinearSolver.h"

template<class T>
std::shared_ptr<spinner::nonlinear_solver::AbstractNonlinearSolver> createConcreteSolver();

#endif  //SPINNER_CREATECONCRETESOLVER_H
