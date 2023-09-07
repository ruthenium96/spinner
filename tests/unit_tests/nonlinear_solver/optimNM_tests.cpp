#include "abstractSolver_tests.h"
#include "tests/tools/concreteSolverConstructors/create_optimNM.h"

INSTANTIATE_TYPED_TEST_SUITE_P(optimNMSolverTests, find_local_minima, optimNM);