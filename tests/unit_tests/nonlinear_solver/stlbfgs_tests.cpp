#include "abstractSolver_tests.h"
#include "tests/tools/concreteSolverConstructors/create_stlbfgs.h"

INSTANTIATE_TYPED_TEST_SUITE_P(stlbfgsSolverTests, find_local_minima, stlbfgs);
