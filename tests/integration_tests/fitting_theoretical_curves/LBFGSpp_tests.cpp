#include "abstractSolver_tests.h"
#include "tests/tools/concreteSolverConstructors/create_LBFGSpp.h"

INSTANTIATE_TYPED_TEST_SUITE_P(LBFGSppSolverTests, fitting_magnetic_susceptibility_simple, LBFGSpp);
