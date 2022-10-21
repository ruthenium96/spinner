#ifndef SPINNER_ABSTRACTNONLINEARSOLVER_H
#define SPINNER_ABSTRACTNONLINEARSOLVER_H

#include <functional>
#include <vector>

namespace nonlinear_solver {

class AbstractNonlinearSolver {
  public:
    virtual void optimize(
        std::function<double(const std::vector<double>&, std::vector<double>&, bool)>
            oneStepFunction,
        std::vector<double>& changeable_values) = 0;
};

}  // namespace nonlinear_solver

#endif  //SPINNER_ABSTRACTNONLINEARSOLVER_H
