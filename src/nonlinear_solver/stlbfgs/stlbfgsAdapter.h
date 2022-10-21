#ifndef SPINNER_STLBFGSADAPTER_H
#define SPINNER_STLBFGSADAPTER_H

#include "src/nonlinear_solver/AbstractNonlinearSolver.h"

namespace nonlinear_solver {

class stlbfgsAdapter: public AbstractNonlinearSolver {
  public:
    void optimize(
        std::function<double(const std::vector<double>&, std::vector<double>&, bool)>
            oneStepFunction,
        std::vector<double>& changeable_values) override;

  private:
    std::function<void(const std::vector<double>&, double&, std::vector<double>&)> adaptSignature(
        const std::function<double(const std::vector<double>&, std::vector<double>&, bool)>&
            oneStepFunction);
};

}  // namespace nonlinear_solver

#endif  //SPINNER_STLBFGSADAPTER_H
