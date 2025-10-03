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
    bool doesGradientsRequired() const override {
        return true;
    };
    std::optional<std::vector<double>> getMainDiagonalOfInverseHessian() const override;
  private:
    std::vector<double> mainDiagonalOfInverseHessian_;
};

}  // namespace nonlinear_solver

#endif  //SPINNER_STLBFGSADAPTER_H
