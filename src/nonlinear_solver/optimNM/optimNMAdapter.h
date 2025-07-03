#ifndef SPINNER_OPTIMNMADAPTER_H
#define SPINNER_OPTIMNMADAPTER_H

#include "src/nonlinear_solver/AbstractNonlinearSolver.h"

namespace nonlinear_solver {

class optimNMAdapter: public AbstractNonlinearSolver {
  public:
    void optimize(
        std::function<double(const std::vector<double>&, std::vector<double>&, bool)>
            oneStepFunction,
        std::vector<double>& changeable_values) override;
    bool doesGradientsRequired() const override {
        return false;
    };
    // we do not have approximation of Hessian in the case of Nelder-Mead method:
    std::optional<std::vector<double>> getMainDiagonalOfInverseHessian() const override {
        return std::nullopt;
    }
};

}  // namespace nonlinear_solver

#endif  //SPINNER_OPTIMNMADAPTER_H
