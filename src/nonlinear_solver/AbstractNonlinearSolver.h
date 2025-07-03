#ifndef SPINNER_ABSTRACTNONLINEARSOLVER_H
#define SPINNER_ABSTRACTNONLINEARSOLVER_H

#include <functional>
#include <optional>
#include <vector>

namespace nonlinear_solver {

class AbstractNonlinearSolver {
  public:
    virtual void optimize(
        std::function<double(const std::vector<double>&, std::vector<double>&, bool)>
            oneStepFunction,
        std::vector<double>& changeable_values) = 0;

    virtual bool doesGradientsRequired() const = 0;

    virtual std::optional<std::vector<double>> getMainDiagonalOfInverseHessian() const = 0;

    virtual ~AbstractNonlinearSolver() = default;
};

}  // namespace nonlinear_solver

#endif  //SPINNER_ABSTRACTNONLINEARSOLVER_H
