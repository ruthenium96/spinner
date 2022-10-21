#include "stlbfgsAdapter.h"

#include <stlbfgs.h>

namespace nonlinear_solver {
void stlbfgsAdapter::optimize(
    std::function<double(const std::vector<double>&, std::vector<double>&, bool)> oneStepFunction,
    std::vector<double>& changeable_values) {
    auto adaptedSignatureFunction = adaptSignature(oneStepFunction);
    STLBFGS::Optimizer optimizer {adaptedSignatureFunction};
    optimizer.verbose = false;
    // Run calculation from initial guess. STLBFGS updates changeable_values every iteration.
    optimizer.run(changeable_values);
}

std::function<void(const std::vector<double>&, double&, std::vector<double>&)>
stlbfgsAdapter::adaptSignature(
    const std::function<double(const std::vector<double>&, std::vector<double>&, bool)>&
        oneStepFunction) {
    auto adaptedSignatureFunction = [oneStepFunction](
                                        const std::vector<double>& changeable_values,
                                        double& residual_error,
                                        std::vector<double>& gradient) {
        residual_error = oneStepFunction(changeable_values, gradient, true);
    };
    return adaptedSignatureFunction;
}

}  // namespace nonlinear_solver