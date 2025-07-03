#include "stlbfgsAdapter.h"

#include <stlbfgs.h>
#include "src/common/Logger.h"

namespace {
std::function<void(const std::vector<double>&, double&, std::vector<double>&)>
adaptSignature(const std::function<double(const std::vector<double>&, std::vector<double>&, bool)>&
                   oneStepFunction) {
    auto adaptedSignatureFunction = [oneStepFunction](
                                        const std::vector<double>& changeable_values,
                                        double& residual_error,
                                        std::vector<double>& gradient) {
        // the last boolean is doesGradientsRequired
        residual_error = oneStepFunction(changeable_values, gradient, true);
    };
    return adaptedSignatureFunction;
}

std::vector<double> constructMainDiagonalOfInverseHessian(
    const STLBFGS::Optimizer& optimizer, size_t size) {
    auto mainDiagonalOfInverseHessian = std::vector<double>(size);
    std::vector<double> row_of_Hessian(size, 0.0);

    for (size_t i = 0; i < size; ++i) {
        std::vector<double> orth(size, 0.0);
        orth[i] = 1.0;
    
        optimizer.invH.mult(orth, row_of_Hessian);
        mainDiagonalOfInverseHessian[i] = row_of_Hessian[i];
    }
    return std::move(mainDiagonalOfInverseHessian);
}

}  // namespace

namespace nonlinear_solver {
void stlbfgsAdapter::optimize(
    std::function<double(const std::vector<double>&, std::vector<double>&, bool)> oneStepFunction,
    std::vector<double>& changeable_values) {
    auto adaptedSignatureFunction = adaptSignature(oneStepFunction);
    STLBFGS::Optimizer optimizer {adaptedSignatureFunction};
    if (common::Logger::get_level() <= common::PrintLevel::debug) {
        optimizer.verbose = true;
    } else {
        optimizer.verbose = false;
    }
    optimizer.ftol = 1e-16;
    optimizer.gtol = 1e-18;

    auto spdlog_buffer = common::SpdlogStreamMsg(common::debug);
    auto coutbuf = std::cout.rdbuf(&spdlog_buffer);
    // Run calculation from initial guess. STLBFGS updates changeable_values every iteration.
    optimizer.run(changeable_values);
    std::cout.rdbuf(coutbuf);

    mainDiagonalOfInverseHessian_ = 
        constructMainDiagonalOfInverseHessian(optimizer, changeable_values.size());
}

std::optional<std::vector<double>> stlbfgsAdapter::getMainDiagonalOfInverseHessian() const {
    return mainDiagonalOfInverseHessian_;
}

}  // namespace nonlinear_solver