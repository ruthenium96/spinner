#include "optimNMAdapter.h"

// Do not use OpenMP in solver: it leads to race condition in Runner.
#define OPTIM_DONT_USE_OPENMP
#define OPTIM_ENABLE_ARMA_WRAPPERS
#include <optim.hpp>

namespace nonlinear_solver {

void optimNMAdapter::optimize(
    std::function<double(const std::vector<double>&, std::vector<double>&, bool)> oneStepFunction,
    std::vector<double>& changeable_values) {
    auto adaptedSignatureFunction = adaptSignature(oneStepFunction);
    auto changeable_values_arma = convertFromSTLToArma(changeable_values);

    optim::nm(changeable_values_arma, adaptedSignatureFunction, nullptr);

    changeable_values = convertFromArmaToSTL(changeable_values_arma);
}

std::function<double(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)>
optimNMAdapter::adaptSignature(
    const std::function<double(const std::vector<double>&, std::vector<double>&, bool)>&
        oneStepFunction) {
    // do nothing with opt_data
    auto adaptedSignatureFunction = [oneStepFunction](
                                        const arma::vec& changeable_values_arma,
                                        arma::vec* gradient_arma,
                                        void* opt_data) {
        std::vector<double> changeable_values_stl = convertFromArmaToSTL(changeable_values_arma);
        // we do not need gradients on NM method:
        std::vector<double> gradient_stl(changeable_values_stl.size(), 0);
        // the last boolean is doesGradientsRequired
        double residual_error = oneStepFunction(changeable_values_stl, gradient_stl, false);
        // do nothing with gradient_stl
        return residual_error;
    };
    return adaptedSignatureFunction;
}

std::vector<double> optimNMAdapter::convertFromArmaToSTL(const arma::vec& arma_vector) {
    return arma::conv_to<std::vector<double>>::from(arma_vector);
}

arma::vec optimNMAdapter::convertFromSTLToArma(const std::vector<double>& std_vector) {
    return arma::conv_to<arma::vec>::from(std_vector);
}
}  // namespace nonlinear_solver