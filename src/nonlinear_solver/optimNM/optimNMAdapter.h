#ifndef SPINNER_OPTIMNMADAPTER_H
#define SPINNER_OPTIMNMADAPTER_H

#include <armadillo>

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

  private:
    std::function<double(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)>
    adaptSignature(
        const std::function<double(const std::vector<double>&, std::vector<double>&, bool)>&
            oneStepFunction);
    static std::vector<double> convertFromArmaToSTL(const arma::vec& arma_vector);
    static arma::vec convertFromSTLToArma(const std::vector<double>& std_vector);
};

}  // namespace nonlinear_solver

#endif  //SPINNER_OPTIMNMADAPTER_H
