#ifndef SPINNER_LBFGSPPADAPTER_H
#define SPINNER_LBFGSPPADAPTER_H

#include "src/nonlinear_solver/AbstractNonlinearSolver.h"

namespace nonlinear_solver {

class LBFGSppAdapter: public AbstractNonlinearSolver {
  public:
    void optimize(
        std::function<double(const std::vector<double>&, std::vector<double>&, bool)>
            oneStepFunction,
        std::vector<double>& changeable_values) override;
    bool doesGradientsRequired() const override {
        return true;
    };
};

}  // namespace nonlinear_solver

#endif  //SPINNER_LBFGSPPADAPTER_H
