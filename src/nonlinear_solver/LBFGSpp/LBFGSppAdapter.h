#ifndef SPINNER_LBFGSPPADAPTER_H
#define SPINNER_LBFGSPPADAPTER_H

#include <Eigen/Core>

#include "src/nonlinear_solver/AbstractNonlinearSolver.h"

namespace nonlinear_solver {

class LBFGSppAdapter: public AbstractNonlinearSolver {
  public:
    void optimize(
        std::function<double(const std::vector<double>&, std::vector<double>&, bool)>
            oneStepFunction,
        std::vector<double>& changeable_values) override;

  private:
    std::function<double(const Eigen::VectorXd&, Eigen::VectorXd&)> adaptSignature(
        const std::function<double(const std::vector<double>&, std::vector<double>&, bool)>&
            oneStepFunction);
    static std::vector<double> convertFromEigenToSTL(const Eigen::VectorXd& eigen_vector);
    static Eigen::VectorXd convertFromSTLToEigen(const std::vector<double>& std_vector);
};

}  // namespace nonlinear_solver

#endif  //SPINNER_LBFGSPPADAPTER_H
