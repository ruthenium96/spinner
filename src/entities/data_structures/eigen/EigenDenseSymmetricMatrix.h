#ifndef SPINNER_EIGENDENSESYMMETRICMATRIX_H
#define SPINNER_EIGENDENSESYMMETRICMATRIX_H

#include <Eigen/Core>

#include "src/entities/data_structures/AbstractSymmetricMatrix.h"

namespace quantum::linear_algebra {
class EigenDenseSymmetricMatrix: public AbstractSymmetricMatrix {
  public:
    void add_to_position(double value, uint32_t i, uint32_t j) override;
    EigenCouple diagonalizeValuesVectors() const override;
    std::unique_ptr<AbstractDenseVector> diagonalizeValues() const override;
    uint32_t size() const override;
    double at(uint32_t i, uint32_t j) const override;
    void print(std::ostream& os) const override;
    void resize(uint32_t matrix_in_space_basis_size_i);

    const Eigen::MatrixXd& getDenseSymmetricMatrix() const;
    Eigen::MatrixXd& modifyDenseSymmetricMatrix();

  private:
    Eigen::MatrixXd denseSymmetricMatrix_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_EIGENDENSESYMMETRICMATRIX_H
