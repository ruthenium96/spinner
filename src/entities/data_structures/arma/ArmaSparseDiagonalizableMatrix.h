#ifndef SPINNER_ARMASPARSEDIAGONALIZABLEMATRIX_H
#define SPINNER_ARMASPARSEDIAGONALIZABLEMATRIX_H

#include <armadillo>

#include "src/entities/data_structures/AbstractDiagonalizableMatrix.h"

namespace quantum::linear_algebra {
class ArmaSparseDiagonalizableMatrix: public AbstractDiagonalizableMatrix {
  public:
    void add_to_position(double value, uint32_t i, uint32_t j) override;
    void resize(uint32_t matrix_in_space_basis_size_i);
    EigenCouple diagonalizeValuesVectors() const override;
    std::unique_ptr<AbstractDenseVector> diagonalizeValues() const override;
    std::unique_ptr<AbstractDiagonalizableMatrix> multiply_by(double multiplier) const override;
    uint32_t size() const override;
    double at(uint32_t i, uint32_t j) const override;
    void print(std::ostream& os) const override;
    const arma::sp_mat& getSparseSymmetricMatrix() const;
    arma::sp_mat& modifySparseSymmetricMatrix();

  private:
    arma::sp_mat sparseDiagonalizableMatrix_;
};

using ArmaSparseSymmetricMatrix = ArmaSparseDiagonalizableMatrix;

}  // namespace quantum::linear_algebra
#endif  //SPINNER_ARMASPARSEDIAGONALIZABLEMATRIX_H
