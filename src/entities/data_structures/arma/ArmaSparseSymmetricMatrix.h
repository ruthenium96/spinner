#ifndef SPINNER_ARMASPARSESYMMETRICMATRIX_H
#define SPINNER_ARMASPARSESYMMETRICMATRIX_H

#include <armadillo>

#include "src/entities/data_structures/AbstractSymmetricMatrix.h"

namespace quantum::linear_algebra {
class ArmaSparseSymmetricMatrix: public AbstractSymmetricMatrix {
  public:
    void add_to_position(double value, uint32_t i, uint32_t j) override;
    void resize(uint32_t matrix_in_space_basis_size_i);
    EigenCouple diagonalizeValuesVectors() const override;
    std::unique_ptr<AbstractDenseVector> diagonalizeValues() const override;
    uint32_t size() const override;
    double at(uint32_t i, uint32_t j) const override;
    void print(std::ostream& os) const override;
    const arma::sp_mat& getSparseSymmetricMatrix() const;
    arma::sp_mat& modifySparseSymmetricMatrix();

  private:
    arma::sp_mat sparseSymmetricMatrix_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ARMASPARSESYMMETRICMATRIX_H
