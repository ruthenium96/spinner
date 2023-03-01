#ifndef SPINNER_EIGENSPARSESYMMETRICMATRIX_H
#define SPINNER_EIGENSPARSESYMMETRICMATRIX_H

#include "Eigen/Sparse"
#include "src/entities/data_structures/AbstractSymmetricMatrix.h"
namespace quantum::linear_algebra {
class EigenSparseSymmetricMatrix: public AbstractSymmetricMatrix {
  public:
    void add_to_position(double value, uint32_t i, uint32_t j) override;
    EigenCouple diagonalizeValuesVectors() const override;
    std::unique_ptr<AbstractDenseVector> diagonalizeValues() const override;
    uint32_t size() const override;
    double at(uint32_t i, uint32_t j) const override;
    void print(std::ostream& os) const override;
    ~EigenSparseSymmetricMatrix() override = default;
    void resize(size_t size);

    const Eigen::SparseMatrix<double>& getSparseSymmetricMatrix() const;
    Eigen::SparseMatrix<double>& modifySparseSymmetricMatrix();

  private:
    Eigen::SparseMatrix<double> sparseSymmetricMatrix_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_EIGENSPARSESYMMETRICMATRIX_H
