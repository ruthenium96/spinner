#ifndef SPINNER_EIGENKRYLOVDENSESEMIUNITARYMATRIX_H
#define SPINNER_EIGENKRYLOVDENSESEMIUNITARYMATRIX_H

#include <Eigen/Core>

#include "src/entities/data_structures/AbstractDenseSemiunitaryMatrix.h"
#include "src/entities/data_structures/AbstractDenseSemiunitaryTransformer.h"

namespace quantum::linear_algebra {
template <typename T>
class EigenKrylovDenseSemiunitaryMatrix: public AbstractDenseSemiunitaryMatrix {
  public:
    EigenKrylovDenseSemiunitaryMatrix();

    uint32_t size_rows() const override;
    uint32_t size_cols() const override;
    void resize(size_t size_rows, size_t size_cols);
    double at(uint32_t i, uint32_t j) const override;

    void print(std::ostream& os) const override;
    std::unique_ptr<AbstractDiagonalizableMatrix> unitaryTransform(
        const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const override;

    const Eigen::Matrix<T, -1, -1>& getKrylovDenseSemiunitaryMatrix() const;
    Eigen::Matrix<T, -1, -1>& modifyKrylovDenseSemiunitaryMatrix();
    const Eigen::Vector<T, -1>& getBackProjectionVector() const;
    Eigen::Vector<T, -1>& modifyBackProjectionVector();
    const Eigen::Vector<T, -1>& getSeedVector() const;
    Eigen::Vector<T, -1>& modifySeedVector();

    const std::unique_ptr<AbstractDenseSemiunitaryTransformer>& getUnitaryTransformer() const override;

  private:
    Eigen::Matrix<T, -1, -1> denseKrylovSemiunitaryMatrix_;
    Eigen::Vector<T, -1> backProjectionVector_;
    Eigen::Vector<T, -1> seed_vector_;

    std::unique_ptr<AbstractDenseSemiunitaryTransformer> transformer_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_EIGENKRYLOVDENSESEMIUNITARYMATRIX_H


