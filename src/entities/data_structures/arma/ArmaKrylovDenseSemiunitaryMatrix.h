#ifndef SPINNER_ARMAKRYLOVDENSESEMIUNITARYMATRIX_H
#define SPINNER_ARMAKRYLOVDENSESEMIUNITARYMATRIX_H

#include <armadillo>

#include "src/entities/data_structures/AbstractDenseSemiunitaryMatrix.h"
#include "src/entities/data_structures/AbstractDenseSemiunitaryTransformer.h"

namespace quantum::linear_algebra {
template <typename T>
class ArmaKrylovDenseSemiunitaryMatrix: public AbstractDenseSemiunitaryMatrix {
  public:
    ArmaKrylovDenseSemiunitaryMatrix();

    uint32_t size_rows() const override;
    uint32_t size_cols() const override;
    double at(uint32_t i, uint32_t j) const override;

    void print(std::ostream& os) const override;
    void resize(size_t size_rows, size_t size_cols);
    const arma::Mat<T>& getKrylovDenseSemiunitaryMatrix() const;
    arma::Mat<T>& modifyKrylovDenseSemiunitaryMatrix();
    const arma::Col<T>& getBackProjectionVector() const;
    arma::Col<T>& modifyBackProjectionVector();
    const arma::Col<T>& getSeedVector() const;
    arma::Col<T>& modifySeedVector();

    const std::unique_ptr<AbstractDenseSemiunitaryTransformer>& getUnitaryTransformer() const override;

  private:
    arma::Mat<T> denseKrylovSemiunitaryMatrix_;
    arma::Col<T> backProjectionVector_;
    arma::Col<T> seed_vector_;

    std::unique_ptr<AbstractDenseSemiunitaryTransformer> transformer_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ARMAKRYLOVDENSESEMIUNITARYMATRIX_H
