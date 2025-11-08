#ifndef SPINNER_EIGENDENSESEMIUNITARYMATRIX_H
#define SPINNER_EIGENDENSESEMIUNITARYMATRIX_H

#include <Eigen/Core>

#include "src/entities/data_structures/AbstractDenseSemiunitaryMatrix.h"
#include "src/entities/data_structures/AbstractDenseSemiunitaryTransformer.h"

namespace quantum::linear_algebra {
template <typename T>
class EigenDenseSemiunitaryMatrix: public AbstractDenseSemiunitaryMatrix {
  public:
    EigenDenseSemiunitaryMatrix();

    uint32_t size_rows() const override;
    uint32_t size_cols() const override;
    void resize(size_t size_rows, size_t size_cols);
    double at(uint32_t i, uint32_t j) const override;

    void print(std::ostream& os) const override;

    const Eigen::Matrix<T, -1, -1>& getDenseSemiunitaryMatrix() const;
    Eigen::Matrix<T, -1, -1>& modifyDenseSemiunitaryMatrix();

    const std::unique_ptr<AbstractDenseSemiunitaryTransformer>& getUnitaryTransformer() const override;

  private:
    Eigen::Matrix<T, -1, -1> denseSemiunitaryMatrix_;

    std::unique_ptr<AbstractDenseSemiunitaryTransformer> transformer_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_EIGENDENSESEMIUNITARYMATRIX_H
