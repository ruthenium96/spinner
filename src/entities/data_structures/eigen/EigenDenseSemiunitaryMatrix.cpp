#include "EigenDenseSemiunitaryMatrix.h"

#include "EigenDenseVector.h"
#include "EigenLogic.h"
#include "EigenDenseSemiunitaryTransformer.h"

namespace quantum::linear_algebra {

template <typename T>
EigenDenseSemiunitaryMatrix<T>::EigenDenseSemiunitaryMatrix() {
    transformer_ = std::make_unique<EigenDenseSemiunitaryTransformer<T>>(this);
}

template <typename T>
const std::unique_ptr<AbstractDenseSemiunitaryTransformer>& EigenDenseSemiunitaryMatrix<T>::getUnitaryTransformer() const {
    return transformer_;
}

template <typename T>
uint32_t EigenDenseSemiunitaryMatrix<T>::size_rows() const {
    return denseSemiunitaryMatrix_.rows();
}

template <typename T>
uint32_t EigenDenseSemiunitaryMatrix<T>::size_cols() const {
    return denseSemiunitaryMatrix_.cols();
}

template <typename T>
void EigenDenseSemiunitaryMatrix<T>::resize(size_t size_rows, size_t size_cols) {
    denseSemiunitaryMatrix_.resize(size_rows, size_cols);
    denseSemiunitaryMatrix_.fill(0);
}

template <typename T>
double EigenDenseSemiunitaryMatrix<T>::at(uint32_t i, uint32_t j) const {
    return denseSemiunitaryMatrix_(j, i);
}

template <typename T>
void EigenDenseSemiunitaryMatrix<T>::print(std::ostream& os) const {
    os << denseSemiunitaryMatrix_ << std::endl;
}

template <typename T>
const Eigen::Matrix<T, -1, -1>& EigenDenseSemiunitaryMatrix<T>::getDenseSemiunitaryMatrix() const {
    return denseSemiunitaryMatrix_;
}

template <typename T>
Eigen::Matrix<T, -1, -1>& EigenDenseSemiunitaryMatrix<T>::modifyDenseSemiunitaryMatrix() {
    return denseSemiunitaryMatrix_;
}

template class EigenDenseSemiunitaryMatrix<double>;
template class EigenDenseSemiunitaryMatrix<float>;
}  // namespace quantum::linear_algebra