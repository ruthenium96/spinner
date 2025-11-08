#include "ArmaDenseSemiunitaryMatrix.h"

#include "ArmaLogic.h"
#include "ArmaDenseSemiunitaryTransformer.h"

namespace quantum::linear_algebra {
template <typename T>
ArmaDenseSemiunitaryMatrix<T>::ArmaDenseSemiunitaryMatrix() {
    transformer_ = std::make_unique<ArmaDenseSemiunitaryTransformer<T>>(this);
}

template <typename T>
const std::unique_ptr<AbstractDenseSemiunitaryTransformer>& ArmaDenseSemiunitaryMatrix<T>::getUnitaryTransformer() const {
    return transformer_;
}

template <typename T>
uint32_t ArmaDenseSemiunitaryMatrix<T>::size_rows() const {
    return denseSemiunitaryMatrix_.n_rows;
}

template <typename T>
uint32_t ArmaDenseSemiunitaryMatrix<T>::size_cols() const {
    return denseSemiunitaryMatrix_.n_cols;
}

template <typename T>
double ArmaDenseSemiunitaryMatrix<T>::at(uint32_t i, uint32_t j) const {
    return denseSemiunitaryMatrix_.at(j, i);
}

template <typename T>
void ArmaDenseSemiunitaryMatrix<T>::print(std::ostream& os) const {
    os << denseSemiunitaryMatrix_ << std::endl;
}

template <typename T>
void ArmaDenseSemiunitaryMatrix<T>::resize(size_t size_rows, size_t size_cols) {
    denseSemiunitaryMatrix_.resize(size_rows, size_cols);
    denseSemiunitaryMatrix_.fill(arma::fill::zeros);
}

template <typename T>
const arma::Mat<T>& ArmaDenseSemiunitaryMatrix<T>::getDenseSemiunitaryMatrix() const {
    return denseSemiunitaryMatrix_;
}

template <typename T>
arma::Mat<T>& ArmaDenseSemiunitaryMatrix<T>::modifyDenseSemiunitaryMatrix() {
    return denseSemiunitaryMatrix_;
}

template class ArmaDenseSemiunitaryMatrix<double>;
template class ArmaDenseSemiunitaryMatrix<float>;
}  // namespace quantum::linear_algebra