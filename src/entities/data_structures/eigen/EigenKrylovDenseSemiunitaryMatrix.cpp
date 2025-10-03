#include "EigenKrylovDenseSemiunitaryMatrix.h"

#include "EigenKrylovDenseSemiunitaryTransformer.h"

namespace quantum::linear_algebra {

template <typename T>
EigenKrylovDenseSemiunitaryMatrix<T>::EigenKrylovDenseSemiunitaryMatrix() {
    transformer_ = std::make_unique<EigenKrylovDenseSemiunitaryTransformer<T>>(this);
}

template <typename T>
const std::unique_ptr<AbstractDenseSemiunitaryTransformer>& EigenKrylovDenseSemiunitaryMatrix<T>::getUnitaryTransformer() const {
    return transformer_;
}

template <typename T>
uint32_t EigenKrylovDenseSemiunitaryMatrix<T>::size_rows() const {
    return denseKrylovSemiunitaryMatrix_.rows();
}

template <typename T>
uint32_t EigenKrylovDenseSemiunitaryMatrix<T>::size_cols() const {
    return denseKrylovSemiunitaryMatrix_.cols();
}

template <typename T>
void EigenKrylovDenseSemiunitaryMatrix<T>::resize(size_t size_rows, size_t size_cols) {
    denseKrylovSemiunitaryMatrix_.resize(size_rows, size_cols);
    denseKrylovSemiunitaryMatrix_.fill(0);
}

template <typename T>
double EigenKrylovDenseSemiunitaryMatrix<T>::at(uint32_t i, uint32_t j) const {
    return denseKrylovSemiunitaryMatrix_(j, i);
}

template <typename T>
void EigenKrylovDenseSemiunitaryMatrix<T>::print(std::ostream& os) const {
    os << denseKrylovSemiunitaryMatrix_ << std::endl;
}

template <typename T>
std::unique_ptr<AbstractDiagonalizableMatrix> EigenKrylovDenseSemiunitaryMatrix<T>::unitaryTransform(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const {
    throw std::invalid_argument("EigenKrylovDenseSemiunitaryMatrix cannot do unitaryTransform");
}

template <typename T>
const Eigen::Matrix<T, -1, -1>& EigenKrylovDenseSemiunitaryMatrix<T>::getKrylovDenseSemiunitaryMatrix() const {
    return denseKrylovSemiunitaryMatrix_;
}

template <typename T>
Eigen::Matrix<T, -1, -1>& EigenKrylovDenseSemiunitaryMatrix<T>::modifyKrylovDenseSemiunitaryMatrix() {
    return denseKrylovSemiunitaryMatrix_;
}

template <typename T>
const Eigen::Vector<T, -1>& EigenKrylovDenseSemiunitaryMatrix<T>::getBackProjectionVector() const {
    return backProjectionVector_;
}

template <typename T>
Eigen::Vector<T, -1>& EigenKrylovDenseSemiunitaryMatrix<T>::modifyBackProjectionVector() {
    return backProjectionVector_;
}

template <typename T>
const Eigen::Vector<T, -1>& EigenKrylovDenseSemiunitaryMatrix<T>::getSeedVector() const {
    return seed_vector_;
}

template <typename T>
Eigen::Vector<T, -1>& EigenKrylovDenseSemiunitaryMatrix<T>::modifySeedVector() {
    return seed_vector_;
}

template <typename T>
void EigenKrylovDenseSemiunitaryMatrix<T>::add_to_position(double value, uint32_t i, uint32_t j) {
    throw std::invalid_argument("EigenKrylovDenseSemiunitaryMatrix cannot do add_to_position");
}

template <typename T>
void EigenKrylovDenseSemiunitaryMatrix<T>::normalize() {
    throw std::invalid_argument("EigenKrylovDenseSemiunitaryMatrix cannot do normalize");
}

template class EigenKrylovDenseSemiunitaryMatrix<double>;
template class EigenKrylovDenseSemiunitaryMatrix<float>;
}  // namespace quantum::linear_algebra