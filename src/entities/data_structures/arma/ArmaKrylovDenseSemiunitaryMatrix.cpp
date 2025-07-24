#include "ArmaKrylovDenseSemiunitaryMatrix.h"
#include "ArmaKrylovDenseSemiunitaryTransformer.h"
#include <stdexcept>

namespace quantum::linear_algebra {

template <typename T>
ArmaKrylovDenseSemiunitaryMatrix<T>::ArmaKrylovDenseSemiunitaryMatrix() {
    transformer_ = std::make_unique<ArmaKrylovDenseSemiunitaryTransformer<T>>(this);
}

template <typename T>
const std::unique_ptr<AbstractDenseSemiunitaryTransformer>& ArmaKrylovDenseSemiunitaryMatrix<T>::getUnitaryTransformer() const {
    return transformer_;
}

template <typename T>
uint32_t ArmaKrylovDenseSemiunitaryMatrix<T>::size_rows() const {
    return denseKrylovSemiunitaryMatrix_.n_rows;
}

template <typename T>
uint32_t ArmaKrylovDenseSemiunitaryMatrix<T>::size_cols() const {
    return denseKrylovSemiunitaryMatrix_.n_cols;
}

template <typename T>
double ArmaKrylovDenseSemiunitaryMatrix<T>::at(uint32_t i, uint32_t j) const {
    return denseKrylovSemiunitaryMatrix_.at(j, i);
}

template <typename T>
void ArmaKrylovDenseSemiunitaryMatrix<T>::print(std::ostream& os) const {
    os << "DenseKrylovSemiunitaryMatrix:\n" << denseKrylovSemiunitaryMatrix_ << std::endl;
    os << "BackProjectionVector:\n" << backProjectionVector_ << std::endl;
    os << "SeedVector:\n" << seed_vector_ << std::endl;
}

template <typename T>
void ArmaKrylovDenseSemiunitaryMatrix<T>::resize(size_t size_rows, size_t size_cols) {
    denseKrylovSemiunitaryMatrix_.resize(size_rows, size_cols);
    denseKrylovSemiunitaryMatrix_.fill(arma::fill::zeros);
}

template <typename T>
std::unique_ptr<AbstractDiagonalizableMatrix> ArmaKrylovDenseSemiunitaryMatrix<T>::unitaryTransform(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const {
    throw std::invalid_argument("ArmaKrylovDenseSemiunitaryMatrix cannot do unitaryTransform");
}

template <typename T>
const arma::Mat<T>& ArmaKrylovDenseSemiunitaryMatrix<T>::getKrylovDenseSemiunitaryMatrix() const {
    return denseKrylovSemiunitaryMatrix_;
}

template <typename T>
arma::Mat<T>& ArmaKrylovDenseSemiunitaryMatrix<T>::modifyKrylovDenseSemiunitaryMatrix() {
    return denseKrylovSemiunitaryMatrix_;
}

template <typename T>
const arma::Col<T>& ArmaKrylovDenseSemiunitaryMatrix<T>::getSeedVector() const {
    return seed_vector_;
}

template <typename T>
arma::Col<T>& ArmaKrylovDenseSemiunitaryMatrix<T>::modifySeedVector() {
    return seed_vector_;
}

template <typename T>
const arma::Col<T>& ArmaKrylovDenseSemiunitaryMatrix<T>::getBackProjectionVector() const {
    return backProjectionVector_;
}

template <typename T>
arma::Col<T>& ArmaKrylovDenseSemiunitaryMatrix<T>::modifyBackProjectionVector() {
    return backProjectionVector_;
}

template <typename T>
void ArmaKrylovDenseSemiunitaryMatrix<T>::add_to_position(double value, uint32_t i, uint32_t j) {
    throw std::invalid_argument("ArmaKrylovDenseSemiunitaryMatrix cannot do add_to_position");
}

template <typename T>
void ArmaKrylovDenseSemiunitaryMatrix<T>::normalize() {
    throw std::invalid_argument("ArmaKrylovDenseSemiunitaryMatrix cannot do normalize");
}

template class ArmaKrylovDenseSemiunitaryMatrix<double>;
template class ArmaKrylovDenseSemiunitaryMatrix<float>;
}  // namespace quantum::linear_algebra