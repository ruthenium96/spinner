#include "DenseMatrix.h"

#include <armadillo>

struct DenseMatrix::Impl {
    arma::dmat RawData;
public:
    Impl() = default;
};

DenseMatrix::DenseMatrix()
: pImpl{std::make_unique<DenseMatrix::Impl>()} {
}

DenseMatrix::~DenseMatrix() = default;
DenseMatrix::DenseMatrix(DenseMatrix&&) noexcept = default;
DenseMatrix& DenseMatrix::operator=(DenseMatrix&&) noexcept = default;

void DenseMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    pImpl->RawData(i, j) += value;
}

std::ostream &operator<<(std::ostream &os, const DenseMatrix &decomposition) {
    os << decomposition.pImpl->RawData << std::endl;
    return os;
}

void DenseMatrix::resize(uint32_t matrix_in_space_basis_size_i, uint32_t matrix_in_space_basis_size_j) {
    // TODO: is it the fastest way to initialize pImpl->RawData?
    pImpl->RawData.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_j);
    pImpl->RawData.fill(arma::fill::zeros);
}

void DenseMatrix::diagonalize(DenseVector& values, DenseMatrix& vectors) const {
    arma::eig_sym(values.pImpl->eigenvalues, vectors.pImpl->RawData, pImpl->RawData);
}

DenseMatrix DenseMatrix::unitary_transform(const DenseMatrix& matrix_to_transform) {
    DenseMatrix transformed_matrix;
    transformed_matrix.pImpl->RawData = pImpl->RawData.t() * matrix_to_transform.pImpl->RawData * pImpl->RawData;
    return std::move(transformed_matrix);
}

DenseVector DenseMatrix::return_main_diagonal() {
    DenseVector main_diagonal;
    main_diagonal.pImpl->eigenvalues = pImpl->RawData.diag();
    return std::move(main_diagonal);
}
