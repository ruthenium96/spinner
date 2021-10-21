#include "RawSubmatrixData.h"

#include <armadillo>

struct RawSubmatrixData::Impl {
    arma::dmat RawData;
public:
    Impl() = default;
};

RawSubmatrixData::RawSubmatrixData()
: pImpl{std::make_unique<Impl>()} {
}

RawSubmatrixData::~RawSubmatrixData() = default;
RawSubmatrixData::RawSubmatrixData(RawSubmatrixData&&) noexcept = default;
RawSubmatrixData& RawSubmatrixData::operator=(RawSubmatrixData&&) noexcept = default;

void RawSubmatrixData::add_to_position(double value, uint32_t i, uint32_t j) {
    pImpl->RawData(i, j) += value;
}

std::ostream &operator<<(std::ostream &os, const RawSubmatrixData &decomposition) {
    os << decomposition.pImpl->RawData << std::endl;
    return os;
}

void RawSubmatrixData::resize(uint32_t matrix_in_space_basis_size_i, uint32_t matrix_in_space_basis_size_j) {
    // TODO: is it the fastest way to initialize pImpl->RawData?
    pImpl->RawData.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_j);
    pImpl->RawData.fill(arma::fill::zeros);
}

void RawSubmatrixData::diagonalize() {
    arma::vec eigenvalues;
    arma::dmat eigenvectors;
    arma::eig_sym(eigenvalues, eigenvectors, pImpl->RawData);
    std::cout << eigenvalues << std::endl;
    std::cout << eigenvectors << std::endl;
}
