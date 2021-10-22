#include "SubmatrixData.h"

#include <armadillo>

struct SubmatrixData::Impl {
    arma::dmat RawData;
public:
    Impl() = default;
};

SubmatrixData::SubmatrixData()
: pImpl{std::make_unique<Impl>()} {
}

SubmatrixData::~SubmatrixData() = default;
SubmatrixData::SubmatrixData(SubmatrixData&&) noexcept = default;
SubmatrixData& SubmatrixData::operator=(SubmatrixData&&) noexcept = default;

void SubmatrixData::add_to_position(double value, uint32_t i, uint32_t j) {
    pImpl->RawData(i, j) += value;
}

std::ostream &operator<<(std::ostream &os, const SubmatrixData &decomposition) {
    os << decomposition.pImpl->RawData << std::endl;
    return os;
}

void SubmatrixData::resize(uint32_t matrix_in_space_basis_size_i, uint32_t matrix_in_space_basis_size_j) {
    // TODO: is it the fastest way to initialize pImpl->RawData?
    pImpl->RawData.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_j);
    pImpl->RawData.fill(arma::fill::zeros);
}

void SubmatrixData::diagonalize() {
    arma::vec eigenvalues;
    arma::dmat eigenvectors;
    arma::eig_sym(eigenvalues, eigenvectors, pImpl->RawData);
    std::cout << eigenvalues << std::endl;
    std::cout << eigenvectors << std::endl;
}
