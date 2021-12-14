#include <armadillo>

#include "SparseMatrix.h"

namespace lexicographic {

struct SparseMatrix::Impl {
    arma::sp_mat sparse_matrix;

  public:
    Impl() = default;
};

SparseMatrix::SparseMatrix(lexicographic::IndexConverter converter) :
    pImpl {std::make_unique<Impl>()},
    converter_(std::move(converter)) {
    resize(converter_.get_total_space_size());
}
SparseMatrix::~SparseMatrix() = default;
SparseMatrix::SparseMatrix(SparseMatrix&&) noexcept = default;
SparseMatrix& SparseMatrix::operator=(SparseMatrix&&) noexcept = default;

std::ostream& operator<<(std::ostream& os, const SparseMatrix& data) {
    os << data.pImpl->sparse_matrix << std::endl;
    return os;
}

uint32_t SparseMatrix::size() const {
    return pImpl->sparse_matrix.n_cols;
}

void SparseMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    pImpl->sparse_matrix(i, j) += value;
}

double SparseMatrix::operator()(uint32_t i, uint32_t j) const {
    return pImpl->sparse_matrix(i, j);
}

void SparseMatrix::resize(uint32_t new_size) {
    pImpl->sparse_matrix.resize(new_size, new_size);
}

void SparseMatrix::add_scalar_product(
    uint32_t index_of_vector,
    uint32_t center_a,
    uint32_t center_b,
    double factor) {
    uint32_t projection_of_center_a =
        converter_.convert_lex_index_to_one_sz_projection(index_of_vector, center_a);
    uint32_t projection_of_center_b =
        converter_.convert_lex_index_to_one_sz_projection(index_of_vector, center_b);

    // (Sa, Sb) = Sax*Sbx + Say*Sby + Saz*Sbz = 0.5 * (Sa+Sb- + Sa-Sb+) + Saz*Sbz

    // Saz Sbz
    double diagonal_value = (projection_of_center_a - converter_.get_spins()[center_a])
        * (projection_of_center_b - converter_.get_spins()[center_b]) * factor;
    add_to_position(diagonal_value, index_of_vector, index_of_vector);

    // Sa+ Sb-
    add_scalar_product_nondiagonal_part(
        index_of_vector,
        center_a,
        center_b,
        projection_of_center_a,
        projection_of_center_b,
        factor);

    // Sa- Sb+
    add_scalar_product_nondiagonal_part(
        index_of_vector,
        center_b,
        center_a,
        projection_of_center_b,
        projection_of_center_a,
        factor);
}

void SparseMatrix::add_scalar_product_nondiagonal_part(
    uint32_t index_of_vector,
    uint32_t plus_center,
    uint32_t minus_center,
    uint32_t projection_of_plus_center,
    uint32_t projection_of_minus_center,
    double factor) {
    if (projection_of_plus_center == converter_.get_mults()[plus_center] - 1
        || projection_of_minus_center == 0) {
        return;
    }
    uint32_t index_of_new_vector = converter_.ladder_projection(
        converter_.ladder_projection(index_of_vector, plus_center, +1),
        minus_center,
        -1);

    // projection m = number n - spin S
    // so S(S+1)-m(m+1) = (2S-n)(n+1)
    // so S(S+1)-m(m-1) = n(2S+1-n)
    double factor_a = (2 * converter_.get_spins()[plus_center] - projection_of_plus_center)
        * (projection_of_plus_center + 1);
    double factor_b = projection_of_minus_center
        * (2 * converter_.get_spins()[minus_center] + 1 - projection_of_minus_center);

    // TODO: fix plus-minus
    add_to_position(0.5 * sqrt(factor_a * factor_b) * factor, index_of_vector, index_of_new_vector);
}
}  // namespace lexicographic