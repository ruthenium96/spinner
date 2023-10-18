#include "ArmaSparseSemiunitaryMatrix.h"

namespace quantum::linear_algebra {

struct IteratorImpl: public AbstractSparseSemiunitaryMatrix::Iterator {
    arma::sp_mat::const_col_iterator iter;
    const arma::sp_mat::const_col_iterator end;
    const size_t size_;

    IteratorImpl(
        arma::sp_mat::const_col_iterator iter1,
        arma::sp_mat::const_col_iterator iter2,
        size_t size) :
        iter(iter1),
        end(iter2),
        size_(size) {}

    bool hasNext() const override {
        return iter != end;
    }

    IndexValueItem getNext() override {
        auto value = *iter;
        auto row = static_cast<uint32_t>(iter.row());
        ++iter;
        return {row, value};
    }

    size_t size() const override {
        return size_;
    }

    ~IteratorImpl() override = default;
};

std::unique_ptr<AbstractSparseSemiunitaryMatrix::Iterator>
ArmaSparseSemiunitaryMatrix::GetNewIterator(size_t index_of_vector) const {
    return std::make_unique<IteratorImpl>(
        sparseSemiunitaryMatrix_.begin_col(index_of_vector),
        sparseSemiunitaryMatrix_.end_col(index_of_vector),
        sparseSemiunitaryMatrix_.col(index_of_vector).n_nonzero);
}

uint32_t ArmaSparseSemiunitaryMatrix::size_rows() const {
    return sparseSemiunitaryMatrix_.n_rows;
}

uint32_t ArmaSparseSemiunitaryMatrix::size_cols() const {
    return sparseSemiunitaryMatrix_.n_cols;
}

bool ArmaSparseSemiunitaryMatrix::empty() const {
    return sparseSemiunitaryMatrix_.n_nonzero == 0;
}

bool ArmaSparseSemiunitaryMatrix::vempty(uint32_t index_of_vector) const {
    return sparseSemiunitaryMatrix_.col(index_of_vector).n_nonzero == 0;
}

void ArmaSparseSemiunitaryMatrix::clear() {
    sparseSemiunitaryMatrix_.clear();
}

void ArmaSparseSemiunitaryMatrix::eraseExplicitZeros() {
    sparseSemiunitaryMatrix_.clean(0.001);
}

bool ArmaSparseSemiunitaryMatrix::is_zero(uint32_t i, uint32_t j) const {
    return (sparseSemiunitaryMatrix_(j, i) == 0);
}

void ArmaSparseSemiunitaryMatrix::move_vector_from(
    uint32_t i,
    std::unique_ptr<AbstractSparseSemiunitaryMatrix>& subspace_from) {
    auto mb_rhs = dynamic_cast<ArmaSparseSemiunitaryMatrix*>(subspace_from.get());
    if (mb_rhs == nullptr) {
        throw std::bad_cast();
    }

    auto old_in_rows = sparseSemiunitaryMatrix_.n_rows;
    auto old_in_cols = sparseSemiunitaryMatrix_.n_cols;
    sparseSemiunitaryMatrix_.resize(old_in_rows, old_in_cols + 1);

    sparseSemiunitaryMatrix_.col(old_in_cols) = mb_rhs->sparseSemiunitaryMatrix_.col(i);
}

void ArmaSparseSemiunitaryMatrix::resize(uint32_t cols, uint32_t rows) {
    sparseSemiunitaryMatrix_.resize(rows, cols);
}

void ArmaSparseSemiunitaryMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    sparseSemiunitaryMatrix_(j, i) += value;
}

double ArmaSparseSemiunitaryMatrix::at(uint32_t i, uint32_t j) const {
    return sparseSemiunitaryMatrix_(j, i);
}

void ArmaSparseSemiunitaryMatrix::normalize() {
    auto tmpSparseSemiunitaryMatrix_ = arma::normalise(sparseSemiunitaryMatrix_);
    sparseSemiunitaryMatrix_ = tmpSparseSemiunitaryMatrix_;
}

void ArmaSparseSemiunitaryMatrix::print(std::ostream& os) const {
    os << sparseSemiunitaryMatrix_ << std::endl;
}

const arma::sp_mat& ArmaSparseSemiunitaryMatrix::getSparseSemiunitaryMatrix() const {
    return sparseSemiunitaryMatrix_;
}

void ArmaSparseSemiunitaryMatrix::unitaryTransform(
    const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrixToTransform,
    std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrixToAdd) const {
    throw std::invalid_argument("unitaryTransform not implemented in ArmaSparseSemiunitaryMatrix");
}
}  // namespace quantum::linear_algebra