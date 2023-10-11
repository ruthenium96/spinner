#include "StdSparseSemiunitaryMatrix.h"

#include <algorithm>
#include <cmath>

namespace quantum::linear_algebra {

struct IteratorStdImpl: public AbstractSparseSemiunitaryMatrix::Iterator {
    StdSparseSemiunitaryMatrix::Map::const_iterator iter;
    const StdSparseSemiunitaryMatrix::Map::const_iterator end;

    IteratorStdImpl(
        StdSparseSemiunitaryMatrix::Map::const_iterator iter1,
        StdSparseSemiunitaryMatrix::Map::const_iterator iter2) :
        iter(iter1),
        end(iter2) {}

    bool hasNext() const override {
        return iter != end;
    }

    IndexValueItem getNext() override {
        auto pair = *iter;
        ++iter;
        return {pair.first, pair.second};
    }

    ~IteratorStdImpl() override = default;
};

std::unique_ptr<AbstractSparseSemiunitaryMatrix::Iterator>
StdSparseSemiunitaryMatrix::GetNewIterator(size_t index_of_vector) const {
    return std::make_unique<IteratorStdImpl>(
        basis_[index_of_vector].cbegin(),
        basis_[index_of_vector].cend());
}

uint32_t StdSparseSemiunitaryMatrix::size_cols() const {
    return basis_.size();
}

uint32_t StdSparseSemiunitaryMatrix::size_rows() const {
    return rows_;
}

bool StdSparseSemiunitaryMatrix::empty() const {
    return basis_.empty();
}

bool StdSparseSemiunitaryMatrix::vempty(uint32_t index_of_vector) const {
    return basis_[index_of_vector].empty();
}

void StdSparseSemiunitaryMatrix::clear() {
    basis_.clear();
}

void StdSparseSemiunitaryMatrix::move_vector_from(
    uint32_t i,
    std::unique_ptr<AbstractSparseSemiunitaryMatrix>& subspace_from) {
    auto rhs_ = downcast_ptr(subspace_from);
    basis_.emplace_back(std::move(rhs_->basis_[i]));
}

void StdSparseSemiunitaryMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    if (!basis_[i].contains(j)) {
        basis_[i].emplace(j, value);
    } else {
        basis_[i][j] += value;
    }
}

double StdSparseSemiunitaryMatrix::at(uint32_t i, uint32_t j) const {
    if (!basis_.at(i).contains(j)) {
        return 0;
    }
    return basis_.at(i).at(j);
}

void StdSparseSemiunitaryMatrix::resize(uint32_t cols, uint32_t rows) {
    basis_.resize(cols);
    rows_ = rows;
}

bool StdSparseSemiunitaryMatrix::is_zero(uint32_t i, uint32_t j) const {
    return !basis_.at(i).contains(j);
}

void StdSparseSemiunitaryMatrix::eraseExplicitZeros() {
    for (auto& mm : basis_) {
        mm.erase_if([](auto& p) { return std::abs(p.second) < 0.001; });
    }
}

void StdSparseSemiunitaryMatrix::normalize() {
    for (auto& v : basis_) {
        double sum_of_squares = 0;
        for (const auto& p : v) {
            sum_of_squares += p.second * p.second;
        }
        double sqrt_of_sum_of_squares = sqrt(sum_of_squares);
        for (auto& p : v) {
            p.second /= sqrt_of_sum_of_squares;
        }
    }
}

const StdSparseSemiunitaryMatrix* StdSparseSemiunitaryMatrix::downcast_ptr(
    const std::unique_ptr<AbstractSparseSemiunitaryMatrix>& ptr) {
    auto answer = dynamic_cast<const StdSparseSemiunitaryMatrix*>(ptr.get());
    if (answer == nullptr) {
        throw std::bad_cast();
    }
    return answer;
}

void StdSparseSemiunitaryMatrix::print(std::ostream& os) const {
    for (uint32_t i = 0; i < size_cols(); ++i) {
        for (const auto& p : basis_[i]) {
            os << p.second << "*[" << p.first << "] ";
        }
        os << std::endl;
    }
    os << std::endl;
}

StdSparseSemiunitaryMatrix*
StdSparseSemiunitaryMatrix::downcast_ptr(std::unique_ptr<AbstractSparseSemiunitaryMatrix>& ptr) {
    auto answer = dynamic_cast<StdSparseSemiunitaryMatrix*>(ptr.get());
    if (answer == nullptr) {
        throw std::bad_cast();
    }
    return answer;
}

const std::vector<StdSparseSemiunitaryMatrix::Map>&
StdSparseSemiunitaryMatrix::getSparseSemiunitaryMatrix() const {
    return basis_;
}

void StdSparseSemiunitaryMatrix::unitaryTransform(
    const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrixToTransform,
    std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrixToAdd) const {
    size_t matrix_in_space_basis_size = this->size_cols();

#pragma omp parallel for shared( \
        matrix_in_space_basis_size, \
            symmetricMatrixToTransform, \
            symmetricMatrixToAdd) default(none)
    for (uint32_t index_of_space_vector_i = 0; index_of_space_vector_i < matrix_in_space_basis_size;
         ++index_of_space_vector_i) {
        for (const auto& i_iter : basis_[index_of_space_vector_i]) {
            uint32_t index_of_lexicographic_vector_k = i_iter.first;
            for (uint32_t index_of_space_vector_j = index_of_space_vector_i;
                 index_of_space_vector_j < matrix_in_space_basis_size;
                 ++index_of_space_vector_j) {
                for (const auto& j_iter : basis_[index_of_space_vector_j]) {
                    uint32_t index_of_lexicographic_vector_l = j_iter.first;
                    double value_in_matrix_in_lexicografical_basis = symmetricMatrixToTransform->at(
                        index_of_lexicographic_vector_k,
                        index_of_lexicographic_vector_l);
                    if (value_in_matrix_in_lexicografical_basis != 0) {
                        // \Delta B_{ij} = \sum_{kl} U_{ki} * \Delta A_{kl} * U_{lj}
                        symmetricMatrixToAdd->add_to_position(
                            i_iter.second * value_in_matrix_in_lexicografical_basis * j_iter.second,
                            index_of_space_vector_i,
                            index_of_space_vector_j);
                    }
                }
            }
        }
    }
}
}  // namespace quantum::linear_algebra
