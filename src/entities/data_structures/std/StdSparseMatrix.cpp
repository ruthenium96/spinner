#include "StdSparseMatrix.h"

#include <algorithm>
#include <cmath>

namespace quantum::linear_algebra {

struct IteratorImpl: public AbstractSparseMatrix::Iterator {
    std::map<uint32_t, double>::const_iterator iter;
    const std::map<uint32_t, double>::const_iterator end;

    IteratorImpl(
        std::map<uint32_t, double>::const_iterator iter1,
        std::map<uint32_t, double>::const_iterator iter2) :
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

    ~IteratorImpl() override = default;
};

std::unique_ptr<AbstractSparseMatrix::Iterator>
StdSparseMatrix::GetNewIterator(size_t index_of_vector) const {
    return std::make_unique<IteratorImpl>(
        basis_[index_of_vector].cbegin(),
        basis_[index_of_vector].cend());
}

std::ostream& operator<<(std::ostream& os, const StdSparseMatrix& decomposition) {}

uint32_t StdSparseMatrix::size() const {
    return basis_.size();
}

bool StdSparseMatrix::empty() const {
    return basis_.empty();
}

bool StdSparseMatrix::vempty(uint32_t index_of_vector) const {
    return basis_[index_of_vector].empty();
}

void StdSparseMatrix::clear() {
    basis_.clear();
}

void StdSparseMatrix::move_vector_from(
    uint32_t i,
    std::unique_ptr<AbstractSparseMatrix>& subspace_from) {
    auto rhs_ = downcast_ptr(subspace_from);
    basis_.emplace_back(std::move(rhs_->basis_[i]));
}

void StdSparseMatrix::move_all_from(std::unique_ptr<AbstractSparseMatrix>& subspace_from) {
    for (uint32_t i = 0; i < subspace_from->size(); ++i) {
        move_vector_from(i, subspace_from);
    }
}

void StdSparseMatrix::copy_vector_from(
    uint32_t i,
    const std::unique_ptr<AbstractSparseMatrix>& subspace_from) {
    auto rhs_ = downcast_ptr(subspace_from);
    basis_.emplace_back(rhs_->basis_[i]);
}

void StdSparseMatrix::copy_all_from(const std::unique_ptr<AbstractSparseMatrix>& subspace_from) {
    for (uint32_t i = 0; i < subspace_from->size(); ++i) {
        copy_vector_from(i, subspace_from);
    }
}

void StdSparseMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    basis_[i][j] += value;
}

double StdSparseMatrix::at(uint32_t i, uint32_t j) const {
    auto mb_iterator = basis_[i].find(j);
    if (mb_iterator == basis_[i].end()) {
        return 0;
    } else {
        return mb_iterator->second;
    }
}

void StdSparseMatrix::resize(uint32_t new_size) {
    basis_.resize(new_size);
}

bool StdSparseMatrix::is_zero(uint32_t i, uint32_t j) const {
    return basis_[i].find(j) == basis_[i].end();
}

void StdSparseMatrix::erase_if_zero() {
    for (std::map<uint32_t, double>& mm : basis_) {
        for (auto i = mm.begin(), last = mm.end(); i != last;) {
            if (std::abs(i->second) < 0.001) {
                i = mm.erase(i);
            } else {
                ++i;
            }
        }
    }
}

void StdSparseMatrix::normalize() {
    for (auto& v : basis_) {
        double sum_of_squares = 0;
        for (const auto p : v) {
            sum_of_squares += p.second * p.second;
        }
        double sqrt_of_sum_of_squares = sqrt(sum_of_squares);
        for (auto& p : v) {
            p.second /= sqrt_of_sum_of_squares;
        }
    }
}

bool StdSparseMatrix::is_equal_up_to_vector_order(
    const std::unique_ptr<AbstractSparseMatrix>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);
    std::vector<std::map<uint32_t, double>> sorted_lhs_basis = basis_;
    std::sort(sorted_lhs_basis.begin(), sorted_lhs_basis.end());

    std::vector<std::map<uint32_t, double>> sorted_rhs_basis = rhs_->basis_;
    std::sort(sorted_rhs_basis.begin(), sorted_rhs_basis.end());

    return (sorted_lhs_basis == sorted_rhs_basis);
}

const StdSparseMatrix*
StdSparseMatrix::downcast_ptr(const std::unique_ptr<AbstractSparseMatrix>& ptr) {
    auto answer = dynamic_cast<const StdSparseMatrix*>(ptr.get());
    if (answer == nullptr) {
        throw std::bad_cast();
    }
    return answer;
}

void StdSparseMatrix::print(std::ostream& os) const {
    for (uint32_t i = 0; i < size(); ++i) {
        for (const auto& p : basis_[i]) {
            os << p.second << "*[" << p.first << "] ";
        }
        os << std::endl;
    }
    os << std::endl;
}

StdSparseMatrix* StdSparseMatrix::downcast_ptr(std::unique_ptr<AbstractSparseMatrix>& ptr) {
    auto answer = dynamic_cast<StdSparseMatrix*>(ptr.get());
    if (answer == nullptr) {
        throw std::bad_cast();
    }
    return answer;
}
}  // namespace quantum::linear_algebra
