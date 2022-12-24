#include "StdSparseSemiunitaryMatrix.h"

#include <algorithm>
#include <cmath>

namespace quantum::linear_algebra {

struct IteratorStdImpl: public AbstractSparseSemiunitaryMatrix::Iterator {
    std::map<uint32_t, double>::const_iterator iter;
    const std::map<uint32_t, double>::const_iterator end;

    IteratorStdImpl(
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

    ~IteratorStdImpl() override = default;
};

std::unique_ptr<AbstractSparseSemiunitaryMatrix::Iterator>
StdSparseSemiunitaryMatrix::GetNewIterator(size_t index_of_vector) const {
    return std::make_unique<IteratorStdImpl>(
        basis_[index_of_vector].cbegin(),
        basis_[index_of_vector].cend());
}

std::ostream& operator<<(std::ostream& os, const StdSparseSemiunitaryMatrix& decomposition) {}

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
    basis_[i][j] += value;
}

double StdSparseSemiunitaryMatrix::at(uint32_t i, uint32_t j) const {
    auto mb_iterator = basis_[i].find(j);
    if (mb_iterator == basis_[i].end()) {
        return 0;
    } else {
        return mb_iterator->second;
    }
}

void StdSparseSemiunitaryMatrix::resize(uint32_t cols, uint32_t rows) {
    basis_.resize(cols);
    rows_ = rows;
}

bool StdSparseSemiunitaryMatrix::is_zero(uint32_t i, uint32_t j) const {
    return basis_[i].find(j) == basis_[i].end();
}

void StdSparseSemiunitaryMatrix::eraseExplicitZeros() {
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

void StdSparseSemiunitaryMatrix::normalize() {
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
}  // namespace quantum::linear_algebra
