#include "EmhashSparseSemiunitaryMatrix.h"

#include <algorithm>
#include <cmath>

#include "EmhashLogic.h"

namespace quantum::linear_algebra {

struct IteratorStdImpl: public AbstractSparseSemiunitaryMatrix::Iterator {
    EmhashSparseSemiunitaryMatrix::Map::const_iterator iter;
    const EmhashSparseSemiunitaryMatrix::Map::const_iterator end;

    IteratorStdImpl(
        EmhashSparseSemiunitaryMatrix::Map::const_iterator iter1,
        EmhashSparseSemiunitaryMatrix::Map::const_iterator iter2) :
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
EmhashSparseSemiunitaryMatrix::GetNewIterator(size_t index_of_vector) const {
    return std::make_unique<IteratorStdImpl>(
        basis_[index_of_vector].cbegin(),
        basis_[index_of_vector].cend());
}

uint32_t EmhashSparseSemiunitaryMatrix::size_cols() const {
    return basis_.size();
}

uint32_t EmhashSparseSemiunitaryMatrix::size_rows() const {
    return rows_;
}

bool EmhashSparseSemiunitaryMatrix::empty() const {
    return basis_.empty();
}

bool EmhashSparseSemiunitaryMatrix::vempty(uint32_t index_of_vector) const {
    return basis_[index_of_vector].empty();
}

void EmhashSparseSemiunitaryMatrix::clear() {
    basis_.clear();
}

void EmhashSparseSemiunitaryMatrix::move_vector_from(
    uint32_t i,
    std::unique_ptr<AbstractSparseSemiunitaryMatrix>& subspace_from) {
    auto rhs_ = downcast_ptr(subspace_from);
    basis_.emplace_back(std::move(rhs_->basis_[i]));
}

void EmhashSparseSemiunitaryMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    if (!basis_[i].contains(j)) {
        basis_[i].emplace(j, value);
    } else {
        basis_[i][j] += value;
    }
}

double EmhashSparseSemiunitaryMatrix::at(uint32_t i, uint32_t j) const {
    if (!basis_.at(i).contains(j)) {
        return 0;
    }
    return basis_.at(i).at(j);
}

void EmhashSparseSemiunitaryMatrix::resize(uint32_t cols, uint32_t rows) {
    basis_.resize(cols);
    rows_ = rows;
}

bool EmhashSparseSemiunitaryMatrix::is_zero(uint32_t i, uint32_t j) const {
    return !basis_.at(i).contains(j);
}

void EmhashSparseSemiunitaryMatrix::eraseExplicitZeros() {
    for (auto& mm : basis_) {
        mm.erase_if([](auto& p) { return std::abs(p.second) < 0.001; });
    }
}

void EmhashSparseSemiunitaryMatrix::normalize() {
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

const EmhashSparseSemiunitaryMatrix* EmhashSparseSemiunitaryMatrix::downcast_ptr(
    const std::unique_ptr<AbstractSparseSemiunitaryMatrix>& ptr) {
    auto answer = dynamic_cast<const EmhashSparseSemiunitaryMatrix*>(ptr.get());
    if (answer == nullptr) {
        throw std::bad_cast();
    }
    return answer;
}

void EmhashSparseSemiunitaryMatrix::print(std::ostream& os) const {
    for (uint32_t i = 0; i < size_cols(); ++i) {
        for (const auto& p : basis_[i]) {
            os << p.second << "*[" << p.first << "] ";
        }
        os << std::endl;
    }
    os << std::endl;
}

EmhashSparseSemiunitaryMatrix*
EmhashSparseSemiunitaryMatrix::downcast_ptr(std::unique_ptr<AbstractSparseSemiunitaryMatrix>& ptr) {
    auto answer = dynamic_cast<EmhashSparseSemiunitaryMatrix*>(ptr.get());
    if (answer == nullptr) {
        throw std::bad_cast();
    }
    return answer;
}

const std::vector<EmhashSparseSemiunitaryMatrix::Map>&
EmhashSparseSemiunitaryMatrix::getSparseSemiunitaryMatrix() const {
    return basis_;
}

void EmhashSparseSemiunitaryMatrix::unitaryTransform(
    const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrixToTransform,
    std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrixToAdd) const {
    EmhashLogic logic;

    logic.unitaryTransform(symmetricMatrixToTransform, symmetricMatrixToAdd, *this);
}
}  // namespace quantum::linear_algebra
