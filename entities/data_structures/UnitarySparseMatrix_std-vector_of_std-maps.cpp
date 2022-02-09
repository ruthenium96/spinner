#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

#include "UnitarySparseMatrix.h"

struct UnitarySparseMatrix::Impl {
    std::vector<std::map<uint32_t, double>> basis;

  public:
    Impl() = default;
};

UnitarySparseMatrix::UnitarySparseMatrix() : pImpl {std::make_unique<Impl>()} {}
UnitarySparseMatrix::~UnitarySparseMatrix() = default;
UnitarySparseMatrix::UnitarySparseMatrix(UnitarySparseMatrix&&) noexcept = default;
UnitarySparseMatrix& UnitarySparseMatrix::operator=(UnitarySparseMatrix&&) noexcept = default;

struct IteratorImpl: public UnitarySparseMatrix::Iterator {
    std::map<uint32_t, double>::iterator iter;
    const std::map<uint32_t, double>::iterator end;

    IteratorImpl(
        std::map<uint32_t, double>::iterator iter1,
        std::map<uint32_t, double>::iterator iter2) :
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

    ~IteratorImpl() = default;
};

std::unique_ptr<UnitarySparseMatrix::Iterator>
UnitarySparseMatrix::GetNewIterator(size_t index_of_vector) const {
    return std::make_unique<IteratorImpl>(
        pImpl->basis[index_of_vector].begin(),
        pImpl->basis[index_of_vector].end());
}

std::ostream& operator<<(std::ostream& os, const UnitarySparseMatrix& decomposition) {
    for (uint32_t i = 0; i < decomposition.size(); ++i) {
        for (const auto& p : decomposition.pImpl->basis[i]) {
            os << p.second << "*[" << p.first << "] ";
        }
        os << std::endl;
    }
    os << std::endl;
    return os;
}

uint32_t UnitarySparseMatrix::size() const {
    return pImpl->basis.size();
}

bool UnitarySparseMatrix::empty() const {
    return pImpl->basis.empty();
}

bool UnitarySparseMatrix::vempty(uint32_t index_of_vector) const {
    return pImpl->basis[index_of_vector].empty();
}

void UnitarySparseMatrix::clear() {
    pImpl->basis.clear();
}

void UnitarySparseMatrix::move_vector_from(uint32_t i, UnitarySparseMatrix& subspace_from) {
    pImpl->basis.emplace_back(std::move(subspace_from.pImpl->basis[i]));
}

void UnitarySparseMatrix::move_all_from(UnitarySparseMatrix& subspace_from) {
    for (uint32_t i = 0; i < subspace_from.size(); ++i) {
        move_vector_from(i, subspace_from);
    }
}

void UnitarySparseMatrix::copy_vector_from(uint32_t i, const UnitarySparseMatrix& subspace_from) {
    pImpl->basis.emplace_back(subspace_from.pImpl->basis[i]);
}

void UnitarySparseMatrix::copy_all_from(const UnitarySparseMatrix& subspace_from) {
    for (uint32_t i = 0; i < subspace_from.size(); ++i) {
        copy_vector_from(i, subspace_from);
    }
}

void UnitarySparseMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    pImpl->basis[i][j] += value;
}

double UnitarySparseMatrix::operator()(uint32_t i, uint32_t j) const {
    return pImpl->basis[i].at(j);
}

void UnitarySparseMatrix::resize(uint32_t new_size) {
    pImpl->basis.resize(new_size);
}

bool UnitarySparseMatrix::is_zero(uint32_t i, uint32_t j) const {
    return pImpl->basis[i].find(j) == pImpl->basis[i].end();
}

void UnitarySparseMatrix::erase_if_zero() {
    for (std::map<uint32_t, double>& mm : pImpl->basis) {
        for (auto i = mm.begin(), last = mm.end(); i != last;) {
            if (std::abs(i->second) < 0.001) {
                i = mm.erase(i);
            } else {
                ++i;
            }
        }
    }
}

void UnitarySparseMatrix::normalize() {
    for (auto& v : pImpl->basis) {
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

bool UnitarySparseMatrix::is_equal_up_to_vector_order(const UnitarySparseMatrix& rhs) const {
    std::vector<std::map<uint32_t, double>> sorted_lhs_basis = pImpl->basis;
    std::sort(sorted_lhs_basis.begin(), sorted_lhs_basis.end());

    std::vector<std::map<uint32_t, double>> sorted_rhs_basis = rhs.pImpl->basis;
    std::sort(sorted_rhs_basis.begin(), sorted_rhs_basis.end());

    return (sorted_lhs_basis == sorted_rhs_basis);
}