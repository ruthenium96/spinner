#include "UnitarySpaseMatrix.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

struct UnitarySpaseMatrix::Impl {
    std::vector<std::map<uint32_t, double>> basis;
public:
    Impl() = default;
};

UnitarySpaseMatrix::UnitarySpaseMatrix() : pImpl{std::make_unique<Impl>()} {}
UnitarySpaseMatrix::~UnitarySpaseMatrix() = default;
UnitarySpaseMatrix::UnitarySpaseMatrix(UnitarySpaseMatrix&&) noexcept = default;
UnitarySpaseMatrix& UnitarySpaseMatrix::operator=(UnitarySpaseMatrix&&) noexcept = default;


struct IteratorImpl : public UnitarySpaseMatrix::Iterator {

    std::map<uint32_t, double>::iterator iter;
    const std::map<uint32_t, double>::iterator end;

    IteratorImpl(std::map<uint32_t, double>::iterator iter1,
                 std::map<uint32_t, double>::iterator iter2): iter(iter1), end(iter2){
    }

    [[nodiscard]] bool hasNext() const override {
        return iter != end;
    }

    IndexValueItem getNext() override {
        auto pair = *iter;
        ++iter;
        return {pair.first, pair.second};
    }

    ~IteratorImpl() = default;
};


std::unique_ptr<UnitarySpaseMatrix::Iterator> UnitarySpaseMatrix::GetNewIterator(size_t index_of_vector) const {
    return std::make_unique<IteratorImpl>(pImpl->basis[index_of_vector].begin(), pImpl->basis[index_of_vector].end());
}


std::ostream &operator<<(std::ostream &os, const UnitarySpaseMatrix &decomposition) {
    for (uint32_t i = 0; i < decomposition.size(); ++i) {
        for (const auto& p : decomposition.pImpl->basis[i]) {
            os << p.second << "*[" << p.first << "] ";
        }
        os << std::endl;
    }
    os << std::endl;
    return os;
}

uint32_t UnitarySpaseMatrix::size() const {
    return pImpl->basis.size();
}

bool UnitarySpaseMatrix::empty() const {
    return pImpl->basis.empty();
}

bool UnitarySpaseMatrix::vempty(uint32_t index_of_vector) const {
    return pImpl->basis[index_of_vector].empty();
}

void UnitarySpaseMatrix::clear() {
    pImpl->basis.clear();
}

void UnitarySpaseMatrix::move_vector_from(uint32_t i, UnitarySpaseMatrix& subspace_from) {
    pImpl->basis.emplace_back(std::move(subspace_from.pImpl->basis[i]));
}

void UnitarySpaseMatrix::move_all_from(UnitarySpaseMatrix& subspace_from) {
    for (uint32_t i = 0; i < subspace_from.size(); ++i) {
        move_vector_from(i, subspace_from);
    }
}

void UnitarySpaseMatrix::copy_vector_from(uint32_t i, const UnitarySpaseMatrix& subspace_from) {
    pImpl->basis.emplace_back(subspace_from.pImpl->basis[i]);
}

void UnitarySpaseMatrix::copy_all_from(const UnitarySpaseMatrix& subspace_from) {
    for (uint32_t i = 0; i < subspace_from.size(); ++i) {
        copy_vector_from(i, subspace_from);
    }
}

void UnitarySpaseMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    pImpl->basis[i][j] += value;
}

double UnitarySpaseMatrix::operator()(uint32_t i, uint32_t j) const {
    return pImpl->basis[i].at(j);
}

void UnitarySpaseMatrix::resize(uint32_t new_size) {
    pImpl->basis.resize(new_size);
}

bool UnitarySpaseMatrix::is_zero(uint32_t i, uint32_t j) const {
    return pImpl->basis[i].find(j) == pImpl->basis[i].end();
}

void UnitarySpaseMatrix::erase_if_zero() {
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

void UnitarySpaseMatrix::normalize() {
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

bool UnitarySpaseMatrix::is_equal_up_to_vector_order(const UnitarySpaseMatrix& rhs) const {
    std::vector<std::map<uint32_t, double>> sorted_lhs_basis = pImpl->basis;
    std::sort(sorted_lhs_basis.begin(), sorted_lhs_basis.end());

    std::vector<std::map<uint32_t, double>> sorted_rhs_basis = rhs.pImpl->basis;
    std::sort(sorted_rhs_basis.begin(), sorted_rhs_basis.end());

    return (sorted_lhs_basis == sorted_rhs_basis);
}