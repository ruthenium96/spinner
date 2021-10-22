#include "SubspaceData.h"

#include <cmath>
#include <map>
#include <vector>

struct SubspaceData::Impl {
    std::vector<std::map<uint32_t, double>> basis;
public:
    Impl() = default;
};

SubspaceData::SubspaceData() : pImpl{std::make_unique<Impl>()} {}
SubspaceData::~SubspaceData() = default;
SubspaceData::SubspaceData(SubspaceData&&) noexcept = default;
SubspaceData& SubspaceData::operator=(SubspaceData&&) noexcept = default;


struct IteratorImpl : public SubspaceData::Iterator {

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


std::unique_ptr<SubspaceData::Iterator> SubspaceData::GetNewIterator(size_t index_of_vector) const {
    return std::make_unique<IteratorImpl>(pImpl->basis[index_of_vector].begin(), pImpl->basis[index_of_vector].end());
}


std::ostream &operator<<(std::ostream &os, const SubspaceData &decomposition) {
    for (uint32_t i = 0; i < decomposition.size(); ++i) {
        for (const auto& p : decomposition.pImpl->basis[i]) {
            os << p.second << "*[" << p.first << "] ";
        }
        os << std::endl;
    }
    os << std::endl;
    return os;
}

uint32_t SubspaceData::size() const {
    return pImpl->basis.size();
}

bool SubspaceData::empty() const {
    return pImpl->basis.empty();
}

bool SubspaceData::vempty(uint32_t index_of_vector) const {
    return pImpl->basis[index_of_vector].empty();
}

void SubspaceData::clear() {
    pImpl->basis.clear();
}

void SubspaceData::move_vector_from(uint32_t i, SubspaceData& subspace_from) {
    pImpl->basis.emplace_back(std::move(subspace_from.pImpl->basis[i]));
}

void SubspaceData::move_all_from(SubspaceData& subspace_from) {
    for (uint32_t i = 0; i < subspace_from.size(); ++i) {
        move_vector_from(i, subspace_from);
    }
}

void SubspaceData::copy_vector_from(uint32_t i, const SubspaceData& subspace_from) {
    pImpl->basis.emplace_back(subspace_from.pImpl->basis[i]);
}

void SubspaceData::copy_all_from(const SubspaceData& subspace_from) {
    for (uint32_t i = 0; i < subspace_from.size(); ++i) {
        copy_vector_from(i, subspace_from);
    }
}

void SubspaceData::add_to_position(double value, uint32_t i, uint32_t j) {
    pImpl->basis[i][j] += value;
}

double SubspaceData::operator()(uint32_t i, uint32_t j) const {
    return pImpl->basis[i].at(j);
}

void SubspaceData::resize(uint32_t new_size) {
    pImpl->basis.resize(new_size);
}

bool SubspaceData::is_zero(uint32_t i, uint32_t j) const {
    return pImpl->basis[i].find(j) == pImpl->basis[i].end();
}

void SubspaceData::erase_if_zero() {
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

void SubspaceData::normalize() {
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

