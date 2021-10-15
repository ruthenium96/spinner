#include "NewBasisDecomposition.h"

#include <map>
#include <vector>

struct NewBasisDecomposition::Impl {
    std::vector<std::map<uint32_t, double>> basis;
public:
    Impl() = default;
};

NewBasisDecomposition::NewBasisDecomposition() : pImpl{std::make_unique<Impl>()} {}
NewBasisDecomposition::~NewBasisDecomposition() = default;
NewBasisDecomposition::NewBasisDecomposition(NewBasisDecomposition&&) noexcept = default;
NewBasisDecomposition& NewBasisDecomposition::operator=(NewBasisDecomposition&&) noexcept = default;


struct IteratorImpl : public NewBasisDecomposition::Iterator {

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


std::unique_ptr<NewBasisDecomposition::Iterator> NewBasisDecomposition::GetNewIterator(size_t index_of_vector) const {
    return std::make_unique<IteratorImpl>(pImpl->basis[index_of_vector].begin(), pImpl->basis[index_of_vector].end());
}


std::ostream &operator<<(std::ostream &os, const NewBasisDecomposition &decomposition) {
    for (uint32_t i = 0; i < decomposition.size(); ++i) {
        for (const auto& p : decomposition.pImpl->basis[i]) {
            os << p.second << "*[" << p.first << "] ";
        }
        os << std::endl;
    }
    os << std::endl;
    return os;
}

uint32_t NewBasisDecomposition::size() const {
    return pImpl->basis.size();
}

bool NewBasisDecomposition::empty() const {
    return pImpl->basis.empty();
}

bool NewBasisDecomposition::vempty(uint32_t index_of_vector) const {
    return pImpl->basis[index_of_vector].empty();
}

void NewBasisDecomposition::clear() {
    pImpl->basis.clear();
}

void NewBasisDecomposition::move_vector_from(uint32_t i, NewBasisDecomposition& subspace_from) {
    pImpl->basis.emplace_back(std::move(subspace_from.pImpl->basis[i]));
}

void NewBasisDecomposition::move_all_from(NewBasisDecomposition& subspace_from) {
    for (uint32_t i = 0; i < subspace_from.size(); ++i) {
        move_vector_from(i, subspace_from);
    }
}

void NewBasisDecomposition::copy_vector_from(uint32_t i, const NewBasisDecomposition& subspace_from) {
    pImpl->basis.emplace_back(subspace_from.pImpl->basis[i]);
}

void NewBasisDecomposition::copy_all_from(const NewBasisDecomposition& subspace_from) {
    for (uint32_t i = 0; i < subspace_from.size(); ++i) {
        copy_vector_from(i, subspace_from);
    }
}

void NewBasisDecomposition::add_to_position(double value, uint32_t i, uint32_t j) {
    pImpl->basis[i][j] += value;
}

double NewBasisDecomposition::operator()(uint32_t i, uint32_t j) const {
    return pImpl->basis[i].at(j);
}

void NewBasisDecomposition::resize(uint32_t new_size) {
    pImpl->basis.resize(new_size);
}

bool NewBasisDecomposition::is_zero(uint32_t i, uint32_t j) const {
    return pImpl->basis[i].find(j) == pImpl->basis[i].end();
}

void NewBasisDecomposition::erase_if_zero() {
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

