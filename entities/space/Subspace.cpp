#include "Subspace.h"


struct Subspace::Impl {
    std::vector<std::map<uint32_t, double>> basis;
public:
    Impl() = default;
};

Subspace::Subspace() : pImpl{std::make_unique<Impl>()} {}
Subspace::~Subspace() = default;
Subspace::Subspace(Subspace&&) noexcept = default;
Subspace& Subspace::operator=(Subspace&&) noexcept = default;


struct SubspaceIteratorImpl : public Subspace::SubspaceIterator {

    std::map<uint32_t, double>::iterator iter;
    const std::map<uint32_t, double>::iterator end;

    SubspaceIteratorImpl(std::map<uint32_t, double>::iterator iter1, 
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

    ~SubspaceIteratorImpl() = default;
};


std::unique_ptr<Subspace::SubspaceIterator> Subspace::GetNewIterator(size_t index_of_vector) const {
    return std::make_unique<SubspaceIteratorImpl>(pImpl->basis[index_of_vector].begin(), pImpl->basis[index_of_vector].end());
}


std::ostream &operator<<(std::ostream &os, const Subspace &subspace) {
    os << subspace.properties;
    for (uint32_t i = 0; i < subspace.size(); ++i) {
        for (const auto& p : subspace.pImpl->basis[i]) {
            os << p.second << "*[" << p.first << "] ";
        }
        os << std::endl;
    }
    os << std::endl;
    return os;
}

uint32_t Subspace::size() const {
    return pImpl->basis.size();
}

bool Subspace::empty() const {
    return pImpl->basis.empty();
}

bool Subspace::vempty(uint32_t index_of_vector) const {
    return pImpl->basis[index_of_vector].empty();
}

void Subspace::clear() {
    pImpl->basis.clear();
}

void Subspace::move_vector_from(uint32_t i, Subspace& subspace_from) {
    pImpl->basis.emplace_back(std::move(subspace_from.pImpl->basis[i]));
}

void Subspace::move_all_from(Subspace& subspace_from) {
    for (uint32_t i = 0; i < subspace_from.size(); ++i) {
        move_vector_from(i, subspace_from);
    }
}

void Subspace::copy_vector_from(uint32_t i, const Subspace& subspace_from) {
    pImpl->basis.emplace_back(subspace_from.pImpl->basis[i]);
}

void Subspace::copy_all_from(const Subspace& subspace_from) {
    for (uint32_t i = 0; i < subspace_from.size(); ++i) {
        copy_vector_from(i, subspace_from);
    }
}

void Subspace::add_to_position(double value, uint32_t i, uint32_t j) {
    pImpl->basis[i][j] += value;
}

double Subspace::operator()(uint32_t i, uint32_t j) const {
    return pImpl->basis[i].at(j);
}

void Subspace::resize(uint32_t new_size) {
    pImpl->basis.resize(new_size);
}

bool Subspace::is_zero(uint32_t i, uint32_t j) const {
    return pImpl->basis[i].find(j) == pImpl->basis[i].end();
}

void Subspace::erase_if_zero() {
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

