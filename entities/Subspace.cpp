#include "Subspace.h"

std::ostream &operator<<(std::ostream &os, const Subspace &subspace) {
    os << subspace.properties;
    for (auto& m: subspace) {
        for (auto& d: m) {
            os << d.second << "*[" << d.first << "] ";
        }
        os << std::endl;
    }
    os << std::endl;
    return os;
}

uint32_t Subspace::size() const {
    return basis.size();
}

bool Subspace::empty() const {
    return basis.empty();
}

void Subspace::clear() {
    basis.clear();
}

Subspace::iterator Subspace::begin() {
    return basis.begin();
}

Subspace::iterator Subspace::end() {
    return basis.end();
}

Subspace::const_iterator Subspace::begin() const {
    return basis.cbegin();
}

Subspace::const_iterator Subspace::end() const {
    return basis.cend();
}

void Subspace::add_new_vector(Decomposition vector) {
    basis.emplace_back(std::move(vector));
}

void Subspace::move_vector_from(uint32_t i, Subspace& subspace_from) {
    basis.emplace_back(std::move(subspace_from.basis[i]));
}

double &Subspace::operator()(uint32_t i, uint32_t j) {
    return basis[i][j];
}


void Subspace::resize(uint32_t new_size) {
    basis.resize(new_size);
}

Subspace::Decomposition &Subspace::back() {
    return basis.back();
}

bool Subspace::is_zero(uint32_t i, uint32_t j) {
    return basis[i].find(j) == basis[i].end();
}