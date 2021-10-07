#include "Subspace.h"

std::ostream &operator<<(std::ostream &os, const Subspace &subspace) {
    os << subspace.properties;
    for (uint32_t i = 0; i < subspace.size(); ++i) {
        for (auto p = subspace.vbegin(i); p != subspace.vend(i); ++p) {
            os << p->second << "*[" << p->first << "] ";
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

bool Subspace::vempty(uint32_t index_of_vector) const {
    return basis[index_of_vector].empty();
}

void Subspace::clear() {
    basis.clear();
}

Subspace::in_col_iterator Subspace::vbegin(uint32_t index_of_vector) {
    return basis[index_of_vector].begin();
}

Subspace::in_col_iterator Subspace::vend(uint32_t index_of_vector) {
    return basis[index_of_vector].end();
}

Subspace::const_in_col_iterator Subspace::vbegin(uint32_t index_of_vector) const {
    return basis[index_of_vector].cbegin();
}

Subspace::const_in_col_iterator Subspace::vend(uint32_t index_of_vector) const {
    return basis[index_of_vector].cend();
}

void Subspace::move_vector_from(uint32_t i, Subspace& subspace_from) {
    basis.emplace_back(std::move(subspace_from.basis[i]));
}

void Subspace::move_all_from(Subspace& subspace_from) {
    for (uint32_t i = 0; i < subspace_from.size(); ++i) {
        move_vector_from(i, subspace_from);
    }
}

void Subspace::copy_vector_from(uint32_t i, const Subspace& subspace_from) {
    basis.emplace_back(subspace_from.basis[i]);
}

void Subspace::copy_all_from(const Subspace& subspace_from) {
    for (uint32_t i = 0; i < subspace_from.size(); ++i) {
        copy_vector_from(i, subspace_from);
    }
}

double &Subspace::operator()(uint32_t i, uint32_t j) {
    return basis[i][j];
}

void Subspace::resize(uint32_t new_size) {
    basis.resize(new_size);
}

bool Subspace::is_zero(uint32_t i, uint32_t j) {
    return basis[i].find(j) == basis[i].end();
}

void Subspace::erase_if_zero() {
    for (std::map<uint32_t, double>& mm : basis) {
        for (auto i = mm.begin(), last = mm.end(); i != last;) {
            if (std::abs(i->second) < 0.001) {
                i = mm.erase(i);
            } else {
                ++i;
            }
        }
    }
}
