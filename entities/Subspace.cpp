#include "Subspace.h"

#if defined(VECTOR_OF_MAPS)

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

in_col_iterator Subspace::vbegin(uint32_t index_of_vector) {
    return basis[index_of_vector].begin();
}

in_col_iterator Subspace::vend(uint32_t index_of_vector) {
    return basis[index_of_vector].end();
}

const_in_col_iterator Subspace::vbegin(uint32_t index_of_vector) const {
    return basis[index_of_vector].cbegin();
}

const_in_col_iterator Subspace::vend(uint32_t index_of_vector) const {
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

void Subspace::add_to_position(double value, uint32_t i, uint32_t j) {
    basis[i][j] += value;
}

double Subspace::operator()(uint32_t i, uint32_t j) const {
    return basis[i].at(j);
}

void Subspace::resize(uint32_t new_size) {
    basis.resize(new_size);
}

bool Subspace::is_zero(uint32_t i, uint32_t j) const {
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

#elif defined(SPARSE_MATRIX)

std::ostream &operator<<(std::ostream &os, const Subspace &subspace) {
    os << subspace.properties;
    os << subspace.sparse_basis << std::endl;
    return os;
}

uint32_t Subspace::size() const {
    return sparse_basis.n_cols;
}

bool Subspace::empty() const {
    return sparse_basis.n_cols == 0;
}

bool Subspace::vempty(uint32_t index_of_vector) const {
    return sparse_basis.col(index_of_vector).is_zero();
}

void Subspace::clear() {
    sparse_basis.reset();
}

in_col_iterator Subspace::vbegin(uint32_t index_of_vector) {
    return sparse_basis.begin_col(index_of_vector);
}

in_col_iterator Subspace::vend(uint32_t index_of_vector) {
    return sparse_basis.end_col(index_of_vector);
}

const_in_col_iterator Subspace::vbegin(uint32_t index_of_vector) const {
    return sparse_basis.begin_col(index_of_vector);
}

const_in_col_iterator Subspace::vend(uint32_t index_of_vector) const {
    return sparse_basis.end_col(index_of_vector);
}

void Subspace::move_vector_from(uint32_t i, Subspace& subspace_from) {
    tensor_size = subspace_from.tensor_size;
    sparse_basis.resize(tensor_size, sparse_basis.n_cols + 1);
    sparse_basis.col(sparse_basis.n_cols - 1) = subspace_from.sparse_basis.col(i);
    subspace_from.sparse_basis.col(i).zeros();
}

void Subspace::move_all_from(Subspace& subspace_from) {
    tensor_size = subspace_from.tensor_size;
    sparse_basis.resize(tensor_size, sparse_basis.n_cols + subspace_from.sparse_basis.n_cols);
    sparse_basis.cols(sparse_basis.n_cols - subspace_from.sparse_basis.n_cols,sparse_basis.n_cols - 1) = subspace_from.sparse_basis;
    subspace_from.clear();
}

void Subspace::copy_vector_from(uint32_t i, const Subspace& subspace_from) {
    tensor_size = subspace_from.tensor_size;
    sparse_basis.resize(tensor_size, sparse_basis.n_cols + 1);
    sparse_basis.col(sparse_basis.n_cols - 1) = subspace_from.sparse_basis.col(i);
}

void Subspace::copy_all_from(const Subspace& subspace_from) {
    tensor_size = subspace_from.tensor_size;
    sparse_basis.resize(tensor_size, sparse_basis.n_cols + subspace_from.sparse_basis.n_cols);
    sparse_basis.cols(sparse_basis.n_cols - subspace_from.sparse_basis.n_cols,sparse_basis.n_cols - 1) = subspace_from.sparse_basis;
}

void Subspace::add_to_position(double value, uint32_t i, uint32_t j) {
    sparse_basis(j, i) += value;
}

double Subspace::operator()(uint32_t i, uint32_t j) const {
    return sparse_basis(j, i);
}

void Subspace::resize(uint32_t new_size) {
    sparse_basis.resize(tensor_size, new_size);
}

bool Subspace::is_zero(uint32_t i, uint32_t j) const {
    // TODO: epsilon
    return std::abs(sparse_basis(j, i)) < 0.001;
}

void Subspace::erase_if_zero() {
    for (auto && el : sparse_basis) {
        // TODO: epsilon
        if (std::abs(el) < 0.001) {
            el = 0.0;
        }
    }
}

#elif defined(VECTOR_OF_SPARSE_VECTORS)

std::ostream &operator<<(std::ostream &os, const Subspace &subspace) {
    os << subspace.properties;
    for (const auto& v : subspace.sparse_vectors) {
        os << v << std::endl;
    }
    return os;
}

uint32_t Subspace::size() const {
    return sparse_vectors.size();
}

bool Subspace::empty() const {
    return sparse_vectors.empty();
}

bool Subspace::vempty(uint32_t index_of_vector) const {
    return sparse_vectors[index_of_vector].n_nonzero == 0;
}

void Subspace::clear() {
    sparse_vectors.clear();
}

in_col_iterator Subspace::vbegin(uint32_t index_of_vector) {
    return sparse_vectors[index_of_vector].begin_col(0);
}

in_col_iterator Subspace::vend(uint32_t index_of_vector) {
    return sparse_vectors[index_of_vector].end_col(0);
}

const_in_col_iterator Subspace::vbegin(uint32_t index_of_vector) const {
    return sparse_vectors[index_of_vector].begin_col(0);
}

const_in_col_iterator Subspace::vend(uint32_t index_of_vector) const {
    return sparse_vectors[index_of_vector].end_col(0);
}

void Subspace::move_vector_from(uint32_t i, Subspace& subspace_from) {
    sparse_vectors.emplace_back(std::move(subspace_from.sparse_vectors[i]));
}

void Subspace::move_all_from(Subspace& subspace_from) {
    for (uint32_t i = 0; i < subspace_from.size(); ++i) {
        move_vector_from(i, subspace_from);
    }
}

void Subspace::copy_vector_from(uint32_t i, const Subspace& subspace_from) {
    sparse_vectors.emplace_back(subspace_from.sparse_vectors[i]);
}

void Subspace::copy_all_from(const Subspace& subspace_from) {
    for (uint32_t i = 0; i < subspace_from.size(); ++i) {
        copy_vector_from(i, subspace_from);
    }
}

void Subspace::add_to_position(double value, uint32_t i, uint32_t j) {
    sparse_vectors[i][j] += value;
}

double Subspace::operator()(uint32_t i, uint32_t j) const {
    return sparse_vectors[i][j];
}

void Subspace::resize(uint32_t new_size) {
    sparse_vectors.resize(new_size, arma::SpCol<double>(tensor_size));
}

bool Subspace::is_zero(uint32_t i, uint32_t j) const {
    // TODO: epsilon
    return std::abs(sparse_vectors[i][j]) < 0.001;
}

void Subspace::erase_if_zero() {
    for (auto & sp_v : sparse_vectors) {
        // TODO: epsilon
        for (auto i = sp_v.begin(), last = sp_v.end(); i != last; ++i) {
            if (std::abs(VALUE(i)) < 0.001) {
                VALUE(i) = 0;
            }
        }
    }
}
#endif
