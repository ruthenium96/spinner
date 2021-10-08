#ifndef JULY_SUBSPACE_H
#define JULY_SUBSPACE_H

#include <armadillo>
#include <cstdint>
#include <map>
#include <vector>
#include <ostream>
#include "BlockProperties.h"

// std::vector<std::map>:
#define INDEX(p) p->first
#define VALUE(p) p->second

// arma::sp_mat:
//#define INDEX(p) p.row()
//#define VALUE(p) *p

// the slowest arma::sp_mat:
//#define INDEX(p) p.col()
//#define VALUE(p) *p

// TODO: create class Decomposition, move all methods to it

struct Subspace {
    BlockProperties properties;
    uint32_t tensor_size = 0;

    [[nodiscard]] uint32_t size() const;
    [[nodiscard]] bool empty() const;
    [[nodiscard]] bool vempty(uint32_t index_of_vector) const;
    void clear();

    void erase_if_zero();

    [[nodiscard]] bool is_zero(uint32_t i, uint32_t j) const;
    void move_vector_from(uint32_t i, Subspace& subspace_from);
    void move_all_from(Subspace& subspace_from);
    void copy_vector_from(uint32_t i, const Subspace& subspace_from);
    void copy_all_from(const Subspace& subspace_from);
    void resize(uint32_t new_size);

    void add_to_position(double value, uint32_t i, uint32_t j);
    double operator()(uint32_t i, uint32_t j) const;

    using in_col_iterator = typename std::map<uint32_t, double>::iterator;
//    using in_col_iterator = typename arma::SpMat<double>::col_iterator;
//    using in_col_iterator = typename arma::SpMat<double>::row_iterator;
    in_col_iterator vbegin(uint32_t index_of_vector);
    in_col_iterator vend(uint32_t index_of_vector);
    using const_in_col_iterator = typename std::map<uint32_t, double>::const_iterator;
//    using const_in_col_iterator = typename arma::SpMat<double>::const_col_iterator;
//    using const_in_col_iterator = typename arma::SpMat<double>::const_row_iterator;
    [[nodiscard]] const_in_col_iterator vbegin(uint32_t index_of_vector) const;
    [[nodiscard]] const_in_col_iterator vend(uint32_t index_of_vector) const;

    friend std::ostream &operator<<(std::ostream &os, const Subspace &subspace);

private:
    std::vector<std::map<uint32_t, double>> basis;
//    arma::sp_mat sparse_basis;
};

#endif // JULY_SUBSPACE_H
