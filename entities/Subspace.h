#ifndef JULY_SUBSPACE_H
#define JULY_SUBSPACE_H

#include <armadillo>
#include <cstdint>
#include <map>
#include <vector>
#include <ostream>
#include "BlockProperties.h"


struct Subspace {
    BlockProperties properties;

    [[nodiscard]] uint32_t size() const;
    [[nodiscard]] bool empty() const;
    [[nodiscard]] bool vempty(uint32_t index_of_vector) const;
    void clear();

    void erase_if_zero();

    bool is_zero(uint32_t i, uint32_t j);
    void move_vector_from(uint32_t i, Subspace& subspace_from);
    void move_all_from(Subspace& subspace_from);
    void copy_vector_from(uint32_t i, const Subspace& subspace_from);
    void copy_all_from(const Subspace& subspace_from);
    void resize(uint32_t new_size);

    double& operator()(uint32_t i, uint32_t j);

    using in_col_iterator = typename std::map<uint32_t, double>::iterator;
    in_col_iterator vbegin(uint32_t index_of_vector);
    in_col_iterator vend(uint32_t index_of_vector);
    using const_in_col_iterator = typename std::map<uint32_t, double>::const_iterator;
    [[nodiscard]] const_in_col_iterator vbegin(uint32_t index_of_vector) const;
    [[nodiscard]] const_in_col_iterator vend(uint32_t index_of_vector) const;

private:
    std::vector<std::map<uint32_t, double>> basis;
    arma::sp_mat sparse_basis;
};

std::ostream &operator<<(std::ostream &os, const Subspace &subspace);


#endif // JULY_SUBSPACE_H
