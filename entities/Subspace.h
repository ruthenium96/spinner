#ifndef JULY_SUBSPACE_H
#define JULY_SUBSPACE_H

#include <armadillo>
#include <cstdint>
#include <map>
#include <vector>
#include <ostream>
#include "BlockProperties.h"

using DecompositionSparse = arma::SpSubview_col<double>;
using DecompositionMap = std::map<uint32_t, double>;

struct Subspace {
    using Decomposition = DecompositionMap;
    BlockProperties properties;

    [[nodiscard]] uint32_t size() const;
    [[nodiscard]] bool empty() const;
    void add_new_vector(Decomposition);
    void clear();
    bool is_zero(uint32_t i, uint32_t j);
    void move_vector_from(uint32_t i, Subspace& subspace_from);
    void resize(uint32_t new_size);

    Decomposition& back();
    double& operator()(uint32_t i, uint32_t j);

    using iterator = typename std::vector<Decomposition>::iterator;
    iterator begin();
    iterator end();
    using const_iterator = typename std::vector<DecompositionMap>::const_iterator;
    [[nodiscard]] const_iterator begin() const;
    [[nodiscard]] const_iterator end() const;

private:
    std::vector<DecompositionMap> basis;
    arma::sp_mat sparse_basis;
};

std::ostream &operator<<(std::ostream &os, const Subspace &subspace);


#endif // JULY_SUBSPACE_H
