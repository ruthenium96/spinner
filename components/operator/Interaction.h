#ifndef JULY_INTERACTION_H
#define JULY_INTERACTION_H

#include <armadillo>
#include <cstdint>
#include "common/lexicographic/IndexConverter.h"
#include "common/lexicographic/SparseMatrix.h"

class ZeroCenterTerm {
public:
    virtual void construct(lexicographic::SparseMatrix& matrix_in_lexicografical_basis, uint32_t index_of_vector) const = 0;
    virtual ~ZeroCenterTerm(){};
};

class OneCenterTerm {
public:
    virtual void construct(lexicographic::SparseMatrix& matrix_in_lexicografical_basis, uint32_t index_of_vector, uint32_t center_a) const = 0;
    virtual ~OneCenterTerm(){};

};

class TwoCenterTerm {
public:
    virtual void construct(lexicographic::SparseMatrix& matrix_in_lexicografical_basis, uint32_t index_of_vector, uint32_t center_a, uint32_t center_b) const = 0;
    [[nodiscard]] virtual arma::dmat get_parameters() const = 0;
    virtual ~TwoCenterTerm(){};

};

#endif //JULY_INTERACTION_H
