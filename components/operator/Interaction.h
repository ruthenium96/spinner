#ifndef JULY_INTERACTION_H
#define JULY_INTERACTION_H

#include <armadillo>
#include <cstdint>
#include "common/LexicographicIndexConverter.h"

// TODO: fix it
typedef arma::sp_mat LexicograficalMatrix;

class ZeroCenterTerm {
public:
    virtual void construct(LexicograficalMatrix& matrix_in_lexicografical_basis, const spaces::LexicographicIndexConverter& converter,
                           uint32_t index_of_vector) const = 0;
};

class OneCenterTerm {
public:
    virtual void construct(LexicograficalMatrix& matrix_in_lexicografical_basis, const spaces::LexicographicIndexConverter& converter,
                           uint32_t index_of_vector, uint32_t center_a) const = 0;
};

class TwoCenterTerm {
public:
    virtual void construct(LexicograficalMatrix& matrix_in_lexicografical_basis, const spaces::LexicographicIndexConverter& converter,
                           uint32_t index_of_vector, uint32_t center_a, uint32_t center_b) const = 0;
};

#endif //JULY_INTERACTION_H
