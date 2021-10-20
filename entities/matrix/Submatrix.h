#ifndef JULY_SUBMATRIX_H
#define JULY_SUBMATRIX_H

#include "entities/BlockProperties.h"

#include <armadillo>

struct Submatrix {
    BlockProperties properties;
    // TODO: create an abstraction
    arma::dmat submatrix;


    friend std::ostream &operator<<(std::ostream &os, const Submatrix &submatrix);
};


#endif //JULY_SUBMATRIX_H
