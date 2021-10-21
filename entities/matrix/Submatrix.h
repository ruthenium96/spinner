#ifndef JULY_SUBMATRIX_H
#define JULY_SUBMATRIX_H

#include "entities/BlockProperties.h"
#include "RawSubmatrixData.h"

#include <armadillo>

struct Submatrix {
    BlockProperties properties;
    RawSubmatrixData raw_data;

    friend std::ostream &operator<<(std::ostream &os, const Submatrix &submatrix);
};


#endif //JULY_SUBMATRIX_H
