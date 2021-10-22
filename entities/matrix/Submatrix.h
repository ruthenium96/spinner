#ifndef JULY_SUBMATRIX_H
#define JULY_SUBMATRIX_H

#include "entities/BlockProperties.h"
#include "SubmatrixData.h"

#include <armadillo>

struct Submatrix {
    BlockProperties properties;
    SubmatrixData raw_data;

    friend std::ostream &operator<<(std::ostream &os, const Submatrix &submatrix);
};


#endif //JULY_SUBMATRIX_H
