#ifndef JULY_SUBMATRIX_H
#define JULY_SUBMATRIX_H

#include "entities/BlockProperties.h"
#include "entities/data_structures/DenseMatrix.h"

struct Submatrix {
    BlockProperties properties;
    DenseMatrix raw_data;

    friend std::ostream& operator<<(std::ostream& os, const Submatrix& submatrix);
};

#endif  //JULY_SUBMATRIX_H
