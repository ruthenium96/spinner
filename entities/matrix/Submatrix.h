#ifndef JULY_SUBMATRIX_H
#define JULY_SUBMATRIX_H

#include "entities/BlockProperties.h"
#include "entities/data_structures/DenseMatrix.h"
#include "entities/operator/Operator.h"
#include "entities/space/Subspace.h"

struct Submatrix {
    BlockProperties properties;
    DenseMatrix raw_data;

    Submatrix() = default;

    Submatrix(
        const Subspace& subspace,
        const Operator& new_operator,
        const lexicographic::IndexConverter& converter);
    friend std::ostream& operator<<(std::ostream& os, const Submatrix& submatrix);
};

#endif  //JULY_SUBMATRIX_H
