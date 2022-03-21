#ifndef JULY_SUBMATRIX_H
#define JULY_SUBMATRIX_H

#include "src/entities/BlockProperties.h"
#include "src/entities/data_structures/DenseMatrix.h"
#include "src/entities/operator/Operator.h"
#include "src/entities/space/Subspace.h"

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
