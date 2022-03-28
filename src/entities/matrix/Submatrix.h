#ifndef JULY_SUBMATRIX_H
#define JULY_SUBMATRIX_H

#include "src/entities/BlockProperties.h"
#include "src/entities/data_structures/DenseMatrix.h"
#include "src/model/operators/Operator.h"
#include "src/space/Subspace.h"

struct Submatrix {
    BlockProperties properties;
    DenseMatrix raw_data;

    Submatrix() = default;

    Submatrix(
        const space::Subspace& subspace,
        const model::operators::Operator& new_operator,
        const lexicographic::IndexConverter& converter);
    friend std::ostream& operator<<(std::ostream& os, const Submatrix& submatrix);
};

#endif  //JULY_SUBMATRIX_H
