#ifndef SPINNER_SUBMATRIX_H
#define SPINNER_SUBMATRIX_H

#include <src/entities/data_structures/FactoriesList.h>

#include "src/entities/BlockProperties.h"
#include "src/model/operators/Operator.h"
#include "src/space/Subspace.h"

struct Submatrix {
    BlockProperties properties;
    std::unique_ptr<quantum::linear_algebra::AbstractSymmetricMatrix> raw_data;

    Submatrix() = delete;

    Submatrix(
        const space::Subspace& subspace,
        const model::operators::Operator& new_operator,
        const lexicographic::IndexConverter& converter,
        const quantum::linear_algebra::FactoriesList& factories);
    friend std::ostream& operator<<(std::ostream& os, const Submatrix& submatrix);
};

#endif  //SPINNER_SUBMATRIX_H
