#ifndef SPINNER_MATRIX_H
#define SPINNER_MATRIX_H

#include <iostream>
#include <vector>

#include "Submatrix.h"
#include "src/entities/data_structures/FactoriesList.h"
#include "src/space/Space.h"

struct Matrix {
    Matrix() = default;
    Matrix(
        const space::Space& space,
        const model::operators::Operator& new_operator,
        std::shared_ptr<const index_converter::AbstractIndexConverter> converter,
        const quantum::linear_algebra::FactoriesList& factories);

    std::vector<Submatrix> blocks;
    explicit Matrix(std::vector<Submatrix>&& m);
};

#endif  //SPINNER_MATRIX_H
