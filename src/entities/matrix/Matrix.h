#ifndef SPINNER_MATRIX_H
#define SPINNER_MATRIX_H

#include <iostream>
#include <vector>

#include "Submatrix.h"
#include "src/space/Space.h"

struct Matrix {
    Matrix() = default;
    Matrix(
        const space::Space& space,
        const model::operators::Operator& new_operator,
        const lexicographic::IndexConverter& converter);

    std::vector<Submatrix> blocks;
    explicit Matrix(std::vector<Submatrix>&& m);
};

std::ostream& operator<<(std::ostream& os, const Matrix& matrix);

#endif  //SPINNER_MATRIX_H
