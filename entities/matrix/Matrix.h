#ifndef JULY_MATRIX_H
#define JULY_MATRIX_H

#include <iostream>
#include <vector>

#include "Submatrix.h"
#include "entities/space/Space.h"

struct Matrix {
    Matrix() = default;
    Matrix(
        const Space& space,
        const Operator& new_operator,
        const lexicographic::IndexConverter& converter);

    std::vector<Submatrix> blocks;
    explicit Matrix(std::vector<Submatrix>&& m);
};

std::ostream& operator<<(std::ostream& os, const Matrix& matrix);

#endif  //JULY_MATRIX_H
