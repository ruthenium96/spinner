#include "Matrix.h"

std::ostream &operator<<(std::ostream &os, const Matrix &matrix) {
    for (const Submatrix& submatrix : matrix.blocks) {
        os << submatrix;
    }
    os << "------" << std::endl;
    return os;
}

Matrix::Matrix(std::vector<Submatrix>&& m) : blocks(std::move(m)) {
}
