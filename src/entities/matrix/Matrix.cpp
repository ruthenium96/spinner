#include "Matrix.h"

std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
    for (const Submatrix& submatrix : matrix.blocks) {
        os << submatrix;
    }
    os << "------" << std::endl;
    return os;
}

Matrix::Matrix(std::vector<Submatrix>&& m) : blocks(std::move(m)) {}

Matrix::Matrix(
    const Space& space,
    const Operator& new_operator,
    const lexicographic::IndexConverter& converter) {
    blocks.resize(space.blocks.size());

    for (size_t i = 0; i < space.blocks.size(); ++i) {
        const Subspace& subspace = space.blocks[i];
        blocks[i] = Submatrix(subspace, new_operator, converter);
    }
}
