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
    const space::Space& space,
    const model::operators::Operator& new_operator,
    const lexicographic::IndexConverter& converter) {
    blocks.resize(space.getBlocks().size());

    for (size_t i = 0; i < space.getBlocks().size(); ++i) {
        const space::Subspace& subspace = space.getBlocks()[i];
        blocks[i] = Submatrix(subspace, new_operator, converter);
    }
}
