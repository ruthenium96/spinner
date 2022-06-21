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
    const lexicographic::IndexConverter& converter,
    const std::unique_ptr<quantum::linear_algebra::AbstractFactory>& factory) {
    blocks.reserve(space.getBlocks().size());
    //    blocks.resize(space.getBlocks().size());

    for (const auto& subspace : space.getBlocks()) {
        blocks.emplace_back(Submatrix(subspace, new_operator, converter, factory));
    }
}
