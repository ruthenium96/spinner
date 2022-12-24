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
    const quantum::linear_algebra::FactoriesList& factories) {
    blocks.reserve(space.getBlocks().size());

    for (const auto& subspace : space.getBlocks()) {
        blocks.emplace_back(Submatrix(subspace, new_operator, converter, factories));
    }
}
