#include "Matrix.h"

Matrix::Matrix(std::vector<Submatrix>&& m) : blocks(std::move(m)) {}

Matrix::Matrix(
    const space::Space& space,
    const model::operators::Operator& new_operator,
    std::shared_ptr<const index_converter::AbstractIndexConverter> converter,
    const quantum::linear_algebra::FactoriesList& factories,
    bool return_sparse_if_possible) {
    blocks.reserve(space.getBlocks().size());

    for (const auto& subspace : space.getBlocks()) {
        blocks.emplace_back(subspace, new_operator, converter, factories, return_sparse_if_possible);
    }
}

MatrixRef::MatrixRef(const Matrix& matrix) {
    for (int i = 0; i < matrix.blocks.size(); ++i) {
        blocks.push_back(matrix.blocks[i]);
    }
}

MatrixRef::MatrixRef(std::vector<std::reference_wrapper<const Submatrix>>&& blocks_) {
    blocks = std::move(blocks_);
}
