#include "Matrix.h"

Matrix::Matrix(std::vector<Submatrix>&& m) : blocks(std::move(m)) {}

Matrix::Matrix(
    const space::Space& space,
    const model::operators::Operator& new_operator,
    std::shared_ptr<const index_converter::AbstractIndexConverter> converter,
    const quantum::linear_algebra::FactoriesList& factories) {
    blocks.reserve(space.getBlocks().size());

    for (const auto& subspace : space.getBlocks()) {
        blocks.emplace_back(subspace, new_operator, converter, factories);
    }
}

MatrixRef::MatrixRef(const Matrix& spectrum) {
    for (int i = 0; i < spectrum.blocks.size(); ++i) {
        blocks.push_back(spectrum.blocks[i]);
    }
}