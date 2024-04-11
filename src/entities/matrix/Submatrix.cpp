#include "Submatrix.h"

#include <set>

Submatrix::Submatrix(
    std::unique_ptr<quantum::linear_algebra::AbstractDiagonalizableMatrix> raw_data_,
    BlockProperties properties_) :
    raw_data(std::move(raw_data_)),
    properties(properties_) {}

Submatrix::Submatrix(
    const space::Subspace& subspace,
    const model::operators::Operator& new_operator,
    const lexicographic::IndexConverter& converter,
    const quantum::linear_algebra::FactoriesList& factories) {
    auto totalSpaceSize = converter.get_total_space_size();
    auto matrix_in_lexicografical_basis = factories.createSparseSymmetricMatrix(totalSpaceSize);

    std::set<unsigned int> lexicografical_vectors_to_built;
    size_t matrix_in_space_basis_size = subspace.decomposition->size_cols();

    for (uint32_t index_of_space_vector_i = 0; index_of_space_vector_i < matrix_in_space_basis_size;
         ++index_of_space_vector_i) {
        auto outer_iterator = subspace.decomposition->GetNewIterator(index_of_space_vector_i);
        while (outer_iterator->hasNext()) {
            auto outer_item = outer_iterator->getNext();
            lexicografical_vectors_to_built.insert(outer_item.index);
        }
    }

    for (auto index_of_lexicographic_vector_k : lexicografical_vectors_to_built) {
        // BUILDING k-th ROW OF INITIAL_MATRIX
        for (auto& term : new_operator.getTerms()) {
            term->construct(*matrix_in_lexicografical_basis, index_of_lexicographic_vector_k);
        }
    }

    properties = subspace.properties;
    if (subspace.dense_semiunitary_matrix.has_value()) {
        raw_data = factories.createSparseDiagonalizableMatrix(matrix_in_space_basis_size);
    } else {
        raw_data = factories.createDenseDiagonalizableMatrix(matrix_in_space_basis_size);
    }

    subspace.decomposition->unitaryTransform(matrix_in_lexicografical_basis, raw_data);
    if (subspace.dense_semiunitary_matrix.has_value()) {
        raw_data = subspace.dense_semiunitary_matrix.value()->unitaryTransform(raw_data);
    }
}
