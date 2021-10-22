#include "MatrixBuilder.h"

#include <unordered_set>

MatrixBuilder::MatrixBuilder(spaces::LexicographicIndexConverter converter) : converter_(std::move(converter)) {
}

Matrix MatrixBuilder::apply(const Space &space, const Operator &new_operator) {

    std::vector<Submatrix> vector_result;
    vector_result.resize(space.blocks.size());

    for (size_t i = 0; i < space.blocks.size(); ++i) {
        const Subspace& subspace = space.blocks[i];
        vector_result[i] = apply_to_subentity(subspace, new_operator);
    }
    return Matrix(std::move(vector_result));
}

Submatrix MatrixBuilder::apply_to_subentity(const Subspace &subspace, const Operator &new_operator) {
    arma::sp_dmat matrix_in_lexicografical_basis(converter_.total_space_size, converter_.total_space_size);
    std::unordered_set<unsigned int> built_lexicografical_vectors;

    uint32_t matrix_in_space_basis_size = subspace.decomposition.size();
    Submatrix submatrix;
    submatrix.properties = subspace.properties;
    submatrix.raw_data.resize(matrix_in_space_basis_size, matrix_in_space_basis_size);

    for (uint32_t index_of_space_vector_i = 0; index_of_space_vector_i < matrix_in_space_basis_size; ++index_of_space_vector_i) {
        auto outer_iterator = subspace.decomposition.GetNewIterator(index_of_space_vector_i);
        while (outer_iterator->hasNext()) {
            auto outer_item = outer_iterator->getNext();
            uint32_t index_of_lexicographic_vector_k = outer_item.index;
            // BUILDING k-th ROW OF INITIAL_MATRIX
            if (built_lexicografical_vectors.count(index_of_lexicographic_vector_k) == 0) {
                for (auto& term : new_operator.zero_center_terms) {
                    term->construct(matrix_in_lexicografical_basis, converter_, index_of_lexicographic_vector_k);
                }
                for (int center_a = 0; center_a < converter_.get_mults().size(); ++center_a) {
                    for (auto& term : new_operator.one_center_terms) {
                        term->construct(matrix_in_lexicografical_basis, converter_, index_of_lexicographic_vector_k, center_a);
                    }
                    for (int center_b = center_a + 1; center_b < converter_.get_mults().size(); ++center_b) {
                        for (auto& term : new_operator.two_center_terms) {
                            term->construct(matrix_in_lexicografical_basis, converter_, index_of_lexicographic_vector_k, center_a, center_b);
                        }
                    }
                }
                built_lexicografical_vectors.insert(index_of_lexicographic_vector_k);
            }
            // TODO: can we start index_of_space_vector_j from index_of_space_vector_i?
            for (uint32_t index_of_space_vector_j = 0; index_of_space_vector_j < matrix_in_space_basis_size; ++index_of_space_vector_j) {
                auto inner_iterator = subspace.decomposition.GetNewIterator(index_of_space_vector_j);
                while (inner_iterator->hasNext()) {
                    auto inner_item = inner_iterator->getNext();
                    uint32_t index_of_lexicographic_vector_l = inner_item.index;
                    double value_in_matrix_in_lexicografical_basis = matrix_in_lexicografical_basis(index_of_lexicographic_vector_k, index_of_lexicographic_vector_l);
                    if (value_in_matrix_in_lexicografical_basis != 0) {
                        submatrix.raw_data.add_to_position(outer_item.value * value_in_matrix_in_lexicografical_basis * inner_item.value,
                                                           index_of_space_vector_i, index_of_space_vector_j);
                    }
                }
            }
        }
    }
    return std::move(submatrix);
}
