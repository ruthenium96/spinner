#include "Submatrix.h"

#include <unordered_set>

std::ostream& operator<<(std::ostream& os, const Submatrix& submatrix) {
    os << submatrix.properties;
    os << submatrix.raw_data;
    os << std::endl;
    return os;
}
Submatrix::Submatrix(
    const Subspace& subspace,
    const model::operators::Operator& new_operator,
    const lexicographic::IndexConverter& converter) {
    lexicographic::SparseMatrix matrix_in_lexicografical_basis(converter);
    std::unordered_set<unsigned int> built_lexicografical_vectors;

    size_t matrix_in_space_basis_size = subspace.decomposition.size();
    properties = subspace.properties;
    raw_data.resize(matrix_in_space_basis_size, matrix_in_space_basis_size);

    for (uint32_t index_of_space_vector_i = 0; index_of_space_vector_i < matrix_in_space_basis_size;
         ++index_of_space_vector_i) {
        auto outer_iterator = subspace.decomposition.GetNewIterator(index_of_space_vector_i);
        while (outer_iterator->hasNext()) {
            auto outer_item = outer_iterator->getNext();
            uint32_t index_of_lexicographic_vector_k = outer_item.index;
            // BUILDING k-th ROW OF INITIAL_MATRIX
            if (built_lexicografical_vectors.count(index_of_lexicographic_vector_k) == 0) {
                for (auto& term : new_operator.getZeroCenterTerms()) {
                    term->construct(
                        matrix_in_lexicografical_basis,
                        index_of_lexicographic_vector_k);
                }
                for (int center_a = 0; center_a < converter.get_mults().size(); ++center_a) {
                    for (auto& term : new_operator.getOneCenterTerms()) {
                        term->construct(
                            matrix_in_lexicografical_basis,
                            index_of_lexicographic_vector_k,
                            center_a);
                    }
                    for (int center_b = center_a + 1; center_b < converter.get_mults().size();
                         ++center_b) {
                        for (auto& term : new_operator.getTwoCenterTerms()) {
                            term->construct(
                                matrix_in_lexicografical_basis,
                                index_of_lexicographic_vector_k,
                                center_a,
                                center_b);
                        }
                    }
                }
                built_lexicografical_vectors.insert(index_of_lexicographic_vector_k);
            }
            // TODO: can we start index_of_space_vector_j from index_of_space_vector_i?
            for (uint32_t index_of_space_vector_j = 0;
                 index_of_space_vector_j < matrix_in_space_basis_size;
                 ++index_of_space_vector_j) {
                auto inner_iterator =
                    subspace.decomposition.GetNewIterator(index_of_space_vector_j);
                while (inner_iterator->hasNext()) {
                    auto inner_item = inner_iterator->getNext();
                    uint32_t index_of_lexicographic_vector_l = inner_item.index;
                    double value_in_matrix_in_lexicografical_basis = matrix_in_lexicografical_basis(
                        index_of_lexicographic_vector_k,
                        index_of_lexicographic_vector_l);
                    if (value_in_matrix_in_lexicografical_basis != 0) {
                        raw_data.add_to_position(
                            outer_item.value * value_in_matrix_in_lexicografical_basis
                                * inner_item.value,
                            index_of_space_vector_i,
                            index_of_space_vector_j);
                    }
                }
            }
        }
    }
}
