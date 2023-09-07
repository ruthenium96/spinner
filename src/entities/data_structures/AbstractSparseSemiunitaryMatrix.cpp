#include "AbstractSparseSemiunitaryMatrix.h"

namespace quantum::linear_algebra {
void AbstractSparseSemiunitaryMatrix::unitaryTransform(
    const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrixToTransform,
    std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrixToAdd) const {
    size_t matrix_in_space_basis_size = this->size_cols();

#pragma omp parallel for shared( \
        matrix_in_space_basis_size, \
            symmetricMatrixToTransform, \
            symmetricMatrixToAdd) default(none)
    for (uint32_t index_of_space_vector_i = 0; index_of_space_vector_i < matrix_in_space_basis_size;
         ++index_of_space_vector_i) {
        auto outer_iterator = this->GetNewIterator(index_of_space_vector_i);
        while (outer_iterator->hasNext()) {
            auto outer_item = outer_iterator->getNext();
            uint32_t index_of_lexicographic_vector_k = outer_item.index;
            for (uint32_t index_of_space_vector_j = index_of_space_vector_i;
                 index_of_space_vector_j < matrix_in_space_basis_size;
                 ++index_of_space_vector_j) {
                auto inner_iterator = this->GetNewIterator(index_of_space_vector_j);
                while (inner_iterator->hasNext()) {
                    auto inner_item = inner_iterator->getNext();
                    uint32_t index_of_lexicographic_vector_l = inner_item.index;
                    double value_in_matrix_in_lexicografical_basis = symmetricMatrixToTransform->at(
                        index_of_lexicographic_vector_k,
                        index_of_lexicographic_vector_l);
                    if (value_in_matrix_in_lexicografical_basis != 0) {
                        // \Delta B_{ij} = \sum_{kl} U_{ki} * \Delta A_{kl} * U_{lj}
                        symmetricMatrixToAdd->add_to_position(
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
}  // namespace quantum::linear_algebra