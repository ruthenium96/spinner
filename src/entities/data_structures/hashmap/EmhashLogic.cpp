#include "EmhashLogic.h"

#include <hash_table8.hpp>

namespace quantum::linear_algebra {
void EmhashLogic::unitaryTransform(
    const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrixToTransform,
    std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrixToAdd,
    const AbstractSparseSemiunitaryMatrix& unitaryMatrix) const {
    const auto maybeSemiunitaryMatrix =
        dynamic_cast<const EmhashSparseSemiunitaryMatrix*>(&unitaryMatrix);
    if (maybeSemiunitaryMatrix == nullptr) {
        throw std::bad_cast();
    }
    const auto& semiunitaryMatrixData = maybeSemiunitaryMatrix->getSparseSemiunitaryMatrix();
    size_t matrix_in_space_basis_size = maybeSemiunitaryMatrix->size_cols();

    auto maybeSymmetricMatrix =
        dynamic_cast<const EmhashSparseSymmetricMatrix*>(symmetricMatrixToTransform.get());
    if (maybeSymmetricMatrix == nullptr) {
        throw std::bad_cast();
    }
    const auto& symmetricMatrixData = maybeSymmetricMatrix->getSparseSymmetricMatrix();

    emhash8::HashMap<uint32_t, EmhashSparseSemiunitaryMatrix::Map> first_mult;

    // cannot easily use pragma omp parallel: race condition inside first_mult hashmap.
    for (uint32_t index_of_space_vector_i = 0; index_of_space_vector_i < matrix_in_space_basis_size;
         ++index_of_space_vector_i) {
        for (const auto& i_iter : semiunitaryMatrixData[index_of_space_vector_i]) {
            uint32_t index_of_lexicographic_vector_k = i_iter.first;
            double coefficient_in_unitary_matrix = i_iter.second;
            if (symmetricMatrixData.contains(index_of_lexicographic_vector_k)) {
                const auto& k_vector = symmetricMatrixData.at(index_of_lexicographic_vector_k);
                for (const auto& k_iter : k_vector) {
                    uint32_t index_of_lexicographic_vector_l = k_iter.first;
                    double value_in_matrix_in_lexicografical_basis = k_iter.second;
                    double product_for_first_mult =
                        value_in_matrix_in_lexicografical_basis * coefficient_in_unitary_matrix;
                    // \Delta X_{il} = \sum_{k} U_{ki} * \Delta A_{kl}
                    first_mult[index_of_lexicographic_vector_l][index_of_space_vector_i] +=
                        product_for_first_mult;
                }
            }
        }
    }

#pragma omp parallel for shared( \
        matrix_in_space_basis_size, \
            first_mult, \
            symmetricMatrixToAdd, \
            semiunitaryMatrixData) default(none)
    for (uint32_t index_of_space_vector_j = 0; index_of_space_vector_j < matrix_in_space_basis_size;
         ++index_of_space_vector_j) {
        for (const auto& j_iter : semiunitaryMatrixData[index_of_space_vector_j]) {
            uint32_t index_of_lexicographic_vector_l = j_iter.first;
            double coefficient_in_unitary_matrix = j_iter.second;
            if (first_mult.contains(index_of_lexicographic_vector_l)) {
                const auto& l_vector = first_mult.at(index_of_lexicographic_vector_l);
                for (const auto& l_iter : l_vector) {
                    uint32_t index_of_space_vector_i = l_iter.first;
                    if (index_of_space_vector_i >= index_of_space_vector_j) {
                        double value_in_matrix_of_first_mult = l_iter.second;
                        double product_for_final_mult =
                            value_in_matrix_of_first_mult * coefficient_in_unitary_matrix;
                        // \Delta B_{ij} = \sum_{l} \Delta X_{il} * U_{lj}
                        symmetricMatrixToAdd->add_to_position(
                            product_for_final_mult,
                            index_of_space_vector_i,
                            index_of_space_vector_j);
                    }
                }
            }
        }
    }
}
}  // namespace quantum::linear_algebra