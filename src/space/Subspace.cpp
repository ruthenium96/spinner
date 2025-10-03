#include "Subspace.h"

namespace space {
Subspace::Subspace(
    std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>&&
        new_basis_decomposition,
    std::optional<std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>&&
        dense_semiunitary_matrix) :
    decomposition(std::move(new_basis_decomposition)),
    dense_semiunitary_matrix(std::move(dense_semiunitary_matrix)) {}

uint32_t Subspace::size() const {
    if (dense_semiunitary_matrix.has_value()) {
        return dense_semiunitary_matrix.value()->size_cols();
    } else {
        return decomposition->size_cols();
    }
}

}  // namespace space
