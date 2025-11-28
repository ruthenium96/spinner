#include "Subspace.h"

namespace space {
Subspace::Subspace(
    std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>&&
        new_basis_decomposition):
    decomposition(std::move(new_basis_decomposition)) {}

uint32_t Subspace::size() const {
    return decomposition->size_cols();
}

}  // namespace space
