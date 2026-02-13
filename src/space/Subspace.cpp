#include "Subspace.h"

namespace spinner::space {
Subspace::Subspace(
    std::unique_ptr<linalg_structures::AbstractSparseSemiunitaryMatrix>&&
        new_basis_decomposition):
    decomposition(std::move(new_basis_decomposition)) {}

uint32_t Subspace::size() const {
    return decomposition->size_cols();
}

} // namespace spinner::space
