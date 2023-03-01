#include "Subspace.h"

namespace space {
std::ostream& operator<<(std::ostream& os, const Subspace& subspace) {
    os << subspace.properties;
    subspace.decomposition->print(os);
    os << std::endl;
    return os;
}

Subspace::Subspace(std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>&&
                       new_basis_decomposition) :
    decomposition(std::move(new_basis_decomposition)) {}
}  // namespace space
