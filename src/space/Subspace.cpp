#include "Subspace.h"

namespace space {
std::ostream& operator<<(std::ostream& os, const Subspace& subspace) {
    os << subspace.properties;
    os << subspace.decomposition;
    os << std::endl;
    return os;
}
}  // namespace space

space::Subspace::Subspace(UnitarySparseMatrix&& new_basis_decomposition) :
    decomposition(std::move(new_basis_decomposition)) {}
