#include "Subspace.h"

std::ostream &operator<<(std::ostream &os, const Subspace &subspace) {
    os << subspace.properties;
    os << subspace.decomposition;
    os << std::endl;
    return os;
}

Subspace::Subspace(UnitarySpaseMatrix && new_basis_decomposition):
decomposition(std::move(new_basis_decomposition)){
}
