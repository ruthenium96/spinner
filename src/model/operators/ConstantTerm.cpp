#include "ConstantTerm.h"

namespace model::operators {
void ConstantTerm::construct(
    UnitarySparseMatrix& matrix_in_lexicografical_basis,
    uint32_t index_of_vector) const {
    matrix_in_lexicografical_basis.add_to_position(constant_, index_of_vector, index_of_vector);
}

ConstantTerm::ConstantTerm(double constant) : constant_(constant) {}

std::unique_ptr<ZeroCenterTerm> ConstantTerm::clone() const {
    return std::make_unique<ConstantTerm>(constant_);
}
}  // namespace model::operators