#include "ConstantTerm.h"

#include <utility>

namespace model::operators {
void ConstantTerm::construct(
    std::unique_ptr<quantum::linear_algebra::AbstractSymmetricMatrix>&
        matrix_in_lexicografical_basis,
    uint32_t index_of_vector) const {
    matrix_in_lexicografical_basis->add_to_position(*constant_, index_of_vector, index_of_vector);
}

ConstantTerm::ConstantTerm(std::shared_ptr<const double> constant) :
    constant_(std::move(constant)) {}

std::unique_ptr<ZeroCenterTerm> ConstantTerm::clone() const {
    return std::make_unique<ConstantTerm>(constant_);
}
}  // namespace model::operators